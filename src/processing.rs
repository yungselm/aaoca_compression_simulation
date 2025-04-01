use rayon::prelude::*;
use std::error::Error;
use std::fs::File;
use std::io::Write;
use std::path::Path;

use crate::contour::Contour;
use crate::io::write_obj_mesh;
use crate::io::ContourPoint;
use crate::texture::{
    compute_displacements, compute_uv_coordinates, create_black_texture,
    create_displacement_texture,
};
use crate::utils::{smooth_contours, trim_to_same_length};

/// Processes a given case by reading diastolic and systolic contours, aligning them,
/// computing displacements and UV coordinates, and finally writing out the OBJ, MTL, and texture files.
/// Additionally it can be specified how many interpolation steps should be used to generate the final meshes
/// used for the animation in blender.
pub fn process_case(
    case_name: &str,
    input_dir: &str,
    output_dir: &str,
    interpolation_steps: usize,
    steps_best_rotation: usize,
    range_rotation_rad: f64,
) -> Result<(), Box<dyn Error>> {
    // Create the output directory if it doesn't exist.
    std::fs::create_dir_all(output_dir)?;

    // === Process Diastolic Contours ===
    println!("--- Processing {} Diastole ---", case_name);
    let diastole_path = Path::new(input_dir).join("diastolic_contours.csv");
    let diastole_points = ContourPoint::read_contour_data(diastole_path.to_str().unwrap())?;
    let diastole_catheter = ContourPoint::create_catheter_points(&diastole_points);

    let mut diastole_contours =
        Contour::create_contours(diastole_points, steps_best_rotation, range_rotation_rad);
    let mut diastole_catheter_contours =
        Contour::create_contours(diastole_catheter, steps_best_rotation, range_rotation_rad);

    // === Process Systolic Contours ===
    println!("--- Processing {} Systole ---", case_name);
    let systole_path = Path::new(input_dir).join("systolic_contours.csv");
    let systole_points = ContourPoint::read_contour_data(systole_path.to_str().unwrap())?;
    let systole_catheter = ContourPoint::create_catheter_points(&systole_points);

    let mut systole_contours =
        Contour::create_contours(systole_points, steps_best_rotation, range_rotation_rad);
    let mut systole_catheter_contours =
        Contour::create_contours(systole_catheter, steps_best_rotation, range_rotation_rad);

    // Align systolic contours to the diastolic reference.
    let diastolic_ref_centroid = Contour::compute_centroid(&diastole_contours[0].1);
    for (_, ref mut contour) in systole_contours
        .iter_mut()
        .chain(systole_catheter_contours.iter_mut())
    {
        let systolic_centroid = Contour::compute_centroid(contour);
        let translation = (
            diastolic_ref_centroid.0 - systolic_centroid.0,
            diastolic_ref_centroid.1 - systolic_centroid.1,
        );
        Contour::translate_contour(contour, translation);
    }

    // Adjust the z-coordinates of systolic contours.
    let z_translation = diastole_contours[0].1[0].z - systole_contours[0].1[0].z;
    for (_, ref mut contour) in systole_contours
        .iter_mut()
        .chain(systole_catheter_contours.iter_mut())
    {
        for point in contour {
            point.z += z_translation;
        }
    }
    
    // Ensure both contour sets have the same number of contours.
    trim_to_same_length(&mut diastole_contours, &mut systole_contours);
    trim_to_same_length(
        &mut diastole_catheter_contours,
        &mut systole_catheter_contours,
    );
    
    // Smooth contours.
    diastole_contours = smooth_contours(&mut diastole_contours);
    systole_contours = smooth_contours(&mut systole_contours);
    
    // Compute a global best rotation for systolic contours (full stack) to best match diastolic contours.
    // This rotates the entire systole stack (all frames) around the z-axis.
    let best_rotation_angle = find_best_rotation_all(
        &diastole_contours,
        &systole_contours,
        steps_best_rotation,   // number of candidate steps (e.g. 200 or 400)
        range_rotation_rad,    // rotation range (e.g. 1.05 for ~±60°)
    );

    println!("Rotating full systole stack by {:.3} rad to best match diastole.", best_rotation_angle);

    // Apply the same rotation to every point in the systolic contours (and catheter contours if desired)
    for (_, contour) in systole_contours.iter_mut().chain(systole_catheter_contours.iter_mut()) {
        for point in contour.iter_mut() {
            let x_new = point.x * best_rotation_angle.cos() - point.y * best_rotation_angle.sin();
            let y_new = point.x * best_rotation_angle.sin() + point.y * best_rotation_angle.cos();
            point.x = x_new;
            point.y = y_new;
        }
    }

    // Interpolate between diastole and systole contours.
    let steps = interpolation_steps; // Number of interpolation steps
    let interpolated_meshes =
        Contour::interpolate_contours(&diastole_contours, &systole_contours, steps)?;
    let interpolated_catheter_meshes = Contour::interpolate_contours(
        &diastole_catheter_contours,
        &systole_catheter_contours,
        steps,
    )?;

    // === Build the Meshes to Process ===
    // Include diastole, systole, and interpolated meshes.
    let all_contour_meshes = std::iter::once(&diastole_contours[..])
        .chain(std::iter::once(&systole_contours[..]))
        .chain(interpolated_meshes.iter().map(|m| &m[..]))
        .collect::<Vec<_>>();

    // STUPID FIX since these two are always wrong: Filter out the ones with indices 1 and 2.
    let all_contour_meshes: Vec<_> = all_contour_meshes
        .into_iter()
        .enumerate()
        .filter_map(|(i, mesh)| if i == 1 || i == 2 { None } else { Some(mesh) })
        .collect();

    // === Compute Maximum Displacement for Normalization ===
    let mut max_disp: f32 = 0.0;
    for mesh in &all_contour_meshes {
        let displacements = compute_displacements(mesh, &diastole_contours);
        println!("max_disp: {:?}", max_disp);
        max_disp = displacements.iter().fold(max_disp, |a, &b| a.max(b));
    }

    let all_catheter_meshes = std::iter::once(&diastole_catheter_contours[..])
        .chain(std::iter::once(&systole_catheter_contours[..]))
        .chain(interpolated_catheter_meshes.iter().map(|m| &m[..]))
        .collect::<Vec<_>>();

    // STUPID FIX since these two are always wrong: Filter out the ones with indices 1 and 2.
    let all_catheter_meshes: Vec<_> = all_catheter_meshes
        .into_iter()
        .enumerate()
        .filter_map(|(i, mesh)| if i == 1 || i == 2 { None } else { Some(mesh) })
        .collect();

    // --- Process Regular (Contour) Meshes ---
    for (i, mesh) in all_contour_meshes.into_iter().enumerate() {
        // Compute UV coordinates.
        let uv_coords = compute_uv_coordinates(mesh);

        // Determine texture dimensions.
        let texture_height = mesh.len() as u32;
        let texture_width = if texture_height > 0 {
            mesh[0].1.len() as u32
        } else {
            0
        };

        // Compute displacement values.
        let displacements = compute_displacements(mesh, &diastole_contours);

        // Save the displacement texture.
        let tex_filename = format!("mesh_{:03}_{}.png", i, case_name);
        let texture_path = Path::new(output_dir).join(&tex_filename);
        create_displacement_texture(
            &displacements,
            texture_width,
            texture_height,
            max_disp,
            texture_path.to_str().unwrap(),
        )?;

        // Write the material file (MTL).
        let mtl_filename = format!("mesh_{:03}_{}.mtl", i, case_name);
        let mtl_path = Path::new(output_dir).join(&mtl_filename);
        let mut mtl_file = File::create(&mtl_path)?;
        writeln!(
            mtl_file,
            "newmtl displacement_material\nKa 1 1 1\nKd 1 1 1\nmap_Kd {}",
            tex_filename
        )?;

        // Write the OBJ file.
        let obj_filename = format!("mesh_{:03}_{}.obj", i, case_name);
        let obj_path = Path::new(output_dir).join(&obj_filename);
        write_obj_mesh(mesh, &uv_coords, obj_path.to_str().unwrap(), &mtl_filename)?;
    }

    // --- Process Catheter Meshes ---
    for (i, mesh) in all_catheter_meshes.into_iter().enumerate() {
        // Determine texture dimensions.
        let texture_height = mesh.len() as u32;
        let texture_width = if texture_height > 0 {
            mesh[0].1.len() as u32
        } else {
            0
        };

        // For catheter meshes we generate a fixed (black) texture.
        let tex_filename = format!("catheter_{:03}_{}.png", i, case_name);
        let texture_path = Path::new(output_dir).join(&tex_filename);
        create_black_texture(
            texture_width,
            texture_height,
            texture_path.to_str().unwrap(),
        )?;

        // Write the material file (MTL) for catheter mesh.
        let mtl_filename = format!("catheter_{:03}_{}.mtl", i, case_name);
        let mtl_path = Path::new(output_dir).join(&mtl_filename);
        let mut mtl_file = File::create(&mtl_path)?;
        // Set both ambient and diffuse to black.
        writeln!(
            mtl_file,
            "newmtl black_material\nKa 0 0 0\nKd 0 0 0\nmap_Kd {}",
            tex_filename
        )?;

        // Compute UV coordinates (or use dummy coordinates if preferred).
        let uv_coords = compute_uv_coordinates(mesh);

        // Write the OBJ file.
        let obj_filename = format!("catheter_{:03}_{}.obj", i, case_name);
        let obj_path = Path::new(output_dir).join(&obj_filename);
        write_obj_mesh(mesh, &uv_coords, obj_path.to_str().unwrap(), &mtl_filename)?;
    }
    Ok(())
}

// pub fn find_best_rotation_all(
//     diastole_contours: &[(u32, Vec<ContourPoint>)],
//     systole_contours: &[(u32, Vec<ContourPoint>)],
//     steps: usize,
//     range: f64,
// ) -> f64 {
//     let steps = steps;
//     let range = range;
//     let increment = (2.0 * range) / steps as f64;

//     (0..=steps)
//         .into_par_iter()
//         .map(|i| {
//             let angle = -range + i as f64 * increment;
//             let total_distance: f64 = diastole_contours
//                 .par_iter()
//                 .zip(systole_contours.par_iter())
//                 .map(|((d_id, d_contour), (s_id, s_contour))| {
//                     assert_eq!(d_id, s_id, "Mismatched contour IDs");
//                     // Rotate each point in systole contour around reference centroid
//                     let rotated_s_contour: Vec<ContourPoint> = s_contour
//                         .iter()
//                         .map(|p| {
//                             let x = p.x * angle.cos() - p.y * angle.sin();
//                             let y = p.x * angle.sin() + p.y * angle.cos();
//                             ContourPoint { x, y, ..*p }
//                         })
//                         .collect();
//                     // Compute Hausdorff distance between corresponding contours
//                     Contour::hausdorff_distance(d_contour, &rotated_s_contour)
//                 })
//                 .sum();
//             let avg_distance = total_distance / diastole_contours.len() as f64;
//             (angle, avg_distance);
//             println!(
//                 "Angle: {:.3} rad, Avg Distance: {:.3}",
//                 angle, avg_distance
//             );
//             (angle, avg_distance)
//         })
//         .min_by(|a, b| {
//             a.1.partial_cmp(&b.1)
//                 .unwrap_or(std::cmp::Ordering::Equal)
//         })
//         .map(|(angle, _)| angle)
//         .unwrap_or(0.0)
// }

/// Best rotation angles per stack, using weights for lower
pub fn find_best_rotation_all(
    diastole_contours: &[(u32, Vec<ContourPoint>)],
    systole_contours: &[(u32, Vec<ContourPoint>)],
    steps: usize,
    range: f64,
) -> f64 {
    let steps = steps;
    let range = range;
    let increment = (2.0 * range) / steps as f64;

    (0..=steps)
        .into_par_iter()
        .map(|i| {
            let angle = -range + i as f64 * increment;
            // For each angle, compute the weighted total distance and the sum of weights.
            let (weighted_total, total_weight): (f64, f64) = diastole_contours
                .par_iter()
                .zip(systole_contours.par_iter())
                .map(|((d_id, d_contour), (s_id, s_contour))| {
                    assert_eq!(d_id, s_id, "Mismatched contour IDs");
                    // Define a weight based on the contour id.
                    // For example, using 1/(id+1) gives lower ids higher weight.
                    let weight = 1.0 / ((*d_id as f64) + 1.0);
                    
                    // Rotate each point in the systole contour.
                    let rotated_s_contour: Vec<ContourPoint> = s_contour
                        .iter()
                        .map(|p| {
                            let x = p.x * angle.cos() - p.y * angle.sin();
                            let y = p.x * angle.sin() + p.y * angle.cos();
                            ContourPoint { x, y, ..*p }
                        })
                        .collect();
                    
                    // Compute the Hausdorff distance between corresponding contours.
                    let distance = Contour::hausdorff_distance(d_contour, &rotated_s_contour);
                    
                    (weight * distance, weight)
                })
                .reduce(
                    || (0.0, 0.0),
                    |(w_total_a, total_weight_a), (w_total_b, total_weight_b)| {
                        (w_total_a + w_total_b, total_weight_a + total_weight_b)
                    },
                );
            
            // Calculate the weighted average distance.
            let avg_distance = if total_weight > 0.0 {
                weighted_total / total_weight
            } else {
                0.0
            };

            println!(
                "Angle: {:.3} rad, Weighted Avg Distance: {:.3}",
                angle, avg_distance
            );
            (angle, avg_distance)
        })
        .min_by(|a, b| {
            a.1.partial_cmp(&b.1)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .map(|(angle, _)| angle)
        .unwrap_or(0.0)
}
