use crate::{
    contour::Contour,
    io::ContourPoint,
    io::{read_obj_mesh, write_obj_mesh},
    texture::{
        compute_displacements, compute_uv_coordinates, create_black_texture,
        create_displacement_texture,
    },
    utils::{smooth_contours, trim_to_same_length},
    processing::find_best_rotation_all,
};
use std::error::Error;
use std::fs::File;
use std::io::Write;
use std::path::Path;

pub fn process_phase_comparison(
    phase_name: &str,
    rest_input_dir: &str,
    stress_input_dir: &str,
    output_dir: &str,
    interpolation_steps: usize,
    steps_best_rotation: usize,
    range_rotation_rad: f64,
) -> Result<(), Box<dyn Error>> {
    let (rest_path, stress_path, rest_catheter_path, stress_catheter_path) =
        if phase_name == "diastolic" {
            (
                format!("{}/mesh_000_rest.obj", rest_input_dir),
                format!("{}/mesh_000_stress.obj", stress_input_dir),
                format!("{}/catheter_000_rest.obj", rest_input_dir),
                format!("{}/catheter_000_stress.obj", stress_input_dir),
            )
        } else {
            (
                format!("{}/mesh_029_rest.obj", rest_input_dir),
                format!("{}/mesh_029_stress.obj", stress_input_dir),
                format!("{}/catheter_029_rest.obj", rest_input_dir),
                format!("{}/catheter_029_stress.obj", stress_input_dir),
            )
        };

    let mut rest_contours = read_obj_mesh(&rest_path)?;
    let mut stress_contours = read_obj_mesh(&stress_path)?;

    let mut rest_catheter_contours = read_obj_mesh(&rest_catheter_path)?;
    let mut stress_catheter_contours = read_obj_mesh(&stress_catheter_path)?;

    // Replace the rest contours with new contours where the x and y coordinates
    // are interpolated between rest and stress, and the z is taken from stress.
    rest_contours = resample_contours_with_reference_z(&rest_contours, &stress_contours);
    rest_catheter_contours =
        resample_contours_with_reference_z(&rest_catheter_contours, &stress_catheter_contours);

    // Align stress contours to rest reference
    let rest_ref_centroid = Contour::compute_centroid(&rest_contours[0].1);
    for (_, contour) in stress_contours
        .iter_mut()
        .chain(stress_catheter_contours.iter_mut())
    {
        let stress_centroid = Contour::compute_centroid(contour);
        let translation = (
            rest_ref_centroid.0 - stress_centroid.0,
            rest_ref_centroid.1 - stress_centroid.1,
        );
        Contour::translate_contour(contour, translation);
    }

    // Z-axis alignment
    let z_translation = rest_contours[0].1[0].z - stress_contours[0].1[0].z;
    for (_, contour) in stress_contours
        .iter_mut()
        .chain(stress_catheter_contours.iter_mut())
    {
        for point in contour {
            point.z += z_translation;
        }
    }

    // Common processing
    trim_to_same_length(&mut rest_contours, &mut stress_contours);
    trim_to_same_length(&mut rest_catheter_contours, &mut stress_catheter_contours);

    rest_contours = smooth_contours(&rest_contours);
    stress_contours = smooth_contours(&stress_contours);

    // Compute a global best rotation for systolic contours (full stack) to best match diastolic contours.
    // This rotates the entire systole stack (all frames) around the z-axis.
    let best_rotation_angle = find_best_rotation_all(
        &rest_contours,
        &stress_contours,
        steps_best_rotation,   // number of candidate steps (e.g. 200 or 400)
        range_rotation_rad,    // rotation range (e.g. 1.05 for ~±60°)
    );

    println!("Rotating full systole stack by {:.3} rad to best match diastole.", best_rotation_angle);

    // Apply the same rotation to every point in the systolic contours (and catheter contours if desired)
    for (_, contour) in stress_contours.iter_mut().chain(stress_catheter_contours.iter_mut()) {
        for point in contour.iter_mut() {
            let x_new = point.x * best_rotation_angle.cos() - point.y * best_rotation_angle.sin();
            let y_new = point.x * best_rotation_angle.sin() + point.y * best_rotation_angle.cos();
            point.x = x_new;
            point.y = y_new;
        }
    }   

    // Interpolation between rest and stress
    let steps = interpolation_steps;
    let interpolated_meshes =
        Contour::interpolate_contours(&rest_contours, &stress_contours, steps)?;
    let interpolated_catheter_meshes =
        Contour::interpolate_contours(&rest_catheter_contours, &stress_catheter_contours, steps)?;

    // === Build the Meshes to Process ===
    // Include diastole, systole, and interpolated meshes.
    let all_contour_meshes = std::iter::once(&rest_contours[..])
        .chain(std::iter::once(&stress_contours[..]))
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
        let displacements = compute_displacements(mesh, &rest_contours);
        println!("max_disp: {:?}", max_disp);
        max_disp = displacements.iter().fold(max_disp, |a, &b| a.max(b));
    }

    let all_catheter_meshes = std::iter::once(&rest_catheter_contours[..])
        .chain(std::iter::once(&stress_catheter_contours[..]))
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
        let displacements = compute_displacements(mesh, &rest_contours);

        // Save the displacement texture.
        let tex_filename = format!("mesh_{:03}_{}.png", i, phase_name);
        let texture_path = Path::new(output_dir).join(&tex_filename);
        create_displacement_texture(
            &displacements,
            texture_width,
            texture_height,
            max_disp,
            texture_path.to_str().unwrap(),
        )?;

        // Write the material file (MTL).
        let mtl_filename = format!("mesh_{:03}_{}.mtl", i, phase_name);
        let mtl_path = Path::new(output_dir).join(&mtl_filename);
        let mut mtl_file = File::create(&mtl_path)?;
        writeln!(
            mtl_file,
            "newmtl displacement_material\nKa 1 1 1\nKd 1 1 1\nmap_Kd {}",
            tex_filename
        )?;

        // Write the OBJ file.
        let obj_filename = format!("mesh_{:03}_{}.obj", i, phase_name);
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
        let tex_filename = format!("catheter_{:03}_{}.png", i, phase_name);
        let texture_path = Path::new(output_dir).join(&tex_filename);
        create_black_texture(
            texture_width,
            texture_height,
            texture_path.to_str().unwrap(),
        )?;

        // Write the material file (MTL) for catheter mesh.
        let mtl_filename = format!("catheter_{:03}_{}.mtl", i, phase_name);
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
        let obj_filename = format!("catheter_{:03}_{}.obj", i, phase_name);
        let obj_path = Path::new(output_dir).join(&obj_filename);
        write_obj_mesh(mesh, &uv_coords, obj_path.to_str().unwrap(), &mtl_filename)?;
    }
    Ok(())
}

/// Resamples a stack of contours so that the new z‑spacing exactly matches the spacing
/// from the reference (stress) contours, while preserving the original z‑direction.
///
/// Assumptions:
///   - `original` contours (e.g. rest) are sorted in order (either increasing or decreasing in z).
///   - Each contour’s z value is constant (all points in a given contour share the same z).
///   - There are at least two contours in both `original` and `reference`.
fn resample_contours_with_reference_z(
    original: &[(u32, Vec<ContourPoint>)],
    reference: &[(u32, Vec<ContourPoint>)],
) -> Vec<(u32, Vec<ContourPoint>)> {
    let n = original.len();
    if n < 2 || reference.len() < 2 {
        return original.to_vec();
    }

    // Determine the z-direction by comparing the first and last contour's z.
    let start_z = original[0].1[0].z;
    let end_z = original[n - 1].1[0].z;
    let direction = if end_z >= start_z { 1.0 } else { -1.0 };

    // Compute the cumulative z distance (absolute differences) along the original stack.
    let mut cum_z = Vec::with_capacity(n);
    let mut total = 0.0;
    cum_z.push(0.0); // starting cumulative distance is 0.
    for i in 1..n {
        let prev_z = original[i - 1].1[0].z;
        let curr_z = original[i].1[0].z;
        total += (curr_z - prev_z).abs();
        cum_z.push(total);
    }

    // Determine target spacing from the reference contours.
    // We use the absolute difference between the first two reference contours.
    let target_spacing = (reference[1].1[0].z - reference[0].1[0].z).abs();
    // Determine the number of new frames to generate.
    let new_frame_count = (total / target_spacing).floor() as usize + 1;

    // Generate new contours.
    let mut new_contours = Vec::with_capacity(new_frame_count);
    for j in 0..new_frame_count {
        // Compute the target cumulative distance and corresponding target z.
        let target_cum = j as f64 * target_spacing;
        let target_z = start_z + direction * target_cum;

        // Find the bracket in the original cum_z where target_cum falls.
        let mut i = 0;
        while i < cum_z.len() - 1 && cum_z[i + 1] < target_cum {
            i += 1;
        }
        if i >= n - 1 {
            // If beyond range, simply copy the last contour and assign target_z.
            let new_points = original[n - 1]
                .1
                .iter()
                .map(|pt| ContourPoint {
                    frame_index: j as u32,
                    x: pt.x,
                    y: pt.y,
                    z: target_z,
                })
                .collect();
            new_contours.push((j as u32, new_points));
        } else {
            // Interpolate between original[i] and original[i+1].
            let z0 = cum_z[i];
            let z1 = cum_z[i + 1];
            let t = if (z1 - z0).abs() < std::f64::EPSILON {
                0.0
            } else {
                (target_cum - z0) / (z1 - z0)
            };

            // For each corresponding point (assuming same number per contour),
            // linearly interpolate x and y, and force z to target_z.
            let new_points: Vec<ContourPoint> = original[i]
                .1
                .iter()
                .zip(original[i + 1].1.iter())
                .map(|(pt0, pt1)| ContourPoint {
                    frame_index: j as u32,
                    x: pt0.x * (1.0 - t) + pt1.x * t,
                    y: pt0.y * (1.0 - t) + pt1.y * t,
                    z: target_z,
                })
                .collect();
            new_contours.push((j as u32, new_points));
        }
    }
    new_contours
}
