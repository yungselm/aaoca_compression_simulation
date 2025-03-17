mod io;
mod contour;
mod processing;
mod utils;
mod texture;
mod comparison;
mod centerline_alignment;

use std::fs::File;
use std::io::Write;
use std::path::Path;
use std::error::Error;

use io::{read_contour_data, create_catheter_points, write_obj_mesh};
use contour::create_contours;
use processing::{compute_centroid, translate_contour};
use utils::{trim_to_same_length, smooth_contours};
use texture::{compute_uv_coordinates, compute_displacements, create_displacement_texture, create_black_texture};
use comparison::process_phase_comparison;
use centerline_alignment::test_function;

fn main() -> Result<(), Box<dyn Error>> {
    // Process both "rest" and "stress" cases.
    test_function()?;
    Ok(())
    // process_case("rest", "input/rest_csv_files", "output/rest")?;
    // process_case("stress", "input/stress_csv_files", "output/stress")?;
    // process_phase_comparison(
    //     "diastolic",
    //     "input/rest_csv_files",
    //     "input/stress_csv_files",
    //     "output/diastole_comparison"
    // )?;

    // process_phase_comparison(
    //     "systolic",
    //     "input/rest_csv_files",
    //     "input/stress_csv_files",
    //     "output/systole_comparison"
    // )?;
    // // align_meshes_to_centerline(
    // //     "resampled_centerline.txt",
    // //     "output/stress",
    // //     "output/aligned_stress"
    // // )?;
    // Ok(())
}

/// Processes a given case by reading diastolic and systolic contours, aligning them,
/// computing displacements and UV coordinates, and finally writing out the OBJ, MTL, and texture files.
fn process_case(case_name: &str, input_dir: &str, output_dir: &str) -> Result<(), Box<dyn Error>> {
    // Create the output directory if it doesn't exist.
    std::fs::create_dir_all(output_dir)?;

    // === Process Diastolic Contours ===
    println!("--- Processing {} Diastole ---", case_name);
    let diastole_path = Path::new(input_dir).join("diastolic_contours.csv");
    let diastole_points = read_contour_data(diastole_path.to_str().unwrap())?;
    let diastole_catheter = create_catheter_points(&diastole_points);

    let mut diastole_contours = create_contours(diastole_points);
    let mut diastole_catheter_contours = create_contours(diastole_catheter);

    // === Process Systolic Contours ===
    println!("--- Processing {} Systole ---", case_name);
    let systole_path = Path::new(input_dir).join("systolic_contours.csv");
    let systole_points = read_contour_data(systole_path.to_str().unwrap())?;
    let systole_catheter = create_catheter_points(&systole_points);

    let mut systole_contours = create_contours(systole_points);
    let mut systole_catheter_contours = create_contours(systole_catheter);

    // Align systolic contours to the diastolic reference.
    let diastolic_ref_centroid = compute_centroid(&diastole_contours[0].1);
    for (_, ref mut contour) in systole_contours
            .iter_mut()
            .chain(systole_catheter_contours.iter_mut()) {
        let systolic_centroid = compute_centroid(contour);
        let translation = (
            diastolic_ref_centroid.0 - systolic_centroid.0,
            diastolic_ref_centroid.1 - systolic_centroid.1,
        );
        translate_contour(contour, translation);
    }

    // Adjust the z-coordinates of systolic contours.
    let z_translation = diastole_contours[0].1[0].z - systole_contours[0].1[0].z;
    for (_, ref mut contour) in systole_contours
        .iter_mut()
        .chain(systole_catheter_contours.iter_mut()) {
        for point in contour {
            point.z += z_translation;
        }
    }

    // Ensure both contour sets have the same number of contours.
    trim_to_same_length(&mut diastole_contours, &mut systole_contours);
    trim_to_same_length(&mut diastole_catheter_contours, &mut systole_catheter_contours);

    // Smooth contours.
    diastole_contours = smooth_contours(&mut diastole_contours);
    systole_contours = smooth_contours(&mut systole_contours);

    // Interpolate between diastole and systole contours.
    let steps = 30; // Number of interpolation steps
    let interpolated_meshes = processing::interpolate_contours(&diastole_contours, &systole_contours, steps)?;
    let interpolated_catheter_meshes = processing::interpolate_contours(&diastole_catheter_contours, &systole_catheter_contours, steps)?;

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
        .filter_map(|(i, mesh)| {
            if i == 1 || i == 2 { None } else { Some(mesh) }
        })
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
        .filter_map(|(i, mesh)| {
            if i == 1 || i == 2 { None } else { Some(mesh) }
        })
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
        write_obj_mesh(
            mesh,
            &uv_coords,
            obj_path.to_str().unwrap(),
            &mtl_filename,
        )?;
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
        create_black_texture(texture_width, texture_height, texture_path.to_str().unwrap())?;

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
        write_obj_mesh(
            mesh,
            &uv_coords,
            obj_path.to_str().unwrap(),
            &mtl_filename,
        )?;
    }
Ok(())
}
