mod io;
mod contour;
mod processing;
mod utils;
mod texture;

use std::fs::File;
use std::io::Write;
use std::path::Path;
use std::error::Error;

use io::{read_contour_data, write_obj_mesh};
use contour::create_contours;
use processing::{compute_centroid, translate_contour};
use utils::trim_to_same_length;
use texture::{compute_uv_coordinates, compute_displacements, create_displacement_texture};

fn main() -> Result<(), Box<dyn Error>> {
    // Process both "rest" and "stress" cases.
    process_case("rest", "input/rest_csv_files", "output/rest")?;
    process_case("stress", "input/stress_csv_files", "output/stress")?;
    Ok(())
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
    let mut diastole_contours = create_contours(diastole_points);

    // === Process Systolic Contours ===
    println!("--- Processing {} Systole ---", case_name);
    let systole_path = Path::new(input_dir).join("systolic_contours.csv");
    let systole_points = read_contour_data(systole_path.to_str().unwrap())?;
    let mut systole_contours = create_contours(systole_points);

    // Align systolic contours to the diastolic reference.
    let diastolic_ref_centroid = compute_centroid(&diastole_contours[0].1);
    for (_, ref mut contour) in systole_contours.iter_mut() {
        let systolic_centroid = compute_centroid(contour);
        let translation = (
            diastolic_ref_centroid.0 - systolic_centroid.0,
            diastolic_ref_centroid.1 - systolic_centroid.1,
        );
        translate_contour(contour, translation);
    }

    // Adjust the z-coordinates of systolic contours.
    let z_translation = diastole_contours[0].1[0].z - systole_contours[0].1[0].z;
    for (_, ref mut contour) in systole_contours.iter_mut() {
        for point in contour {
            point.z += z_translation;
        }
    }

    // Ensure both contour sets have the same number of contours.
    trim_to_same_length(&mut diastole_contours, &mut systole_contours);

    // Interpolate between diastole and systole contours.
    let steps = 30; // Number of interpolation steps
    let interpolated_meshes = processing::interpolate_contours(&diastole_contours, &systole_contours, steps)?;

    // === Build the Meshes to Process ===
    // Include diastole, systole, and interpolated meshes.
    let all_meshes = std::iter::once(&diastole_contours[..])
        .chain(std::iter::once(&systole_contours[..]))
        .chain(interpolated_meshes.iter().map(|m| &m[..]))
        .collect::<Vec<_>>();

    // === Compute Maximum Displacement for Normalization ===
    let mut max_disp: f32 = 0.0;
    for mesh in &all_meshes {
        let displacements = compute_displacements(mesh, &diastole_contours);
        println!("max_disp: {:?}", max_disp);
        max_disp = displacements.iter().fold(max_disp, |a, &b| a.max(b));
    }

    // === Process Each Mesh ===
    for (i, mesh) in all_meshes.into_iter().enumerate() {
        // Compute UV coordinates for the current mesh.
        let uv_coords = compute_uv_coordinates(mesh);

        // Determine texture dimensions:
        // - The texture height equals the number of contours.
        // - The texture width equals the number of points per contour.
        let texture_height = mesh.len() as u32;
        let texture_width = if texture_height > 0 {
            mesh[0].1.len() as u32
        } else {
            0
        };

        // Compute the displacement values.
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

        // Write the OBJ file with vertices, normals, and UV coordinates.
        let obj_filename = format!("mesh_{:03}_{}.obj", i, case_name);
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
