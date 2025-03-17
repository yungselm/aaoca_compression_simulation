// src/comparison.rs
use crate::{
    io::{read_contour_data, create_catheter_points, write_obj_mesh},
    contour::create_contours,
    processing::{self, compute_centroid, translate_contour, interpolate_contours},
    utils::{trim_to_same_length, smooth_contours},
    texture::{compute_uv_coordinates, compute_displacements, create_displacement_texture, create_black_texture},
};
use std::path::Path;
use std::error::Error;
use std::fs::File;
use std::io::Write;

/// Processes a comparison between rest and stress for a specific cardiac phase
pub fn process_phase_comparison(
    phase_name: &str,
    rest_input_dir: &str,
    stress_input_dir: &str,
    output_dir: &str,
) -> Result<(), Box<dyn Error>> {
    std::fs::create_dir_all(output_dir)?;

    // === Process Rest Phase ===
    println!("--- Processing Rest {} ---", phase_name);
    let rest_path = Path::new(rest_input_dir).join(format!("{}_contours.csv", phase_name));
    let rest_points = read_contour_data(rest_path.to_str().unwrap())?;
    let rest_catheter = create_catheter_points(&rest_points);
    let mut rest_contours = create_contours(rest_points);
    let mut rest_catheter_contours = create_contours(rest_catheter);

    // === Process Stress Phase ===
    println!("--- Processing Stress {} ---", phase_name);
    let stress_path = Path::new(stress_input_dir).join(format!("{}_contours.csv", phase_name));
    let stress_points = read_contour_data(stress_path.to_str().unwrap())?;
    let stress_catheter = create_catheter_points(&stress_points);
    let mut stress_contours = create_contours(stress_points);
    let mut stress_catheter_contours = create_contours(stress_catheter);

    // Align stress contours to rest reference
    let rest_ref_centroid = compute_centroid(&rest_contours[0].1);
    for (_, contour) in stress_contours.iter_mut().chain(stress_catheter_contours.iter_mut()) {
        let stress_centroid = compute_centroid(contour);
        let translation = (
            rest_ref_centroid.0 - stress_centroid.0,
            rest_ref_centroid.1 - stress_centroid.1,
        );
        translate_contour(contour, translation);
    }

    // Z-axis alignment
    let z_translation = rest_contours[0].1[0].z - stress_contours[0].1[0].z;
    for (_, contour) in stress_contours.iter_mut().chain(stress_catheter_contours.iter_mut()) {
        for point in contour {
            point.z += z_translation;
        }
    }

    // Common processing
    trim_to_same_length(&mut rest_contours, &mut stress_contours);
    trim_to_same_length(&mut rest_catheter_contours, &mut stress_catheter_contours);
    
    rest_contours = smooth_contours(&rest_contours);
    stress_contours = smooth_contours(&stress_contours);

    // Interpolation between rest and stress
    let steps = 30;
    let interpolated_meshes = interpolate_contours(&rest_contours, &stress_contours, steps)?;
    let interpolated_catheter_meshes = interpolate_contours(&rest_catheter_contours, &stress_catheter_contours, steps)?;

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
        .filter_map(|(i, mesh)| {
            if i == 1 || i == 2 { None } else { Some(mesh) }
        })
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
        let tex_filename = format!("catheter_{:03}_{}.png", i, phase_name);
        let texture_path = Path::new(output_dir).join(&tex_filename);
        create_black_texture(texture_width, texture_height, texture_path.to_str().unwrap())?;

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
        write_obj_mesh(
            mesh,
            &uv_coords,
            obj_path.to_str().unwrap(),
            &mtl_filename,
        )?;
    }
Ok(())
}
