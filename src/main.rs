mod data_read;
mod mesh_builder;
mod utils;

use data_read::read_contour_data;
use mesh_builder::{
    align_contours, compute_centroid, interpolate_contours, translate_contour, write_obj_mesh,
};
use std::error::Error;
use std::path::Path;
use utils::trim_to_same_length;

fn main() -> Result<(), Box<dyn Error>> {
    // Process both rest and stress cases
    process_case("rest", "input/rest_csv_files", "output/rest")?;
    process_case("stress", "input/stress_csv_files", "output/stress")?;
    Ok(())
}

fn process_case(
    case_name: &str,
    input_dir: &str,
    output_dir: &str,
) -> Result<(), Box<dyn Error>> {
    // Create output directory if it doesn't exist
    std::fs::create_dir_all(output_dir)?;

    // === DIASTOLE PROCESSING ===
    println!("--- Processing {} Diastole ---", case_name);
    let diastole_path = Path::new(input_dir).join("diastolic_contours.csv");
    let diastole_points = read_contour_data(diastole_path.to_str().unwrap())?;
    let mut diastole_contours = creating_contours(diastole_points);

    // === SYSTOLE PROCESSING ===
    println!("--- Processing {} Systole ---", case_name);
    let systole_path = Path::new(input_dir).join("systolic_contours.csv");
    let systole_points = read_contour_data(systole_path.to_str().unwrap())?;
    let mut systole_contours = creating_contours(systole_points);

    // Translate systolic contours to diastolic centroid
    let diastolic_ref_centroid = compute_centroid(&diastole_contours[0].1);
    for (_, ref mut contour) in systole_contours.iter_mut() {
        let systolic_centroid = compute_centroid(contour);
        let translation = (
            diastolic_ref_centroid.0 - systolic_centroid.0,
            diastolic_ref_centroid.1 - systolic_centroid.1,
        );
        translate_contour(contour, translation);
    }

    // Adjust z-axis coordinates
    let z_translation = diastole_contours[0].1[0].z - systole_contours[0].1[0].z;
    for (_, ref mut contour) in systole_contours.iter_mut() {
        for p in contour {
            p.z += z_translation;
        }
    }

    trim_to_same_length(&mut diastole_contours, &mut systole_contours);

    // Write base meshes
    let diastole_output = Path::new(output_dir).join(format!("diastole_{}.obj", case_name));
    let systole_output = Path::new(output_dir).join(format!("systole_{}.obj", case_name));
    write_obj_mesh(&diastole_contours, diastole_output.to_str().unwrap())?;
    write_obj_mesh(&systole_contours, systole_output.to_str().unwrap())?;

    // === INTERPOLATION ===
    println!("--- Interpolating {} Meshes ---", case_name);
    let steps = 30;
    let interpolated_meshes = interpolate_contours(&diastole_contours, &systole_contours, steps)?;

    // Write interpolated meshes
    for (i, intermediate) in interpolated_meshes.iter().enumerate() {
        let filename = Path::new(output_dir).join(format!("mesh_{:03}_{}.obj", i, case_name));
        write_obj_mesh(intermediate, filename.to_str().unwrap())?;
    }

    Ok(())
}

fn creating_contours(
    points: Vec<data_read::ContourPoint>,
) -> Vec<(u32, Vec<data_read::ContourPoint>)> {
    // Group points by their frame index
    use std::collections::HashMap;
    let mut groups = HashMap::new();
    for p in points {
        groups.entry(p.frame_index).or_insert(vec![]).push(p);
    }

    let mut contours: Vec<_> = groups.into_iter().collect();
    
    // Align and sort contours
    contours = align_contours(contours);
    contours.sort_by_key(|(frame, _)| std::cmp::Reverse(*frame));
    
    // Re-index contours
    contours.into_iter()
        .enumerate()
        .map(|(i, (_, points))| (i as u32, points))
        .collect()
}