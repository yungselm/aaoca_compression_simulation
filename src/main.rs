mod data_read;
mod mesh_builder;
mod utils;

use data_read::read_contour_data;
use mesh_builder::{
    align_contours, compute_centroid, interpolate_contours, translate_contour, write_obj_mesh,
};
use std::collections::HashMap;
use std::error::Error;
use utils::trim_to_same_length;

fn main() -> Result<(), Box<dyn Error>> {
    // === DIASTOLE PROCESSING ===
    println!("--- Processing Diastole ---");
    let diastole_points = read_contour_data("input/rest_csv_files/diastolic_contours.csv")?;

    let mut diastole_contours = creating_contours(diastole_points);

    // === SYSTOLE PROCESSING ===
    println!("--- Processing Systole ---");
    let systole_points = read_contour_data("input/rest_csv_files/systolic_contours.csv")?;

    let mut systole_contours = creating_contours(systole_points);

    // *** Translate systolic contours to diastolic (ostium) centroid:
    let diastolic_ref_centroid = compute_centroid(&diastole_contours[0].1);
    println!(
        "Diastolic (ostium) centroid: ({:.3}, {:.3})",
        diastolic_ref_centroid.0, diastolic_ref_centroid.1
    );
    for (_, ref mut contour) in systole_contours.iter_mut() {
        let systolic_centroid = compute_centroid(contour);
        let translation = (
            diastolic_ref_centroid.0 - systolic_centroid.0,
            diastolic_ref_centroid.1 - systolic_centroid.1,
        );
        translate_contour(contour, translation);
    }

    // move whole stack on z-axis so both stacks have same z-coordinates
    let z_translation = diastole_contours[0].1[0].z - systole_contours[0].1[0].z;
    for (_, ref mut contour) in systole_contours.iter_mut() {
        for p in contour {
            p.z += z_translation;
        }
    }

    let _min_len = trim_to_same_length(&mut diastole_contours, &mut systole_contours);
    // Write the meshes
    write_obj_mesh(&diastole_contours, "output/diastole_rest.obj")?;
    write_obj_mesh(&systole_contours, "output/systole_rest.obj")?;

    // === INTERPOLATION BETWEEN DIASTOLE AND SYSTOLE ===
    println!("--- Interpolating Meshes ---");
    let steps = 30;
    let interpolated_meshes = interpolate_contours(&diastole_contours, &systole_contours, steps)?;

    // Write each interpolated mesh.
    for (i, intermediate) in interpolated_meshes.iter().enumerate() {
        let filename = format!("output/mesh_{:03}_rest.obj", i);
        write_obj_mesh(intermediate, &filename)?;
    }

    Ok(())
}

fn creating_contours(
    points: Vec<data_read::ContourPoint>,
) -> Vec<(u32, Vec<data_read::ContourPoint>)> {
    // Group points by their (arbitrary) frame id.
    let mut groups: HashMap<u32, Vec<_>> = HashMap::new();
    for p in points {
        groups.entry(p.frame_index).or_default().push(p);
    }
    let mut contours: Vec<(u32, Vec<_>)> = groups.into_iter().collect();

    // Align contours (each contour’s internal order is sorted, but we leave the outer vector order unchanged)
    contours = align_contours(contours);

    // *** Sort the aligned diastolic contours once here:
    // Sort in descending order so that the contour with the highest original frame index (the ostium) comes first.
    contours.sort_by_key(|(frame, _)| std::cmp::Reverse(*frame));

    // *** Re-index so that index 0 is the ostium.
    contours = contours
        .into_iter()
        .enumerate()
        .map(|(i, (_, contour))| (i as u32, contour))
        .collect();
    contours
}
