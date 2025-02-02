mod data_read;
mod mesh_builder;

use data_read::read_contour_data;
use mesh_builder::{align_contours, interpolate_contours, write_obj_mesh};
use std::collections::HashMap;
use std::error::Error;

fn main() -> Result<(), Box<dyn Error>> {
    // === DIASTOLE PROCESSING ===
    println!("--- Processing Diastole ---");
    let diastole_points = read_contour_data("input/rest_csv_files/diastolic_contours.csv")?;

    // Group points by contour id (using frame_index).
    let mut diastole_groups: HashMap<u32, Vec<_>> = HashMap::new();
    for p in diastole_points {
        diastole_groups.entry(p.frame_index).or_default().push(p);
    }
    let mut diastole_contours: Vec<(u32, Vec<_>)> = diastole_groups.into_iter().collect();
    // Align diastole contours.
    diastole_contours = align_contours(diastole_contours);
    // After aligning diastole contours:
    diastole_contours = align_contours(diastole_contours);

    // Reindex diastole contours:
    diastole_contours = diastole_contours
        .into_iter()
        .enumerate()
        .map(|(i, (_, contour))| (i as u32, contour))
        .collect();

    // Write the diastole mesh.
    write_obj_mesh(&diastole_contours, "diastole.obj")?;

    // === SYSTOLE PROCESSING ===
    let systole_points = read_contour_data("input/rest_csv_files/systolic_contours.csv")?;
    let mut systole_groups: HashMap<u32, Vec<_>> = HashMap::new();
    for p in systole_points {
        systole_groups.entry(p.frame_index).or_default().push(p);
    }
    let mut systole_contours: Vec<(u32, Vec<_>)> = systole_groups.into_iter().collect();
    systole_contours = align_contours(systole_contours);

    // Reindex systole contours:
    systole_contours = systole_contours
        .into_iter()
        .enumerate()
        .map(|(i, (_, contour))| (i as u32, contour))
        .collect();

    // Write the systole mesh.
    write_obj_mesh(&systole_contours, "systole.obj")?;

    // === INTERPOLATION BETWEEN DIASTOLE AND SYSTOLE ===
    println!("--- Interpolating Meshes ---");
    // Generate, for example, 30 intermediate steps.
    let steps = 30;
    let interpolated_meshes = interpolate_contours(&diastole_contours, &systole_contours, steps)?;

    // Write each interpolated mesh to its own OBJ file.
    for (i, intermediate) in interpolated_meshes.iter().enumerate() {
        let filename = format!("mesh_{:03}.obj", i);
        write_obj_mesh(intermediate, &filename)?;
    }

    Ok(())
}
