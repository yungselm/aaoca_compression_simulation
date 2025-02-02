// mod data_read;
// use data_read::{read_contour_data, ContourPoint};

// use std::collections::HashMap;
// use std::error::Error;
// use std::f64::consts::PI;
// use std::fs::File;
// use std::io::{BufWriter, Write};

// /// Compute the centroid (average x,y) of a contour.
// fn compute_centroid(contour: &[ContourPoint]) -> (f64, f64) {
//     let (sum_x, sum_y) = contour.iter().fold((0.0, 0.0), |(sx, sy), p| (sx + p.x, sy + p.y));
//     let n = contour.len() as f64;
//     (sum_x / n, sum_y / n)
// }

// /// Rotate a point about a given center by an angle (in radians).
// fn rotate_point(p: &ContourPoint, angle: f64, center: (f64, f64)) -> ContourPoint {
//     let (cx, cy) = center;
//     // Translate point to origin.
//     let x = p.x - cx;
//     let y = p.y - cy;
//     // Rotate.
//     let cos_a = angle.cos();
//     let sin_a = angle.sin();
//     let new_x = x * cos_a - y * sin_a;
//     let new_y = x * sin_a + y * cos_a;
//     // Translate back.
//     ContourPoint {
//         frame_index: p.frame_index,
//         x: new_x + cx,
//         y: new_y + cy,
//         z: p.z,
//     }
// }

// /// Rotate all points in a contour about a given center.
// fn rotate_contour(contour: &mut [ContourPoint], angle: f64, center: (f64, f64)) {
//     for p in contour.iter_mut() {
//         let rotated = rotate_point(p, angle, center);
//         p.x = rotated.x;
//         p.y = rotated.y;
//         // We assume z remains unchanged.
//     }
// }

// /// Translate a contour by a (dx,dy) vector.
// fn translate_contour(contour: &mut [ContourPoint], translation: (f64, f64)) {
//     let (dx, dy) = translation;
//     for p in contour.iter_mut() {
//         p.x += dx;
//         p.y += dy;
//     }
// }

// /// Given a contour, find the pair of points with the maximum distance between them.
// /// Returns ((p1, p2), distance). (This is O(n^2); assume n is small.)
// fn find_farthest_points(contour: &[ContourPoint]) -> ((&ContourPoint, &ContourPoint), f64) {
//     let mut max_dist = 0.0;
//     let mut farthest_pair = (&contour[0], &contour[0]);
//     for i in 0..contour.len() {
//         for j in i + 1..contour.len() {
//             let dx = contour[i].x - contour[j].x;
//             let dy = contour[i].y - contour[j].y;
//             let dist = (dx * dx + dy * dy).sqrt();
//             if dist > max_dist {
//                 max_dist = dist;
//                 farthest_pair = (&contour[i], &contour[j]);
//             }
//         }
//     }
//     (farthest_pair, max_dist)
// }

// /// Given two contours (assumed to have the same number of points), 
// /// find a "best" rotation angle (in radians) for the target contour that minimizes
// /// the sum of squared distances to the reference contour. 
// /// Both contours are assumed already centered (i.e. centroids subtracted).
// fn find_best_rotation(reference: &[ContourPoint], target: &[ContourPoint]) -> f64 {
//     // Search over a small range of angles (say, ±0.2 radians) in fine increments.
//     let mut best_angle = 0.0;
//     let mut best_error = std::f64::MAX;
//     let steps = 400;
//     let range = 0.2;
//     let start = -range;
//     let end = range;
//     let increment = (end - start) / (steps as f64);
    
//     // We assume the contours are ordered (same number of points).
//     for i in 0..=steps {
//         let angle = start + (i as f64) * increment;
//         let mut error = 0.0;
//         // For each point, rotate target by the candidate angle and compute error.
//         for (p_ref, p_target) in reference.iter().zip(target.iter()) {
//             // Rotate target point around the origin.
//             let x = p_target.x * angle.cos() - p_target.y * angle.sin();
//             let y = p_target.x * angle.sin() + p_target.y * angle.cos();
//             let dx = x - p_ref.x;
//             let dy = y - p_ref.y;
//             error += dx * dx + dy * dy;
//         }
//         if error < best_error {
//             best_error = error;
//             best_angle = angle;
//         }
//     }
//     best_angle
// }

// /// Write an OBJ file connecting each adjacent pair of contours by creating two triangles per quad.
// fn write_obj_mesh(contours: Vec<(u32, Vec<ContourPoint>)>) -> Result<(), Box<dyn Error>> {
//     // Ensure the contours are sorted by frame_index.
//     let mut sorted_contours = contours;
//     sorted_contours.sort_by_key(|(frame_index, _)| *frame_index);

//     if sorted_contours.len() < 2 {
//         return Err("Need at least two contours to create a mesh.".into());
//     }

//     let points_per_contour = sorted_contours[0].1.len();
//     for (_, contour) in &sorted_contours {
//         if contour.len() != points_per_contour {
//             return Err("All contours must have the same number of points.".into());
//         }
//     }

//     let file = File::create("output.obj")?;
//     let mut writer = BufWriter::new(file);
//     let mut vertex_offsets = Vec::new();
//     let mut current_offset = 1;

//     // Write vertices.
//     for (_, contour) in &sorted_contours {
//         vertex_offsets.push(current_offset);
//         for point in contour {
//             writeln!(writer, "v {} {} {}", point.x, point.y, point.z)?;
//             current_offset += 1;
//         }
//     }

//     // Write faces by stitching adjacent contours.
//     for c in 0..(sorted_contours.len() - 1) {
//         let offset1 = vertex_offsets[c];
//         let offset2 = vertex_offsets[c + 1];
//         for j in 0..points_per_contour {
//             let j_next = (j + 1) % points_per_contour;
//             // First triangle.
//             let v1 = offset1 + j;
//             let v2 = offset1 + j_next;
//             let v3 = offset2 + j;
//             writeln!(writer, "f {} {} {}", v1, v2, v3)?;
//             // Second triangle.
//             let v1 = offset2 + j;
//             let v2 = offset1 + j_next;
//             let v3 = offset2 + j_next;
//             writeln!(writer, "f {} {} {}", v1, v2, v3)?;
//         }
//     }

//     println!("OBJ mesh successfully written to output.obj");
//     Ok(())
// }

// fn main() -> Result<(), Box<dyn Error>> {
//     // Read contours from file.
//     let contour_points = read_contour_data("input/rest_csv_files/diastolic_contours.csv")?;

//     // Group points by contour (using frame_index).
//     let mut groups: HashMap<u32, Vec<ContourPoint>> = HashMap::new();
//     for point in contour_points {
//         groups.entry(point.frame_index).or_default().push(point);
//     }

//     // Sort each contour's points (counterclockwise starting from highest y).
//     for contour in groups.values_mut() {
//         // (Reusing our previous sort function; here we assume it exists or inline it.)
//         // For simplicity, we use the same sorting as before.
//         // Compute the centroid.
//         let (sum_x, sum_y) = contour.iter().fold((0.0, 0.0), |(sx, sy), p| (sx + p.x, sy + p.y));
//         let count = contour.len() as f64;
//         let center = (sum_x / count, sum_y / count);
//         contour.sort_by(|a, b| {
//             let angle_a = (a.y - center.1).atan2(a.x - center.0);
//             let angle_b = (b.y - center.1).atan2(b.x - center.0);
//             angle_a.partial_cmp(&angle_b).unwrap()
//         });
//         let start_idx = contour
//             .iter()
//             .enumerate()
//             .max_by(|(_, a), (_, b)| a.y.partial_cmp(&b.y).unwrap())
//             .map(|(i, _)| i)
//             .unwrap_or(0);
//         contour.rotate_left(start_idx);
//     }

//     // Build a vector of contours (frame_index, contour).
//     let mut contours: Vec<(u32, Vec<ContourPoint>)> = groups.into_iter().collect();
//     // Sort contours by frame_index.
//     contours.sort_by_key(|(frame_index, _)| *frame_index);

//     // --- Alignment Phase ---

//     // Use the contour with the highest frame_index as the reference.
//     let reference_index = contours.iter().map(|(id, _)| *id).max().unwrap();
//     let reference_pos = contours
//         .iter()
//         .position(|(id, _)| *id == reference_index)
//         .expect("Reference contour not found");
//     let (ref_frame, ref_contour) = &contours[reference_pos];
//     println!("Using contour {} as reference.", ref_frame);

//     // Find farthest points in the reference contour.
//     let ((p1, p2), _dist) = find_farthest_points(ref_contour);
//     // Compute the angle (from the x-axis) of the line connecting these points.
//     let dx = p2.x - p1.x;
//     let dy = p2.y - p1.y;
//     let line_angle = dy.atan2(dx);
//     // We want this line to be vertical (aligned with the y-axis), i.e. at pi/2.
//     let rotation_to_y = (PI / 2.0) - line_angle;
//     println!(
//         "Reference line angle: {:.3} rad; rotating reference by {:.3} rad",
//         line_angle, rotation_to_y
//     );

//     // Rotate the reference contour about its centroid.
//     let ref_centroid = compute_centroid(ref_contour);
//     let mut aligned_reference = ref_contour.clone();
//     rotate_contour(&mut aligned_reference, rotation_to_y, ref_centroid);

//     // Replace the reference contour in our contours vector.
//     contours[reference_pos].1 = aligned_reference.clone();

//     // For the other contours, align them to the reference.
//     // The approach:
//     //   1. Subtract each contour's centroid so that rotation is around (0,0).
//     //   2. Do a search for the best rotation (minimizing sum of squared differences to reference).
//     //   3. Rotate and then add back the reference centroid.
//     for (id, contour) in contours.iter_mut() {
//         // Skip the reference contour.
//         if *id == reference_index {
//             continue;
//         }
//         let orig_centroid = compute_centroid(contour);
//         // Center the contour.
//         let centered: Vec<ContourPoint> = contour
//             .iter()
//             .map(|p| ContourPoint {
//                 frame_index: p.frame_index,
//                 x: p.x - orig_centroid.0,
//                 y: p.y - orig_centroid.1,
//                 z: p.z,
//             })
//             .collect();
//         // Also center the reference.
//         let centered_reference: Vec<ContourPoint> = aligned_reference
//             .iter()
//             .map(|p| ContourPoint {
//                 frame_index: p.frame_index,
//                 x: p.x - ref_centroid.0,
//                 y: p.y - ref_centroid.1,
//                 z: p.z,
//             })
//             .collect();

//         // Find best rotation for the target contour.
//         let best_rot = find_best_rotation(&centered_reference, &centered);
//         println!("Rotating contour {} by {:.3} rad for best correlation", id, best_rot);

//         // Rotate the original (non-centered) contour about its own centroid.
//         rotate_contour(contour, best_rot, orig_centroid);
//         // After rotation, translate the contour so its centroid matches the reference centroid.
//         let new_centroid = compute_centroid(contour);
//         let translation = (ref_centroid.0 - new_centroid.0, ref_centroid.1 - new_centroid.1);
//         translate_contour(contour, translation);
//     }

//     // --- End Alignment Phase ---

//     // (Optional) Print final aligned contours for debugging.
//     for (id, contour) in &contours {
//         println!("Aligned Contour {}:", id);
//         for (i, p) in contour.iter().enumerate() {
//             println!("  Point {}: ({:.3}, {:.3}, {:.3})", i, p.x, p.y, p.z);
//         }
//         println!("");
//     }

//     // Create the mesh from the aligned contours.
//     write_obj_mesh(contours)?;

//     Ok(())
// }

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
    // (You may wish to sort each contour’s points here if needed.)
    // Align diastole contours.
    diastole_contours = align_contours(diastole_contours);
    // Write the diastole mesh.
    write_obj_mesh(&diastole_contours, "diastole.obj")?;

    // === SYSTOLE PROCESSING ===
    println!("--- Processing Systole ---");
    let systole_points = read_contour_data("input/rest_csv_files/systolic_contours.csv")?;
    let mut systole_groups: HashMap<u32, Vec<_>> = HashMap::new();
    for p in systole_points {
        systole_groups.entry(p.frame_index).or_default().push(p);
    }
    let mut systole_contours: Vec<(u32, Vec<_>)> = systole_groups.into_iter().collect();
    // Align systole contours.
    systole_contours = align_contours(systole_contours);
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
