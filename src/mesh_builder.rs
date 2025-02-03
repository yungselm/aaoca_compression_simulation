use crate::data_read::ContourPoint;
use std::error::Error;
use std::f64::consts::PI;
use std::fs::File;
use std::io::{BufWriter, Write};

/// Compute the centroid (average x,y) of a contour.
pub fn compute_centroid(contour: &[ContourPoint]) -> (f64, f64) {
    let (sum_x, sum_y) = contour
        .iter()
        .fold((0.0, 0.0), |(sx, sy), p| (sx + p.x, sy + p.y));
    let n = contour.len() as f64;
    (sum_x / n, sum_y / n)
}

/// Sort a contour’s points in counterclockwise order around its centroid,
/// and then rotate the vector so that the highest y-value is first.
pub fn sort_contour_points(contour: &mut Vec<ContourPoint>) {
    let center = compute_centroid(contour);
    contour.sort_by(|a, b| {
        let angle_a = (a.y - center.1).atan2(a.x - center.0);
        let angle_b = (b.y - center.1).atan2(b.x - center.0);
        angle_a.partial_cmp(&angle_b).unwrap()
    });
    let start_idx = contour
        .iter()
        .enumerate()
        .max_by(|(_, a), (_, b)| a.y.partial_cmp(&b.y).unwrap())
        .map(|(i, _)| i)
        .unwrap_or(0);
    contour.rotate_left(start_idx);
}

/// Rotate a point about a given center by an angle (in radians).
pub fn rotate_point(p: &ContourPoint, angle: f64, center: (f64, f64)) -> ContourPoint {
    let (cx, cy) = center;
    let x = p.x - cx;
    let y = p.y - cy;
    let cos_a = angle.cos();
    let sin_a = angle.sin();
    ContourPoint {
        frame_index: p.frame_index,
        x: x * cos_a - y * sin_a + cx,
        y: x * sin_a + y * cos_a + cy,
        z: p.z,
    }
}

/// Rotate all points in a contour about a given center.
pub fn rotate_contour(contour: &mut [ContourPoint], angle: f64, center: (f64, f64)) {
    for p in contour.iter_mut() {
        let rotated = rotate_point(p, angle, center);
        p.x = rotated.x;
        p.y = rotated.y;
    }
}

/// Translate a contour by a (dx,dy) vector.
pub fn translate_contour(contour: &mut [ContourPoint], translation: (f64, f64)) {
    let (dx, dy) = translation;
    for p in contour.iter_mut() {
        p.x += dx;
        p.y += dy;
    }
}

/// Given a contour, find the pair of points with the maximum distance between them.
/// Returns ((p1, p2), distance).
pub fn find_farthest_points(contour: &[ContourPoint]) -> ((&ContourPoint, &ContourPoint), f64) {
    let mut max_dist = 0.0;
    let mut farthest_pair = (&contour[0], &contour[0]);
    for i in 0..contour.len() {
        for j in i + 1..contour.len() { // currently O(n²), could be optimized
            let dx = contour[i].x - contour[j].x;
            let dy = contour[i].y - contour[j].y;
            let dist = (dx * dx + dy * dy).sqrt();
            if dist > max_dist {
                max_dist = dist;
                farthest_pair = (&contour[i], &contour[j]);
            }
        }
    }
    (farthest_pair, max_dist)
}

/// Given two contours (assumed to have the same number of points),
/// find a "best" rotation angle (in radians) for the target contour that minimizes
/// the sum of squared distances to the reference contour. 
/// Both contours are assumed to be centered (i.e. their centroids subtracted).
pub fn find_best_rotation(reference: &[ContourPoint], target: &[ContourPoint]) -> f64 {
    let mut best_angle = 0.0;
    let mut best_error = std::f64::MAX;
    let steps = 400;
    let range = 1.05; // around +/- 60 degrees
    let start = -range;
    let end = range;
    let increment = (end - start) / (steps as f64);

    for i in 0..=steps {
        let angle = start + (i as f64) * increment;
        let mut error = 0.0;
        for (p_ref, p_target) in reference.iter().zip(target.iter()) {
            let x = p_target.x * angle.cos() - p_target.y * angle.sin();
            let y = p_target.x * angle.sin() + p_target.y * angle.cos();
            let dx = x - p_ref.x;
            let dy = y - p_ref.y;
            error += dx * dx + dy * dy;
        }
        if error < best_error {
            best_error = error;
            best_angle = angle;
        }
    }
    best_angle
}

/// Writes an OBJ file connecting each neighboring pair of contours (each as a slice of ContourPoint)
/// by creating two triangles per quad.
pub fn write_obj_mesh(
    contours: &[(u32, Vec<ContourPoint>)],
    filename: &str,
) -> Result<(), Box<dyn Error>> {
    // Make a copy and sort by frame_index.
    let mut sorted_contours = contours.to_owned();
    sorted_contours.sort_by_key(|(frame_index, _)| *frame_index);

    if sorted_contours.len() < 2 {
        return Err("Need at least two contours to create a mesh.".into());
    }

    let points_per_contour = sorted_contours[0].1.len();
    for (_, contour) in &sorted_contours {
        if contour.len() != points_per_contour {
            return Err("All contours must have the same number of points.".into());
        }
    }

    let file = File::create(filename)?;
    let mut writer = BufWriter::new(file);
    let mut vertex_offsets = Vec::new();
    let mut current_offset = 1;

    // Write vertices.
    for (_, contour) in &sorted_contours {
        vertex_offsets.push(current_offset);
        for point in contour {
            writeln!(writer, "v {} {} {}", point.x, point.y, point.z)?;
            current_offset += 1;
        }
    }

    // Write faces.
    for c in 0..(sorted_contours.len() - 1) {
        let offset1 = vertex_offsets[c];
        let offset2 = vertex_offsets[c + 1];
        for j in 0..points_per_contour {
            let j_next = (j + 1) % points_per_contour;
            // Triangle 1.
            let v1 = offset1 + j;
            let v2 = offset1 + j_next;
            let v3 = offset2 + j;
            writeln!(writer, "f {} {} {}", v1, v2, v3)?;
            // Triangle 2.
            let v1 = offset2 + j;
            let v2 = offset1 + j_next;
            let v3 = offset2 + j_next;
            writeln!(writer, "f {} {} {}", v1, v2, v3)?;
        }
    }

    println!("OBJ mesh successfully written to {}", filename);
    Ok(())
}

/// Align a vector of contours by choosing the one with the highest frame index as reference (always ostium).
/// The alignment rotates the reference contour so that its farthest-points line is vertical,
/// then rotates each other contour to best match the reference, and finally translates all
/// contours so that their centroids coincide with the reference centroid.
/// 
/// The input is a vector of tuples: (contour_id, contour_points).
/// Returns the aligned contours.
pub fn align_contours(mut contours: Vec<(u32, Vec<ContourPoint>)>) -> Vec<(u32, Vec<ContourPoint>)> {
    // First, ensure every contour is sorted in the same manner.
    for (_, contour) in contours.iter_mut() {
        sort_contour_points(contour);
    }

    // Sort contours by frame_index.
    contours.sort_by_key(|(frame_index, _)| *frame_index);

    // Identify the reference contour: the one with the highest frame_index.
    let reference_index = contours.iter().map(|(id, _)| *id).max().unwrap();
    let reference_pos = contours
        .iter()
        .position(|(id, _)| *id == reference_index)
        .expect("Reference contour not found");
    let (ref_frame, ref_contour) = &contours[reference_pos];
    println!("Using contour {} as reference.", ref_frame);

    // Find the farthest points in the reference contour.
    let ((p1, p2), _dist) = find_farthest_points(ref_contour);
    let dx = p2.x - p1.x;
    let dy = p2.y - p1.y;
    let line_angle = dy.atan2(dx);
    // We want the line to be vertical (i.e. at PI/2).
    let rotation_to_y = (PI / 2.0) - line_angle;
    println!(
        "Reference line angle: {:.3} rad; rotating reference by {:.3} rad",
        line_angle, rotation_to_y
    );

    // Rotate the reference contour about its centroid.
    let ref_centroid = compute_centroid(ref_contour);
    let mut aligned_reference = (*ref_contour).clone();
    rotate_contour(&mut aligned_reference, rotation_to_y, ref_centroid);
    // Re-sort the reference after rotation to ensure the highest-y point is first.
    sort_contour_points(&mut aligned_reference);
    contours[reference_pos].1 = aligned_reference.clone();

    // Extract a clone of the reference contour so we can use it immutably.
    let reference_clone = contours[reference_pos].1.clone();

    // Align each non-reference contour.
    for (id, contour) in contours.iter_mut() {
        if *id == reference_index {
            continue;
        }
        let orig_centroid = compute_centroid(contour);
        // Center the contour.
        let centered: Vec<ContourPoint> = contour
            .iter()
            .map(|p| ContourPoint {
                frame_index: p.frame_index,
                x: p.x - orig_centroid.0,
                y: p.y - orig_centroid.1,
                z: p.z,
            })
            .collect();
        // Center the reference using our pre-cloned reference.
        let centered_reference: Vec<ContourPoint> = reference_clone
            .iter()
            .map(|p| ContourPoint {
                frame_index: p.frame_index,
                x: p.x - ref_centroid.0,
                y: p.y - ref_centroid.1,
                z: p.z,
            })
            .collect();

        let best_rot = find_best_rotation(&centered_reference, &centered);
        println!("Rotating contour {} by {:.3} rad for best correlation", id, best_rot);
        rotate_contour(contour, best_rot, orig_centroid);
        // Re-sort each contour after rotation.
        sort_contour_points(contour);
        let new_centroid = compute_centroid(contour);
        let translation = (ref_centroid.0 - new_centroid.0, ref_centroid.1 - new_centroid.1);
        translate_contour(contour, translation);
    }

    contours
}

// Interpolate between two aligned sets of contours.
// For each corresponding point in start and end, a linear interpolation is performed.
// If the number of contours differs, the longer set is trimmed to the shorter one.
pub fn interpolate_contours(
    contours_start: &[(u32, Vec<ContourPoint>)],
    contours_end: &[(u32, Vec<ContourPoint>)],
    steps: usize,
) -> Result<Vec<Vec<(u32, Vec<ContourPoint>)>>, Box<dyn Error>> {
    use std::cmp::min;

    // Trim the sets to the same number.
    let n = min(contours_start.len(), contours_end.len());
    let start = &contours_start[0..n];
    let end = &contours_end[0..n];

    let mut interpolated = Vec::with_capacity(steps);
    for step in 0..steps {
        let t = step as f64 / (steps - 1) as f64;
        let mut intermediate = Vec::with_capacity(n);
        for ((id_start, contour_start), (id_end, contour_end)) in start.iter().zip(end.iter()) {
            if id_start != id_end {
                return Err("Contour IDs do not match between start and end.".into());
            }
            if contour_start.len() != contour_end.len() {
                return Err("Contour point counts do not match between start and end.".into());
            }
            let interp_contour: Vec<ContourPoint> = contour_start
                .iter()
                .zip(contour_end.iter())
                .map(|(p_start, p_end)| ContourPoint {
                    frame_index: *id_start,
                    x: p_start.x * (1.0 - t) + p_end.x * t,
                    y: p_start.y * (1.0 - t) + p_end.y * t,
                    z: p_start.z * (1.0 - t) + p_end.z * t,
                })
                .collect();
            intermediate.push((*id_start, interp_contour));
        }
        interpolated.push(intermediate);
    }
    Ok(interpolated)
}
