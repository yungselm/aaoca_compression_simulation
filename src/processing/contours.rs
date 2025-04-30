use rayon::prelude::*;
use std::f64::consts::PI;

use crate::io::input::ContourPoint;
use crate::io::Geometry;

const MIDDLE_IMAGE: f64 = 4.5; // this is for an IVUS image with 512x512 pixels and 0.01759 pixelsspacing

pub fn align_frames_in_geometry(mut geometry: Geometry, steps: usize, range: f64) -> Geometry {
    // Sort contours by frame index.
    geometry.contours.sort_by_key(|contour| contour.id);

    // Use sort_contour_points on all catheter points to ensure they are counterclockwise
    for catheter in &mut geometry.catheter {
        catheter.sort_contour_points();
    }

    // Use the contour with the highest frame index as reference.
    let reference_index = geometry
        .contours
        .iter()
        .map(|contour| contour.id)
        .max()
        .unwrap();
    let reference_pos = geometry
        .contours
        .iter()
        .position(|contour| contour.id == reference_index)
        .expect("Reference contour not found");
    let ref_contour = &mut geometry.contours.remove(reference_pos);
    println!("Using contour {} as reference.", reference_index);

    let ((p1, p2), _dist) = ref_contour.find_farthest_points();
    let dx = p2.x - p1.x;
    let dy = p2.y - p1.y;
    let line_angle = dy.atan2(dx);
    let mut rotation_to_y = (PI / 2.0) - line_angle;
    println!(
        "Reference line angle: {:.3} rad; rotating reference by {:.3} rad",
        &line_angle, &rotation_to_y
    );

    // Clone the reference contour to apply initial rotation without affecting the original yet
    let mut rotated_ref = ref_contour.clone(); //
    rotated_ref.rotate_contour(rotation_to_y);
    rotated_ref.sort_contour_points();

    let n = rotated_ref.points.len() / 2;

    // --- Begin aortic determination ---
    // We assume there are exactly 500 points; split into two halves: indices 0..249 and 250..499.
    // Compute the cumulative distance for each half from the rotated reference.
    let first_half_distance: f64 = rotated_ref
        .points
        .iter()
        .take(n)
        .map(|pt| pt.distance_to(&geometry.reference_point))
        .sum();
    let second_half_distance: f64 = rotated_ref
        .points
        .iter()
        .skip(n)
        .map(|pt| pt.distance_to(&geometry.reference_point))
        .sum();

    // Mark the half that is closer as 'aortic' (i.e. set to true).
    if first_half_distance < second_half_distance {
        for pt in rotated_ref.points.iter_mut().take(n) {
            pt.aortic = true;
        }
        println!("First half (points 0 to 249) marked as aortic.");
        if first_half_distance < MIDDLE_IMAGE {
            rotation_to_y += PI;
            rotated_ref.rotate_contour(PI);
        }
    } else {
        for pt in rotated_ref.points.iter_mut().skip(n) {
            pt.aortic = true;
        }
        println!("Second half (points 250 to 499) marked as aortic.");
        if second_half_distance < MIDDLE_IMAGE {
            rotation_to_y += PI;
            rotated_ref.rotate_contour(PI);
        }
    }

    // Update the reference contour in geometry with rotated version
    geometry.contours.insert(reference_pos, rotated_ref.clone());

    // Align each contour to the reference
    // Store processed contours' points and centroids as reference for each following contour
    let mut processed_refs = std::collections::HashMap::new();
    let mut id_translation = Vec::new();

    // Process contours in reverse order (highest ID first)
    for contour in geometry.contours.iter_mut().rev() {
        // Skip the reference contour
        if contour.id == reference_index {
            continue;
        }
        // Determine reference points and centroid
        let (ref_points, ref_centroid) = match processed_refs.get(&(contour.id + 1)) {
            Some((points, centroid)) => (points, centroid),
            None => (&rotated_ref.points, &rotated_ref.centroid),
        };

        // Initial rotation alignment (same as original)
        contour.rotate_contour(rotation_to_y);

        // Translate to reference centroid
        let tx = ref_centroid.0 - contour.centroid.0;
        let ty = ref_centroid.1 - contour.centroid.1;

        id_translation.push((contour.id, (tx, ty, 0.0)));

        contour.translate_contour((tx, ty, 0.0));

        // Optimize rotation using dynamic reference
        let best_rot = find_best_rotation(
            ref_points,
            &contour.points,
            steps,
            range,
            &contour.centroid,
        );
        println!("Matching Contour {:?} -> Contour {:?}, Best rotation: {:?}", &contour.id, contour.id + 1, &best_rot);
        contour.rotate_contour(best_rot);
        contour.sort_contour_points();

        // Store this processed contour's state for future reference
        processed_refs.insert(contour.id, (contour.points.clone(), contour.centroid));

        // Mark second half as aortic (original logic)
        let half_len = contour.points.len() / 2;
        for pt in contour.points.iter_mut().skip(half_len) {
            pt.aortic = true;
        }
    }

    for catheter in geometry.catheter.iter_mut() {
        // Translate catheter points to match the reference contour
        for (id, translation) in &id_translation {
            if catheter.id == *id {
                catheter.translate_contour(*translation);
            }
        }
    }

    geometry
}

/// Finds the best rotation angle (in radians) that minimizes the Hausdorff distance
/// leveraging parallel computation for performance.
pub fn find_best_rotation(
    reference: &[ContourPoint],
    target: &[ContourPoint],
    steps: usize,
    range: f64,
    centroid: &(f64, f64, f64),
) -> f64 {
    let increment = (2.0 * range) / (steps as f64);

    (0..=steps)
        .into_par_iter() // Parallel iteration (requires Rayon)
        .map(|i| {
            let angle = -range + (i as f64) * increment;
            // Rotate target contour inline without extra allocation
            let rotated: Vec<ContourPoint> = target
                .iter()
                .map(|p| p.rotate_point(angle, (centroid.0, centroid.1)))
                .collect();

            let hausdorff_dist = hausdorff_distance(reference, &rotated);
            (angle, hausdorff_dist)
        })
        .reduce(
            || (std::f64::NEG_INFINITY, std::f64::MAX),
            |(best_a, best_d), (angle, dist)| {
                if dist < best_d {
                    (angle, dist)
                } else {
                    (best_a, best_d)
                }
            },
        )
        .0
}

/// Computes the Hausdorff distance between two point sets.
pub fn hausdorff_distance(set1: &[ContourPoint], set2: &[ContourPoint]) -> f64 {
    let forward = directed_hausdorff(set1, set2);
    let backward = directed_hausdorff(set2, set1);
    forward.max(backward) // Hausdorff distance is max of both directed distances
}

/// Computes directed Hausdorff distance from A to B
fn directed_hausdorff(contour_a: &[ContourPoint], contour_b: &[ContourPoint]) -> f64 {
    contour_a
        .par_iter() // Use parallel iteration
        .map(|pa| {
            contour_b
                .iter()
                .map(|pb| {
                    let dx = pa.x - pb.x;
                    let dy = pa.y - pb.y;
                    (dx * dx + dy * dy).sqrt()
                })
                .fold(std::f64::MAX, f64::min) // Directly find min without storing a Vec
        })
        .reduce(|| 0.0, f64::max) // Directly find max without extra allocation
}
