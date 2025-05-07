use rayon::prelude::*;
use std::f64::consts::PI;
use std::collections::HashSet;

use crate::io::input::ContourPoint;
use crate::io::Geometry;

pub fn align_frames_in_geometry(mut geometry: Geometry, steps: usize, range: f64) -> Geometry {
    // Sort contours by frame index.
    geometry.contours.sort_by_key(|contour| contour.id);

    // Use sort_contour_points on all contour points to ensure they are counterclockwise
    for contour in &mut geometry.contours {
        contour.sort_contour_points();
    }

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

    let ref_contour_copy = ref_contour.clone();

    let ((p1, p2), _dist) = ref_contour_copy.find_farthest_points();
    
    // Calculate indices in sorted contour
    let p1_pos = ref_contour_copy.points.iter().position(|pt| pt == p1).unwrap();
    let p2_pos = ref_contour_copy.points.iter().position(|pt| pt == p2).unwrap();
    
    // Create index sets for aortic regions
    let (first_half_indices, second_half_indices) = if p1_pos < p2_pos {
        (
            (p1_pos..=p2_pos).collect::<HashSet<_>>(),
            (0..p1_pos).chain(p2_pos+1..ref_contour_copy.points.len()).collect::<HashSet<_>>()
        )
    } else {
        (
            (p1_pos..ref_contour_copy.points.len()).chain(0..=p2_pos).collect(),
            (p2_pos+1..p1_pos).collect()
        )
    };
    
    // Calculate distances using actual positions in sorted contour
    let dist_first: f64 = ref_contour_copy.points.iter().enumerate()
    .filter(|(i, _)| first_half_indices.contains(i))
    .map(|(_, p)| p.distance_to(&geometry.reference_point))
    .sum();

    let dist_second: f64 = ref_contour_copy.points.iter().enumerate()
    .filter(|(i, _)| second_half_indices.contains(i))
    .map(|(_, p)| p.distance_to(&geometry.reference_point))
    .sum();

    // Assign aortic flags using current indices
    for pt in ref_contour.points.iter_mut() {
        pt.aortic = if dist_first < dist_second {
                first_half_indices.contains(&(pt.point_index as usize))
            } else {
                second_half_indices.contains(&(pt.point_index as usize))
            };
        }
    
    // Find a rotation that makes contour vertical and ensures aortic points are on the right side
    let dx = p2.x - p1.x;
    let dy = p2.y - p1.y;
    let line_angle = dy.atan2(dx);
    let mut rotation_to_y = (PI / 2.0) - line_angle;
    
    let ref_contour_copy = ref_contour.clone();

    let ((p3, p4), _dist) = ref_contour_copy.find_closest_opposite();
    p3.rotate_point(rotation_to_y, (ref_contour_copy.centroid.0, ref_contour_copy.centroid.1));
    p3.rotate_point(rotation_to_y, (ref_contour_copy.centroid.0, ref_contour_copy.centroid.1));
    
    // Adjust rotation if necessary to ensure aortic points are on the right side
    if p4.aortic && p3.x < p4.x {
        rotation_to_y += PI;
    }

    println!("----------------------Aligning frames----------------------");
    println!(
        "Reference line angle: {:.3} rad; rotating reference by {:.3} rad",
        &line_angle, &rotation_to_y
    );

    // Clone the reference contour to apply initial rotation without affecting the original yet
    let mut rotated_ref = ref_contour.clone(); //
    rotated_ref.rotate_contour(rotation_to_y);
    rotated_ref.sort_contour_points();

    // Update the reference contour in geometry with rotated version
    geometry.contours.insert(reference_pos, rotated_ref.clone());

    let ref_catheter = &mut geometry.catheter.remove(reference_pos);
    ref_catheter.rotate_contour_around_point(rotation_to_y, (ref_contour.centroid.0, ref_contour.centroid.1));
    ref_catheter.sort_contour_points();
    geometry.catheter.insert(reference_pos, ref_catheter.clone());

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

        
        contour.translate_contour((tx, ty, 0.0));
        
        // Optimize rotation using dynamic reference
        let best_rot = find_best_rotation(
            ref_points,
            &contour.points,
            steps,
            range,
            &contour.centroid,
        );

        let rotation = rotation_to_y + best_rot;

        id_translation.push((contour.id, (tx, ty, 0.0), rotation, (contour.centroid.0, contour.centroid.1)));

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
        for (id, translation, best_rot, center) in &id_translation {
            if catheter.id == *id {
                catheter.translate_contour(*translation);
                catheter.rotate_contour_around_point(*best_rot, *center);
                catheter.sort_contour_points();
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
