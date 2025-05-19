use rayon::prelude::*;
use std::f64::consts::PI;
use std::collections::HashSet;

use crate::io::input::{ContourPoint, Contour};
use crate::io::Geometry;

pub fn align_frames_in_geometry(geometry: Geometry, steps: usize, range: f64) -> Geometry {
    let (mut geometry, reference_index, reference_pos, ref_contour) = 
        prep_data_geometry(geometry);

    let ref_contour_copy = ref_contour.clone();
    
    let (p1, p2, updated_ref) = 
        assign_aortic(ref_contour_copy, &geometry);
    let ref_contour = updated_ref.clone();

    let (_line_angle, rotation_to_y, rotated_ref) = 
        rotate_reference_contour(p1, p2, ref_contour.clone());

    // Update the reference contour in geometry with rotated version
    geometry.contours.insert(reference_pos, rotated_ref.clone());

    // prepare reference catheter
    let ref_catheter = &mut geometry.catheter.remove(reference_pos);
    ref_catheter.rotate_contour_around_point(
        rotation_to_y, 
        (ref_contour.centroid.0, ref_contour.centroid.1)
    );
    ref_catheter.sort_contour_points();
    geometry.catheter.insert(reference_pos, ref_catheter.clone());

    let (mut geometry, id_translation) = align_remaining_contours(
        geometry,
        reference_index,
        rotated_ref,
        rotation_to_y,
        steps,
        range,
    );

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

fn prep_data_geometry(mut geometry: Geometry) -> (Geometry, u32, usize, Contour) {
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

    (geometry, reference_index, reference_pos, ref_contour.clone())    
}

fn assign_aortic(
    contour: Contour,
    geometry: &Geometry, 
) -> (ContourPoint, ContourPoint, Contour) {
    let ((p1, p2), _dist) = contour.find_farthest_points();

    let p1_pos = contour.points.iter().position(|pt| pt == p1).unwrap();
    let p2_pos = contour.points.iter().position(|pt| pt == p2).unwrap();

    let (first_half_indices, second_half_indices) = if p1_pos < p2_pos {
        (
            (p1_pos..=p2_pos).collect::<HashSet<_>>(),
            (0..p1_pos).chain(p2_pos+1..contour.points.len()).collect::<HashSet<_>>()
        )
    } else {
        (
            (p1_pos..contour.points.len()).chain(0..=p2_pos).collect::<HashSet<_>>(),
            (p2_pos+1..p1_pos).collect::<HashSet<_>>()
        )
    };

    // Compute distances first — no borrows beyond this point
    let dist_first = first_half_indices.iter()
        .map(|&i| contour.points[i].distance_to(&geometry.reference_point))
        .sum::<f64>();

    let dist_second = second_half_indices.iter()
        .map(|&i| contour.points[i].distance_to(&geometry.reference_point))
        .sum::<f64>();

    let use_first = dist_first < dist_second;

    // borrow checker complained, maybe better way
    let mut new_contour = contour.clone();

    for (i, pt) in new_contour.points.iter_mut().enumerate() {
        pt.aortic = if use_first {
            first_half_indices.contains(&i)
        } else {
            second_half_indices.contains(&i)
        };
    }

    (p1.clone(), p2.clone(), new_contour)
}

/// takes a contour (should be contour with highest index)
/// aligns it vertically, and ensures aortic is to the right
/// returns the angle, rotation and the new contour
fn rotate_reference_contour(
    p1: ContourPoint,
    p2: ContourPoint,
    contour: Contour,
) -> (f64, f64, Contour) {
    // Find a rotation that makes contour vertical and ensures aortic points are on the right side
    let dx = p2.x - p1.x;
    let dy = p2.y - p1.y;
    let line_angle = dy.atan2(dx);
    let mut rotation_to_y = (PI / 2.0) - line_angle;
    
    let ref_contour_copy = contour.clone();

    let ((p3, p4), _dist) = ref_contour_copy.find_closest_opposite();
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
    let mut rotated_ref = contour.clone(); //
    rotated_ref.rotate_contour(rotation_to_y);
    rotated_ref.sort_contour_points();

    (line_angle, rotation_to_y, rotated_ref)
}

fn align_remaining_contours(
    mut geometry: Geometry,
    ref_idx: u32,
    ref_contour: Contour,
    rot: f64,
    steps: usize,
    range: f64,
) -> (Geometry, Vec<(u32, (f64, f64, f64), f64, (f64, f64))>) {
    let mut processed_refs = std::collections::HashMap::new();
    let mut id_translation = Vec::new();

    // Process contours in reverse order (highest ID first)
    for contour in geometry.contours.iter_mut().rev() {
        if contour.id == ref_idx {
            continue;
        }

        // Determine reference points and centroid
        let (ref_points, ref_centroid) = match processed_refs.get(&(contour.id + 1)) {
            Some((points, centroid)) => (points, centroid),
            None => (&ref_contour.points, &ref_contour.centroid),
        };

        // Initial rotation alignment
        contour.rotate_contour(rot);

        // Calculate translation
        let tx = contour.centroid.0 - ref_centroid.0;
        let ty = contour.centroid.1 - ref_centroid.1;
        contour.translate_contour((tx, ty, 0.0));

        // Optimize rotation
        let best_rot = find_best_rotation(
            ref_points,
            &contour.points,
            steps,
            range,
            &contour.centroid,
        );

        // Apply final rotation
        let total_rotation = rot + best_rot;
        contour.rotate_contour(best_rot);
        contour.sort_contour_points();

        // Store transformation data
        id_translation.push((
            contour.id,
            (tx, ty, 0.0),
            total_rotation,
            (contour.centroid.0, contour.centroid.1),
        ));
        println!("Matching Contour {:?} -> Contour {:?}, Best rotation: {:?}", &contour.id, contour.id + 1, &best_rot);

        // Update reference tracking
        processed_refs.insert(
            contour.id,
            (contour.points.clone(), contour.centroid),
        );

        // Mark aortic points
        let half_len = contour.points.len() / 2;
        for pt in contour.points.iter_mut().skip(half_len) {
            pt.aortic = true;
        }
    }

    (geometry, id_translation)
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

#[cfg(test)]
mod alignment_tests {
    use super::*;
    use approx::assert_relative_eq;
    use serde_json::Value;
    use std::fs::File;

    fn load_test_manifest(mode: &str) -> Value {
        let manifest_path = format!(
            "python_src/test_geometries/output/{}_csv_files/test_manifest.json", 
            mode
        );
        let file = File::open(manifest_path).expect("Failed to open manifest");
        serde_json::from_reader(file).expect("Failed to parse manifest")
    }

    #[test]
    fn test_aortic_assignment() {
        let geometry = Geometry::new(
            "python_src/test_geometries/output/rest_csv_files",
            "test".to_string(),
            true
        ).expect("Failed to load geometry");

        let n = geometry.contours.len();
        let ref_contour = &geometry.contours[n -1];

        let (_p1, _p2, contour) = assign_aortic(ref_contour.clone(), &geometry);

        // Find the point closest to (3.7, 5.3) -> Numbers from plot in python script
        let aortic_pt = contour.points.iter()
            .min_by(|a, b| {
                let da = ((a.x - 3.7).powi(2) + (a.y - 5.3).powi(2)).partial_cmp(
                    &((b.x - 3.7).powi(2) + (b.y - 5.3).powi(2))
                ).unwrap();
                da
            })
            .expect("No points found near (3.7, 5.3)");
        assert!(aortic_pt.aortic, "Point near (3.7, 5.3) should be aortic");

        // Find the point closest to (5.01, 3.6) -> Numbers from plot in python script
        let nonaortic_pt = contour.points.iter()
            .min_by(|a, b| {
                let da = ((a.x - 5.01).powi(2) + (a.y - 3.6).powi(2)).partial_cmp(
                    &((b.x - 5.01).powi(2) + (b.y - 3.6).powi(2))
                ).unwrap();
                da
            })
            .expect("No points found near (5.01, 3.6)");
        assert!(!nonaortic_pt.aortic, "Point near (5.01, 3.6) should NOT be aortic");
    }

    #[test]
    fn test_rotate_reference_contour() {
        let manifest = load_test_manifest("rest");
        let dia_expected: Vec<u32> = manifest["dia"]["rotations"]
            .as_array().unwrap()
            .iter()
            .map(|v| v.as_u64().unwrap() as u32)
            .collect();
        
        let geometry = Geometry::new(
            "python_src/test_geometries/output/rest_csv_files",
            "test".to_string(),
            true
        ).expect("Failed to load geometry");
        
        let n = geometry.contours.len();
        let ref_contour = &geometry.contours[n -1];
        // let rotation_exp = (dia_expected[n - 1] as f64).to_radians();
        let rotation_exp = dia_expected[n - 1] as f64;

        let (p1, p2, contour) = assign_aortic(ref_contour.clone(), &geometry);

        let (_line_angle, rotation_to_y, rotated_contour) = rotate_reference_contour(p1, p2, contour);

        // Python creates the rotation counterclockwise, Rust has it clockwise.
        let rotation_to_y = 360.0 - rotation_to_y.to_degrees();
        assert_relative_eq!(rotation_to_y, rotation_exp, epsilon = 1e-2);

        // Additionally assert aortic is on the "right" (higher cummulative x-values)
        // split points by their `aortic` flag
        let (sum_aortic, count_aortic) = rotated_contour.points.iter()
            .filter(|pt| pt.aortic)
            .fold((0.0, 0), |(sum, cnt), pt| (sum + pt.x, cnt + 1));
        let (sum_non, count_non) = rotated_contour.points.iter()
            .filter(|pt| !pt.aortic)
            .fold((0.0, 0), |(sum, cnt), pt| (sum + pt.x, cnt + 1));
        
        // compute means (or you could just compare sums directly)
        let mean_aortic = sum_aortic / (count_aortic as f64);
        let mean_non   = sum_non   / (count_non   as f64);

        assert!(
            mean_aortic > mean_non,
            "Expected aortic half to lie to the right: mean_aortic={} ≤ mean_non={}",
            mean_aortic,
            mean_non
        );
    }
    
    #[test]
    fn test_align_remaining_contours() {
        let mode = "dia";
        let manifest = load_test_manifest("rest");
        let dia_config = &manifest[mode];

        // Load expected values from manifest
        let exp_trans: Vec<[f64; 2]> = dia_config["translations"]
            .as_array().unwrap()
            .iter()
            .map(|pair| {
                let arr = pair.as_array().unwrap();
                [arr[0].as_f64().unwrap(), arr[1].as_f64().unwrap()]
            })
            .collect();
        
        let exp_rots_deg: Vec<f64> = dia_config["rotations"]
            .as_array().unwrap()
            .iter()
            .map(|r| r.as_f64().unwrap())
            .collect();

        // Load and align geometry
        let input_dir = "python_src/test_geometries/output/rest_csv_files";
        let geometry = Geometry::new(input_dir, "test".into(), true).unwrap();
        
        let (mut geometry, reference_index, reference_pos, ref_contour) = 
            prep_data_geometry(geometry);
        
        // copy paste from above fix later
        let ref_contour_copy = ref_contour.clone();
        
        let (p1, p2, updated_ref) = 
            assign_aortic(ref_contour_copy, &geometry);
        let ref_contour = updated_ref.clone();

        let (_line_angle, rotation_to_y, rotated_ref) = 
            rotate_reference_contour(p1, p2, ref_contour.clone());

        // Update the reference contour in geometry with rotated version
        geometry.contours.insert(reference_pos, rotated_ref.clone());

        // prepare reference catheter
        let ref_catheter = &mut geometry.catheter.remove(reference_pos);
        ref_catheter.rotate_contour_around_point(
            rotation_to_y, 
            (ref_contour.centroid.0, ref_contour.centroid.1)
        );
        ref_catheter.sort_contour_points();
        geometry.catheter.insert(reference_pos, ref_catheter.clone());

        let (_, id_translation) = align_remaining_contours(
            geometry,
            reference_index,
            rotated_ref,
            rotation_to_y,
            300, // from config
            1.57, // from config
        );
        let mut ids = Vec::new();
        let mut rots = Vec::new();
        for t in id_translation.iter() {
            ids.push(t.0);
            let rot = t.2;
            let rot = rot.to_degrees();
            rots.push(rot);
        }
        println!("Rotations per frame indices: {:?} {:?}", ids, rots);
        // Assertions

        for &(id, (tx, ty, _), rot_rad, _) in &id_translation {
            let index = id as usize;
            let [exp_tx, exp_ty] = exp_trans[index];
            let exp_rot = exp_rots_deg[index];

            // single-line asserts to avoid indentation errors
            assert_relative_eq!(tx, exp_tx, epsilon = 1e-2);
            assert_relative_eq!(ty, exp_ty, epsilon = 1e-2);
            let rot_deg = 360.0 - rot_rad.to_degrees();
            assert_relative_eq!(rot_deg, exp_rot, epsilon = 1e-2);
            }
    }
}