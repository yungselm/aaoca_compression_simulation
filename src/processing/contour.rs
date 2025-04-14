use rayon::prelude::*;
use std::collections::HashMap;
use std::error::Error;
use std::f64::consts::PI;

use crate::io::input::ContourPoint;
use crate::io::Geometry;

pub fn create_and_align_geometries(
    mut geometry: Geometry,
    steps: usize,
    range: f64,
) -> Geometry {
    // Sort contours by frame index.
    geometry.contours.sort_by_key(|contour| contour.id);

    // Use sort_contour_points on all catheter points to ensure they are counterclockwise
    for catheter in &mut geometry.catheter {
        catheter.sort_contour_points();
    }

    // Use the contour with the highest frame index as reference.
    let reference_index = geometry.contours.iter().map(|contour| contour.id).max().unwrap();
    let reference_pos = geometry.contours
        .iter()
        .position(|contour| contour.id == reference_index)
        .expect("Reference contour not found");
    let ref_contour = &mut geometry.contours[reference_pos];
    println!("Using contour {} as reference.", reference_index);

    let ((p1, p2), _dist) = ref_contour.find_farthest_points();
    let dx = p2.x - p1.x;
    let dy = p2.y - p1.y;
    let line_angle = dy.atan2(dx);
    let mut rotation_to_y = (PI / 2.0) - line_angle;
    println!(
        "Reference line angle: {:.3} rad; rotating reference by {:.3} rad",
        line_angle, rotation_to_y
    );
    ref_contour.rotate_contour(rotation_to_y);
    ref_contour.sort_contour_points();

    let n = ref_contour.points.len() / 2;

    // --- Begin aortic determination ---
    // We assume there are exactly 500 points; split into two halves: indices 0..249 and 250..499.
    // Compute the cumulative distance for each half from the rotated reference.
    let first_half_distance: f64 = ref_contour.points
        .iter()
        .take(n)
        .map(|pt| pt.distance_to(&geometry.reference_point))
        .sum();
    let second_half_distance: f64 = ref_contour.points
        .iter()
        .skip(n)
        .map(|pt| pt.distance_to(&geometry.reference_point))
        .sum();

    // Debug output: print average distances for each half.
    println!(
        "First half avg distance: {:.3}, Second half avg distance: {:.3}",
        first_half_distance / (n as f64),
        second_half_distance / (n as f64)
    );

    // Mark the half that is closer as 'aortic' (i.e. set to true).
    if first_half_distance < second_half_distance {
        for pt in ref_contour.points.iter_mut().take(n) {
            pt.aortic = true;
        }
        println!("First half (points 0 to 249) marked as aortic.");
        if first_half_distance < 4.5 {
            rotation_to_y += PI;
            ref_contour.rotate_contour(PI);
        }
    } else {
        for pt in ref_contour.points.iter_mut().skip(n) {
            pt.aortic = true;
        }
        println!("Second half (points 250 to 499) marked as aortic.");
        if second_half_distance < 4.5 {
            rotation_to_y += PI;
            ref_contour.rotate_contour(PI);
        }
    }

    todo!();
}