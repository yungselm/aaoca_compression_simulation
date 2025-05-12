use std::ops::RangeInclusive;

use crate::io::input::{Contour, ContourPoint};
use crate::io::Geometry;

pub fn create_wall_geometry(
    geometry: &Geometry,
    with_pulmonary: bool,
) -> Geometry {
    let mut new_contours = Vec::new();

    for contour in &geometry.contours {
        let new_contour = if with_pulmonary {
            create_wall_contour_with_pulmonary(contour)
        } else {
            create_wall_contour_aortic_only(contour)
        };

        new_contours.push(new_contour);
    }

    let new_geometry = Geometry {
        contours: new_contours,
        catheter: geometry.catheter.clone(),
        reference_point: geometry.reference_point.clone(),
        label: geometry.label.clone(),
    };

    new_geometry
}

fn create_wall_contour_aortic_only(
    contour: &Contour,
) -> Contour {
    if contour.aortic_thickness.is_none() {
        let new_contour = enlarge_contour(contour, 1.5, None);
        new_contour
    } else {
        let new_contour = create_aortic_wall(contour);
        new_contour
    }
}

fn create_wall_contour_with_pulmonary(
    contour: &Contour,
) -> Contour {
    todo!("not yet implemented");
}

pub fn enlarge_contour(
    contour: &Contour,
    factor: f64,
    point_range: Option<RangeInclusive<u32>>,
) -> Contour {
    let (cx, cy, cz) = contour.centroid;

    let new_points = contour.points.iter().map(|point| {
        let mut new_point = point.clone();

        // Only scale if point_index is in the given range (or if no range was passed)
        let do_scale = point_range
            .as_ref()
            .map(|r| r.contains(&point.point_index))
            .unwrap_or(true);

        if do_scale {
            let dx = point.x - cx;
            let dy = point.y - cy;
            let dz = point.z - cz;

            new_point.x = cx + factor * dx;
            new_point.y = cy + factor * dy;
            new_point.z = cz + factor * dz;
        }

        new_point
    });

    Contour {
        id: contour.id,
        points: new_points.collect(),
        centroid: (cx, cy, cz),
        aortic_thickness: contour.aortic_thickness,
        pulmonary_thickness: contour.pulmonary_thickness,
    }
}

/// Create aortic wall contour with the thickness from IVUS images
/// assumption: using preexisting point indices since the contour is already aligned
/// 0 is the highest point, 250 is the lowest point, point 375 used to build the thickness
/// by adding the distance to point x-coordinates. y-coordinates aortic side are +1.5
/// to the y-coodinate of point 0 and -1.5 to the y-coordinates of point 250. The side
/// with indices 0-250 will use the enlarge contour function to create a coronary wall.
/// For easier geometry creation again 500 points are used.
fn create_aortic_wall(contour: &Contour) -> Contour {
    // Key measurements
    let ref_pt    = contour.points[375];
    let thickness = contour.aortic_thickness
        .expect("aortic_thickness[375] must be Some");
    let outer_x   = ref_pt.x + thickness;
    let y_up      = contour.points[0].y   + 1.5;
    let y_low     = contour.points[250].y - 1.5;

    // Left half: indices 0..=250 inclusive
    let left_wall = enlarge_contour(contour, 1.5, Some(0..=250)).points;
    let left_len  = left_wall.len();            // 251

    // Total length of your new contour = original contour length
    let total_pts = contour.points.len();       // 501
    let right_len = total_pts - left_len;       // 250

    // Right half: interpolate between (outer_x, y_up) → (outer_x, y_low)
    let mut right_wall = Vec::with_capacity(right_len);
    for i in 0..right_len {
        // t = 0.0 .. 1.0 over 250 points
        let t = i as f64 / (right_len - 1) as f64;
        let y = y_up * (1.0 - t) + y_low * t;

        // original source point to copy frame_index, point_index, z, aortic flag
        let src = &contour.points[left_len + i];

        right_wall.push( ContourPoint {
            frame_index: src.frame_index,
            point_index: src.point_index,
            x: outer_x,
            y,
            z: src.z,
            aortic: src.aortic,
        });
    }

    // Stitch them together
    let mut new_points = Vec::with_capacity(total_pts);
    new_points.extend(left_wall.into_iter());
    new_points.extend(right_wall.into_iter());

    Contour {
        id: contour.id,
        points: new_points,
        centroid: contour.centroid,
        aortic_thickness: contour.aortic_thickness,
        pulmonary_thickness: contour.pulmonary_thickness,
    }
}