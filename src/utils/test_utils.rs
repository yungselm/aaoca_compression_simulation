use std::f64::consts::PI;
use crate::io::input::{Contour, ContourPoint, Centerline, CenterlinePoint, Record};
use nalgebra::Vector3;
use crate::io::Geometry;
use crate::processing::geometries::GeometryPair;

/// Generates ellipse contour points for testing
pub fn generate_ellipse_points(
    major: f64,
    minor: f64,
    num_points: usize,
    rotation: f64,
    translation: (f64, f64),
    frame_idx: u32,
) -> Vec<ContourPoint> {
    let mut points = Vec::with_capacity(num_points);
    for i in 0..num_points {
        let theta = 2.0 * PI * (i as f64) / (num_points as f64);
        let x = major * theta.cos();
        let y = minor * theta.sin();
        let (x_rot, y_rot) = rotate_point((x, y), rotation);
        points.push(ContourPoint {
            frame_index: frame_idx,
            point_index: i as u32,
            x: x_rot + translation.0,
            y: y_rot + translation.1,
            z: 0.0,
            aortic: false,
        });
    }
    points
}

/// Rotates a point around origin
pub fn rotate_point(point: (f64, f64), angle: f64) -> (f64, f64) {
    let (x, y) = point;
    let cos = angle.cos();
    let sin = angle.sin();
    (x * cos - y * sin, x * sin + y * cos)
}

/// Creates a dummy Contour with given id
pub fn new_dummy_contour(id: u32) -> Contour {
    Contour {
        id,
        points: Vec::new(),
        centroid: (0.0, 0.0, 0.0),
        aortic_thickness: None,
        pulmonary_thickness: None,
    }
}

// /// Creates a Geometry with aligned ellipse contours
// pub fn generate_aligned_geometry(
//     majors: &[f64],
//     minors: &[f64],
//     num_points: usize,
//     rotations: &[f64],
//     translations: &[(f64, f64)],
//     label: &str,
// ) -> Geometry {
//     let mut contours = Vec::new();
//     for (i, ((&maj, &min), (&rot, &trans))) in majors.iter()
//         .zip(minors.iter())
//         .zip(rotations.iter().zip(translations.iter()))
//         .enumerate() {
//         let pts = generate_ellipse_points(maj, min, num_points, rot, trans, i as u32);
//         contours.push(Contour { id: i as u32, points: pts, centroid: (0.0,0.0,0.0), aortic_thickness: None, pulmonary_thickness: None });
//     }
//     // reference_point: first point of first contour or dummy
//     let reference_point = contours.get(0)
//         .and_then(|c| c.points.get(0)).cloned()
//         .unwrap_or(ContourPoint { frame_index: 0, point_index: 0, x:0.0, y:0.0, z:0.0, aortic:false });

//     Geometry {
//         contours,
//         catheter: Vec::new(),
//         reference_point,
//         label: label.to_string(),
//     }
// }

// /// Creates a dummy Geometry with no contours
// pub fn new_dummy_geometry() -> Geometry {
//     let reference_point = ContourPoint { frame_index: 0, point_index: 0, x: 0.0, y: 0.0, z: 0.0, aortic: false };
//     Geometry { contours: Vec::new(), catheter: Vec::new(), reference_point, label: "dummy".into() }
// }

// /// Creates a pair of geometries for systolic and diastolic phases
// pub fn new_dummy_geometry_pair() -> GeometryPair {
//     GeometryPair {
//         dia_geom: new_dummy_geometry(),
//         sys_geom: new_dummy_geometry(),
//     }
// }

// /// Creates a Centerline with simple points along z-axis
// pub fn new_dummy_centerline(num_points: usize) -> Centerline {
//     let pts = (0..num_points).map(|i| point_on_centerline(i as f64)).collect();
//     Centerline { points: pts }
// }

// fn point_on_centerline(z: f64) -> CenterlinePoint {
//     CenterlinePoint { contour_point: ContourPoint { frame_index: 0, point_index: 0, x: 0.0, y: 0.0, z, aortic: false }, normal: Vector3::new(0.0,0.0,1.0) }
// }

// /// Creates a Record with specified values
// pub fn new_dummy_record(
//     frame: u32,
//     phase: &str,
//     m1: Option<f64>,
//     m2: Option<f64>
// ) -> Record {
//     Record { frame, phase: phase.to_string(), measurement_1: m1, measurement_2: m2 }
// }
