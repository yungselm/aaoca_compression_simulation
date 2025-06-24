// use std::f64::consts::PI;
// use crate::io::input::{Contour, ContourPoint};
// use crate::io::Geometry;

// pub fn generate_aligned_geometry(
//     majors: &[f64],
//     minors: &[f64],
//     num_points: usize,
//     rotations: &[f64],
//     translations: &[(f64, f64)],
// ) -> Geometry {
//     let mut contours = Vec::new();

//     for (i, ((&major, &minor), (&rotation, &translation))) in majors.iter()
//         .zip(minors.iter())
//         .zip(rotations.iter().zip(translations.iter()))
//         .enumerate()
//     {
//         let points = generate_ellipse_points(major, minor, num_points, rotation, translation, i as u32);
//         let contour = Contour {
//             id: i as u32,
//             points,
//             centroid: (0.0, 0.0, 0.0),
//             aortic_thickness: None,
//             pulmonary_thickness: None,
//         };
//         contours.push(contour);
//     }

//     Geometry {
//         contours,
//         catheter: None,
//         reference_point: None,
//         label: None,
//     }
// }

// /// Generates ellipse contour points for testing
// pub fn generate_ellipse_points(
//     major: f64,
//     minor: f64,
//     num_points: usize,
//     rotation: f64,
//     translation: (f64, f64),
//     frame_idx: u32,
// ) -> Vec<ContourPoint> {
//     let mut points = Vec::with_capacity(num_points);
//     for i in 0..num_points {
//         let theta = 2.0 * PI * (i as f64) / (num_points as f64);
//         let x = major * theta.cos();
//         let y = minor * theta.sin();
//         let (x_rot, y_rot) = rotate_point((x, y), rotation);
//         let x_trans = x_rot + translation.0;
//         let y_trans = y_rot + translation.1;
//         points.push(ContourPoint {
//             frame_index: frame_idx,
//             point_index: i as u32,
//             x: x_trans,
//             y: y_trans,
//             z: 0.0,
//             aortic: false,
//         });
//     }
//     points
// }

// /// Rotates a point around origin
// pub fn rotate_point(point: (f64, f64), angle: f64) -> (f64, f64) {
//     let (x, y) = point;
//     let cos = angle.cos();
//     let sin = angle.sin();
//     (x * cos - y * sin, x * sin + y * cos)
// }

// /// Creates a dummy contour for testing
// pub fn new_dummy_contour(id: u32) -> Contour {
//     Contour {
//         id,
//         points: Vec::new(),
//         centroid: (0.0, 0.0, 0.0),
//         aortic_thickness: None,
//         pulmonary_thickness: None,
//     }
// }