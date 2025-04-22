// const TEST_IDS_ORIGINAL: [u32; 10] = [50, 123, 245, 331, 401, 498, 567, 607, 772, 803];
// const TEST_IDS_SORTED: [u32; 10] = [50, 123, 331, 567, 607, 245, 401, 772, 498, 803];
// const TEST_X_CENTER: [f64; 10] = [4.4, 4.6, 4.3, 4.5, 4.5, 4.4, 4.3, 4.7, 4.5, 4.5];
// const TEST_Y_CENTER: [f64; 10] = [4.6, 4.5, 4.4, 4.5, 4.6, 4.5, 4.7, 4.5, 4.5, 4.5];
// const TEST_Z_CENTER: [f64; 10] = [0.0, 0.9, 1.9, 2.7, 3.7, 4.5, 5.3, 6.2, 7.2, 8.1];

// #[cfg(test)]
// mod test_utils {
//     use super::*;
//     use rand::prelude::*;
//     use crate::io::{input::{ContourPoint, Contour}, Geometry}; // Import missing types

//     pub fn dummy_point() -> ContourPoint {
//         ContourPoint {
//             frame_index: 0,
//             point_index: 0,
//             x: 1.0,
//             y: 2.0,
//             z: 3.0,
//             aortic: false,
//         }
//     }

//     pub fn dummy_contour(id: u32) -> Contour {
//         Contour {
//             id,
//             points: vec![dummy_point()],
//             centroid: (1.0, 2.0, 3.0),
//             aortic_thickness: vec![Some(1.2)],
//             pulmonary_thickness: vec![None],
//         }
//     }

//     pub fn dummy_geometry() -> Geometry {
//         Geometry {
//             contours: vec![dummy_contour(1)],
//             catheter: vec![],
//             reference_point: dummy_point(),
//             label: "TestGeometry".to_string(),
//         }
//     }

//     pub fn dummy_vec_contour_dia(
//         TEST_IDS_ORIGINAL: &[u32; 10],
//         TEST_X_CENTER: &[f64; 10],
//         TEST_Y_CENTER: &[f64; 10],
//         TEST_Z_CENTER: &[f64; 10],
//     ) -> Vec<Contour> {
//         let centroids: Vec<(f64, f64, f64)> = TEST_X_CENTER
//             .iter()
//             .zip(TEST_Y_CENTER.iter())
//             .zip(TEST_Z_CENTER.iter())
//             .map(|((&x, &y), &z)| (x, y, z))
//             .collect();
//         todo!();
//     }

//     pub fn dummy_vec_contour_sys() -> Vec<Contour> {
//         todo!();
//     }

//     pub fn generate_ellipse_contour(
//         id: u32,
//         num_points: usize,
//         origin: (f64, f64, f64),
//         radius_x: f64,
//         radius_y: f64,
//         z: f64,
//     ) -> Contour {
//         let mut points = Vec::with_capacity(num_points);
//         for i in 0..num_points {
//             let angle = 2.0 * std::f64::consts::PI * i as f64 / num_points as f64;
//             let x = origin.0 + radius_x * angle.cos();
//             let y = origin.1 + radius_y * angle.sin();
//             points.push(ContourPoint {
//                 frame_index: 0,
//                 point_index: i as u32,
//                 x,
//                 y,
//                 z: origin.2,
//                 aortic: false,
//             });
//         }

//         Contour {
//             id,
//             points,
//             centroid: origin.0 + radius_x / 2.0, origin.1 + radius_y / 2.0, // rough center
//             aortic_thickness: vec![Some(1.2); num_points],
//             pulmonary_thickness: vec![None; num_points],
//         }
//     }

//     pub fn dummy_geometry_from_test_ids() -> Geometry {
//         let num_points = 501;
//         let mut contours = Vec::new();

//         for (i, id) in TEST_IDS_SORTED.iter().enumerate() {
//             // progressively more elliptical
//             let ratio = 1.0 + (i as f64) * 0.3; // from 1.0 (circle) to ~4.0
//             let origin = (4.5 + i as f64 * 0.05, 4.5 + i as f64 * 0.05);
//             let contour = generate_ellipse_contour(
//                 *id,
//                 num_points,
//                 origin,
//                 1.5 * ratio,       // radius_x
//                 1.5 / ratio,       // radius_y
//                 0.0,
//             );
//             contours.push(contour);
//         }

//         Geometry {
//             contours,
//             catheter: vec![],
//             reference_point: dummy_point(),
//             label: "EllipticalTestGeometry".to_string(),
//         }
//     }
// }
