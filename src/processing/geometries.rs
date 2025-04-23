use rayon::prelude::*;

use crate::io::Geometry;
use crate::processing::contours::hausdorff_distance;
use std::error::Error;

use super::contours::align_frames_in_geometry;
use crate::io::input::ContourPoint;

pub struct GeometryPair {
    pub dia_geom: Geometry,
    pub sys_geom: Geometry,
}

impl GeometryPair {
    pub fn new(input_dir: &str, label: String) -> Result<GeometryPair, Box<dyn Error>> {
        let dia_geom = Geometry::new(input_dir, label.clone(), true)?;
        println!("geometry pair: diastolic geometry generated");
        let sys_geom = Geometry::new(input_dir, label, false)?;
        println!("geometry pair: systolic geometry generated");
        Ok(GeometryPair { dia_geom, sys_geom })
    }

    /// aligns the frames in each geomtery by rotating them based on best Hausdorff distance
    /// then translates systolic contours to the diastolic contours, aligns z-coordinates and
    /// trims them to same length.
    pub fn process_geometry_pair(
        self,
        steps_best_rotation: usize,
        range_rotation_rad: f64,
    ) -> GeometryPair {
        let diastole =
            align_frames_in_geometry(self.dia_geom, steps_best_rotation, range_rotation_rad);
        let mut systole =
            align_frames_in_geometry(self.sys_geom, steps_best_rotation, range_rotation_rad);

        let diastolic_ref_centroid = diastole.contours[0].centroid;

        for ref mut contour in systole
            .contours
            .iter_mut()
            .chain(systole.catheter.iter_mut())
        {
            let systolic_centroid = contour.centroid;
            let t = (
                diastolic_ref_centroid.0 - systolic_centroid.0,
                diastolic_ref_centroid.1 - systolic_centroid.1,
            );
            contour.translate_contour((t.0, t.1, 0.0))
        }

        // Adjust the z-coordinates of systolic contours. (later replaceed by adjust_z_coordinates)
        let z_translation = diastole.contours[0].centroid.2 - systole.contours[0].centroid.2;
        for ref mut contour in systole
            .contours
            .iter_mut()
            .chain(systole.catheter.iter_mut())
        {
            for point in contour.points.iter_mut() {
                point.z += z_translation;
            }
            contour.centroid.2 += z_translation;
        }

        let best_rotation_angle = find_best_rotation_all(
            &diastole,
            &systole,
            steps_best_rotation, // number of candidate steps (e.g. 200 or 400)
            range_rotation_rad,  // rotation range (e.g. 1.05 for ~±60°)
        );

        for ref mut contour in systole
            .contours
            .iter_mut()
            .chain(systole.catheter.iter_mut())
        {
            contour.rotate_contour(best_rotation_angle);
        }
        GeometryPair {
            dia_geom: diastole,
            sys_geom: systole,
        }
    }

    pub fn adjust_z_coordinates(mut self) -> GeometryPair {
        let mut z_coords: Vec<f64> = self
            .dia_geom
            .contours
            .iter()
            .chain(self.sys_geom.contours.iter())
            .map(|contour| contour.centroid.2)
            .collect();

        for i in (0..z_coords.len()).rev() {
            z_coords[i] /= (i + 1) as f64;
        }

        let mean_z_coords = z_coords.iter().sum::<f64>() / z_coords.len() as f64;

        self.dia_geom
            .contours
            .iter_mut()
            .chain(self.sys_geom.contours.iter_mut())
            .chain(self.dia_geom.catheter.iter_mut())
            .chain(self.sys_geom.catheter.iter_mut())
            .for_each(|contour| {
                contour.centroid.2 = (contour.centroid.2 / mean_z_coords).floor() * mean_z_coords;
                for point in contour.points.iter_mut() {
                    point.z = (point.z / mean_z_coords).floor() * mean_z_coords;
                }
            });
        self
    }

    pub fn trim_geometries_same_length(mut self) -> GeometryPair {
        let dia_len = self.dia_geom.contours.len();
        let sys_len = self.sys_geom.contours.len();
        let min_len = std::cmp::min(dia_len, sys_len);
    
        if dia_len > min_len {
            self.dia_geom.contours.drain(0..(dia_len - min_len));
        }
        if sys_len > min_len {
            self.sys_geom.contours.drain(0..(sys_len - min_len));
        }
    
        let dia_catheter_len = self.dia_geom.catheter.len();
        let sys_catheter_len = self.sys_geom.catheter.len();
        let min_catheter_len = std::cmp::min(dia_catheter_len, sys_catheter_len);
    
        if dia_catheter_len > min_catheter_len {
            self.dia_geom.catheter.drain(0..(dia_catheter_len - min_catheter_len));
        }
        if sys_catheter_len > min_catheter_len {
            self.sys_geom.catheter.drain(0..(sys_catheter_len - min_catheter_len));
        }
    
        self
    }
}

pub fn find_best_rotation_all(
    diastole: &Geometry,
    systole: &Geometry,
    steps: usize,
    range: f64,
) -> f64 {
    let increment = (2.0 * range) / steps as f64;

    (0..=steps)
        .into_par_iter()
        .map(|i| {
            let angle = -range + i as f64 * increment;
            let total_distance: f64 = diastole
                .contours
                .par_iter()
                .zip(systole.contours.par_iter())
                .map(|(d_contour, s_contour)| {
                    assert_eq!(d_contour.id, s_contour.id, "Mismatched contour IDs");

                    // Rotate each point in systole contour
                    let rotated_points: Vec<ContourPoint> = s_contour
                        .points
                        .iter()
                        .map(|p| {
                            let x = p.x * angle.cos() - p.y * angle.sin();
                            let y = p.x * angle.sin() + p.y * angle.cos();
                            ContourPoint { x, y, ..*p }
                        })
                        .collect();

                    // Compute Hausdorff distance between corresponding contours
                    hausdorff_distance(&d_contour.points, &rotated_points)
                })
                .sum();

            let avg_distance = total_distance / diastole.contours.len() as f64;
            println!("Angle: {:.3} rad, Avg Distance: {:.3}", angle, avg_distance);
            (angle, avg_distance)
        })
        .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
        .map(|(angle, _)| angle)
        .unwrap_or(0.0)
}
