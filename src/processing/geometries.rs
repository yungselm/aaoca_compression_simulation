use rayon::prelude::*;

use crate::io::Geometry;
use crate::processing::contours::hausdorff_distance;

use super::contours::align_frames_in_geometry;
use crate::io::input::{ContourPoint, Contour};

#[derive(Clone)]
pub struct GeometryPair {
    pub dia_geom: Geometry,
    pub sys_geom: Geometry,
}

impl GeometryPair {
    pub fn new(input_dir: &str, label: String) -> anyhow::Result<GeometryPair> {
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

        Self::translate_contours_to_match(&diastole, &mut systole);

        // Adjust the z-coordinates of systolic contours. (later replaceed by adjust_z_coordinates)
        Self::apply_z_transformation(&diastole, &mut systole);

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

    fn translate_contours_to_match(dia: &Geometry, sys: &mut Geometry) {
        let dia_ref = dia.contours.last().unwrap().centroid;
        let sys_ref = sys.contours.last().unwrap().centroid;
        let offset = (dia_ref.0 - sys_ref.0, dia_ref.1 - sys_ref.1);

        for contour in &mut sys.contours {
            contour.translate_contour((offset.0, offset.1, 0.0));
        }

        for catheter in &mut sys.catheter {
            catheter.translate_contour((offset.0, offset.1, 0.0));
        }
    }

    fn apply_z_transformation(dia: &Geometry, sys: &mut Geometry) {
        let dia_last_z = dia.contours.last().unwrap().centroid.2;
        let sys_last_z = sys.contours.last().unwrap().centroid.2;
        let z_offset = dia_last_z - sys_last_z;

        for contour in &mut sys.contours {
            contour.points.iter_mut().for_each(|p| p.z += z_offset);
            contour.centroid.2 += z_offset;
        }

        for catheter in &mut sys.catheter {
            catheter.points.iter_mut().for_each(|p| p.z += z_offset);
            catheter.centroid.2 += z_offset;
        }
    }

    pub fn adjust_z_coordinates(mut self) -> GeometryPair {
        let mut z_coords_dia: Vec<f64> = self
            .dia_geom
            .contours
            .iter()
            .skip(1) // Skip the first entry since 0.0
            .map(|contour| contour.centroid.2)
            .collect();

        let mut z_coords_sys: Vec<f64> = self
            .sys_geom
            .contours
            .iter()
            .skip(1) // Skip the first entry sicne 0.0
            .map(|contour| contour.centroid.2)
            .collect();
        
        for i in (0..z_coords_dia.len()).rev() {
            z_coords_dia[i] /= (i + 1) as f64;
        }

        for i in (0..z_coords_sys.len()).rev() {
            z_coords_sys[i] /= (i + 1) as f64;
        }

        let mut z_coords = z_coords_sys;
        z_coords.extend(z_coords_dia);

        let mean_z_coords = z_coords.iter().sum::<f64>() / z_coords.len() as f64;

        // If there are missing frames in between this will create false results, but probably
        // still more accurate then taking the actual frame position due to breathing artefacts
        // and the resampling performed in combined_sorted_manual to counter this. 
        // 4) Assign evenly spaced Z in parallel across DIA & SYS contours and catheters
        let n_slices = 
            self.dia_geom.contours.len()
            .max(self.sys_geom.contours.len())
            .max(self.dia_geom.catheter.len())
            .max(self.sys_geom.catheter.len());

        let mut current_z = 0.0;
        for i in 0..n_slices {
            // helper to set z for a mutable slice element
            let assign_z = |cont_opt: Option<&mut Contour>| {
                if let Some(cont) = cont_opt {
                    cont.centroid.2 = current_z;
                    for pt in &mut cont.points {
                        pt.z = current_z;
                    }
                }
            };

            assign_z(self.dia_geom.contours.get_mut(i));
            assign_z(self.sys_geom.contours.get_mut(i));
            assign_z(self.dia_geom.catheter.get_mut(i));
            assign_z(self.sys_geom.catheter.get_mut(i));

            println!("layer {} => z = {}", i, current_z);
            current_z += mean_z_coords;
        }
      
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

    /// Adjusts the aortic and pulmonary thicknesses of the contours in both geometries
    /// to be the average of the two. This is done for each contour in both geometries.
    /// The function ensures that the lengths of the thickness vectors are equal by resizing
    /// them to the maximum length found in either geometry. The average is calculated
    /// for each corresponding element in the vectors.
    /// If one of the elements is None, it takes the value from the other element.
    /// If both are None, it remains None.
    /// This function is called after the geometries have been aligned and translated.
    /// It is important to ensure that the geometries are aligned before calling this function.
    /// The function assumes that the contours in both geometries are in the same order.
    /// It does not check for matching IDs, so it is the caller's responsibility to ensure
    /// that the contours correspond to the same anatomical structures.
    pub fn thickness_adjustment(mut self) -> GeometryPair {
        let min_contours = std::cmp::min(
            self.dia_geom.contours.len(),
            self.sys_geom.contours.len(),
        );
        for i in 0..min_contours {
            let dia = &mut self.dia_geom.contours[i];
            let sys = &mut self.sys_geom.contours[i];

            // Combine aortic thickness
            let combined_aortic = match (dia.aortic_thickness, sys.aortic_thickness) {
                (Some(a), Some(b)) => Some((a + b) / 2.0),
                (Some(a), None   ) => Some(a),
                (None   , Some(b)) => Some(b),
                (None   , None   ) => None,
            };
            dia.aortic_thickness = combined_aortic;
            sys.aortic_thickness = combined_aortic;

            // Combine pulmonary thickness
            let combined_pulmonary = match (dia.pulmonary_thickness, sys.pulmonary_thickness) {
                (Some(a), Some(b)) => Some((a + b) / 2.0),
                (Some(a), None   ) => Some(a),
                (None   , Some(b)) => Some(b),
                (None   , None   ) => None,
            };
            dia.pulmonary_thickness = combined_pulmonary;
            sys.pulmonary_thickness = combined_pulmonary;
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
    println!("---------------------Finding optimal rotation Diastole/Systole---------------------");
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
