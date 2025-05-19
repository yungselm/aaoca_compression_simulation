pub mod input;
pub mod load_geometry;
pub mod output;

use input::{read_records, Contour, ContourPoint, Record};
use std::path::Path;

#[derive(Debug, Clone)]
pub struct Geometry {
    pub contours: Vec<Contour>,
    pub catheter: Vec<Contour>,
    pub reference_point: ContourPoint, // needs to be set on aortic wall ostium!
    pub label: String,
}

impl Geometry {
    /// Creates a new Geometry instance by loading all required data files
    pub fn new(input_dir: &str, label: String, diastole: bool) -> anyhow::Result<Self> {
        let label = if diastole {
            format!("{}_diastole", label)
        } else {
            format!("{}_systole", label)
        };

        let base_path = Path::new(input_dir);
        let records_path = Path::new(base_path).join("combined_sorted_manual.csv");
        let (contour_path, reference_path) = if diastole {
            (
                Path::new(base_path).join("diastolic_contours.csv"),
                Path::new(base_path).join("diastolic_reference_points.csv"),
            )
        } else {
            (
                Path::new(base_path).join("systolic_contours.csv"),
                Path::new(base_path).join("systolic_reference_points.csv"),
            )
        };

        // Load core components
        let mut contours = Self::load_contours(&contour_path, &records_path)?;
        println!("Loaded contours");
        let reference_point = Self::load_reference_point(&reference_path)?;
        println!("Loaded reference_point");
        let records = Self::load_results(&records_path)?;
        println!("Loaded results");
        
        // since reordeing the frames, destroys the z-coordinates of everyframe they need to be stored here
        // and then be reused after reordering them
        let mut z_coords = Vec::new();

        for contour in &contours {
            z_coords.push(contour.centroid.2)
        }
        // order z_coords by f64 ascending
        z_coords.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

        Self::reorder_contours(&mut contours, &records, diastole, &z_coords);

        let mut catheter = Contour::create_catheter_contours(
            &contours
                .iter()
                .flat_map(|c| c.points.clone())
                .collect::<Vec<_>>(),
        )?;
        
        //sort catheter in ascending order
        catheter.sort_by_key(|c| c.id);

        println!("Created catheter contours");
        println!("Generating geometry for {:?}", input_dir);

        Ok(Self {
            contours,
            catheter,
            reference_point,
            label,
        })
    }

    fn load_contours(
        contour_path: &Path,
        records_path: &Path,
    ) -> anyhow::Result<Vec<Contour>> {
        let raw_points = ContourPoint::read_contour_data(contour_path)?;
        let results = read_records(records_path)?;

        // Create validated contours with measurements
        Contour::create_contours(raw_points, results)
    }

    fn load_reference_point(reference_path: &Path) -> anyhow::Result<ContourPoint> {
        ContourPoint::read_reference_point(reference_path)
    }

    fn load_results(records_path: &Path) -> anyhow::Result<Vec<Record>> {
        read_records(records_path)
    }

    /// Reorders contours by record frame order, updates z-coordinates and ids
    fn reorder_contours(
        contours: &mut Vec<Contour>,
        records: &[Record],
        diastole: bool,
        z_coords: &[f64],
    ) {
        // reorder contours by records frame order, first filter only phase == 'D' if diastole true
        // otherwise only phase == 'S'
        let phase = if diastole { "D" } else { "S" };
        let filtered: Vec<u32> = records
            .iter()
            .filter(|r| r.phase == phase)
            .map(|r| r.frame)
            .collect();

        println!("Records id order: {:?}", &filtered);
        let mut original_contour_ids = Vec::new();
        for c in contours.iter() {
            original_contour_ids.push(c.id)
        }
        println!("Original contour id order: {:?}", original_contour_ids);

        // Sort contours to match filtered record frame order
        contours.sort_by_key(|c| {
            filtered
                .iter()
                .position(|&f| f == c.id)
                .unwrap_or(usize::MAX)
        });

        let mut post_contour_ids = Vec::new();
        for c in contours.iter() {
            post_contour_ids.push(c.id)
        }
        println!("Post contour id order: {:?}", post_contour_ids);

        // Update the z-coordinates of contours and their points using z_coords
        for (i, contour) in contours.iter_mut().enumerate() {
            // Replace the centroid's z-coordinate
            contour.centroid.2 = z_coords[i];

            // Replace the z-coordinate for every point in the contour
            for pt in contour.points.iter_mut() {
                pt.z = z_coords[i];
            }
        }

        // Reassign indices for contours and update their points' frame_index accordingly.
        for (new_id, contour) in contours.iter_mut().enumerate() {
            contour.id = new_id as u32;
            for pt in contour.points.iter_mut() {
                pt.frame_index = new_id as u32;
            }
        } // new order has now highest index for the ostium
    }

    /// Smooths the x and y coordinates of the contours using a 3‐point moving average.
    ///
    /// For each point i in contour j, the new x and y values are computed as:
    ///     new_x = (prev_contour[i].x + current_contour[i].x + next_contour[i].x) / 3.0
    ///     new_y = (prev_contour[i].y + current_contour[i].y + next_contour[i].y) / 3.0
    /// while the z coordinate remains unchanged (taken from the current contour).
    ///
    /// For the first and last contours, the current contour is used twice to simulate a mirror effect.
    pub fn smooth_contours(mut self) -> Geometry {
        let n = self.contours.len();
        if n == 0 {
            return self;
        }

        // Ensure all contours have the same number of points.
        let point_count = self.contours[0].points.len();
        if !self.contours.iter().all(|c| c.points.len() == point_count) {
            panic!("All contours must have the same number of points for smoothing.");
        }

        let mut smoothed_contours = Vec::with_capacity(n);

        for j in 0..n {
            let current_contour = &self.contours[j];
            let mut new_points = Vec::with_capacity(current_contour.points.len());

            for i in 0..current_contour.points.len() {
                let (prev_contour, next_contour) = if j == 0 {
                    // First contour: use current for previous
                    (&self.contours[j], &self.contours[j + 1])
                } else if j == n - 1 {
                    // Last contour: use current for next
                    (&self.contours[j - 1], &self.contours[j])
                } else {
                    (&self.contours[j - 1], &self.contours[j + 1])
                };

                let prev_point = &prev_contour.points[i];
                let curr_point = &current_contour.points[i];
                let next_point = &next_contour.points[i];

                let avg_x = (prev_point.x + curr_point.x + next_point.x) / 3.0;
                let avg_y = (prev_point.y + curr_point.y + next_point.y) / 3.0;

                let new_point = ContourPoint {
                    frame_index: curr_point.frame_index,
                    point_index: curr_point.point_index,
                    x: avg_x,
                    y: avg_y,
                    z: curr_point.z,
                    aortic: curr_point.aortic,
                };

                new_points.push(new_point);
            }

            // Create new Contour with smoothed points and updated centroid
            let centroid = Contour::compute_centroid(&new_points);
            let new_contour = Contour {
                id: current_contour.id,
                points: new_points,
                centroid,
                aortic_thickness: current_contour.aortic_thickness.clone(),
                pulmonary_thickness: current_contour.pulmonary_thickness.clone(),
            };

            smoothed_contours.push(new_contour);
        }
        // Replace the original contours with smoothed_contours
        self.contours = smoothed_contours;

        self
    }
}

#[cfg(test)]
mod geometry_tests {
    use super::*;
    use approx::assert_relative_eq;
    use serde_json::Value;
    use std::fs::File;

    const NUM_POINTS_CATHETER: usize = 20;

    fn load_test_manifest(mode: &str) -> Value {
        let manifest_path = format!(
            "python_src/test_geometries/output/{}_csv_files/test_manifest.json", 
            mode
        );
        let file = File::open(manifest_path).expect("Failed to open manifest");
        serde_json::from_reader(file).expect("Failed to parse manifest")
    }

    impl Contour {
        pub fn elliptic_ratio(&self) -> f64 {
            let major_length = self.find_farthest_points().1;
            let minor_length = self.find_closest_opposite().1;
            // Ensure major is always larger than minor
            if major_length < minor_length {
                minor_length / major_length
            } else {
                major_length / minor_length
            }
        }

        pub fn area(&self) -> f64 {
            let major_length = self.find_farthest_points().1;
            let minor_length = self.find_closest_opposite().1;
            let a = major_length / 2.0;
            let b = minor_length / 2.0;
            std::f64::consts::PI * a * b
        }
    }

    #[test]
    fn test_reorder_matches_manifest_indices() {
        let mode = "rest";
        let manifest = load_test_manifest(mode);
        let dia_expected: Vec<u32> = manifest["dia"]["expected_indices"]
            .as_array().unwrap()
            .iter()
            .map(|v| v.as_u64().unwrap() as u32)
            .collect();

        // Load raw records and geometry
        let input_dir = format!("python_src/test_geometries/output/{0}_csv_files", mode);
        let geometry = Geometry::new(&input_dir, "test".into(), true).unwrap();
        let records = Geometry::load_results(&Path::new(&input_dir).join("combined_sorted_manual.csv")).unwrap();
        let filtered: Vec<u32> = records.into_iter()
            .filter(|r| r.phase == "D")
            .map(|r| r.frame)
            .collect();

        // Map reordered contours back to original frame indices
        let actual_sequence: Vec<u32> = geometry.contours.iter()
            .map(|c| filtered[c.id as usize])
            .collect();

        for (i, (got, want)) in actual_sequence.iter().zip(&dia_expected).enumerate() {
            assert_eq!(
                got, want,
                "Mismatch at position {}: got frame {} but expected {}",
                i, got, want
            );
        }
    }

    #[test]
    fn test_rest_diastolic_config_match() {
        let geometry = Geometry::new(
            "python_src/test_geometries/output/rest_csv_files",
            "test".to_string(),
            true
        ).expect("Failed to load geometry");
        
        let manifest = load_test_manifest("rest");
        let dia_config = &manifest["dia"];

        // Test contour count
        assert_eq!(
            geometry.contours.len(), 
            dia_config["num_contours"].as_u64().unwrap() as usize,
            "Contour count mismatch"
        );
        let n = geometry.contours.len() as u32;
        // Test frame indices ordering
        let expected_indices: Vec<u32> = (0..=(n - 1)).collect();
        
        let actual_indices: Vec<u32> = geometry.contours.iter()
            .map(|c| c.id)
            .collect();

        assert_eq!(
            actual_indices, expected_indices,
            "Frame indices ordering mismatch"
        );
    }

    #[test]
    fn test_contour_property_consistency() {
        let geometry = Geometry::new(
            "python_src/test_geometries/output/rest_csv_files",
            "test".to_string(),
            true
        ).expect("Failed to load geometry");
        
        let manifest = load_test_manifest("rest");
        let dia_config = &manifest["dia"];

        for (i, contour) in geometry.contours.iter().enumerate().rev() {
            // Verify elliptic ratio
            let expected_ratio = dia_config["elliptic_ratios"][i].as_f64().unwrap();
            assert_relative_eq!(
                contour.elliptic_ratio(),
                expected_ratio,
                epsilon = 0.1  // Allow some tolerance
            );

            // Verify area
            let expected_area = dia_config["areas"][i].as_f64().unwrap();
            assert_relative_eq!(
                contour.area(),
                expected_area,
                epsilon = 0.1
            );

            // Verify aortic thickness
            let expected_thickness = match dia_config["aortic_thickness"][i].as_f64() {
                Some(v) => Some(v),
                None => None,
            };
            assert_eq!(
                contour.aortic_thickness, 
                expected_thickness,
                "Aortic thickness mismatch at index {}",
                i
            );
        }
    }

    #[test]
    fn test_catheter_contour_properties() {
        let geometry = Geometry::new(
            "python_src/test_geometries/output/rest_csv_files",
            "test".to_string(),
            true
        ).expect("Failed to load geometry");

        // Verify number of catheter points per contour
        for catheter_contour in &geometry.catheter {
            assert_eq!(
                catheter_contour.points.len(),
                NUM_POINTS_CATHETER,
                "Incorrect number of catheter points"
            );
        }

        // Verify z-coordinate consistency
        for (contour, catheter) in geometry.contours.iter().zip(&geometry.catheter) {
            assert_relative_eq!(
                catheter.centroid.2,
                contour.centroid.2,
                epsilon = 1e-6
            );
        }
    }
}