pub mod input;
pub mod load_geometry;
pub mod output;

use input::{read_records, Contour, ContourPoint, Record};
use std::error::Error;
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
    pub fn new(input_dir: &str, label: String, diastole: bool) -> Result<Self, Box<dyn Error>> {
        let label = if diastole {
            format!("{}_diastole", label)
        } else {
            format!("{}_systole", label)
        };

        let base_path = Path::new(input_dir);
        let records_path = Path::new(base_path).join("combined_original_manual.csv");
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

        // reorder contours by records frame order, first filter only phase == 'D' if diastole true
        // otherwise only phase == 'S'
        let desired_phase = if diastole { "D" } else { "S" };
        let filtered_records: Vec<_> = records
            .into_iter()
            .filter(|r| r.phase == desired_phase)
            .collect();

        // Sort contours to match filtered record frame order
        let desired_order: Vec<u32> = filtered_records.iter().map(|r| r.frame).collect();
        contours.sort_by_key(|c| {
            let id = c.id as u32;
            desired_order
                .iter()
                .position(|&frame| frame == id)
                .unwrap_or(usize::MAX)
        });

        // Reassign indices for contours and update their points' frame_index accordingly.
        for (new_idx, contour) in contours.iter_mut().enumerate() {
            // Set the new contour id
            contour.id = new_idx as u32;
            // Update frame_index for each point in the contour
            for point in contour.points.iter_mut() {
                point.frame_index = new_idx as u32;
            }
        } // new order has now highest index for the ostium

        let catheter = Contour::create_catheter_contours(
            &contours
                .iter()
                .flat_map(|c| c.points.clone())
                .collect::<Vec<_>>(),
        )?;
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
    ) -> Result<Vec<Contour>, Box<dyn Error>> {
        let raw_points = ContourPoint::read_contour_data(contour_path)?;
        let results = read_records(records_path)?;

        // Create validated contours with measurements
        Contour::create_contours(raw_points, results)
    }

    fn load_reference_point(reference_path: &Path) -> Result<ContourPoint, Box<dyn Error>> {
        ContourPoint::read_reference_point(reference_path)
    }

    fn load_results(records_path: &Path) -> Result<Vec<Record>, Box<dyn Error>> {
        read_records(records_path)
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
