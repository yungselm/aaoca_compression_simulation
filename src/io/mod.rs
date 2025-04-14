pub mod input;
pub mod output;

use std::error::Error;
use std::path::{Path, PathBuf};
use input::{Contour, ContourPoint, Record, read_centerline_txt, read_records};

#[derive(Debug)]
pub struct Geometry<'a> {
    pub contours: Vec<Contour>,
    pub catheter: Vec<Contour>,
    pub reference_point: ContourPoint, // needs to be set on aortic wall ostium!
    pub label: &'a str,
}

impl<'a> Geometry<'a> {
    /// Creates a new Geometry instance by loading all required data files
    pub fn new(
        input_dir: &str,
        label: &'a str,
        diastole: bool,
    ) -> Result<Self, Box<dyn Error>> {
        let base_path = Path::new(input_dir);
        let records_path = Path::new(base_path).join("combined_original_manual.csv");
        let (contour_path, reference_path) = if diastole {
            (Path::new(base_path).join("diastolic_contours.csv"), 
            Path::new(base_path).join("diastolic_reference_points.csv"))
        } else {
            (Path::new(base_path).join("systolic_contours.csv"), 
            Path::new(base_path).join("systolic_reference_points.csv"))
        };
        
        // Load core components
        let mut contours = Self::load_contours(&contour_path, &records_path)?;
        let reference_point = Self::load_reference_point(&reference_path)?;
        let records = Self::load_results(&records_path)?;

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
            desired_order.iter().position(|&frame| frame == id).unwrap_or(usize::MAX)
        });

        // Reassign indices for contours and update their points' frame_index accordingly.
        for (new_idx, contour) in contours.iter_mut().enumerate() {
            // Set the new contour id
            contour.id = new_idx as u32;
            // Update frame_index for each point in the contour
            for point in contour.points.iter_mut() {
                point.frame_index = new_idx as u32;
            }
        }

        let catheter = Contour::create_catheter_contours(
            &contours.iter().flat_map(|c| c.points.clone()).collect::<Vec<_>>()
        )?;

        Ok(Self {
            contours,
            catheter,
            reference_point,
            label,
        })
    }

    fn load_contours(contour_path: &Path, records_path: &Path) -> Result<Vec<Contour>, Box<dyn Error>> {
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
}
