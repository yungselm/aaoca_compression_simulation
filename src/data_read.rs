use serde::Deserialize;
use std::error::Error;
use std::fs::File;
use std::path::Path;

#[derive(Debug, Deserialize, Clone)]
#[allow(dead_code)]
pub struct ContourPoint {
    pub frame_index: u32,
    pub x: f64,
    pub y: f64,
    pub z: f64,
}


pub fn read_contour_data<P: AsRef<Path>>(path: P) -> Result<Vec<ContourPoint>, Box<dyn Error>> {
    let file = File::open(path)?;
    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(false) // No headers in the file
        .delimiter(b'\t')   // Use tab as the delimiter
        .from_reader(file);

    let mut points = Vec::new();

    for result in rdr.records() {
        match result {
            Ok(record) => {
                match record.deserialize(None) {
                    Ok(point) => points.push(point),
                    Err(e) => eprintln!("Skipping invalid record: {:?}", e),
                }
            }
            Err(e) => eprintln!("Skipping invalid row: {:?}", e),
        }
    }

    Ok(points)
}