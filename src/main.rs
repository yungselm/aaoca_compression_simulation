mod data_read;
use data_read::read_contour_data;
use std::error::Error;

fn main() -> Result<(), Box<dyn Error>> {
    let contour_points = read_contour_data("input/rest_csv_files/diastolic_contours.csv")?;

    Ok(())
}