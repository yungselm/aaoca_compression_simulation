use serde::Deserialize;
use std::error::Error;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

#[derive(Debug, Deserialize, Clone)]
pub struct ContourPoint {
    pub frame_index: u32,
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

/// Reads contour points from a CSV file.
pub fn read_contour_data<P: AsRef<Path>>(path: P) -> Result<Vec<ContourPoint>, Box<dyn Error>> {
    let file = File::open(path)?;
    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(false)
        .delimiter(b'\t')
        .from_reader(file);

    let mut points = Vec::new();
    for result in rdr.records() {
        match result {
            Ok(record) => match record.deserialize(None) {
                Ok(point) => points.push(point),
                Err(e) => eprintln!("Skipping invalid record: {:?}", e),
            },
            Err(e) => eprintln!("Skipping invalid row: {:?}", e),
        }
    }
    Ok(points)
}

/// Writes an OBJ mesh file for the given contours.
pub fn write_obj_mesh(
    contours: &[(u32, Vec<ContourPoint>)],
    filename: &str,
) -> Result<(), Box<dyn Error>> {
    // (We reuse your existing logic here; you might consider refactoring parts further.)
    let sorted_contours = contours.to_owned();

    if sorted_contours.len() < 2 {
        return Err("Need at least two contours to create a mesh.".into());
    }

    let points_per_contour = sorted_contours[0].1.len();
    for (_, contour) in &sorted_contours {
        if contour.len() != points_per_contour {
            return Err("All contours must have the same number of points.".into());
        }
    }

    let file = File::create(filename)?;
    let mut writer = BufWriter::new(file);
    let mut vertex_offsets = Vec::new();
    let mut current_offset = 1;
    let mut normals = Vec::new();

    // Write vertices and compute normals
    for (_, contour) in &sorted_contours {
        vertex_offsets.push(current_offset);
        let centroid = crate::processing::compute_centroid(contour);
        for point in contour {
            writeln!(writer, "v {} {} {}", point.x, point.y, point.z)?;
            let dx = point.x - centroid.0;
            let dy = point.y - centroid.1;
            let length = (dx * dx + dy * dy).sqrt();
            let (nx, ny, nz) = if length > 0.0 {
                (dx / length, dy / length, 0.0)
            } else {
                (0.0, 0.0, 0.0)
            };
            normals.push((nx, ny, nz));
            current_offset += 1;
        }
    }

    // Write normals
    for (nx, ny, nz) in &normals {
        writeln!(writer, "vn {} {} {}", nx, ny, nz)?;
    }

    // Write faces with normals
    for c in 0..(sorted_contours.len() - 1) {
        let offset1 = vertex_offsets[c];
        let offset2 = vertex_offsets[c + 1];
        for j in 0..points_per_contour {
            let j_next = (j + 1) % points_per_contour;
            let v1 = offset1 + j;
            let v2 = offset1 + j_next;
            let v3 = offset2 + j;
            writeln!(writer, "f {}//{} {}//{} {}//{}", v1, v1, v2, v2, v3, v3)?;
            let v1_t2 = offset2 + j;
            let v2_t2 = offset1 + j_next;
            let v3_t2 = offset2 + j_next;
            writeln!(writer, "f {}//{} {}//{} {}//{}", v1_t2, v1_t2, v2_t2, v2_t2, v3_t2, v3_t2)?;
        }
    }

    println!("OBJ mesh with normals written to {}", filename);
    Ok(())
}
