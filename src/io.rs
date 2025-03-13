use serde::Deserialize;
use std::error::Error;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::collections::HashMap;

use std::f64::consts::PI;

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

pub fn create_catheter_points(points: &Vec<ContourPoint>) -> Vec<ContourPoint> {
    // Map to store unique frame indices and one associated z coordinate per frame.
    let mut frame_z: HashMap<u32, f64> = HashMap::new();
    for point in points {
        // Use the first encountered z-coordinate for each frame index.
        frame_z.entry(point.frame_index).or_insert(point.z);
    }

    let mut catheter_points = Vec::new();
    // Sort the frame indices to ensure a predictable order.
    let mut frames: Vec<u32> = frame_z.keys().cloned().collect();
    frames.sort();

    // Parameters for the catheter circle.
    let center_x = 4.5;
    let center_y = 4.5;
    let radius = 0.5;
    let num_points = 20;

    // For each unique frame, generate 10 catheter points around a circle.
    for frame in frames {
        let z = frame_z[&frame];
        for i in 0..num_points {
            let angle = 2.0 * PI * (i as f64) / (num_points as f64);
            let x = center_x + radius * angle.cos();
            let y = center_y + radius * angle.sin();
            catheter_points.push(ContourPoint {
                frame_index: frame,
                x,
                y,
                z,
            });
        }
    }
    catheter_points
}

pub fn write_obj_mesh(
    contours: &[(u32, Vec<ContourPoint>)],
    uv_coords: &[(f32, f32)],
    filename: &str,
    mtl_filename: &str,
) -> Result<(), Box<dyn Error>> {
    let sorted_contours = contours.to_owned();

    // Validation
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

    // Validate UV coordinates
    if uv_coords.len() != current_offset - 1 {
        return Err(format!(
            "UV coordinates must match the number of vertices. Expected {}, got {}.",
            current_offset - 1,
            uv_coords.len()
        ).into());
    }

    // Write material reference
    writeln!(writer, "mtllib {}", mtl_filename)?;
    writeln!(writer, "usemtl displacement_material")?;
    

    // Write UV coordinates
    for (u, v) in uv_coords {
        writeln!(writer, "vt {} {}", u, v)?;
    }

    // Write normals
    for (nx, ny, nz) in &normals {
        writeln!(writer, "vn {} {} {}", nx, ny, nz)?;
    }

    // Write faces with normals and UVs
    for c in 0..(sorted_contours.len() - 1) {
        let offset1 = vertex_offsets[c];
        let offset2 = vertex_offsets[c + 1];
        for j in 0..points_per_contour {
            let j_next = (j + 1) % points_per_contour;
            let v1 = offset1 + j;
            let v2 = offset1 + j_next;
            let v3 = offset2 + j;
            writeln!(
                writer,
                "f {0}/{0}/{0} {1}/{1}/{1} {2}/{2}/{2}",
                v1, v2, v3
            )?;
            let v1_t2 = offset2 + j;
            let v2_t2 = offset1 + j_next;
            let v3_t2 = offset2 + j_next;
            writeln!(
                writer,
                "f {0}/{0}/{0} {1}/{1}/{1} {2}/{2}/{2}",
                v1_t2, v2_t2, v3_t2
            )?;
        }
    }

    println!("OBJ mesh with normals written to {}", filename);
    Ok(())
}