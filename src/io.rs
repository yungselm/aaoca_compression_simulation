use nalgebra::Vector3;
use serde::Deserialize;
use std::collections::HashMap;
use std::error::Error;
use std::f64::consts::PI;
use std::fs::File;
use std::io::BufRead;
use std::io::{BufReader, BufWriter, Write};
use std::path::Path;

#[derive(Debug, Deserialize, Clone, PartialEq)]
pub struct ContourPoint {
    pub frame_index: u32,
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl ContourPoint {
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
}

#[derive(Debug, Clone, PartialEq)]
pub struct Centerline {
    pub points: Vec<CenterlinePoint>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct CenterlinePoint {
    pub contour_point: ContourPoint,
    pub normal: Vector3<f64>,
}

impl Centerline {
    pub fn from_contour_points(contour_points: Vec<ContourPoint>) -> Self {
        let mut points: Vec<CenterlinePoint> = Vec::with_capacity(contour_points.len());

        // Calculate normals for all but the last point.
        for i in 0..contour_points.len() {
            let current = &contour_points[i];
            let normal = if i < contour_points.len() - 1 {
                let next = &contour_points[i + 1];
                Vector3::new(next.x - current.x, next.y - current.y, next.z - current.z).normalize()
            } else if !contour_points.is_empty() {
                points[i - 1].normal
            } else {
                Vector3::zeros()
            };

            points.push(CenterlinePoint {
                contour_point: current.clone(),
                normal,
            });
        }

        Centerline { points }
    }

    /// Retrieves a centerline point by matching frame index.
    pub fn get_by_frame(&self, frame_index: u32) -> Option<&CenterlinePoint> {
        self.points
            .iter()
            .find(|p| p.contour_point.frame_index == frame_index)
    }
}

pub fn read_centerline_txt(path: &str) -> Result<Vec<ContourPoint>, Box<dyn Error>> {
    let file = File::open(path)?;
    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(false)
        .delimiter(b' ')
        .from_reader(file);

    let mut points = Vec::new();
    for (i, result) in rdr.records().enumerate() {
        match result {
            Ok(record) => {
                let mut iter = record.iter();
                let x = iter.next().unwrap().parse::<f64>()?;
                let y = iter.next().unwrap().parse::<f64>()?;
                let z = iter.next().unwrap().parse::<f64>()?;
                points.push(ContourPoint {
                    frame_index: i as u32,
                    x,
                    y,
                    z,
                });
            }
            Err(e) => eprintln!("Skipping invalid row: {:?}", e),
        }
    }
    Ok(points)
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
        let centroid = crate::contour::Contour::compute_centroid(contour);
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
        )
        .into());
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
            writeln!(writer, "f {0}/{0}/{0} {1}/{1}/{1} {2}/{2}/{2}", v1, v2, v3)?;
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

/// Writes a new OBJ file with updated vertices, normals, UVs, and faces.
/// This function is used in centerline_alignment.rs to write mesh fitted to CCTA centerline.
pub fn write_updated_obj_mesh(
    filename: &str,
    mtl_filename: &str,
    vertices: &[(f32, f32, f32)],
    uv_coords: &[(f32, f32)],
    normals: &[(f32, f32, f32)],
    faces: &[(usize, usize, usize)],
) -> Result<(), Box<dyn Error>> {
    let file = File::create(filename)?;
    let mut writer = BufWriter::new(file);

    // Write material reference.
    writeln!(writer, "mtllib {}", mtl_filename)?;
    writeln!(writer, "usemtl displacement_material")?;

    // Write vertices.
    for &(x, y, z) in vertices {
        writeln!(writer, "v {} {} {}", x, y, z)?;
    }

    // Write UV coordinates.
    for &(u, v) in uv_coords {
        writeln!(writer, "vt {} {}", u, v)?;
    }

    // Write normals.
    for &(nx, ny, nz) in normals {
        writeln!(writer, "vn {} {} {}", nx, ny, nz)?;
    }

    // Write faces (OBJ indices start at 1).
    for &(v1, v2, v3) in faces {
        writeln!(
            writer,
            "f {}/{}/{} {}/{}/{} {}/{}/{}",
            v1 + 1,
            v1 + 1,
            v1 + 1,
            v2 + 1,
            v2 + 1,
            v2 + 1,
            v3 + 1,
            v3 + 1,
            v3 + 1
        )?;
    }
    println!("Updated OBJ mesh written to {}", filename);
    Ok(())
}

pub fn read_obj_mesh(filename: &str) -> Result<Vec<(u32, Vec<ContourPoint>)>, Box<dyn Error>> {
    let file = File::open(filename)?;
    let reader = BufReader::new(file);

    let mut vertices: Vec<ContourPoint> = Vec::new();
    // This will be set based on the first face encountered.
    let mut points_per_contour: Option<usize> = None;

    // Loop over every line in the file.
    for line_result in reader.lines() {
        let line = line_result?;
        let trimmed = line.trim();
        if trimmed.starts_with("v ") {
            // This is a vertex line. Split on whitespace.
            let parts: Vec<&str> = trimmed.split_whitespace().collect();
            if parts.len() >= 4 {
                let x = parts[1].parse::<f64>()?;
                let y = parts[2].parse::<f64>()?;
                let z = parts[3].parse::<f64>()?;
                // We temporarily set frame_index to 0.
                vertices.push(ContourPoint {
                    frame_index: 0,
                    x,
                    y,
                    z,
                });
            }
        } else if trimmed.starts_with("f ") && points_per_contour.is_none() {
            // Use the first face line to deduce the number of vertices per contour.
            // A face line written by write_obj_mesh has the form:
            // f {v1}/{v1}/{v1} {v2}/{v2}/{v2} {v3}/{v3}/{v3}
            // where v1 comes from the first contour and v3 from the next.
            let parts: Vec<&str> = trimmed.split_whitespace().collect();
            if parts.len() >= 4 {
                // Extract the first vertex index from parts[1] and parts[3].
                let v1_str = parts[1].split('/').next().unwrap();
                let v3_str = parts[3].split('/').next().unwrap();
                let v1: usize = v1_str.parse()?;
                let v3: usize = v3_str.parse()?;
                if v3 > v1 {
                    points_per_contour = Some(v3 - v1);
                }
            }
        }
        // Other lines (vt, vn, mtllib, usemtl, etc.) are ignored.
    }

    // If no face line was found, assume all vertices form one single contour.
    let points_per_contour = points_per_contour.unwrap_or(vertices.len());

    // Ensure the total number of vertices divides evenly into contours.
    if vertices.len() % points_per_contour != 0 {
        return Err(format!(
            "Vertex count {} is not divisible by points per contour {}.",
            vertices.len(),
            points_per_contour
        )
        .into());
    }

    let num_contours = vertices.len() / points_per_contour;
    let mut contours: Vec<(u32, Vec<ContourPoint>)> = Vec::with_capacity(num_contours);

    for i in 0..num_contours {
        let start = i * points_per_contour;
        let end = start + points_per_contour;
        let mut contour = vertices[start..end].to_vec();
        // Set the contour's frame_index to the contour index.
        for point in &mut contour {
            point.frame_index = i as u32;
        }
        contours.push((i as u32, contour));
    }

    Ok(contours)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::EPSILON;
    use std::fs;

    // If you don't already have a compute_centroid function in your crate,
    // uncomment the following dummy version for testing purposes.
    #[allow(dead_code)]
    fn compute_centroid(contour: &[ContourPoint]) -> (f64, f64, f64) {
        let (sum_x, sum_y, sum_z) = contour.iter().fold((0.0, 0.0, 0.0), |(sx, sy, sz), p| {
            (sx + p.x, sy + p.y, sz + p.z)
        });
        let n = contour.len() as f64;
        (sum_x / n, sum_y / n, sum_z / n)
    }

    #[test]
    fn test_read_write_obj_mesh() {
        // Create sample contours: two contours, each with 2 vertices.
        let input_contours = vec![
            (
                0,
                vec![
                    ContourPoint {
                        frame_index: 0,
                        x: 0.0,
                        y: 0.0,
                        z: 0.0,
                    },
                    ContourPoint {
                        frame_index: 0,
                        x: 1.0,
                        y: 0.0,
                        z: 0.0,
                    },
                ],
            ),
            (
                1,
                vec![
                    ContourPoint {
                        frame_index: 1,
                        x: 0.0,
                        y: 0.0,
                        z: 1.0,
                    },
                    ContourPoint {
                        frame_index: 1,
                        x: 1.0,
                        y: 0.0,
                        z: 1.0,
                    },
                ],
            ),
        ];

        // Total vertices: 2 contours × 2 points = 4.
        // Provide a UV coordinate for each vertex.
        let total_vertices: usize = input_contours.iter().map(|(_, pts)| pts.len()).sum();
        let uv_coords = vec![(0.0_f32, 0.0_f32); total_vertices];

        // Define temporary filenames.
        let obj_filename = "test_roundtrip.obj";
        let mtl_filename = "test_roundtrip.mtl";

        // Write the OBJ mesh.
        write_obj_mesh(&input_contours, &uv_coords, obj_filename, mtl_filename)
            .expect("Failed to write OBJ mesh");

        // Read the mesh back.
        let output_contours = read_obj_mesh(obj_filename).expect("Failed to read OBJ mesh");

        println!("Input {:?}", input_contours);
        println!("Output {:?}", output_contours);

        // Check that the number of contours is the same.
        assert_eq!(input_contours.len(), output_contours.len());

        // Verify that each contour has the expected number of points and matching vertex values.
        for (i, (frame_idx, contour)) in output_contours.iter().enumerate() {
            assert_eq!(contour.len(), input_contours[i].1.len());
            for (p_written, p_read) in input_contours[i].1.iter().zip(contour.iter()) {
                assert!((p_written.x - p_read.x).abs() < EPSILON);
                assert!((p_written.y - p_read.y).abs() < EPSILON);
                assert!((p_written.z - p_read.z).abs() < EPSILON);
            }
            // Confirm that the new frame indices are sequential starting at 0.
            assert_eq!(*frame_idx, i as u32);
        }

        // Clean up temporary files.
        fs::remove_file(obj_filename).unwrap();
    }
}
