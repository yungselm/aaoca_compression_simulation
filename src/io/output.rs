use std::error::Error;
use std::fs::File;
use std::io::{BufWriter, Write};
use crate::io::input::Contour;


pub fn write_obj_mesh(
    contours: &Vec<Contour>,
    uv_coords: &[(f64, f64)],
    filename: &str,
    mtl_filename: &str,
) -> Result<(), Box<dyn Error>> {
    let sorted_contours = contours.to_owned();

    // Validation
    if sorted_contours.len() < 2 {
        return Err("Need at least two contours to create a mesh.".into());
    }

    let points_per_contour = sorted_contours[0].points.len();
    for contour in &sorted_contours {
        if contour.points.len() != points_per_contour {
            return Err("All contours must have the same number of points.".into());
        }
    }

    let file = File::create(filename)?;
    let mut writer = BufWriter::new(file);
    let mut vertex_offsets = Vec::new();
    let mut current_offset = 1;
    let mut normals = Vec::new();

    // Write vertices and compute normals
    for contour in &sorted_contours {
        vertex_offsets.push(current_offset);
        let centroid = &contour.centroid;
        for point in &contour.points {
            writeln!(writer, "v {} {} {}", point.x, point.y, point.z)?;
            let dx = point.x - centroid.0;
            let dy = point.y - centroid.1;
            let length = (dx * dx + dy * dy).sqrt();
            let (nx, ny, nz) = if length > 0.0 {
                (dx / length, dy / length, 0.0)
            } else {
                (0.0, 0.0, 0.0)
            };
            normals.push((nx * -1.0, ny * -1.0, nz * -1.0));
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


