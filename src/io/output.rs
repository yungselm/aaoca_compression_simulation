use anyhow::{bail, anyhow};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use crate::io::input::Contour;
use crate::io::Geometry;

pub fn write_obj_mesh(
    contours: &Vec<Contour>,
    uv_coords: &[(f64, f64)],
    filename: &str,
    mtl_filename: &str,
) -> anyhow::Result<()> {
    let sorted_contours = contours.to_owned();

    // Validation
    if sorted_contours.len() < 2 {
        bail!("Need at least two contours to create a mesh.");
    }

    let points_per_contour = sorted_contours[0].points.len();
    for contour in &sorted_contours {
        if contour.points.len() != points_per_contour {
            bail!("All contours must have the same number of points.");
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
        return Err(anyhow!(
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

#[derive(Copy, Clone)]
pub enum GeometryType {
    Contour,
    Catheter,
    Wall,
}

impl GeometryType {
    // Get the contour data from a Geometry based on the enum variant
    pub fn get_contours<'a>(&self, geometry: &'a Geometry) -> &'a Vec<Contour> {
        match self {
            GeometryType::Contour => &geometry.contours,
            GeometryType::Catheter => &geometry.catheter,
            GeometryType::Wall => &geometry.contours,
        }
    }

    // Get the object string for filenames
    pub fn as_str(&self) -> &'static str {
        match self {
            GeometryType::Contour => "mesh",
            GeometryType::Catheter => "catheter",
            GeometryType::Wall => "wall",
        }
    }
}

pub fn write_geometry_vec_to_obj(
    geometry_type: GeometryType,
    case_name: &str,
    output_dir: &str,
    geometries: &[Geometry],
    uv_coords: &[Vec<(f64, f64)>],
) -> anyhow::Result<()> {
    for (i, (geometry, uv_coords_mesh)) in geometries
        .iter()
        .zip(uv_coords.iter())
        .enumerate()
    {
        let obj_filename = format!(
            "{}_{:03}_{}.obj",
            geometry_type.as_str(),
            i,
            case_name
        );
        let mtl_filename = format!(
            "{}_{:03}_{}.mtl",
            geometry_type.as_str(),
            i,
            case_name
        );
        let obj_path = Path::new(output_dir).join(&obj_filename);
        let obj_path_str = obj_path.to_str()
            .ok_or_else(|| anyhow!("Invalid path for OBJ file"))?;

        // Get the appropriate contour data (contours or catheter)
        let contours = geometry_type.get_contours(geometry);

        write_obj_mesh(
            contours,
            uv_coords_mesh,
            obj_path_str,
            &mtl_filename,
        ).map_err(|e| anyhow!("Failed to write OBJ mesh: {}", e))?;
    }
    Ok(())
}