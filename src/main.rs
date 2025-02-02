mod data_read;
use data_read::read_contour_data;
use data_read::ContourPoint;
use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::{BufWriter, Write};

/// Sorts the points of a contour in counterclockwise order starting at the point with the highest y value. Easier to create a mesh like this.
fn sort_contour_points(points: &mut [ContourPoint]) {
    // Compute the center of the contour (the average of x and y).
    let (sum_x, sum_y) = points.iter().fold((0.0, 0.0), |(sx, sy), p| (sx + p.x, sy + p.y));
    let count = points.len() as f64;
    let center = (sum_x / count, sum_y / count);

    // Sort the points in counterclockwise order using the angle relative to the center.
    points.sort_by(|a, b| {
        let angle_a = (a.y - center.1).atan2(a.x - center.0);
        let angle_b = (b.y - center.1).atan2(b.x - center.0);
        angle_a.partial_cmp(&angle_b).unwrap()
    });

    // Rotate the vector so that the point with the highest y is at the beginning.
    let start_idx = points
        .iter()
        .enumerate()
        .max_by(|(_, a), (_, b)| a.y.partial_cmp(&b.y).unwrap())
        .map(|(i, _)| i)
        .unwrap_or(0);
    points.rotate_left(start_idx);
}

/// Writes an OBJ file connecting each adjacent pair of contours by creating two triangles per quad.
fn write_obj_mesh(contours: Vec<(u32, Vec<ContourPoint>)>) -> Result<(), Box<dyn Error>> {
    let mut sorted_contours = contours;
    sorted_contours.sort_by_key(|(frame_index, _)| *frame_index); //why need dereference here?

    // Check that we have at least two contours.
    if sorted_contours.len() < 2 {
        return Err("Need at least two contours to create a mesh.".into());
    }

    // Assume all contours have the same number of points.
    let points_per_contour = sorted_contours[0].1.len();
    for (_, contour) in &sorted_contours {
        if contour.len() != points_per_contour {
            return Err("All contours must have the same number of points.".into());
        }
    }

    // Create (or overwrite) the OBJ file.
    let file = File::create("output.obj")?;
    let mut writer = BufWriter::new(file);

    // We'll need to track the starting index of each contour’s vertices.
    // OBJ indices start at 1.
    let mut vertex_offsets = Vec::new();
    let mut current_offset = 1;

    for (_, contour) in &sorted_contours {
        vertex_offsets.push(current_offset);
        for point in contour {
            // Write each vertex line: "v x y z"
            writeln!(writer, "v {} {} {}", point.x, point.y, point.z)?;
            current_offset += 1;
        }
    }

    // For each adjacent pair of contours, create faces (two triangles per quad).
    for c in 0..(sorted_contours.len() - 1) {
        let offset1 = vertex_offsets[c];
        let offset2 = vertex_offsets[c + 1];
        for j in 0..points_per_contour {
            let j_next = (j + 1) % points_per_contour;

            // First triangle: contour i point j, contour i point j+1, contour i+1 point j
            let v1 = offset1 + j;
            let v2 = offset1 + j_next;
            let v3 = offset2 + j;
            writeln!(writer, "f {} {} {}", v1, v2, v3)?;

            // Second triangle: contour i+1 point j, contour i point j+1, contour i+1 point j+1
            let v1 = offset2 + j;
            let v2 = offset1 + j_next;
            let v3 = offset2 + j_next;
            writeln!(writer, "f {} {} {}", v1, v2, v3)?;
        }
    }

    println!("OBJ mesh successfully written to output.obj"); // courld print with error handling
    Ok(())
}

fn main() -> Result<(), Box<dyn Error>> {
    let contour_points = read_contour_data("input/rest_csv_files/diastolic_contours.csv")?;

    let mut groups: HashMap<u32, Vec<ContourPoint>> = HashMap::new();
    for point in contour_points {
        groups.entry(point.frame_index).or_default().push(point);
    }

    let mut contours = Vec::new();
    for (frame_index, mut contour) in groups {
        sort_contour_points(&mut contour);
        contours.push((frame_index, contour));
    }

    // (Optional) Print sorted contours for debugging.
    for (frame_index, contour) in &contours {
        println!("Contour {}:", frame_index);
        for (i, point) in contour.iter().enumerate() {
            println!("  Point {}: {:?}", i, point);
        }
        println!("");
    }

    write_obj_mesh(contours)?;

    Ok(())
}
