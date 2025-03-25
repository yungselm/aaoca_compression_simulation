use std::error::Error;
use nalgebra::{Point3, Rotation3, Unit, Vector3};
use crate::io::{Centerline, CenterlinePoint};
use crate::io::{read_centerline_txt, read_obj_mesh, write_updated_obj_mesh, ContourPoint};
use crate::contour::Contour;
 
// const FIXED_ROTATION_DEG: f64 = 235.0;
const FIXED_ROTATION_DEG: f64 = 220.0; // needs to be replaced with function that finds automatically.

/// Applies a fixed rotation around the z-axis to preserve the initial in-plane orientation.
fn rotate_contours_around_z(mesh: &mut Vec<(u32, Vec<ContourPoint>)>, degrees: f64) {
    let radians = degrees.to_radians();
    let cos_theta = radians.cos();
    let sin_theta = radians.sin();

    for (_, contours) in mesh.iter_mut() {
        for point in contours.iter_mut() {
            let x_new = point.x * cos_theta - point.y * sin_theta;
            let y_new = point.x * sin_theta + point.y * cos_theta;
            point.x = x_new;
            point.y = y_new;
            // z remains unchanged as we're rotating around the z-axis
        }
    }
}

/// Rotates a single contour around the z-axis by `degrees`.
fn rotate_single_contour_around_z(contour: &[ContourPoint], degrees: f64) -> Vec<ContourPoint> {
    let radians = degrees.to_radians();
    let cos_theta = radians.cos();
    let sin_theta = radians.sin();

    contour
        .iter()
        .map(|point| ContourPoint {
            frame_index: point.frame_index,
            x: point.x * cos_theta - point.y * sin_theta,
            y: point.x * sin_theta + point.y * cos_theta,
            z: point.z,
        })
        .collect()
}

/// Finds the optimal rotation angle by minimizing the distance between the closest opposite point 
/// and the reference coordinate.
fn find_optimal_rotation(
    contour: &[ContourPoint],
    target_x: f64,
    target_y: f64,
    angle_step: f64,
) -> f64 {
    let mut best_angle = 0.0;
    let mut min_distance = f64::MAX;

    let mut angle_deg = 0.0;
    while angle_deg < 360.0 {
        let rotated = rotate_single_contour_around_z(contour, angle_deg);
        let (closest_pair, _) = Contour::find_closest_opposite(&rotated);
        let p = &closest_pair.0;
        let distance = ((p.x - target_x).powi(2) + (p.y - target_y).powi(2)).sqrt();

        if distance < min_distance {
            min_distance = distance;
            best_angle = angle_deg;
        }
        angle_deg += angle_step;
    }

    best_angle
}

#[derive(Debug, Clone, PartialEq)]
pub struct ContourFrame {
    pub frame_index: u32,
    pub points: Vec<ContourPoint>,  // 500 points per frame
    pub centroid_3d: ContourPoint,
    pub normal: Vector3<f64>,       // Frame's normal vector
    pub translation: Point3<f64>,   // Translation from centerline
}

impl ContourFrame {
    /// Creates a new ContourFrame from raw contour points, calculating the centroid and normal.
    fn from_contour(frame_index: u32, points: Vec<ContourPoint>) -> Self {
        assert!(!points.is_empty(), "Cannot create ContourFrame from empty points");

        let centroid_3d = Self::calculate_centroid_3d(&points);
        let normal = Self::calculate_normal(&points, &centroid_3d);
        
        ContourFrame {
            frame_index,
            points,
            centroid_3d,
            normal,
            translation: Point3::origin(), // Will be updated during alignment
        }
    }

    /// Calculates the geometric centroid of the contour points.
    fn calculate_centroid_3d(points: &[ContourPoint]) -> ContourPoint {
        let count = points.len() as f64;
        let (sum_x, sum_y, sum_z) = points.iter().fold((0.0, 0.0, 0.0), |(sx, sy, sz), p| {
            (sx + p.x, sy + p.y, sz + p.z)
        });

        ContourPoint {
            frame_index: points.first().map(|p| p.frame_index).unwrap_or(0),
            x: sum_x / count,
            y: sum_y / count,
            z: sum_z / count,
        }
    }

    /// Calculates the normal vector using a stable cross product method.
    fn calculate_normal(points: &[ContourPoint], centroid: &ContourPoint) -> Vector3<f64> {
        // Select three non-collinear points for robust normal calculation.
        let p1 = &points[0];
        let p2 = &points[points.len() / 3];
        let p3 = &points[2 * points.len() / 3];
        
        let v1 = Vector3::new(p1.x - centroid.x, p1.y - centroid.y, p1.z - centroid.z);
        let v2 = Vector3::new(p2.x - centroid.x, p2.y - centroid.y, p2.z - centroid.z);
        let v3 = Vector3::new(p3.x - centroid.x, p3.y - centroid.y, p3.z - centroid.z);

        (v1.cross(&v2) + v2.cross(&v3)).normalize()
    }
}

pub struct FrameTransformation {
    pub frame_index: u32,
    pub translation: Vector3<f64>,
    pub rotation: Rotation3<f64>,
    pub pivot: Point3<f64>,
}

/// Modified align_frame returns the transformation applied.
fn align_frame(frame: &mut ContourFrame, cl_point: &CenterlinePoint) -> FrameTransformation {
    // === Translation Step ===
    // Compute the translation vector to bring the frame's centroid to the centerline point.
    let translation_vec = Vector3::new(
        cl_point.contour_point.x - frame.centroid_3d.x,
        cl_point.contour_point.y - frame.centroid_3d.y,
        cl_point.contour_point.z - frame.centroid_3d.z,
    );
    for point in frame.points.iter_mut() {
        point.x += translation_vec.x;
        point.y += translation_vec.y;
        point.z += translation_vec.z;
    }
    // Update the frame's centroid and translation.
    frame.centroid_3d.x += translation_vec.x;
    frame.centroid_3d.y += translation_vec.y;
    frame.centroid_3d.z += translation_vec.z;
    frame.translation = Point3::new(
        cl_point.contour_point.x,
        cl_point.contour_point.y,
        cl_point.contour_point.z,
    );

    // === Rotation Step ===
    // Compute the rotation needed to align the frame's normal with the centerline normal.
    let current_normal = frame.normal;
    let desired_normal = cl_point.normal;
    let angle = current_normal.angle(&desired_normal);
    let rotation: Rotation3<f64> = if angle.abs() < 1e-6 {
        Rotation3::identity()
    } else {
        let rotation_axis = current_normal.cross(&desired_normal);
        if rotation_axis.norm() < 1e-6 {
            Rotation3::identity()
        } else {
            let rotation_axis_unit = Unit::new_normalize(rotation_axis);
            Rotation3::from_axis_angle(&rotation_axis_unit, angle)
        }
    };

    // Define the pivot as the centerline point.
    let pivot = Point3::new(
        cl_point.contour_point.x,
        cl_point.contour_point.y,
        cl_point.contour_point.z,
    );
    // Apply rotation to every point if needed.
    if angle.abs() >= 1e-6 {
        for point in frame.points.iter_mut() {
            let current_point = Point3::new(point.x, point.y, point.z);
            let relative_vector = current_point - pivot;
            let rotated_relative = rotation * relative_vector;
            let rotated_point = pivot + rotated_relative;
            point.x = rotated_point.x;
            point.y = rotated_point.y;
            point.z = rotated_point.z;
        }
        // Also update the frame's normal.
        frame.normal = rotation * frame.normal;
    }

    // Return the transformation details for later use.
    FrameTransformation {
        frame_index: frame.frame_index,
        translation: translation_vec,
        rotation,
        pivot,
    }
}

/// Updated function that processes mesh frames and collects transformations.
pub fn process_mesh_frames(
    mesh: Vec<(u32, Vec<ContourPoint>)>,
    centerline: &Centerline,
) -> (Vec<ContourFrame>, Vec<FrameTransformation>) {
    let mut aligned_frames = Vec::new();
    let mut transformations = Vec::new();

    for (frame_idx, points) in mesh.into_iter() {
        let mut frame = ContourFrame::from_contour(frame_idx, points);
        if let Some(cl_point) = centerline.get_by_frame(frame_idx) {
            let transformation = align_frame(&mut frame, cl_point);
            transformations.push(transformation);
        }
        aligned_frames.push(frame);
    }
    (aligned_frames, transformations)
}

/// Generates new face connectivity and UV coordinates from aligned frames.
/// Assumes every frame has the same number of points.
fn generate_faces_and_uv(
    frames: &[ContourFrame]
) -> (Vec<(usize, usize, usize)>, Vec<(f32, f32)>) {
    let num_frames = frames.len();
    let points_per_frame = frames[0].points.len();
    let mut faces = Vec::new();
    let mut uv_coords = Vec::with_capacity(num_frames * points_per_frame);
    
    // Generate new UV coordinates.
    for i in 0..num_frames {
        let v = (i as f32 + 0.5) / num_frames as f32;
        for j in 0..points_per_frame {
            let u = (j as f32 + 0.5) / points_per_frame as f32;
            uv_coords.push((u, v));
        }
    }
    
    // Generate new faces: For each adjacent pair of frames,
    // create two triangles per quad.
    for i in 0..(num_frames - 1) {
        let base1 = i * points_per_frame;
        let base2 = (i + 1) * points_per_frame;
        for j in 0..points_per_frame {
            let j_next = (j + 1) % points_per_frame;
            // First triangle.
            faces.push((base1 + j, base1 + j_next, base2 + j));
            // Second triangle.
            faces.push((base2 + j, base1 + j_next, base2 + j_next));
        }
    }
    (faces, uv_coords)
}

pub fn create_centerline_aligned_meshes(
    state: &str,
    centerline_path: &str,
    input_dir: &str,
    output_dir: &str,
    interpolation_steps: usize,
    x_coord_ref: f64,
    y_coord_ref: f64,
) -> Result<(), Box<dyn Error>> {
    // ----- Build the common centerline -----
    let raw_centerline = read_centerline_txt(centerline_path)?;
    let centerline = Centerline::from_contour_points(raw_centerline);

    // ----- Process the reference mesh: mesh_000_rest.obj -----
    let ref_mesh_path = format!("{}/mesh_000_{}.obj", input_dir, state);
    let mut reference_mesh = read_obj_mesh(&ref_mesh_path)?;

    // Compute optimal rotation using the first frame's contour
    let first_frame_contour = &reference_mesh[0].1;
    let optimal_angle = find_optimal_rotation(
        first_frame_contour,
        x_coord_ref,
        y_coord_ref,
        1.0, // Step size in degrees
    );
    println!("Optimal rotation angle: {:.2} degrees", optimal_angle);
    rotate_contours_around_z(&mut reference_mesh, optimal_angle);
    
    // rotate_contours_around_z(&mut reference_mesh, FIXED_ROTATION_DEG);
    
    // Compute aligned frames and capture transformation parameters.
    let (aligned_frames, transformations) = process_mesh_frames(reference_mesh, &centerline);

    // Generate new connectivity and UVs for the reference mesh.
    let (new_faces, new_uvs) = generate_faces_and_uv(&aligned_frames);
    let new_vertices: Vec<(f32, f32, f32)> = aligned_frames
        .iter()
        .flat_map(|frame| frame.points.iter().map(|pt| (pt.x as f32, pt.y as f32, pt.z as f32)))
        .collect();
    let new_normals: Vec<(f32, f32, f32)> = aligned_frames
        .iter()
        .flat_map(|frame| {
            let n = frame.points.len();
            std::iter::repeat((
                frame.normal.x as f32,
                frame.normal.y as f32,
                frame.normal.z as f32,
            ))
            .take(n)
        })
        .collect();

    // Ensure output directory exists.
    std::fs::create_dir_all(output_dir)?;
    let ref_mtl_src = format!("{}/mesh_000_{}.mtl", input_dir, state);
    let ref_png_src = format!("{}/mesh_000_{}.png", input_dir, state);
    let ref_mtl_dest = format!("{}/mesh_000_{}.mtl", output_dir, state);
    let ref_png_dest = format!("{}/mesh_000_{}.png", output_dir, state);
    std::fs::copy(ref_mtl_src, ref_mtl_dest)?;
    std::fs::copy(ref_png_src, ref_png_dest)?;

    let ref_obj_dest = format!("{}/mesh_000_{}.obj", output_dir, state);
    let ref_mtl_dest = format!("mesh_000_{}.mtl", state);
    write_updated_obj_mesh(
        &ref_obj_dest,
        &ref_mtl_dest,
        &new_vertices,
        &new_uvs,
        &new_normals,
        &new_faces,
    )?;
    println!(
        "Reference mesh (mesh_000_{}.obj) processed and transformation parameters stored.", state
    );

    // ----- Process subsequent meshes using stored transformation parameters -----
    // Loop over prefixes "mesh" and "catheter" for indices 1 to 31.
    for prefix in ["mesh", "catheter"].iter() {
        let start_index = if *prefix == "mesh" { 1 } else { 0 };
        for i in start_index..=(interpolation_steps - 1) {
            let obj_filename = format!("{}/{}_{:03}_{}.obj", input_dir, prefix, i, state);
            let mut mesh = read_obj_mesh(&obj_filename)?;
            rotate_contours_around_z(&mut mesh, FIXED_ROTATION_DEG);

            // Apply transformation for each frame in the mesh.
            for (frame_index, contours) in mesh.iter_mut() {
                if let Some(transformation) =
                    transformations.iter().find(|t| t.frame_index == *frame_index)
                {
                    for point in contours.iter_mut() {
                        // Translation.
                        point.x += transformation.translation.x;
                        point.y += transformation.translation.y;
                        point.z += transformation.translation.z;
                        // Rotation about stored pivot.
                        let current_point = Point3::new(point.x, point.y, point.z);
                        let relative_vector = current_point - transformation.pivot;
                        let rotated_relative = transformation.rotation * relative_vector;
                        let rotated_point = transformation.pivot + rotated_relative;
                        point.x = rotated_point.x;
                        point.y = rotated_point.y;
                        point.z = rotated_point.z;
                    }
                } else {
                    println!("Warning: No transformation found for frame index {}", frame_index);
                }
            }

            // Rebuild vertex list.
            let new_vertices: Vec<(f32, f32, f32)> = mesh
                .iter()
                .flat_map(|(_frame_index, contours)| {
                    contours.iter().map(|pt| (pt.x as f32, pt.y as f32, pt.z as f32))
                })
                .collect();
            let num_frames = mesh.len();
            let points_per_frame = if num_frames > 0 { mesh[0].1.len() } else { 0 };

            // Generate new UV coordinates.
            let mut new_uvs = Vec::with_capacity(num_frames * points_per_frame);
            for i in 0..num_frames {
                let v = (i as f32 + 0.5) / num_frames as f32;
                for j in 0..points_per_frame {
                    let u = (j as f32 + 0.5) / points_per_frame as f32;
                    new_uvs.push((u, v));
                }
            }

            // Generate new face connectivity.
            let mut new_faces = Vec::new();
            for i in 0..(num_frames - 1) {
                let base1 = i * points_per_frame;
                let base2 = (i + 1) * points_per_frame;
                for j in 0..points_per_frame {
                    let j_next = (j + 1) % points_per_frame;
                    new_faces.push((base1 + j, base1 + j_next, base2 + j));
                    new_faces.push((base2 + j, base1 + j_next, base2 + j_next));
                }
            }

            // For normals we assume a default (or compute per-frame normals if available).
            // Here we simply assign (0.0, 0.0, 1.0) to every vertex.
            let new_normals: Vec<(f32, f32, f32)> = (0..(num_frames * points_per_frame))
                .map(|_| (0.0, 0.0, 1.0))
                .collect();

            // Copy corresponding .mtl and texture files.
            let mtl_src = format!("{}/{}_{:03}_{}.mtl", input_dir, prefix, i, state);
            let texture_src = format!("{}/{}_{:03}_{}.png", input_dir, prefix, i, state);
            let mtl_dest = format!("{}/{}_{:03}_{}.mtl", output_dir, prefix, i, state);
            let texture_dest = format!("{}/{}_{:03}_{}.png", output_dir, prefix, i, state);
            std::fs::copy(mtl_src, mtl_dest)?;
            std::fs::copy(texture_src, texture_dest)?;

            // Write the updated OBJ file.
            let obj_dest = format!("{}/{}_{:03}_{}.obj", output_dir, prefix, i, state);
            let mtl_basename = format!("{}_{:03}_{}.mtl", prefix, i, state);
            write_updated_obj_mesh(
                &obj_dest,
                &mtl_basename,
                &new_vertices,
                &new_uvs,
                &new_normals,
                &new_faces,
            )?;
            println!(
                "OBJ {}_{:03}_{}.obj processed using stored transformation parameters.",
                prefix, i, state
            );
        }
    }
    Ok(())
}