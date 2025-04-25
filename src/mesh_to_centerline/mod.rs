pub mod preprocessing;
pub mod operations;

use std::error::Error;

use preprocessing::{prepare_data_3d_alignment, read_interpolated_meshes};
use operations::{find_optimal_rotation, get_transformations};
use crate::io::input::Contour;
use crate::texture::write_mtl_geometry;

use std::path::Path;
use crate::io::output::write_obj_mesh;
use crate::utils::utils::{write_centerline_to_csv, write_geometry_to_csv, write_frame_transformation_to_csv};

pub fn create_centerline_aligned_meshes(
    state: &str,
    centerline_path: &str,
    input_dir: &str,
    output_dir: &str,
    interpolation_steps: usize,
    x_coord_ref: f64,
    y_coord_ref: f64,
    z_coord_ref: f64,
    x_coord_upper: f64,
    y_coord_upper: f64,
    z_coord_upper: f64,
    x_coord_lower: f64,
    y_coord_lower: f64,
    z_coord_lower: f64,
) -> Result<(), Box<dyn Error>> {
    println!("the input_dir is: {:?}", &input_dir);
    println!("the output_dir is: {:?}", &output_dir);
    let (centerline, ref_mesh_dia, ref_mesh_sys) = prepare_data_3d_alignment(state, centerline_path, input_dir, interpolation_steps)?;
    let interpolated_meshes = read_interpolated_meshes(state, input_dir, interpolation_steps);

    let mut geometries = vec![ref_mesh_dia.clone()];
    geometries.extend(interpolated_meshes);
    geometries.extend(vec![ref_mesh_sys.clone()]);

    let n = ref_mesh_dia.contours.len();
    let n_cl = centerline.points.len();

    println!("Contour 0 has this z: {:?} and index: {:?}", &ref_mesh_dia.contours[0].centroid.2, ref_mesh_dia.contours[0].id);
    println!("Contour n -1 has this z: {:?} and index: {:?}", &ref_mesh_dia.contours[n - 1].centroid.2, ref_mesh_dia.contours[n - 1].id);
    println!("Centerline 0 has this z: {:?} and index: {:?}", &centerline.points[0].contour_point.z, &centerline.points[0].contour_point.frame_index);
    println!("Centerline n - 1 has this z: {:?} and index: {:?}", &centerline.points[n_cl - 1].contour_point.z, &centerline.points[n_cl - 1].contour_point.frame_index);
    let optimal_angle = find_optimal_rotation(
        &ref_mesh_dia.contours[0],
        x_coord_ref,
        y_coord_ref,
        z_coord_ref,
        x_coord_upper,
        y_coord_upper,
        z_coord_upper,
        x_coord_lower,
        y_coord_lower,
        z_coord_lower,
        0.01745329, // Step size in degrees
        &centerline.points[0],
    );

    println!("find optimal rotation worked");

    // Rotate every contour in every geometry in geometries by the optimal angle
    for geometry in &mut geometries {
        for contour in &mut geometry.contours {
            contour.rotate_contour(optimal_angle);
        }
    }

    // print elliptic ratio of the last and first frame
    let n = ref_mesh_dia.contours.len();
    let ((_, _), dist) = ref_mesh_dia.contours[0].find_farthest_points();
    let ((_, _), dist2) = ref_mesh_dia.contours[0].find_closest_opposite();
    // calculate distance between closest opposite points
    println!("Elliptic ratio Contour 0: {:?}", (dist / dist2));
    let ((_, _), dist) = ref_mesh_dia.contours[n - 1].find_farthest_points();
    let((_, _), dist2) = ref_mesh_dia.contours[n - 1].find_closest_opposite();
    println!("Elliptic ratio Contour n - 1: {:?}", (dist / dist2));

    let transformations = get_transformations(ref_mesh_dia, &centerline);

    for geometry in &mut geometries {
        for contour in &mut geometry.contours {
            if let Some(transformation) = transformations.iter().find(|t| t.frame_index == contour.id) {
                for pt in &mut contour.points {
                    *pt = transformation.apply_to_point(pt);
                }
                // Recompute the centroid after transformation
                contour.centroid = Contour::compute_centroid(&contour.points);
            } else {
                eprintln!("No transformation found for contour {}", contour.id);
            }
        }
        // Apply transformations to catheter points if applicable
        for catheter in &mut geometry.catheter {
            for pt in &mut catheter.points {
                // Assuming catheter points use the same frame index as their corresponding contour
                // Find the transformation based on catheter's frame index (if applicable)
                // This part depends on how catheter points are associated with frames
                // Example (adjust as needed):
                if let Some(transformation) = transformations.iter().find(|t| t.frame_index == catheter.id) {
                    *pt = transformation.apply_to_point(pt);
                }
            }
        }
    }
    
    // let n_transf = transformations.len();
    
    // write_centerline_to_csv("output/centerline.csv", &centerline)?;
    // write_frame_transformation_to_csv("output/transformation.csv", &transformations[0])?;
    // write_frame_transformation_to_csv("output/transformation_2.csv", &transformations[n_transf - 1])?;
    // write_geometry_to_csv("output/geometry.csv", &geometries[0])?;

    let (uv_coords_contours, uv_coords_catheter) = 
        write_mtl_geometry(&geometries, output_dir, state);

    // Write the interpolated geometries contours to OBJ files
    for (i, (mesh, uv_coords)) in geometries
        .iter()
        .zip(uv_coords_contours.iter()) // iterate UVs by reference
        .enumerate()
    {
        let obj_filename = format!("mesh_{:03}_{}.obj", i, state);
        let mtl_filename = format!("mesh_{:03}_{}.mtl", i, state);
        let obj_path = Path::new(output_dir).join(&obj_filename);
        let obj_path_str = obj_path.to_str().unwrap();

        // borrow mesh.contours and pass your &Vec<(f64,f64)> -> &[…] coerces automatically:
        write_obj_mesh(
            &mesh.contours, // ← borrow, don’t move
            uv_coords,      // ← &Vec<(f64,f64)> coerces to &[(f64,f64)]
            obj_path_str,
            &mtl_filename,
        )?;
    }

    // Write the interpolated geometries catheter to OBJ files
    for (i, (mesh, uv_coords)) in geometries
        .iter()
        .zip(uv_coords_catheter.iter()) // iterate UVs by reference
        .enumerate()
    {
        let obj_filename = format!("catheter_{:03}_{}.obj", i, state);
        let mtl_filename = format!("catheter_{:03}_{}.mtl", i, state);
        let obj_path = Path::new(output_dir).join(&obj_filename);
        let obj_path_str = obj_path.to_str().unwrap();

        // borrow mesh.contours and pass your &Vec<(f64,f64)> -> &[…] coerces automatically:
        write_obj_mesh(
            &mesh.catheter, // ← borrow, don’t move
            uv_coords,      // ← &Vec<(f64,f64)> coerces to &[(f64,f64)]
            obj_path_str,
            &mtl_filename,
        )?;
    }

    Ok(())
}
