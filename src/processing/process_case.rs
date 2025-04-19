use std::error::Error;
use std::path::Path;

use crate::io::Geometry;
use crate::io::output::write_obj_mesh;
use crate::io::input::{Contour, ContourPoint};
use crate::texture::write_mtl_geometry;
use crate::processing::geometries::GeometryPair;

/// Processes a given case by reading diastolic and systolic contours, aligning them,
/// computing displacements and UV coordinates, and finally writing out the OBJ, MTL, and texture files.
/// Additionally it can be specified how many interpolation steps should be used to generate the final meshes
/// used for the animation in blender.
pub fn process_case(
    case_name: String,
    input_dir: &str,
    output_dir: &str,
    interpolation_steps: usize,
    steps_best_rotation: usize,
    range_rotation_rad: f64,
) -> Result<GeometryPair, Box<dyn Error>> {
    std::fs::create_dir_all(output_dir)?;
    
    let geometries = GeometryPair::new(input_dir, case_name.clone())?;
    let mut geometries = geometries.adjust_z_coordinates();
    geometries = geometries.process_geometry_pair(steps_best_rotation, range_rotation_rad);
    geometries = geometries.trim_geometries_same_length();
    
    let dia_geom = geometries.dia_geom;
    let dia_geom = dia_geom.smooth_contours();
    let sys_geom = geometries.sys_geom;
    let sys_geom = sys_geom.smooth_contours();
    
    let case_name_str = case_name.as_str();

    // Interpolate between diastolic and systolic geometries
    let interpolated_geometries = interpolate_contours(
        &dia_geom,
        &sys_geom,
        interpolation_steps.clone())?;

    let (uv_coords_contours, uv_coords_catheter) = write_mtl_geometry(
        &interpolated_geometries, 
        output_dir,
        case_name_str);

    // Write the interpolated geometries contours to OBJ files
    for (i, (mesh, uv_coords)) in interpolated_geometries
    .iter()
    .zip(uv_coords_contours.iter())    // iterate UVs by reference
    .enumerate()
    {
        let obj_filename = format!("mesh_{:03}_{}.obj", i, case_name);
        let mtl_filename = format!("mesh_{:03}_{}.mtl", i, case_name);
        let obj_path    = Path::new(output_dir).join(&obj_filename);
        let obj_path_str = obj_path.to_str().unwrap();
    
        // borrow mesh.contours and pass your &Vec<(f64,f64)> -> &[…] coerces automatically:
        write_obj_mesh(
            &mesh.contours,    // ← borrow, don’t move
            uv_coords,         // ← &Vec<(f64,f64)> coerces to &[(f64,f64)]
            obj_path_str,
            &mtl_filename,
        )?;
    }

    // Write the interpolated geometries catheter to OBJ files
    for (i, (mesh, uv_coords)) in interpolated_geometries
    .iter()
    .zip(uv_coords_catheter.iter())    // iterate UVs by reference
    .enumerate()
    {
        let obj_filename = format!("catheter_{:03}_{}.obj", i, case_name_str);
        let mtl_filename = format!("catheter_{:03}_{}.mtl", i, case_name_str);
        let obj_path    = Path::new(output_dir).join(&obj_filename);
        let obj_path_str = obj_path.to_str().unwrap();
    
        // borrow mesh.contours and pass your &Vec<(f64,f64)> -> &[…] coerces automatically:
        write_obj_mesh(
            &mesh.catheter,    // ← borrow, don’t move
            uv_coords,         // ← &Vec<(f64,f64)> coerces to &[(f64,f64)]
            obj_path_str,
            &mtl_filename,
        )?;
    }
    
    Ok(GeometryPair{
        dia_geom, 
        sys_geom,
    })
}

/// Interpolates between two aligned Geometry configurations with number of steps
/// used to visualize deformation over a cardiac cycle.
pub fn interpolate_contours(
    contours_start: &Geometry,
    contours_end: &Geometry,
    steps: usize,
) -> Result<Vec<Geometry>, Box<dyn Error>> {
    use std::cmp::min;

    // Ensure equal number of contours and matching IDs
    let n = min(contours_start.contours.len(), contours_end.contours.len());
    let start_contours = &contours_start.contours[0..n];
    let end_contours = &contours_end.contours[0..n];

    let start_catheter = &contours_start.catheter[0..n];
    let end_catheter = &contours_end.catheter[0..n];

    let lumen_contours = interpolate_points(start_contours, end_contours, steps, n);
    let catheter_contours = interpolate_points(start_catheter, end_catheter, steps, n);

    let mut interpolated_geometries = Vec::with_capacity(steps);

    for (step, (contour, catheter)) in lumen_contours.iter().zip(catheter_contours).enumerate() {
        interpolated_geometries.push(Geometry {
            contours: contour.clone(),
            catheter: catheter.clone(),
            reference_point: contours_start.reference_point.clone(), // reference point is fix
            label: format!("{}_inter_{}", contours_start.label, step),
        })
    }
    // Ensure diastole and systole are included in the interpolated_geometries
    interpolated_geometries.insert(0, contours_start.clone());
    interpolated_geometries.push(contours_end.clone());

    Ok(interpolated_geometries)
}

/// Helper function to interpolate Vec<Contour>
fn interpolate_points(
    start_contours: &[Contour],
    end_contours: &[Contour],
    steps: usize,
    n: usize,
) -> Result<Vec<Contour>, Box<dyn Error>> {
    let mut intermediate_contours = Vec::with_capacity(n);

    for step in 0..steps {
        let t = step as f64 / (steps - 1) as f64;

        // Pair matching contours between start and end geometries
        for (start_contour, end_contour) in start_contours.iter().zip(end_contours.iter()) {
            if start_contour.id != end_contour.id {
                return Err("Contour IDs do not match between start and end geometries".into());
            }
            if start_contour.points.len() != end_contour.points.len() {
                return Err("Contour point counts do not match between start and end".into());
            }

            // Interpolate points between contours
            let interp_points: Vec<ContourPoint> = start_contour.points
                .iter()
                .zip(end_contour.points.iter())
                .map(|(p_start, p_end)| {
                    // Linear interpolation for each coordinate
                    ContourPoint {
                        frame_index: p_start.frame_index,
                        point_index: p_start.point_index,
                        x: p_start.x * (1.0 - t) + p_end.x * t,
                        y: p_start.y * (1.0 - t) + p_end.y * t,
                        z: p_start.z * (1.0 - t) + p_end.z * t,
                        aortic: p_start.aortic,
                    }
                })
                .collect();

            // Create new Contour with interpolated values
            let interp_contour = Contour {
                id: start_contour.id,
                points: interp_points,
                centroid: (
                    start_contour.centroid.0 * (1.0 - t) + end_contour.centroid.0 * t,
                    start_contour.centroid.1 * (1.0 - t) + end_contour.centroid.1 * t,
                    start_contour.centroid.2 * (1.0 - t) + end_contour.centroid.2 * t,
                ),
                aortic_thickness: interpolate_thickness(&start_contour.aortic_thickness, &end_contour.aortic_thickness, t),
                pulmonary_thickness: interpolate_thickness(&start_contour.pulmonary_thickness, &end_contour.pulmonary_thickness, t),
            };

            intermediate_contours.push(interp_contour);
        }
    }

    Ok(intermediate_contours)    
}

/// Helper function to interpolate thickness values
fn interpolate_thickness(
    start: &[Option<f64>],
    end: &[Option<f64>],
    t: f64
) -> Vec<Option<f64>> {
    start.iter()
        .zip(end.iter())
        .map(|(s, e)| match (s, e) {
            (Some(s_val), Some(e_val)) => Some(s_val * (1.0 - t) + e_val * t),
            _ => None
        })
        .collect()
}