use rayon::prelude::*;
use std::error::Error;
use std::fs::File;
use std::io::Write;
use std::path::Path;

use crate::io::Geometry;
use crate::io::output::write_obj_mesh;
use crate::io::input::{Contour, ContourPoint};
use crate::texture::{
    compute_displacements, compute_uv_coordinates, create_black_texture,
    create_displacement_texture,
};
use crate::utils::trim_to_same_length;

use crate::processing::{geometries::GeometryPair};

use super::contours::align_frames_in_geometry;

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
) -> Result<(), Box<dyn Error>> {
    std::fs::create_dir_all(output_dir)?;

    let geometries = GeometryPair::new(input_dir, case_name)?;
    let mut geometries = geometries.adjust_z_coordinates();
    geometries = geometries.process_geometry_pair(steps_best_rotation, range_rotation_rad);
    geometries = geometries.trim_geometries_same_length();

    let dia_geom = geometries.dia_geom;
    dia_geom.smooth_contours();
    let sys_geom = geometries.sys_geom;
    sys_geom.smooth_contours();

    // Interpolate between diastolic and systolic geometries
    let interpolated_geometries: Vec<Geometry> = (0..=interpolation_steps)
        .into_par_iter()
        .map(|step| {
            let t = step as f64 / interpolation_steps as f64;
            interpolate_contours(&dia_geom, &sys_geom, interpolation_steps)
        })
        .collect();

    // Write the interpolated geometries to OBJ files
    for (i, geometry) in interpolated_geometries.iter().enumerate() {
        let file_name = format!("{}/interpolated_{:03}.obj", output_dir, i);
        let path = Path::new(&file_name);
        let mut file = File::create(path)?;
        write_obj_mesh(&mut file, geometry)?;
    }

    Ok(())
}

/// Interpolates between two aligned Geometry configurations
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

    let mut interpolated_geometries = Vec::with_capacity(steps);
    
    for step in 0..steps {
        let t = step as f64 / (steps - 1) as f64;
        let mut intermediate_contours = Vec::with_capacity(n);

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

        // Create new Geometry for this interpolation step
        interpolated_geometries.push(Geometry {
            contours: intermediate_contours,
            catheter: Vec::new(), // You'll need to implement catheter interpolation separately
            reference_point: interpolate_reference_point(
                &contours_start.reference_point,
                &contours_end.reference_point,
                t
            ),
            label: format!("{}_interp_{}", contours_start.label, step),
        });
    }

    Ok(interpolated_geometries)
}

/// Helper function to interpolate thickness values
pub fn interpolate_thickness(
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

/// Helper function to interpolate reference points
pub fn interpolate_reference_point(
    start: &ContourPoint,
    end: &ContourPoint,
    t: f64
) -> ContourPoint {
    ContourPoint {
        frame_index: start.frame_index,
        point_index: start.point_index,
        x: start.x * (1.0 - t) + end.x * t,
        y: start.y * (1.0 - t) + end.y * t,
        z: start.z * (1.0 - t) + end.z * t,
        aortic: start.aortic,
    }
}