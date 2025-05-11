use std::path::Path;
use anyhow::bail;

use crate::io::input::{Contour, ContourPoint};
use crate::io::output::write_obj_mesh;
use crate::io::Geometry;
use crate::processing::geometries::GeometryPair;
use crate::texture::write_mtl_geometry;

pub fn create_geometry_pair(
    case_name: String,
    input_dir: &str,
    steps_best_rotation: usize,
    range_rotation_rad: f64,
) -> anyhow::Result<GeometryPair> {
    let geometries = GeometryPair::new(input_dir, case_name.clone())?;
    println!("-------------------------Z-coordinates before adjustment-------------------------");
    let mut geometries = geometries.adjust_z_coordinates();

    geometries = geometries.process_geometry_pair(steps_best_rotation, range_rotation_rad);
    geometries = geometries.trim_geometries_same_length();
    geometries = geometries.thickness_adjustment();

    let dia_geom = geometries.dia_geom;
    let dia_geom = dia_geom.smooth_contours();
    let sys_geom = geometries.sys_geom;
    let sys_geom = sys_geom.smooth_contours(); 

    Ok(GeometryPair {
        dia_geom: dia_geom,
        sys_geom: sys_geom,
    })
}

/// Processes a given case by reading diastolic and systolic contours, aligning them,
/// computing displacements and UV coordinates, and finally writing out the OBJ, MTL, and texture files.
/// Additionally it can be specified how many interpolation steps should be used to generate the final meshes
/// used for the animation in blender.
pub fn process_case(
    case_name: &str,
    geometries: GeometryPair,
    output_dir: &str,
    interpolation_steps: usize,
) -> anyhow::Result<GeometryPair> {
    std::fs::create_dir_all(output_dir)?;

    let dia_geom = geometries.dia_geom;
    let sys_geom = geometries.sys_geom;

    // Interpolate between diastolic and systolic geometries
    let interpolated_geometries =
        interpolate_contours(&dia_geom, &sys_geom, interpolation_steps.clone())?;

    let (uv_coords_contours, uv_coords_catheter) =
        write_mtl_geometry(&interpolated_geometries, output_dir, case_name);

    // Write the interpolated geometries contours to OBJ files
    for (i, (mesh, uv_coords)) in interpolated_geometries
        .iter()
        .zip(uv_coords_contours.iter()) // iterate UVs by reference
        .enumerate()
    {
        let obj_filename = format!("mesh_{:03}_{}.obj", i, case_name);
        let mtl_filename = format!("mesh_{:03}_{}.mtl", i, case_name);
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
    for (i, (mesh, uv_coords)) in interpolated_geometries
        .iter()
        .zip(uv_coords_catheter.iter()) // iterate UVs by reference
        .enumerate()
    {
        let obj_filename = format!("catheter_{:03}_{}.obj", i, case_name);
        let mtl_filename = format!("catheter_{:03}_{}.mtl", i, case_name);
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

    Ok(GeometryPair { dia_geom, sys_geom })
}

/// Interpolates between two aligned Geometry configurations with number of steps
/// used to visualize deformation over a cardiac cycle.
pub fn interpolate_contours(
    start: &Geometry,
    end:   &Geometry,
    steps: usize,
) -> anyhow::Result<Vec<Geometry>> {
    use std::cmp::min;
    let n = min(start.contours.len(), end.contours.len());
    let sc = &start.contours[..n];
    let ec = &end.contours[..n];
    let sat = &start.catheter[..n];
    let eat = &end.catheter[..n];

    let mut geoms = Vec::with_capacity(steps + 2);
    // First frame
    geoms.push(start.clone());

    for step in 0..steps {
        let t = step as f64 / (steps - 1) as f64;

        let contours = sc.iter().zip(ec).map(|(s, e)| Contour {
            id: s.id,
            points: s.points.iter().zip(&e.points)
                .map(|(ps, pe)| ContourPoint {
                    frame_index: ps.frame_index,
                    point_index: ps.point_index, // or reassign later
                    x: ps.x * (1.0 - t) + pe.x * t,
                    y: ps.y * (1.0 - t) + pe.y * t,
                    z: ps.z * (1.0 - t) + pe.z * t,
                    aortic: ps.aortic,
                })
                .collect(),
            centroid: (
                s.centroid.0 * (1.0 - t) + e.centroid.0 * t,
                s.centroid.1 * (1.0 - t) + e.centroid.1 * t,
                s.centroid.2 * (1.0 - t) + e.centroid.2 * t,
            ),
            aortic_thickness: interpolate_thickness(&s.aortic_thickness, &e.aortic_thickness, t),
            pulmonary_thickness: interpolate_thickness(&s.pulmonary_thickness, &e.pulmonary_thickness, t),
        }).collect();

        let catheter = sat.iter().zip(eat).map(|(s, e)| {
            // same inner logic, or extract into a helper…
            Contour { 
                id: s.id,
                points: s.points.iter().zip(&e.points)
                .map(|(ps, pe)| ContourPoint {
                    frame_index: ps.frame_index,
                    point_index: ps.point_index,
                    x: ps.x * (1.0 - t) + pe.x * t,
                    y: ps.y * (1.0 - t) + pe.y * t,
                    z: ps.z * (1.0 - t) + pe.z * t,
                    aortic: ps.aortic,
                })
                .collect(),
                centroid: (
                    s.centroid.0 * (1.0 - t) + e.centroid.0 * t,
                    s.centroid.1 * (1.0 - t) + e.centroid.1 * t,
                    s.centroid.2 * (1.0 - t) + e.centroid.2 * t,
                ),
                aortic_thickness: Vec::new(),
                pulmonary_thickness: Vec::new(),
            }
        }).collect();

        geoms.push( Geometry {
            contours,
            catheter,
            reference_point: start.reference_point.clone(),
            label: format!("{}_inter_{}", start.label, step),
        } );
    }

    // Last frame
    geoms.push(end.clone());
    Ok(geoms)
}

/// Helper function to interpolate Vec<Contour>
#[allow(dead_code)]
fn interpolate_points(
    start_contours: &[Contour],
    end_contours: &[Contour],
    steps: usize,
    n: usize,
) -> anyhow::Result<Vec<Contour>> {
    let mut intermediate_contours = Vec::with_capacity(n);

    for step in 0..steps {
        let t = step as f64 / (steps - 1) as f64;

        // Pair matching contours between start and end geometries
        for (start_contour, end_contour) in start_contours.iter().zip(end_contours.iter()) {
            if start_contour.id != end_contour.id {
                bail!("Contour IDs do not match between start and end geometries");
            }
            if start_contour.points.len() != end_contour.points.len() {
                bail!("Contour point counts do not match between start and end");
            }

            // Interpolate points between contours
            let interp_points: Vec<ContourPoint> = start_contour
                .points
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
                aortic_thickness: interpolate_thickness(
                    &start_contour.aortic_thickness,
                    &end_contour.aortic_thickness,
                    t,
                ),
                pulmonary_thickness: interpolate_thickness(
                    &start_contour.pulmonary_thickness,
                    &end_contour.pulmonary_thickness,
                    t,
                ),
            };

            intermediate_contours.push(interp_contour);
        }
    }

    Ok(intermediate_contours)
}

/// Helper function to interpolate thickness values
fn interpolate_thickness(start: &[Option<f64>], end: &[Option<f64>], t: f64) -> Vec<Option<f64>> {
    start
        .iter()
        .zip(end.iter())
        .map(|(s, e)| match (s, e) {
            (Some(s_val), Some(e_val)) => Some(s_val * (1.0 - t) + e_val * t),
            _ => None,
        })
        .collect()
}