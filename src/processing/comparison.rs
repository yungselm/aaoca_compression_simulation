use super::geometries::GeometryPair;
use crate::io::input::{Contour, ContourPoint};
use crate::io::Geometry;
use std::error::Error;

pub fn prepare_geometries_comparison(
    geometries_rest: GeometryPair,
    geometries_stress: GeometryPair,
) -> (GeometryPair, GeometryPair) {
    let dia_rest = geometries_rest.dia_geom;
    let sys_rest = geometries_rest.sys_geom;
    let dia_stress = geometries_stress.dia_geom;
    let sys_stress = geometries_stress.sys_geom;

    let dia_rest = resample_contours_with_reference_z(&dia_rest, &dia_stress).unwrap();
    let sys_rest = resample_contours_with_reference_z(&sys_rest, &sys_stress).unwrap();

    let dia_rest = align_geometries(&dia_stress, dia_rest);
    let sys_rest = align_geometries(&sys_stress, sys_rest);

    let dia_pair = GeometryPair {
        dia_geom: dia_rest,
        sys_geom: dia_stress,
    };
    let sys_pair = GeometryPair {
        dia_geom: sys_rest,
        sys_geom: sys_stress,
    };
    
    let dia_pair = dia_pair.trim_geometries_same_length();
    let sys_pair = sys_pair.trim_geometries_same_length();
    
    let dia_pair = dia_pair.adjust_z_coordinates();
    let sys_pair = sys_pair.adjust_z_coordinates();

    (dia_pair, sys_pair)
}

/// Resamples a stack of contours so that the new z‑spacing exactly matches the spacing
/// from the reference (stress) contours, while preserving the original z‑direction.
///
/// Assumptions:
///   - `original` contours (e.g. rest) are sorted in order (either increasing or decreasing in z).
///   - Each contour’s z value is constant (all points in a given contour share the same z).
///   - There are at least two contours in both `original` and `reference`.
fn resample_contours_with_reference_z(
    original: &Geometry,
    reference: &Geometry,
) -> Result<Geometry, Box<dyn Error>> {
    let original_geometry = original.clone();
    let n = original.contours.len();

    // Determine the z-direction by comparing the first and last contour's z.
    let start_z = original.contours[0].points[0].z;
    let end_z = original.contours[n - 1].points[0].z;
    let direction = if end_z >= start_z { 1.0 } else { -1.0 };
    
    // Compute the cumulative z distance (absolute differences) along the original stack.
    let mut cum_z = Vec::with_capacity(n);
    let mut total = 0.0;
    cum_z.push(0.0); // starting cumulative distance is 0.
    for i in 1..n {
        let prev_z = original.contours[i - 1].points[0].z;
        let curr_z = original.contours[i].points[0].z;
        total += (curr_z - prev_z).abs();
        cum_z.push(total);
    }
    
    // Determine target spacing from the reference contours.
    // We use the absolute difference between the first two reference contours.
    let total_spacing: f64 = reference
        .contours
        .windows(2)
        .map(|pair| (pair[1].centroid.2 - pair[0].centroid.2).abs())
        .sum();
    let target_spacing = total_spacing / (reference.contours.len() - 1) as f64;
    // Determine the number of new frames to generate.
    let new_frame_count = (total / target_spacing).floor() as usize + 1;
    
    // Generate new contours.
    let new_contours: Vec<Vec<ContourPoint>> = contours_new_z_spacing(
        n.clone(),
        new_frame_count.clone(),
        target_spacing.clone(),
        start_z.clone(),
        direction.clone(),
        cum_z.clone(),
        &original.contours,
    );
    
    let new_catheter: Vec<Vec<ContourPoint>> = contours_new_z_spacing(
        n,
        new_frame_count,
        target_spacing,
        start_z,
        direction,
        cum_z,
        &original.catheter,
    );
    
    let new_contours = vec_contour_points_to_contour(new_contours);
    let new_catheter = vec_contour_points_to_contour(new_catheter);

    let new_geometry = Geometry {
        contours: new_contours,
        catheter: new_catheter,
        reference_point: original_geometry.reference_point,
        label: original_geometry.label,
    };
    Ok(new_geometry)
}

fn contours_new_z_spacing(
    n: usize,
    new_frame_count: usize,
    target_spacing: f64,
    start_z: f64,
    direction: f64,
    cum_z: Vec<f64>,
    contours: &Vec<Contour>,
) -> Vec<Vec<ContourPoint>> {
    let mut new_contours = Vec::with_capacity(new_frame_count);
    
    for j in 0..new_frame_count {
        // Compute the target cumulative distance and corresponding target z.
        let target_cum = j as f64 * target_spacing;
        let target_z = start_z + direction * target_cum;

        // Find the bracket in the original cum_z where target_cum falls.
        let mut i = 0;
        while i < cum_z.len() - 1 && cum_z[i + 1] < target_cum {
            i += 1;
        }
        if i >= n - 1 {
            // If beyond range, simply copy the last contour and assign target_z.
            let new_points = contours[n - 1]
            .points
            .iter()
            .map(|pt| ContourPoint {
                frame_index: j as u32,
                point_index: pt.point_index,
                x: pt.x,
                y: pt.y,
                z: target_z,
                aortic: pt.aortic,
            })
            .collect();
        new_contours.push(new_points);
    } else {
            // Interpolate between original[i] and original[i+1].
            let z0 = cum_z[i];
            let z1 = cum_z[i + 1];
            let t = if (z1 - z0).abs() < std::f64::EPSILON {
                0.0
            } else {
                (target_cum - z0) / (z1 - z0)
            };

            // For each corresponding point (assuming same number per contour),
            // linearly interpolate x and y, and force z to target_z.
            let new_points: Vec<ContourPoint> = contours[i]
            .points
            .iter()
            .zip(contours[i + 1].points.iter())
            .map(|(pt0, pt1)| ContourPoint {
                frame_index: j as u32,
                point_index: pt0.point_index,
                    x: pt0.x * (1.0 - t) + pt1.x * t,
                    y: pt0.y * (1.0 - t) + pt1.y * t,
                    z: target_z,
                    aortic: pt0.aortic,
                })
                .collect();
            new_contours.push(new_points);
        }
    }
    new_contours
}

fn vec_contour_points_to_contour(vec_contour: Vec<Vec<ContourPoint>>) -> Vec<Contour> {
    let mut contours = Vec::new();
    
    for vec in vec_contour {
        let centroid = Contour::compute_centroid(&vec);
        
        let contour = Contour {
            id: vec[0].frame_index,
            points: vec,
            centroid: centroid,
            aortic_thickness: Vec::new(),
            pulmonary_thickness: Vec::new(),
        };

        contours.push(contour);
    }
    contours
}

fn align_geometries(
    ref_geom: &Geometry,
    mut geom: Geometry,
) -> Geometry {
    let len_ref_geom = ref_geom.contours.len();
    let n_geom = geom.contours.len();
    let ref_centroid = ref_geom.contours[len_ref_geom - 1].centroid;

    for ref mut contour in geom
        .contours
        .iter_mut()
        .chain(geom.catheter.iter_mut())
        {
            let centroid = contour.centroid;
            let t = (
                ref_centroid.0 - centroid.0,
                ref_centroid.1 - centroid.1,
            );
            contour.translate_contour((t.0, t.1, 0.0));
        }
    
    // Adjust the z-coordinates of systolic contours. (later replaceed by adjust_z_coordinates)
    let z_translation = ref_geom.contours[len_ref_geom - 1].centroid.2 - geom.contours[n_geom - 1].centroid.2;
    for ref mut contour in geom
        .contours
        .iter_mut()
        .chain(geom.catheter.iter_mut())
    {
        for point in contour.points.iter_mut() {
            point.z += z_translation;
        }
        contour.centroid.2 += z_translation;
    }
    geom
}