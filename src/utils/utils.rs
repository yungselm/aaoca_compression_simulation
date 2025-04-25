use std::error::Error;
use std::path::Path;
use csv::Writer;

use crate::io::Geometry;
use crate::io::input::{Contour, Centerline};
use crate::mesh_to_centerline::operations::FrameTransformation;

// debugging functions
pub fn write_frame_transformation_to_csv<P: AsRef<Path>>(
    path: P,
    frame: &FrameTransformation,
) -> Result<(), Box<dyn Error>> {
    let mut wtr = Writer::from_path(path)?;

    // Write header
    wtr.write_record(&[
        "frame_index",
        "translation_x", "translation_y", "translation_z",
        "rotation_00", "rotation_01", "rotation_02",
        "rotation_10", "rotation_11", "rotation_12",
        "rotation_20", "rotation_21", "rotation_22",
        "pivot_x", "pivot_y", "pivot_z",
    ])?;

    let rot = frame.rotation.matrix(); // 3x3 rotation matrix

    // Write data row
    wtr.write_record(&[
        frame.frame_index.to_string(),
        frame.translation.x.to_string(),
        frame.translation.y.to_string(),
        frame.translation.z.to_string(),
        rot[(0, 0)].to_string(), rot[(0, 1)].to_string(), rot[(0, 2)].to_string(),
        rot[(1, 0)].to_string(), rot[(1, 1)].to_string(), rot[(1, 2)].to_string(),
        rot[(2, 0)].to_string(), rot[(2, 1)].to_string(), rot[(2, 2)].to_string(),
        frame.pivot.x.to_string(),
        frame.pivot.y.to_string(),
        frame.pivot.z.to_string(),
    ])?;

    wtr.flush()?;
    Ok(())
}

pub fn write_centerline_to_csv<P: AsRef<Path>>(
    path: P,
    centerline: &Centerline,
) -> Result<(), Box<dyn Error>> {
    let mut wtr = Writer::from_path(path)?;

    // Write header
    wtr.write_record(&[
        "frame_index", "point_index",
        "x", "y", "z",
        "aortic",
        "normal_x", "normal_y", "normal_z",
    ])?;

    for point in &centerline.points {
        let cp = &point.contour_point;
        wtr.write_record(&[
            cp.frame_index.to_string(),
            cp.point_index.to_string(),
            cp.x.to_string(),
            cp.y.to_string(),
            cp.z.to_string(),
            cp.aortic.to_string(),
            point.normal.x.to_string(),
            point.normal.y.to_string(),
            point.normal.z.to_string(),
        ])?;
    }

    wtr.flush()?;
    Ok(())
}

pub fn write_geometry_to_csv<P: AsRef<Path>>(
    path: P,
    geometry: &Geometry,
) -> Result<(), Box<dyn Error>> {
    let mut wtr = Writer::from_path(path)?;

    // Write header
    wtr.write_record(&[
        "label",
        "source",
        "contour_id",
        "point_index",
        "frame_index",
        "x", "y", "z",
        "aortic",
        "aortic_thickness",
        "pulmonary_thickness",
    ])?;

    // Helper to write from a source (contours or catheter)
    let write_contours = |source: &str, contours: &Vec<Contour>, wtr: &mut Writer<std::fs::File>| -> Result<(), Box<dyn Error>> {
        for contour in contours {
            for (i, point) in contour.points.iter().enumerate() {
                let a_thickness = contour.aortic_thickness.get(i).and_then(|x| *x);
                let p_thickness = contour.pulmonary_thickness.get(i).and_then(|x| *x);

                let record = vec![
                    geometry.label.clone(),
                    source.to_string(),
                    contour.id.to_string(),
                    point.point_index.to_string(),
                    point.frame_index.to_string(),
                    point.x.to_string(),
                    point.y.to_string(),
                    point.z.to_string(),
                    point.aortic.to_string(),
                    a_thickness.map_or("".into(), |v| v.to_string()),
                    p_thickness.map_or("".into(), |v| v.to_string()),
                ];

                wtr.write_record(&record)?;
            }
        }
        Ok(())
    };

    write_contours("contour", &geometry.contours, &mut wtr)?;
    write_contours("catheter", &geometry.catheter, &mut wtr)?;

    wtr.flush()?;
    Ok(())
}