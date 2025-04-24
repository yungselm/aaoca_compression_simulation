use crate::io::input::{read_centerline_txt, Contour, ContourPoint};
use crate::io::Geometry;
use crate::io::load_geometry::rebuild_geometry;
use crate::io::input::Centerline;
use std::error::Error;

pub fn prepare_data_3d_alignment(
    state: &str,
    centerline_path: &str,
    input_dir: &str,
) -> Result<(Centerline, Geometry), Box<dyn Error>> {
    // ----- Build the common centerline -----
    let centerline_path = format!("{}/resampled_centerline.txt", centerline_path);
    let raw_centerline = read_centerline_txt(&centerline_path)?;
    let centerline = Centerline::from_contour_points(raw_centerline);

    // ----- Process the reference mesh: mesh_000_rest.obj or mesh_000_stress.obj -----
    let ref_mesh_path = format!("{}/mesh_000_{}.obj", input_dir, state);
    let catheter_path = format!("{}/catheter_000_{}.obj", input_dir, state);

    let mut reference_mesh = rebuild_geometry(&ref_mesh_path, &catheter_path);

    let n = reference_mesh.contours.len();
    let ((pt1, pt2), _) = reference_mesh.contours[n - 1].find_closest_opposite();
    
    let reference_point = if pt1.aortic {
        pt1.clone()
    } else {
        pt2.clone()
    };
    
    reference_mesh.label = state.to_string();
    reference_mesh.reference_point = reference_point;

    Ok((centerline, reference_mesh))
}