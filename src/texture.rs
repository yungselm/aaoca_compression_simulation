use std::error::Error;
use image::{ImageBuffer, Rgb};
use crate::io::ContourPoint;

pub fn compute_uv_coordinates(contours: &[(u32, Vec<ContourPoint>)]) -> Vec<(f32, f32)> {
    let points_per_contour = contours[0].1.len();
    let num_contours = contours.len();
    let mut uvs = Vec::with_capacity(num_contours * points_per_contour);

    for c in 0..num_contours {
        // Use 0.5_f32 so that the literal is f32.
        let v = (c as f32 + 0.5_f32) / num_contours as f32;
        for p in 0..points_per_contour {
            let u = (p as f32 + 0.5_f32) / points_per_contour as f32;
            uvs.push((u, v));
        }
    }

    uvs
}

pub fn compute_displacements(
    mesh: &[(u32, Vec<ContourPoint>)],
    diastole: &[(u32, Vec<ContourPoint>)],
) -> Vec<f32> {
    let mut displacements = Vec::new();
    for (c, (_, contour)) in mesh.iter().enumerate() {
        let diastole_contour = &diastole[c].1;
        for (p, point) in contour.iter().enumerate() {
            let diastole_point = &diastole_contour[p];
            let dx = point.x - diastole_point.x;
            let dy = point.y - diastole_point.y;
            let dz = point.z - diastole_point.z;
            // Compute the displacement and cast to f32.
            displacements.push(((dx * dx + dy * dy + dz * dz).sqrt()) as f32);
        }
    }
    displacements
}

pub fn create_displacement_texture(
    displacements: &[f32],
    width: u32,
    height: u32,
    max_displacement: f32,
    filename: &str,
) -> Result<(), Box<dyn Error>> {
    let mut img = ImageBuffer::new(width, height);
    for (i, &disp) in displacements.iter().enumerate() {
        let x = (i % width as usize) as u32;
        let y = (i / width as usize) as u32;
        let normalized = (disp / max_displacement).clamp(0.0, 1.0);
        let color = Rgb([
            (normalized * 255.0) as u8,
            0,
            ((1.0 - normalized) * 255.0) as u8,
        ]);
        img.put_pixel(x, y, color);
    }
    img.save(filename)?; // Save as PNG
    Ok(())
}

pub fn create_black_texture(width: u32, height: u32, filename: &str) -> Result<(), Box<dyn Error>> {
    let black = Rgb([0u8, 0u8, 0u8]); // Ensure pixel values are u8
    let img = ImageBuffer::from_pixel(width, height, black);
    img.save(filename)?; // Save as PNG
    Ok(())
}