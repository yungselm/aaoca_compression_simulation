
use std::error::Error;
use std::fs::{create_dir_all, read_dir, File};
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::f64::consts::PI;

use nalgebra::{center, Point3, Rotation3, Unit, Vector3};
use rand::Rng;

use crate::contour;
use crate::io::{read_centerline_txt, read_obj_mesh, write_obj_mesh, ContourPoint};

const FIXED_ROTATION_DEG: f64 = 30.0;

#[derive(Debug, Clone, PartialEq)]
pub struct ContourFrame {
    pub frame_index: u32,
    pub points: Vec<ContourPoint>,  // 500 points
    pub centroid_3d: ContourPoint,
    pub normal: Vector3<f64>,       // Normal vector for orientation
    pub translation: Point3<f64>,   // Translation from centerline
}

impl ContourFrame {
    pub fn from_contour(frame_index: u32, points: Vec<ContourPoint>) -> Self {
        let centroid_3d = Self::calculate_centroid(&points);
        let normal = Self::calculate_normal(&points, &centroid_3d);
        
        ContourFrame {
            frame_index,
            points,
            centroid_3d,
            normal,
            translation: Point3::origin(), // Will be set later using centerline
        }
    }

    fn calculate_centroid(points: &[ContourPoint]) -> ContourPoint {
        let count = points.len() as f64;
        let (sum_x, sum_y, sum_z) = points.iter()
            .fold((0.0, 0.0, 0.0), |(sx, sy, sz), p| {
                (sx + p.x, sy + p.y, sz + p.z)
            });

        ContourPoint {
            frame_index: points.first().map(|p| p.frame_index).unwrap_or(0),
            x: sum_x / count,
            y: sum_y / count,
            z: sum_z / count,
        }
    }

    fn calculate_normal(points: &[ContourPoint], centroid: &ContourPoint) -> Vector3<f64> {
        // Use first point and midpoint to create vectors
        let p1 = &points[0];
        let p2 = &points[points.len() / 2];
        
        // Create vectors relative to centroid
        let v1 = Vector3::new(
            p1.x - centroid.x,
            p1.y - centroid.y,
            p1.z - centroid.z
        );
        
        let v2 = Vector3::new(
            p2.x - centroid.x,
            p2.y - centroid.y,
            p2.z - centroid.z
        );

        // Calculate cross product and normalize
        v1.cross(&v2).normalize()
    }
}

// Conversion function for your mesh data
pub fn process_mesh_frames(
    mesh: Vec<(u32, Vec<ContourPoint>)>,
    centerline: &Centerline
) -> Vec<ContourFrame> {
    mesh.into_iter()
        .map(|(frame_idx, points)| {
            let mut frame = ContourFrame::from_contour(frame_idx, points);
            
            // Set translation from centerline
            if let Some(cl_point) = centerline.get_by_frame(frame_idx) {
                frame.translation = Point3::new(
                    cl_point.contour_point.x,
                    cl_point.contour_point.y,
                    cl_point.contour_point.z
                );
            }
            
            frame
        })
        .collect()
}

#[derive(Debug, Clone, PartialEq)]
pub struct Centerline {
    pub points: Vec<CenterlinePoint>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct CenterlinePoint {
    pub contour_point: ContourPoint,
    pub normal: Vector3<f64>,
}

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

// fn calculate_centroid_3d(mesh: &mut Vec<(u32, Vec<ContourPoint>)>) {

// }

impl Centerline {
    pub fn from_contour_points(contour_points: Vec<ContourPoint>) -> Self {
        let mut points: Vec<CenterlinePoint> = Vec::with_capacity(contour_points.len());
        
        // Calculate normals for all but last point
        for i in 0..contour_points.len() {
            let current = &contour_points[i];
            
            // Calculate direction to next point (or previous if last)
            let normal = if i < contour_points.len() - 1 {
                let next = &contour_points[i + 1];
                Vector3::new(
                    next.x - current.x,
                    next.y - current.y,
                    next.z - current.z
                ).normalize()
            } else if !contour_points.is_empty() {
                // For last point, use previous normal
                points[i - 1].normal
            } else {
                Vector3::zeros()
            };

            points.push(CenterlinePoint {
                contour_point: current.clone(),
                normal,
            });
        }

        Centerline { points }
    }

    // Helper to find centerline point by frame index
    pub fn get_by_frame(&self, frame_index: u32) -> Option<&CenterlinePoint> {
        self.points.iter().find(|p| p.contour_point.frame_index == frame_index)
    }
}

pub fn test_function() -> Result<(), Box<dyn Error>> {
    let raw_centerline = read_centerline_txt("resampled_centerline.txt")?;
    let mut mesh = read_obj_mesh("output/stress/mesh_000_stress.obj")?;
    
    // Rotate first
    rotate_contours_around_z(&mut mesh, FIXED_ROTATION_DEG);
    
    // Process into frames
    let centerline = Centerline::from_contour_points(raw_centerline);
    let frames = process_mesh_frames(mesh, &centerline);
    
    // Print first frame's data
    if let Some(first_frame) = frames.first() {
        println!("First Contour Frame:");
        println!("Frame Index: {}", first_frame.frame_index);
        println!("Centroid: ({:.4}, {:.4}, {:.4})", 
            first_frame.centroid_3d.x,
            first_frame.centroid_3d.y,
            first_frame.centroid_3d.z
        );
        println!("Normal: [{:.4}, {:.4}, {:.4}]", 
            first_frame.normal.x,
            first_frame.normal.y,
            first_frame.normal.z
        );
        println!("Translation: {}", first_frame.translation);
    }
    
    Ok(())
}
