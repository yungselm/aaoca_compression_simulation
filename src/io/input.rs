use csv::ReaderBuilder;
use nalgebra::Vector3;
use serde::Deserialize;
use std::collections::HashMap;
use std::error::Error;
use std::f64::consts::PI;
use std::fs::File;
use std::path::Path;

#[derive(Debug, Deserialize, Clone, PartialEq)]
pub struct Contour {
    pub id: u32,
    pub points: Vec<ContourPoint>,
    pub centroid: (f64, f64, f64),
    pub aortic_thickness: Vec<Option<f64>>,
    pub pulmonary_thickness: Vec<Option<f64>>,
}

impl Contour {
    pub fn create_contours(
        points: Vec<ContourPoint>,
        result: Vec<Record>,
    ) -> Result<Vec<Contour>, Box<dyn Error>> {
        let mut groups: HashMap<u32, Vec<ContourPoint>> = HashMap::new();
        for p in points {
            groups.entry(p.frame_index).or_default().push(p);
        }

        let mut contours = Vec::new();
        for (frame_index, group_points) in groups {
            let centroid = Self::compute_centroid(&group_points);
            let aortic_thickness = result
                .iter()
                .find(|r| r.frame == frame_index as u32)
                .and_then(|r| r.measurement_1);

            let pulmonary_thickness = result
                .iter()
                .find(|r| r.frame == frame_index as u32)
                .and_then(|r| r.measurement_2);

            contours.push(Contour {
                id: frame_index,
                points: group_points,
                centroid,
                aortic_thickness: vec![aortic_thickness],
                pulmonary_thickness: vec![pulmonary_thickness],
            });

        for contour in &mut contours {
            for (i, point) in contour.points.iter_mut().enumerate() {
                point.point_index = i as u32;
            }
        }
        }
        Ok(contours)
    }

    pub fn create_catheter_contours(
        points: &Vec<ContourPoint>,
    ) -> Result<Vec<Contour>, Box<dyn Error>> {
        let catheter_points = ContourPoint::create_catheter_points(&points);

        let mut groups: HashMap<u32, Vec<ContourPoint>> = HashMap::new();
        for p in catheter_points {
            groups.entry(p.frame_index).or_default().push(p);
        }

        let mut contours = Vec::new();
        for (frame_index, group_points) in groups {
            let centroid = Self::compute_centroid(&group_points);
            let aortic_thickness = None;
            let pulmonary_thickness = None;

            contours.push(Contour {
                id: frame_index,
                points: group_points,
                centroid,
                aortic_thickness: vec![aortic_thickness],
                pulmonary_thickness: vec![pulmonary_thickness],
            });
        }
        Ok(contours)
    }

    pub fn compute_centroid(points: &Vec<ContourPoint>) -> (f64, f64, f64) {
        let (sum_x, sum_y, sum_z) = points.iter().fold((0.0, 0.0, 0.0), |(sx, sy, sz), p| {
            (sx + p.x, sy + p.y, sz + p.z)
        });
        let n = points.len() as f64;
        (sum_x / n, sum_y / n, sum_z / n)
    }

    /// Finds the pair of farthest points in the current contour.
    pub fn find_farthest_points(&self) -> ((&ContourPoint, &ContourPoint), f64) {
        let mut max_dist = 0.0;
        let mut farthest_pair = (&self.points[0], &self.points[0]);

        for i in 0..self.points.len() {
            for j in i + 1..self.points.len() {
                let dx = self.points[i].x - self.points[j].x;
                let dy = self.points[i].y - self.points[j].y;
                let dist = (dx * dx + dy * dy).sqrt();
                if dist > max_dist {
                    max_dist = dist;
                    farthest_pair = (&self.points[i], &self.points[j]);
                }
            }
        }

        (farthest_pair, max_dist)
    }

    /// Find the closest opposite points on a contour, where "opposite" means a point
    /// and the point halfway around the contour vector (e.g., 0 with n/2, 1 with n/2+1, etc.).
    /// If the number of points in the contour is odd, the function will ignore the last point.
    pub fn find_closest_opposite(&self) -> ((&ContourPoint, &ContourPoint), f64) {
        let n = self.points.len();
        // If odd, ignore the last point
        let effective_n = if n % 2 == 0 { n } else { n - 1 };
        if effective_n < 2 {
            panic!("Not enough points to evaluate");
        }

        let half = effective_n / 2;
        let mut min_dist = f64::MAX;
        let mut closest_pair = (&self.points[0], &self.points[half]);

        for i in 0..half {
            let j = i + half; // Directly opposite index (ignoring the extra odd element)
            let dx = self.points[i].x - self.points[j].x;
            let dy = self.points[i].y - self.points[j].y;
            let dist = (dx * dx + dy * dy).sqrt(); // Euclidean distance

            if dist < min_dist {
                min_dist = dist;
                closest_pair = (&self.points[i], &self.points[j]);
            }
        }
        (closest_pair, min_dist)
    }

    /// Rotates all points in a contour about a center.
    pub fn rotate_contour(&mut self, angle: f64) {
        for p in self.points.iter_mut() {
            let rotated = ContourPoint::rotate_point(p, angle, (self.centroid.0, self.centroid.1));
            p.x = rotated.x;
            p.y = rotated.y;
        }
    }

    pub fn rotate_contour_around_point(&mut self, angle: f64, center: (f64, f64)) {
        for p in self.points.iter_mut() {
            let rotated = ContourPoint::rotate_point(p, angle, center);
            p.x = rotated.x;
            p.y = rotated.y;
        }
    }

    /// Reorders the point indices so that the point with the highest y-value is 0,
    /// and the others are numbered counterclockwise around the centroid.
    pub fn sort_contour_points(&mut self) {
        // Find the index of the point with the highest y-value
        let start_idx = self
            .points
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.y.partial_cmp(&b.y).unwrap())
            .map(|(i, _)| i)
            .unwrap_or(0);

        // Rotate the points so that the highest y-value point is first
        self.points.rotate_left(start_idx);

        // Reassign point indices in counterclockwise order
        for (i, point) in self.points.iter_mut().enumerate() {
            point.point_index = i as u32;
        }
    }

    /// Translates a contour by a given (dx, dy, dz) offset and recalculates the centroid.
    pub fn translate_contour(&mut self, translation: (f64, f64, f64)) {
        let (dx, dy, dz) = translation;
        for p in self.points.iter_mut() {
            p.x += dx;
            p.y += dy;
            p.z += dz;
        }
        // Recalculate the centroid
        self.centroid = Self::compute_centroid(&self.points);
    }
}

#[derive(Debug, Deserialize, Clone, PartialEq)]
pub struct ContourPoint {
    pub frame_index: u32,

    #[serde(default, skip_deserializing)]
    pub point_index: u32,
    
    pub x: f64,
    pub y: f64,
    pub z: f64,

    #[serde(default)]
    pub aortic: bool,
}

impl ContourPoint {
    /// Reads contour points from a CSV file.
    pub fn read_contour_data<P: AsRef<Path> + std::fmt::Debug + Clone>(path: P) -> Result<Vec<ContourPoint>, Box<dyn Error>> {
        let debug_path = path.clone();
        let file = File::open(path)?;
        let mut rdr = csv::ReaderBuilder::new()
            .has_headers(false)
            .delimiter(b'\t')
            .from_reader(file);

        let mut points = Vec::new();
        for result in rdr.records() {
            match result {
                Ok(record) => match record.deserialize(None) {
                    Ok(point) => points.push(point),
                    Err(e) => eprintln!("Skipping invalid record: {:?}", e),
                },
                Err(e) => eprintln!("Skipping invalid row: {:?}", e),
            }
        }
        println!("Reading from {:?}", debug_path);
        Ok(points)
    }

    pub fn read_reference_point<P: AsRef<Path>>(path: P) -> Result<ContourPoint, Box<dyn Error>> {
        let file = File::open(path)?;
        let mut rdr = csv::ReaderBuilder::new()
            .has_headers(false)
            .delimiter(b'\t')
            .from_reader(file);

        for result in rdr.records() {
            match result {
                Ok(record) => match record.deserialize(None) {
                    Ok(point) => return Ok(point),
                    Err(e) => return Err(Box::new(e)),
                },
                Err(e) => return Err(Box::new(e)),
            }
        }

        Err("No valid reference point found in file".into())
    }

    pub fn create_catheter_points(points: &Vec<ContourPoint>) -> Vec<ContourPoint> {
        // Map to store unique frame indices and one associated z coordinate per frame.
        let mut frame_z: HashMap<u32, f64> = HashMap::new();
        for point in points {
            // Use the first encountered z-coordinate for each frame index.
            frame_z.entry(point.frame_index).or_insert(point.z);
        }

        let mut catheter_points = Vec::new();
        // Sort the frame indices to ensure a predictable order.
        let mut frames: Vec<u32> = frame_z.keys().cloned().collect();
        frames.sort();

        // Parameters for the catheter circle.
        let center_x = 4.5;
        let center_y = 4.5;
        let radius = 0.5;
        let num_points = 20;

        // For each unique frame, generate 20 catheter points around a circle.
        for frame in frames {
            let z = frame_z[&frame];
            for i in 0..num_points {
                let angle = 2.0 * PI * (i as f64) / (num_points as f64);
                let x = center_x + radius * angle.cos();
                let y = center_y + radius * angle.sin();
                catheter_points.push(ContourPoint {
                    frame_index: frame,
                    point_index: i,
                    x,
                    y,
                    z,
                    aortic: false,
                });
            }
        }
        catheter_points
    }

    /// Computes the Euclidean distance between two contour points.
    pub fn distance_to(&self, other: &ContourPoint) -> f64 {
        let dx = self.x - other.x;
        let dy = self.y - other.y;
        let dz = self.z - other.z;
        (dx * dx + dy * dy + dz * dz).sqrt()
    }

    /// Rotates a single point about a given center (cx, cy) by a specified angle (in radians).
    pub fn rotate_point(&self, angle: f64, center: (f64, f64)) -> ContourPoint {
        let (cx, cy) = center;
        let x = self.x - cx;
        let y = self.y - cy;
        let cos_a = angle.cos();
        let sin_a = angle.sin();
        ContourPoint {
            frame_index: self.frame_index,
            point_index: self.point_index,
            x: x * cos_a - y * sin_a + cx,
            y: x * sin_a + y * cos_a + cy,
            z: self.z,
            aortic: self.aortic,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Deserialize)]
pub struct Record {
    pub frame: u32,
    pub phase: String,
    #[serde(deserialize_with = "csv::invalid_option")]
    pub measurement_1: Option<f64>,
    #[serde(deserialize_with = "csv::invalid_option")]
    pub measurement_2: Option<f64>,
}

pub fn read_records<P: AsRef<Path>>(path: P) -> Result<Vec<Record>, Box<dyn Error>> {
    let file = File::open(path)?;
    let mut reader = ReaderBuilder::new()
        .delimiter(b',')
        .has_headers(true)
        .from_reader(file);

    let mut records = Vec::new();
    for result in reader.deserialize() {
        let record: Record = result?;
        records.push(record);
    }
    Ok(records)
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

impl Centerline {
    pub fn from_contour_points(contour_points: Vec<ContourPoint>) -> Self {
        let mut points: Vec<CenterlinePoint> = Vec::with_capacity(contour_points.len());

        // Calculate normals for all but the last point.
        for i in 0..contour_points.len() {
            let current = &contour_points[i];
            let normal = if i < contour_points.len() - 1 {
                let next = &contour_points[i + 1];
                Vector3::new(next.x - current.x, next.y - current.y, next.z - current.z).normalize()
            } else if !contour_points.is_empty() {
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

    /// Retrieves a centerline point by matching frame index.
    pub fn get_by_frame(&self, frame_index: u32) -> Option<&CenterlinePoint> {
        self.points
            .iter()
            .find(|p| p.contour_point.frame_index == frame_index)
    }
}

pub fn read_centerline_txt(path: &str) -> Result<Vec<ContourPoint>, Box<dyn Error>> {
    let file = File::open(path)?;
    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(false)
        .delimiter(b' ')
        .from_reader(file);

    let mut points = Vec::new();
    for (i, result) in rdr.records().enumerate() {
        match result {
            Ok(record) => {
                let mut iter = record.iter();
                let x = iter.next().unwrap().parse::<f64>()?;
                let y = iter.next().unwrap().parse::<f64>()?;
                let z = iter.next().unwrap().parse::<f64>()?;
                points.push(ContourPoint {
                    frame_index: i as u32,
                    point_index: i as u32,
                    x,
                    y,
                    z,
                    aortic: false,
                });
            }
            Err(e) => eprintln!("Skipping invalid row: {:?}", e),
        }
    }
    Ok(points)
}
