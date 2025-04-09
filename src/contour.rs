use rayon::prelude::*;
use std::collections::HashMap;
use std::error::Error;
use std::f64::consts::PI;

use crate::io::ContourPoint;

#[allow(dead_code)] // used in processing.rs, but warning unused here otherwise
pub struct Contour {
    pub id: u32,
    pub points: Vec<ContourPoint>,
    pub centroid: (f64, f64, f64),
}

impl Contour {
    /// Creates contours from a list of contour points.
    ///
    /// Groups points by their frame index, aligns them, and re-indexes.
    pub fn create_contours(
        points: Vec<ContourPoint>,
        steps: usize,
        range: f64,
        reference_point: &ContourPoint, // needs to be set on aortic wall ostium!
    ) -> Vec<(u32, Vec<ContourPoint>)> {
        // Group points by frame index.
        let mut groups: HashMap<u32, Vec<ContourPoint>> = HashMap::new();
        for p in points {
            groups.entry(p.frame_index).or_default().push(p);
        }

        let mut contours: Vec<_> = groups.into_iter().collect();

        // Align and sort contours using processing routines.
        contours = Self::align_contours(contours, steps, range, reference_point);
        
        // Optionally, sort them in descending order by frame (or re-index).
        // This ensures that ostium has index 0, because usually last frame in pullback
        contours.sort_by_key(|(frame, _)| std::cmp::Reverse(*frame));

        contours
            .into_iter()
            .enumerate()
            .map(|(i, (_, points))| (i as u32, points))
            .collect()
    }

    /// Aligns contours by rotating and translating them so that a designated
    /// reference contour (e.g., with the highest frame index) is used as the basis.
    fn align_contours(
        mut contours: Vec<(u32, Vec<ContourPoint>)>,
        steps: usize,
        range: f64,
        reference_point: &ContourPoint,
    ) -> Vec<(u32, Vec<ContourPoint>)> {
        // Sort contours by frame index.
        contours.sort_by_key(|(frame_index, _)| *frame_index);

        // Use the contour with the highest frame index as reference.
        let reference_index = contours.iter().map(|(id, _)| *id).max().unwrap();
        let reference_pos = contours
            .iter()
            .position(|(id, _)| *id == reference_index)
            .expect("Reference contour not found");
        let (ref_frame, ref_contour) = &contours[reference_pos];
        println!("Using contour {} as reference.", ref_frame);

        // Rotate the reference contour so that its farthest-points line is vertical.
        let ((p1, p2), _dist) = Self::find_farthest_points(ref_contour);
        let dx = p2.x - p1.x;
        let dy = p2.y - p1.y;
        let line_angle = dy.atan2(dx);
        let mut rotation_to_y = (PI / 2.0) - line_angle;
        println!(
            "Reference line angle: {:.3} rad; rotating reference by {:.3} rad",
            line_angle, rotation_to_y
        );
        let ref_centroid = Self::compute_centroid(ref_contour);
        let mut aligned_reference = ref_contour.clone();
        Self::rotate_contour(&mut aligned_reference, rotation_to_y, ref_centroid);
        Self::sort_contour_points(&mut aligned_reference);
        let mut aligned_reference_point = reference_point.clone();
        Self::rotate_point(&mut aligned_reference_point, rotation_to_y, (5.0, 5.0));

        let n = aligned_reference.len() / 2;

        // --- Begin aortic determination ---
        // We assume there are exactly 500 points; split into two halves: indices 0..249 and 250..499.
        // Compute the cumulative distance for each half from the rotated reference.
        let first_half_distance: f64 = aligned_reference
            .iter()
            .take(n)
            .map(|pt| pt.distance_to(&aligned_reference_point))
            .sum();
        let second_half_distance: f64 = aligned_reference
            .iter()
            .skip(n)
            .map(|pt| pt.distance_to(&aligned_reference_point))
            .sum();

        // Debug output: print average distances for each half.
        println!(
            "First half avg distance: {:.3}, Second half avg distance: {:.3}",
            first_half_distance / (n as f64),
            second_half_distance / (n as f64)
        );

        // Mark the half that is closer as 'aortic' (i.e. set to true).
        if first_half_distance < second_half_distance {
            for pt in aligned_reference.iter_mut().take(n) {
                pt.aortic = true;
            }
            println!("First half (points 0 to 249) marked as aortic.");
            if first_half_distance < 5.0 {
                rotation_to_y += PI;
                Self::rotate_contour(&mut aligned_reference, PI, ref_centroid);
            }
        } else {
            for pt in aligned_reference.iter_mut().skip(n) {
                pt.aortic = true;
            }
            println!("Second half (points 250 to 499) marked as aortic.");
            if second_half_distance < 5.0 {
                rotation_to_y += PI;
                Self::rotate_contour(&mut aligned_reference, PI, ref_centroid);
            }
        }

        contours[reference_pos].1 = aligned_reference.clone();

        let reference_clone = contours[reference_pos].1.clone();

        // Align every non-reference contour.
        for (id, contour) in contours.iter_mut() {
            if *id == reference_index {
                continue;
            }
            let orig_centroid = Self::compute_centroid(contour);
            Self::rotate_contour(contour, rotation_to_y, orig_centroid);
            // Center the contour and the reference.
            let centered: Vec<ContourPoint> = contour
                .iter()
                .map(|p| ContourPoint {
                    frame_index: p.frame_index,
                    x: p.x - orig_centroid.0,
                    y: p.y - orig_centroid.1,
                    z: p.z,
                    aortic: p.aortic,
                })
                .collect();
            let centered_reference: Vec<ContourPoint> = reference_clone
                .iter()
                .map(|p| ContourPoint {
                    frame_index: p.frame_index,
                    x: p.x - ref_centroid.0,
                    y: p.y - ref_centroid.1,
                    z: p.z,
                    aortic: p.aortic,
                })
                .collect();

            // Translate contour to match reference.
            let translation = (
                ref_centroid.0 - orig_centroid.0,
                ref_centroid.1 - orig_centroid.1,
            );
            Self::translate_contour(contour, translation);
            let best_rot = Self::find_best_rotation(&centered_reference, &centered, steps, range);
            println!(
                "Rotating contour {} by {:.3} rad for best correlation",
                id, best_rot
            );
            Self::rotate_contour(contour, best_rot, ref_centroid);
            Self::sort_contour_points(contour);
            
            // set first half to aortic false and second to true
            for pt in contour.iter_mut().skip(n) {
                pt.aortic = true;
            }
        }
        contours
    }

    /// Interpolates between two aligned sets of contours.
    pub fn interpolate_contours(
        contours_start: &[(u32, Vec<ContourPoint>)],
        contours_end: &[(u32, Vec<ContourPoint>)],
        steps: usize,
    ) -> Result<Vec<Vec<(u32, Vec<ContourPoint>)>>, Box<dyn Error>> {
        use std::cmp::min;

        let n = min(contours_start.len(), contours_end.len());
        let start = &contours_start[0..n];
        let end = &contours_end[0..n];

        let mut interpolated = Vec::with_capacity(steps);
        for step in 0..steps {
            let t = step as f64 / (steps - 1) as f64;
            let mut intermediate = Vec::with_capacity(n);
            for ((id_start, contour_start), (id_end, contour_end)) in start.iter().zip(end.iter()) {
                if id_start != id_end {
                    return Err("Contour IDs do not match between start and end.".into());
                }
                if contour_start.len() != contour_end.len() {
                    return Err("Contour point counts do not match between start and end.".into());
                }
                let interp_contour: Vec<ContourPoint> = contour_start
                    .iter()
                    .zip(contour_end.iter())
                    .map(|(p_start, p_end)| ContourPoint {
                        frame_index: step as u32, // new sequential id for this interpolation step
                        x: p_start.x * (1.0 - t) + p_end.x * t,
                        y: p_start.y * (1.0 - t) + p_end.y * t,
                        z: p_start.z * (1.0 - t) + p_end.z * t,
                        aortic: p_start.aortic,
                    })
                    .collect();
                intermediate.push((*id_start, interp_contour));
            }
            interpolated.push(intermediate);
        }
        Ok(interpolated)
    }

    /// Finds the best rotation angle (in radians) that minimizes the Hausdorff distance
    /// leveraging parallel computation for performance.
    pub fn find_best_rotation(
        reference: &[ContourPoint],
        target: &[ContourPoint],
        steps: usize,
        range: f64,
    ) -> f64 {
        let steps = steps; // Reduce to optimize performance
        let range = range; // +/- 60 degrees
        let increment = (2.0 * range) / (steps as f64);

        (0..=steps)
            .into_par_iter() // Parallel iteration (requires Rayon)
            .map(|i| {
                let angle = -range + (i as f64) * increment;

                // Rotate target contour inline without extra allocation
                let rotated: Vec<ContourPoint> = target
                    .iter()
                    .map(|p| {
                        let x = p.x * angle.cos() - p.y * angle.sin();
                        let y = p.x * angle.sin() + p.y * angle.cos();
                        ContourPoint {
                            frame_index: p.frame_index,
                            x,
                            y,
                            z: p.z,
                            aortic: p.aortic,
                        }
                    })
                    .collect();

                let hausdorff_dist = Self::hausdorff_distance(reference, &rotated);
                (angle, hausdorff_dist)
            })
            .reduce(
                || (std::f64::NEG_INFINITY, std::f64::MAX),
                |(best_a, best_d), (angle, dist)| {
                    if dist < best_d {
                        (angle, dist)
                    } else {
                        (best_a, best_d)
                    }
                },
            )
            .0
    }

    /// Computes the Hausdorff distance between two point sets.
    pub fn hausdorff_distance(set1: &[ContourPoint], set2: &[ContourPoint]) -> f64 {
        let forward = Self::directed_hausdorff(set1, set2);
        let backward = Self::directed_hausdorff(set2, set1);
        forward.max(backward) // Hausdorff distance is max of both directed distances
    }

    /// Computes directed Hausdorff distance from A to B
    fn directed_hausdorff(contour_a: &[ContourPoint], contour_b: &[ContourPoint]) -> f64 {
        contour_a
            .par_iter() // Use parallel iteration
            .map(|pa| {
                contour_b
                    .iter()
                    .map(|pb| {
                        let dx = pa.x - pb.x;
                        let dy = pa.y - pb.y;
                        (dx * dx + dy * dy).sqrt()
                    })
                    .fold(std::f64::MAX, f64::min) // Directly find min without storing a Vec
            })
            .reduce(|| 0.0, f64::max) // Directly find max without extra allocation
    }

    /// Computes the centroid (average x, y) of a contour.
    pub fn compute_centroid(contour: &[ContourPoint]) -> (f64, f64) {
        let (sum_x, sum_y) = contour
            .iter()
            .fold((0.0, 0.0), |(sx, sy), p| (sx + p.x, sy + p.y));
        let n = contour.len() as f64;
        (sum_x / n, sum_y / n)
    }

    /// Sorts contour points in counterclockwise order around the centroid
    /// and rotates so that the highest y-value is first.
    fn sort_contour_points(contour: &mut Vec<ContourPoint>) {
        let center = Self::compute_centroid(contour);
        contour.sort_by(|a, b| {
            let angle_a = (a.y - center.1).atan2(a.x - center.0);
            let angle_b = (b.y - center.1).atan2(b.x - center.0);
            angle_a.partial_cmp(&angle_b).unwrap()
        });
        let start_idx = contour
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.y.partial_cmp(&b.y).unwrap())
            .map(|(i, _)| i)
            .unwrap_or(0);
        contour.rotate_left(start_idx);
    }

    /// Rotates a single point about a center by a given angle (in radians).
    fn rotate_point(p: &ContourPoint, angle: f64, center: (f64, f64)) -> ContourPoint {
        let (cx, cy) = center;
        let x = p.x - cx;
        let y = p.y - cy;
        let cos_a = angle.cos();
        let sin_a = angle.sin();
        ContourPoint {
            frame_index: p.frame_index,
            x: x * cos_a - y * sin_a + cx,
            y: x * sin_a + y * cos_a + cy,
            z: p.z,
            aortic: p.aortic,
        }
    }

    /// Rotates all points in a contour about a center.
    fn rotate_contour(contour: &mut [ContourPoint], angle: f64, center: (f64, f64)) {
        for p in contour.iter_mut() {
            let rotated = Self::rotate_point(p, angle, center);
            p.x = rotated.x;
            p.y = rotated.y;
        }
    }

    /// Translates a contour by a given (dx, dy) offset.
    pub fn translate_contour(contour: &mut [ContourPoint], translation: (f64, f64)) {
        let (dx, dy) = translation;
        for p in contour.iter_mut() {
            p.x += dx;
            p.y += dy;
        }
    }

    /// Finds the pair of farthest points in a contour.
    fn find_farthest_points(contour: &[ContourPoint]) -> ((&ContourPoint, &ContourPoint), f64) {
        let mut max_dist = 0.0;
        let mut farthest_pair = (&contour[0], &contour[0]);
        for i in 0..contour.len() {
            for j in i + 1..contour.len() {
                let dx = contour[i].x - contour[j].x;
                let dy = contour[i].y - contour[j].y;
                let dist = (dx * dx + dy * dy).sqrt();
                if dist > max_dist {
                    max_dist = dist;
                    farthest_pair = (&contour[i], &contour[j]);
                }
            }
        }
        (farthest_pair, max_dist)
    }

    /// Find the closest opposite points
    pub fn find_closest_opposite(
        contour: &[ContourPoint],
    ) -> ((&ContourPoint, &ContourPoint), f64) {
        let mut min_dist = f64::MAX;
        let mut closest_pair = (&contour[0], &contour[0]);

        let n = contour.len();
        for i in 0..n / 2 {
            let j = n - 1 - i; // Opposite index of i
            let dx = contour[i].x - contour[j].x;
            let dy = contour[i].y - contour[j].y;
            let dist = (dx * dx + dy * dy).sqrt(); // Euclidean distance

            if dist < min_dist {
                min_dist = dist;
                closest_pair = (&contour[i], &contour[j]);
            }
        }

        (closest_pair, min_dist)
    }
}
