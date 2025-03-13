use std::collections::HashMap;
use crate::io::ContourPoint;

/// Trims two vectors so that both have the same length.
/// Returns the new (minimum) length.
pub fn trim_to_same_length<T>(v1: &mut Vec<T>, v2: &mut Vec<T>) -> usize {
    let min_len = std::cmp::min(v1.len(), v2.len());
    v1.truncate(min_len);
    v2.truncate(min_len);
    min_len
}

/// Smooths the x and y coordinates of the contours using a 3‐point moving average.
///
/// For each point i in contour j, the new x and y values are computed as:
///     new_x = (prev_contour[i].x + current_contour[i].x + next_contour[i].x) / 3.0
///     new_y = (prev_contour[i].y + current_contour[i].y + next_contour[i].y) / 3.0
/// while the z coordinate remains unchanged (taken from the current contour).
///
/// For the first and last contours, the current contour is used twice to simulate a mirror effect.
pub fn smooth_contours(
    contours: &[(u32, Vec<ContourPoint>)]
) -> Vec<(u32, Vec<ContourPoint>)> {
    let n = contours.len();
    if n == 0 {
        return Vec::new();
    }
    // Create a new vector to hold the smoothed contours.
    let mut smoothed = Vec::with_capacity(n);

    for j in 0..n {
        let (frame, points) = &contours[j];
        let mut new_points = Vec::with_capacity(points.len());
        for i in 0..points.len() {
            // For each point index i, get the x and y from the previous, current, and next contours.
            let (prev_pt, curr_pt, next_pt) = if j == 0 {
                // For the first contour, use the current contour twice.
                (&contours[j].1[i], &contours[j].1[i], &contours[j+1].1[i])
            } else if j == n - 1 {
                // For the last contour, use the current contour twice.
                (&contours[j-1].1[i], &contours[j].1[i], &contours[j].1[i])
            } else {
                (&contours[j-1].1[i], &contours[j].1[i], &contours[j+1].1[i])
            };

            let avg_x = (prev_pt.x + curr_pt.x + next_pt.x) / 3.0;
            let avg_y = (prev_pt.y + curr_pt.y + next_pt.y) / 3.0;
            // Leave z unchanged (from the current point).
            let smoothed_point = ContourPoint {
                frame_index: *frame,
                x: avg_x,
                y: avg_y,
                z: curr_pt.z,
            };
            new_points.push(smoothed_point);
        }
        smoothed.push((*frame, new_points));
    }
    smoothed
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_trim_to_same_length() {
        let mut v1 = vec![1, 2, 3, 4, 5];
        let mut v2 = vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
        assert_eq!(trim_to_same_length(&mut v1, &mut v2), 5);
        assert_eq!(v1, vec![1, 2, 3, 4, 5]);
        assert_eq!(v2, vec![1, 2, 3, 4, 5]);
    }
}
