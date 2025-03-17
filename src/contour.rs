use crate::io::ContourPoint;
use crate::processing::align_contours;
use std::collections::HashMap;

/// Creates contours from a list of contour points.
///
/// Groups points by their frame index, aligns them, and re-indexes.
pub fn create_contours(points: Vec<ContourPoint>) -> Vec<(u32, Vec<ContourPoint>)> {
    // Group points by frame index.
    let mut groups: HashMap<u32, Vec<ContourPoint>> = HashMap::new();
    for p in points {
        groups.entry(p.frame_index).or_default().push(p);
    }

    let mut contours: Vec<_> = groups.into_iter().collect();

    // Align and sort contours using processing routines.
    contours = align_contours(contours);
    // Optionally, sort them in descending order by frame (or re-index).
    contours.sort_by_key(|(frame, _)| std::cmp::Reverse(*frame));
    contours
        .into_iter()
        .enumerate()
        .map(|(i, (_, points))| (i as u32, points))
        .collect()
}