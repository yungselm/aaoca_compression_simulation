/// Trims two vectors of contours so that both have the same length.
/// Returns the new (minimum) length.
pub fn trim_to_same_length<T>(v1: &mut Vec<T>, v2: &mut Vec<T>) -> usize {
    let min_len = std::cmp::min(v1.len(), v2.len());
    v1.truncate(min_len);
    v2.truncate(min_len);
    min_len
}
