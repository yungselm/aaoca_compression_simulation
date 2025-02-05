/// Trims two vectors so that both have the same length.
/// Returns the new (minimum) length.
pub fn trim_to_same_length<T>(v1: &mut Vec<T>, v2: &mut Vec<T>) -> usize {
    let min_len = std::cmp::min(v1.len(), v2.len());
    v1.truncate(min_len);
    v2.truncate(min_len);
    min_len
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
