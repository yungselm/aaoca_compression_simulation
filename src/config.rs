use serde::Deserialize;
use std::fs;

#[derive(Debug, Deserialize)]
pub struct GeneralConfig {
    pub rest_input_path: String,
    pub stress_input_path: String,
    pub rest_output_path: String,
    pub stress_output_path: String,
    pub diastole_comparison_path: String,
    pub systole_comparison_path: String,
    pub aligned_rest_path: String,
    pub aligned_stress_path: String,
}

#[derive(Debug, Deserialize)]
pub struct ProcessingConfig {
    pub run_process_case: bool,
    pub run_centerline_alignment: bool,
}

#[derive(Debug, Deserialize)]
pub struct SettingsConfig {
    pub interpolation_steps: usize,
    pub x_coord_ref: f64,
    pub y_coord_ref: f64,
    pub z_coord_ref: f64,
    pub x_coord_upper: f64,
    pub y_coord_upper: f64,
    pub z_coord_upper: f64,
    pub x_coord_lower: f64,
    pub y_coord_lower: f64,
    pub z_coord_lower: f64,
    pub steps_best_rotation: usize,
    pub range_rotation_rad: f64,
}

#[derive(Debug, Deserialize)]
pub struct AppConfig {
    pub general: GeneralConfig,
    pub processing: ProcessingConfig,
    pub settings: SettingsConfig,
}

pub fn load_config(path: &str) -> anyhow::Result<AppConfig> {
    let config_str = fs::read_to_string(path)?;
    let config: AppConfig = toml::from_str(&config_str)?;
    Ok(config)
}
