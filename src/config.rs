use serde::Deserialize;
use std::fs;
use std::error::Error;

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
    pub run_phase_comparison: bool,
    pub run_centerline_alignment: bool,
}

#[derive(Debug, Deserialize)]
pub struct SettingsConfig {
    pub interpolation_steps: usize,
}

#[derive(Debug, Deserialize)]
pub struct AppConfig {
    pub general: GeneralConfig,
    pub processing: ProcessingConfig,
    pub settings: SettingsConfig,
}

pub fn load_config(path: &str) -> Result<AppConfig, Box<dyn Error>> {
    let config_str = fs::read_to_string(path)?;
    let config: AppConfig = toml::from_str(&config_str)?;
    Ok(config)
}