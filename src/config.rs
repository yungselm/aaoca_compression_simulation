use serde::Deserialize;
use std::fs;
use std::error::Error;

#[derive(Debug, Deserialize)]
pub struct GeneralConfig {
    pub data_dir: String,
    pub output_dir: String,
}

#[derive(Debug, Deserialize)]
pub struct ProcessingConfig {
    pub run_process_case: bool,
    pub run_phase_comparison: bool,
    pub run_centerline_alignment: bool,
    pub cases: Vec<String>,
}

#[derive(Debug, Deserialize)]
pub struct AppConfig {
    pub general: GeneralConfig,
    pub processing: ProcessingConfig,
}

pub fn load_config() -> Result<AppConfig, Box<dyn Error>> {
    let config_str = fs::read_to_string("config.toml")?;
    let config: AppConfig = toml::from_str(&config_str)?;
    Ok(config)
}