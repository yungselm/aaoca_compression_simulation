mod io;
mod contour;
mod processing;
mod utils;
mod texture;
mod comparison;
mod centerline_alignment;
mod config;

use std::error::Error;

use processing::process_case;
use comparison::process_phase_comparison;
use centerline_alignment::create_centerline_aligned_meshes;
use config::load_config;

fn main() -> Result<(), Box<dyn Error>> {
    let config = load_config("config.toml")?;

    // Run process_case if enabled.
    if config.processing.run_process_case {
        process_case("rest", 
        &config.general.rest_input_path, 
        &config.general.rest_output_path, 
        config.settings.interpolation_steps)?;
        process_case("stress", 
        &config.general.stress_input_path, 
        &config.general.stress_output_path, 
        config.settings.interpolation_steps)?;
    }

    // Run phase comparison if enabled.
    if config.processing.run_phase_comparison {
        process_phase_comparison(
            "diastolic",
            &config.general.rest_input_path,
            &config.general.stress_input_path,
            &config.general.diastole_comparison_path,
            config.settings.interpolation_steps,
        )?;
        process_phase_comparison(
            "systolic",
            &config.general.rest_input_path,
            &config.general.stress_input_path,
            &config.general.systole_comparison_path,
            config.settings.interpolation_steps,
        )?;
    }

    // Run centerline alignment if enabled.
    if config.processing.run_centerline_alignment {
        create_centerline_aligned_meshes(
            "rest",
            "resampled_centerline.txt",
            &config.general.rest_output_path,
            &config.general.aligned_rest_path,
            config.settings.interpolation_steps,
        )?;
        create_centerline_aligned_meshes(
            "stress",
            "resampled_centerline_stress.txt",
            &config.general.stress_output_path,
            &config.general.aligned_stress_path,
            config.settings.interpolation_steps,
        )?;
    }
    Ok(())
}
