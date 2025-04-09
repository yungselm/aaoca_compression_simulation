mod centerline_alignment;
mod comparison;
mod config;
mod contour;
mod io;
mod processing;
mod texture;
mod utils;

use std::error::Error;

use centerline_alignment::create_centerline_aligned_meshes;
use comparison::process_phase_comparison;
use config::load_config;
use processing::process_case;

fn main() -> Result<(), Box<dyn Error>> {
    let config = load_config("config.toml")?;

    // Run process_case if enabled.
    if config.processing.run_process_case {
        process_case(
            "rest",
            &config.general.rest_input_path,
            &config.general.rest_output_path,
            config.settings.interpolation_steps,
            config.settings.steps_best_rotation,
            config.settings.range_rotation_rad,
        )?;
        process_case(
            "stress",
            &config.general.stress_input_path,
            &config.general.stress_output_path,
            config.settings.interpolation_steps,
            config.settings.steps_best_rotation,
            config.settings.range_rotation_rad,
        )?;
    }

    // Run phase comparison if enabled.
    if config.processing.run_phase_comparison {
        process_phase_comparison(
            "diastolic",
            &config.general.rest_output_path,
            &config.general.stress_output_path,
            &config.general.diastole_comparison_path,
            config.settings.interpolation_steps,
            config.settings.steps_best_rotation,
            config.settings.range_rotation_rad,
        )?;
        process_phase_comparison(
            "systolic",
            &config.general.rest_output_path,
            &config.general.stress_output_path,
            &config.general.systole_comparison_path,
            config.settings.interpolation_steps,
            config.settings.steps_best_rotation,
            config.settings.range_rotation_rad,
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
            config.settings.x_coord_ref,
            config.settings.y_coord_ref,
            config.settings.z_coord_ref,
        )?;
        create_centerline_aligned_meshes(
            "stress",
            "resampled_centerline_stress.txt",
            &config.general.stress_output_path,
            &config.general.aligned_stress_path,
            config.settings.interpolation_steps,
            config.settings.x_coord_ref,
            config.settings.y_coord_ref,
            config.settings.z_coord_ref,
        )?;
    }
    Ok(())
}
