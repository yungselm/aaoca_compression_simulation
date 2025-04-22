mod config;
mod io;
mod mesh_to_centerline;
mod processing;
mod texture;
mod utils;

use std::error::Error;

// use mesh_to_centerline::create_centerline_aligned_meshes;
use config::load_config;
use processing::comparison::prepare_geometries_comparison;
use processing::process_case::{create_geometry_pair, process_case};

fn main() -> Result<(), Box<dyn Error>> {
    let config = load_config("config.toml")?;

    // Run process_case if enabled.
    if config.processing.run_process_case {
        let geometries_rest = create_geometry_pair(
            "rest".to_string(),
            &config.general.rest_input_path,
            config.settings.steps_best_rotation,
            config.settings.range_rotation_rad,
        )?;
        let geometries_rest = process_case(
            "rest",
            geometries_rest,
            &config.general.rest_output_path,
            config.settings.interpolation_steps,
        )?;
        let geometries_stress = create_geometry_pair(
            "stress".to_string(),
            &config.general.stress_input_path,
            config.settings.steps_best_rotation,
            config.settings.range_rotation_rad,
        )?;
        let geometries_stress = process_case(
            "stress",
            geometries_stress,
            &config.general.stress_output_path,
            config.settings.interpolation_steps,
        )?;

        let (dia_geom_pair, sys_geom_pair) =
            prepare_geometries_comparison(geometries_rest, geometries_stress);
        process_case(
            "diastolic",
            dia_geom_pair,
            &config.general.diastole_comparison_path,
            config.settings.interpolation_steps,
        )?;
        process_case(
            "systolic",
            sys_geom_pair,
            &config.general.systole_comparison_path,
            config.settings.interpolation_steps,
        )?;
    }
    Ok(())

    // // Run centerline alignment if enabled.
    // if config.processing.run_centerline_alignment {
    //     create_centerline_aligned_meshes(
    //         "rest",
    //         "resampled_centerline.txt",
    //         &config.general.rest_output_path,
    //         &config.general.aligned_rest_path,
    //         config.settings.interpolation_steps,
    //         config.settings.x_coord_ref,
    //         config.settings.y_coord_ref,
    //         config.settings.z_coord_ref,
    //         config.settings.x_coord_upper,
    //         config.settings.y_coord_upper,
    //         config.settings.z_coord_upper,
    //         config.settings.x_coord_lower,
    //         config.settings.y_coord_lower,
    //         config.settings.z_coord_lower,
    //     )?;
    //     create_centerline_aligned_meshes(
    //         "stress",
    //         "resampled_centerline_stress.txt",
    //         &config.general.stress_output_path,
    //         &config.general.aligned_stress_path,
    //         config.settings.interpolation_steps,
    //         config.settings.x_coord_ref,
    //         config.settings.y_coord_ref,
    //         config.settings.z_coord_ref,
    //         config.settings.x_coord_upper,
    //         config.settings.y_coord_upper,
    //         config.settings.z_coord_upper,
    //         config.settings.x_coord_lower,
    //         config.settings.y_coord_lower,
    //         config.settings.z_coord_lower,
    //     )?;
    // }
    // Ok(())
}
