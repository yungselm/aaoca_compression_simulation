use anyhow::{Result, Context, anyhow};
use crossbeam::thread;

mod config;
mod io;
mod mesh_to_centerline;
mod processing;
mod texture;
mod utils;

use config::load_config;
use processing::comparison::prepare_geometries_comparison;
use processing::process_case::{create_geometry_pair, process_case};
use mesh_to_centerline::create_centerline_aligned_meshes;

fn main() -> Result<()> {
    // 1. Load configuration
    let config = load_config("config.toml")
        .context("Failed to load configuration from config.toml")?;

    // 2. Parallel rest / stress processing
    if config.processing.run_process_case {
        // Chain directly on thread::scope without a stray semicolon
        let _result = thread::scope(|s| -> Result<()> {
            // REST thread
            let rest_handle = s.spawn(|_| -> Result<_> {
                let geom_rest = create_geometry_pair(
                    "rest".to_string(),
                    &config.general.rest_input_path,
                    config.settings.steps_best_rotation,
                    config.settings.range_rotation_rad,
                )
                .context("create_geometry_pair(rest) failed")?;
                
                let processed_rest = process_case(
                    "rest",
                    geom_rest,
                    &config.general.rest_output_path,
                    config.settings.interpolation_steps,
                )
                .context("process_case(rest) failed")?;
                
                Ok(processed_rest)
            });

            // STRESS thread
            let stress_handle = s.spawn(|_| -> Result<_> {
                let geom_stress = create_geometry_pair(
                    "stress".to_string(),
                    &config.general.stress_input_path,
                    config.settings.steps_best_rotation,
                    config.settings.range_rotation_rad,
                )
                .context("create_geometry_pair(stress) failed")?;
                
                let processed_stress = process_case(
                    "stress",
                    geom_stress,
                    &config.general.stress_output_path,
                    config.settings.interpolation_steps,
                )
                .context("process_case(stress) failed")?;
                
                Ok(processed_stress)
            });

            // Join threads & propagate any processing errors
            let geometries_rest   = rest_handle.join().unwrap()?;
            let geometries_stress = stress_handle.join().unwrap()?;

            // Diastole / systole comparison
            let (dia_geom_pair, sys_geom_pair) =
                prepare_geometries_comparison(geometries_rest, geometries_stress);

            process_case(
                "diastolic",
                dia_geom_pair,
                &config.general.diastole_comparison_path,
                config.settings.interpolation_steps,
            )
            .context("process_case(diastolic) failed")?;

            process_case(
                "systolic",
                sys_geom_pair,
                &config.general.systole_comparison_path,
                config.settings.interpolation_steps,
            )
            .context("process_case(systolic) failed")?;

            Ok(())
        })
        .map_err(|panic_payload| {
            anyhow!("Parallel processing threads panicked: {:?}", panic_payload)
        })?;  
    }

    // 3. Optional centerline alignment
    if config.processing.run_centerline_alignment {
        create_centerline_aligned_meshes(
            "rest",
            "input/rest_csv_files/resampled_centerline.txt",
            &config.general.rest_output_path,
            &config.general.aligned_rest_path,
            config.settings.interpolation_steps,
            config.settings.x_coord_ref,
            config.settings.y_coord_ref,
            config.settings.z_coord_ref,
            config.settings.x_coord_upper,
            config.settings.y_coord_upper,
            config.settings.z_coord_upper,
            config.settings.x_coord_lower,
            config.settings.y_coord_lower,
            config.settings.z_coord_lower,
        )
        .context("create_centerline_aligned_meshes(rest) failed")?;

        create_centerline_aligned_meshes(
            "stress",
            "input/stress_csv_files/resampled_centerline_stress.txt",
            &config.general.stress_output_path,
            &config.general.aligned_stress_path,
            config.settings.interpolation_steps,
            config.settings.x_coord_ref,
            config.settings.y_coord_ref,
            config.settings.z_coord_ref,
            config.settings.x_coord_upper,
            config.settings.y_coord_upper,
            config.settings.z_coord_upper,
            config.settings.x_coord_lower,
            config.settings.y_coord_lower,
            config.settings.z_coord_lower,
        )
        .context("create_centerline_aligned_meshes(stress) failed")?;
    }

    Ok(())
}

