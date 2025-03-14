
# AAOCA compression simulation
## Description
This code takes the output of the ["AIVUS-CAA"](https://github.com/AI-in-Cardiovascular-Medicine/AIVUS-CAA) and visualizes either pulsatile lumen deformation or stress-induced lumen deformation.

![Dynamic lumen changes](media/dynamic_lumen_changes.png)

<!-- An example for rest pulsatile lumen deformation:

![Phasic Compression](media/phasic_compression.gif)

And with additional uv texture map, depicting the change in distance in red scale:

![Phasic Compression UV](media/uv_map.gif) -->
An example for rest pulsatile lumen deformation with additionally the IVUS catheter and UV-mapping to depict displacements.

![Rest Pulsatile Lumen Deformation](media/animation_pulsatile_lumen_deformation_rest.gif)

An example for stress pulsatile lumen deformation with additionally the IVUS catheter and UV-mapping to depict displacements.

![Stress Pulsatile Lumen Deformation](animation_pulsatile_lumen_deformation_stress.gif)

## Installation and Running
```bash
    cargo build
    cargo run
```
