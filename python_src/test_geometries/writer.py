import os
import math
import pandas as pd
import numpy as np

from typing import Tuple, List
from geometry import Contour, GeometryAssembler

from config import (
    FRAME_RATE, PULLBACK_SPEED,
    IDX_DIA_REST_SORTED, IDX_SYS_REST_SORTED,
    ELLIP_DIA_REST, ELLIP_SYS_REST,
    AREA_DIA_REST, AREA_SYS_REST,
    AORTIC_DIA_REST, AORTIC_SYS_REST,
    PULMONARY_DIA_REST, PULMONARY_SYS_REST,
    TRANSLATION_DIA_REST, TRANSLATION_SYS_REST,
    ROTATION_DIA_REST, ROTATION_SYS_REST,
    IDX_DIA_STRESS_SORTED, IDX_SYS_STRESS_SORTED,
    ELLIP_DIA_STRESS, ELLIP_SYS_STRESS,
    AREA_DIA_STRESS, AREA_SYS_STRESS,
    AORTIC_DIA_STRESS, AORTIC_SYS_STRESS,
    PULMONARY_DIA_STRESS, PULMONARY_SYS_STRESS,
    TRANSLATION_DIA_STRESS, TRANSLATION_SYS_STRESS,
    ROTATION_DIA_STRESS, ROTATION_SYS_STRESS,
)

def write_contours(geom: List[Contour],
                   out_dir: str        = "test_geometries/output",
                   out_name: str       = "contours.csv",
                   sep: str            = "\t",
                   include_header: bool= False) -> None:
    """
    Write a list of Contour objects to a tab-delimited CSV with columns:
      idx    x    y    z

    Parameters
    ----------
    geom : List[Contour]
        Your contours (with .idx, .points_x, .points_y, .points_z).
    out_dir : str
        Directory where the file will be written.
    out_name : str
        Name of the CSV file.
    sep : str
        Field separator (default is tab).
    include_header : bool
        Whether to write a header row (default is False).

    Example row:
      126\t4.055600479244161\t1.751602097089961\t0.0
    """
    # make sure output dir exists
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, out_name)

    # collect all rows
    rows = []
    for c in geom:
        print(c.idx)
        # ensure z is iterable
        if isinstance(c.points_z, (list, np.ndarray)):
            zs = c.points_z
        else:
            zs = [c.points_z] * len(c.points_x)

        for x, y, z in zip(c.points_x, c.points_y, zs):
            rows.append((c.idx, x, y, z))

    # turn into DataFrame
    df = pd.DataFrame(rows, columns=["idx", "x", "y", "z"])
    # sort by idx
    df = df.sort_values("idx")

    # write to disk
    df.to_csv(out_path,
              sep=sep,
              header=include_header,
              index=False,
              float_format="%.12g")  # tweak float formatting if you like

    print(f"Wrote {len(df)} rows to {out_path!r}")


def write_reference_point(
    ref_pt: Tuple[float, float, float, float],
    out_dir: str        = "test_geometries/output",
    out_name: str       = "reference_point.csv",
    sep: str            = "\t",
    include_header: bool= False
) -> None:
    """
    Writes a single reference point (idx, x, y, z) to a CSV.

    Parameters
    ----------
    ref_pt : tuple of (idx, x, y, z)
        e.g. (803, 3.5506694412222464, 2.595584306716262, 22.77458)
    out_dir : str
        Directory to write into.
    out_name : str
        Name of the CSV file.
    sep : str
        Field separator (default is tab).
    include_header : bool
        Whether to write a header row (default False).

    Produces a file containing, for example:

      803\t3.5506694412222464\t2.595584306716262\t22.77458
    """
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, out_name)

    print(ref_pt)

    # if the user wants a header, write it first
    with open(path, "w", newline="") as f:
        if include_header:
            f.write(sep.join(("idx", "x", "y", "z")) + "\n")
        # format floats with full precision but strip trailing zeros
        line = sep.join(f"{v:.12g}" if isinstance(v, float) else str(v)
                        for v in ref_pt)
        f.write(line + "\n")

    print(f"Wrote reference point to {path!r}")


def write_full_csv(
    out_dir: str        = "test_geometries/output",
    out_name: str       = "full_measurements.csv",
    sep: str            = ",",
    include_header: bool= True,
    mode: str           = 'rest'
) -> None:
    """
    Build a combined DIA+SYS table with all requested columns.
    """
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, out_name)

    asm = GeometryAssembler(mode=mode)
    z_coords_dia, z_coords_sys = asm.calculate_z_coords()

    rows = []

    def process_phase(
        idx_list:        List[int],
        ellips:          List[float],
        areas:           List[float],
        aortic:          List[float],
        pulmonary:       List[float],
        trans:           List[tuple],
        rots:            List[float],
        phase_label:     str
    ):
        z_coords = z_coords_dia if phase_label == "D" else z_coords_sys
        start_frame = idx_list[0]

        for pos, frame in enumerate(idx_list):
            a_ratio = ellips[pos]
            area   = areas[pos]
            # recover a,b from area = π a b
            # and elliptic_ratio = b/a  ⇒  b = a * ratio
            a = math.sqrt(area / (math.pi * a_ratio))
            b = a * a_ratio

            # circumference approximation
            C = 2 * math.pi * math.sqrt((a*a + b*b) / 2)

            dx, dy = trans[pos]
            vlen = math.hypot(dx, dy)
            vang = rots[pos]

            rows.append({
                "frame": frame,
                "position": float(z_coords[pos]),
                "phase": phase_label,
                "lumen_area": area,
                "lumen_circumf": C,
                "longest_distance": 2*a,
                "shortest_distance": 2*b,
                "elliptic_ratio": a_ratio,
                "vector_length": vlen,
                "vector_angle": vang,
                "measurement_1": aortic[pos] if aortic[pos] is not None else "",
                "measurement_2": pulmonary[pos] if pulmonary[pos] is not None else "",
                "pullback_speed": PULLBACK_SPEED,
                "pullback_start_frame": start_frame,
                "frame_rate": FRAME_RATE
            })

    if mode == 'rest':
        # DIA
        process_phase(
            IDX_DIA_REST_SORTED,
            ELLIP_DIA_REST,
            AREA_DIA_REST,
            AORTIC_DIA_REST,
            PULMONARY_DIA_REST,
            TRANSLATION_DIA_REST,
            ROTATION_DIA_REST,
            phase_label="D"
        )

        # SYS
        process_phase(
            IDX_SYS_REST_SORTED,
            ELLIP_SYS_REST,
            AREA_SYS_REST,
            AORTIC_SYS_REST,
            PULMONARY_SYS_REST,
            TRANSLATION_SYS_REST,
            ROTATION_SYS_REST,
            phase_label="S"
        )
    elif mode == 'stress':
        # DIA
        process_phase(
            IDX_DIA_STRESS_SORTED,
            ELLIP_DIA_STRESS,
            AREA_DIA_STRESS,
            AORTIC_DIA_STRESS,
            PULMONARY_DIA_STRESS,
            TRANSLATION_DIA_STRESS,
            ROTATION_DIA_STRESS,
            phase_label="D"
        )

        # SYS
        process_phase(
            IDX_SYS_STRESS_SORTED,
            ELLIP_SYS_STRESS,
            AREA_SYS_STRESS,
            AORTIC_SYS_STRESS,
            PULMONARY_SYS_STRESS,
            TRANSLATION_SYS_STRESS,
            ROTATION_SYS_STRESS,
            phase_label="S"
        )
    else:
        print("Either mode == 'rest' or mode == 'stress'")

    # dump to DataFrame
    df = pd.DataFrame(rows)
    df = df.sort_values(["phase", "frame"])

    df.to_csv(
        out_path,
        sep=sep,
        index=False,
        header=include_header,
        float_format="%.12g"
    )
    print(f"Wrote {len(df)} rows to {out_path!r}")