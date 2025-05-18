import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from geometry import GeometryAssembler

from writer import write_contours, write_reference_point, write_full_csv


def main():
    ref_dia_rest, geom_dia_rest, ref_sys_rest, geom_sys_rest = GeometryAssembler(
        "rest"
    )()

    write_contours(
        geom_dia_rest,
        out_dir="test_geometries/output/rest_csv_files",
        out_name="diastolic_contours.csv",
    )
    write_contours(
        geom_sys_rest,
        out_dir="test_geometries/output/rest_csv_files",
        out_name="systolic_contours.csv",
    )
    write_reference_point(
        ref_dia_rest,
        out_dir="test_geometries/output/rest_csv_files",
        out_name="diastolic_reference_points.csv",
    )
    write_reference_point(
        ref_sys_rest,
        out_dir="test_geometries/output/rest_csv_files",
        out_name="systolic_reference_points.csv",
    )
    write_full_csv(
        out_dir="test_geometries/output/rest_csv_files",
        out_name="combined_sorted_manual.csv",
    )

    ref_dia_stress, geom_dia_stress, ref_sys_stress, geom_sys_stress = (
        GeometryAssembler("stress")()
    )

    write_contours(
        geom_dia_stress,
        out_dir="test_geometries/output/stress_csv_files",
        out_name="diastolic_contours.csv",
    )
    write_contours(
        geom_sys_stress,
        out_dir="test_geometries/output/stress_csv_files",
        out_name="systolic_contours.csv",
    )
    write_reference_point(
        ref_dia_stress,
        out_dir="test_geometries/output/stress_csv_files",
        out_name="diastolic_reference_points.csv",
    )
    write_reference_point(
        ref_sys_stress,
        out_dir="test_geometries/output/stress_csv_files",
        out_name="systolic_reference_points.csv",
    )
    write_full_csv(
        out_dir="test_geometries/output/stress_csv_files",
        out_name="combined_sorted_manual.csv",
    )


if __name__ == "__main__":
    main()
