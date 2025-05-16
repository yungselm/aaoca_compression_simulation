import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from config import (
    FRAME_RATE,
    PULLBACK_SPEED,
    START_DIA,
    START_SYS,
    NUM_POINTS_CONTOUR,
    NUM_POINTS_CATHETER,
    IDX_DIA_REST_SORTED,
    IDX_SYS_REST_SORTED,
    ELLIP_DIA_REST,
    ELLIP_SYS_REST,
    AREA_DIA_REST,
    AREA_SYS_REST,
    AORTIC_DIA_REST,
    AORTIC_SYS_REST,
    PULMONARY_DIA_REST,
    PULMONARY_SYS_REST,
)


class Contour:
    def __init__(self):
        self.idx = None
        self.points_x = []
        self.points_y = []
        self.points_z = []
        self.aortic_thickness = None
        self.pulmonary_thickness = None


class GeometryAssembler:
    def __init__(self):
        self.frame_rate = FRAME_RATE
        self.pullback_speed = PULLBACK_SPEED
        self.start_dia = START_DIA
        self.start_sys = START_SYS
        self.num_points_contour = NUM_POINTS_CONTOUR
        self.num_points_catheter = NUM_POINTS_CATHETER
        self.idx_dia_rest_sorted = IDX_DIA_REST_SORTED
        self.idx_sys_rest_sorted = IDX_SYS_REST_SORTED
        self.ellip_dia_rest = ELLIP_DIA_REST
        self.ellip_sys_rest = ELLIP_SYS_REST
        self.area_dia_rest = AREA_DIA_REST
        self.area_sys_rest = AREA_SYS_REST
        self.aortic_dia_rest = AORTIC_DIA_REST
        self.aortic_sys_rest = AORTIC_SYS_REST
        self.pulmonary_dia_rest = PULMONARY_DIA_REST
        self.pulmonary_sys_rest = PULMONARY_SYS_REST
        # initialize empty lists for contours
        self.z_coords_dia, self.z_coords_sys = self.calculate_z_coords()
        self.geom_dia = []  # list to hold Contour instances for DIA
        self.geom_sys = []  # list to hold Contour instances for SYS

    def __call__(self):
        self.create_geoms()
        self._plot_contour()

    def create_geoms(self):
        for (idx, ellip, area, aortic, pulmonary, z) in zip(
                self.idx_dia_rest_sorted,
                self.ellip_dia_rest,
                self.area_dia_rest,
                self.aortic_dia_rest,
                self.pulmonary_dia_rest,
                self.z_coords_dia,
            ):
            new_contour = Contour()
            x, y = self.create_ellipse(area, ellip)

            new_contour.idx = idx
            new_contour.points_x = x
            new_contour.points_y = y
            new_contour.points_z = z  # scalar or list as needed
            new_contour.aortic_thickness = aortic
            new_contour.pulmonary_thickness = pulmonary

            self.geom_dia.append(new_contour)

    def calculate_z_coords(self):
        """
        Returns two NumPy arrays (z_coords_dia, z_coords_sys) giving the cumulative
        pullback distance (in mm) at each contour index. The first entry in each array is 0.0.
        """
        def _compute_z(idxs):
            # frame-to-frame differences, with first diff = 0
            frame_diffs = np.diff(idxs, prepend=idxs[0])
            # convert frame diffs → time (s) → distance (mm)
            dist_mm = frame_diffs / self.frame_rate * self.pullback_speed
            # cumulative sum so z[0] == 0.0
            return np.cumsum(dist_mm)

        z_coords_dia = _compute_z(self.idx_dia_rest_sorted)
        z_coords_sys = _compute_z(self.idx_sys_rest_sorted)
        return z_coords_dia, z_coords_sys

    def create_ellipse(self, area, elliptic_ratio, num_points=501):
        """
        Create an ellipse given the area and elliptic ratio, centered at (4.5, 4.5).
        """
        a = np.sqrt(area / (np.pi * elliptic_ratio))
        b = elliptic_ratio * a
        t = np.linspace(0, 2 * np.pi, num_points)
        x = a * np.cos(t) + 4.5
        y = b * np.sin(t) + 4.5
        return x, y

    def _plot_contour(self):
        # plot the first DIA contour
        if not self.geom_dia:
            return
        contour = self.geom_dia[20]
        plt.scatter(contour.points_x, contour.points_y)
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.xlim(0, 9)
        plt.ylim(0, 9)
        plt.title(f'Contour at Z={contour.points_z:.2f} mm')
        plt.gca().set_aspect('equal', adjustable='box')
        plt.show()
