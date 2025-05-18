import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from typing import Tuple, List

from config import (
    FRAME_RATE,
    PULLBACK_SPEED,
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
    TRANSLATION_DIA_REST,
    TRANSLATION_SYS_REST,
    ROTATION_DIA_REST,
    ROTATION_SYS_REST,
    IDX_DIA_STRESS_SORTED,
    IDX_SYS_STRESS_SORTED,
    ELLIP_DIA_STRESS,
    ELLIP_SYS_STRESS,
    AREA_DIA_STRESS,
    AREA_SYS_STRESS,
    AORTIC_DIA_STRESS,
    AORTIC_SYS_STRESS,
    PULMONARY_DIA_STRESS,
    PULMONARY_SYS_STRESS,
    TRANSLATION_DIA_STRESS,
    TRANSLATION_SYS_STRESS,
    ROTATION_DIA_STRESS,
    ROTATION_SYS_STRESS,
)


class Contour:
    def __init__(self):
        self.idx = None
        self.points_x = []
        self.points_y = []
        self.points_z = []
        self.aortic_thickness = None
        self.pulmonary_thickness = None
        self.centroid = None

    def calculate_centroid(self):
        if len(self.points_x) == 0 or len(self.points_y) == 0:
            self.centroid = None
            return
        x = np.array(self.points_x)
        y = np.array(self.points_y)
        # z may be scalar or array
        z_arr = (np.array(self.points_z)
                 if isinstance(self.points_z, (list, np.ndarray)) else
                 np.full_like(x, self.points_z))
        self.centroid = (x.mean(), y.mean(), z_arr.mean())


class GeometryAssembler:
    def __init__(self, mode='rest'):
        self.frame_rate = FRAME_RATE
        self.pullback_speed = PULLBACK_SPEED
        self.num_points_contour = NUM_POINTS_CONTOUR
        self.num_points_catheter = NUM_POINTS_CATHETER
        if mode == 'rest':
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
            self.translations_dia_rest = TRANSLATION_DIA_REST
            self.translations_sys_rest = TRANSLATION_SYS_REST
            self.rotations_dia_rest = ROTATION_DIA_REST
            self.rotations_sys_rest = ROTATION_SYS_REST
        if mode == 'stress':
            self.idx_dia_rest_sorted = IDX_DIA_STRESS_SORTED
            self.idx_sys_rest_sorted = IDX_SYS_STRESS_SORTED
            self.ellip_dia_rest = ELLIP_DIA_STRESS
            self.ellip_sys_rest = ELLIP_SYS_STRESS
            self.area_dia_rest = AREA_DIA_STRESS
            self.area_sys_rest = AREA_SYS_STRESS
            self.aortic_dia_rest = AORTIC_DIA_STRESS
            self.aortic_sys_rest = AORTIC_SYS_STRESS
            self.pulmonary_dia_rest = PULMONARY_DIA_STRESS
            self.pulmonary_sys_rest = PULMONARY_SYS_STRESS
            self.translations_dia_rest = TRANSLATION_DIA_STRESS
            self.translations_sys_rest = TRANSLATION_SYS_STRESS
            self.rotations_dia_rest = ROTATION_DIA_STRESS
            self.rotations_sys_rest = ROTATION_SYS_STRESS
        # initialize empty lists for contours
        self.z_coords_dia, self.z_coords_sys = self.calculate_z_coords()
        self.geom_dia: List[Contour] = []  # list to hold Contour instances for DIA
        self.geom_sys: List[Contour] = []  # list to hold Contour instances for SYS
        self.reference_dia = None
        self.reference_sys = None

    def __call__(self):
        self.create_geoms()
        self.create_reference_points()
        self._plot_contour()
        return self.reference_dia, self.geom_dia, self.reference_sys, self.geom_sys

    def create_geoms(self):
        """
        Build both DIA and SYS contour lists in one loop.
        """
        # clear any existing contours
        self.geom_dia.clear()
        self.geom_sys.clear()

        # zip all DIA params together, and all SYS params together
        dia_iter = zip(
            self.idx_dia_rest_sorted,
            self.ellip_dia_rest,
            self.area_dia_rest,
            self.aortic_dia_rest,
            self.pulmonary_dia_rest,
            self.z_coords_dia,
            self.translations_dia_rest,
            self.rotations_dia_rest
        )
        sys_iter = zip(
            self.idx_sys_rest_sorted,
            self.ellip_sys_rest,
            self.area_sys_rest,
            self.aortic_sys_rest,
            self.pulmonary_sys_rest,
            self.z_coords_sys,
            self.translations_sys_rest,
            self.rotations_sys_rest
        )

        for dia_params, sys_params in zip(dia_iter, sys_iter):
            # build one DIA contour
            contour_dia = self._make_contour(*dia_params)
            self.geom_dia.append(contour_dia)

            # build one SYS contour
            contour_sys = self._make_contour(*sys_params)
            self.geom_sys.append(contour_sys)

    def create_reference_points(self):
        """
        Compute reference points for DIA and SYS by translating and then rotating
        the base point (7.5, 4.5) around each contour's centroid.
        """
        idx_dia = len(self.geom_dia) - 1
        idx_sys = len(self.geom_sys) - 1
        base = np.array([6.5, 4.5])

        max_idx_dia = max(c.idx for c in self.geom_dia)
        max_idx_sys = max(c.idx for c in self.geom_sys)

        def _compute_ref(idx, max_idx, geom_list, translations, rotations):
            z = geom_list[idx].points_z
            pt = base.copy()
            dx, dy = translations[idx]
            pt += np.array([dx, dy])
            pivot = np.array(geom_list[idx].centroid[:2])
            θ = np.deg2rad(rotations[idx])
            c, s = np.cos(θ), np.sin(θ)
            rel = pt - pivot
            rot = np.array([rel[0] * c - rel[1] * s,
                            rel[0] * s + rel[1] * c])
            final = pivot + rot
            return (max_idx, final[0], final[1], float(z))

        self.reference_dia = _compute_ref(idx_dia, max_idx_dia, self.geom_dia, self.translations_dia_rest, self.rotations_dia_rest)
        self.reference_sys = _compute_ref(idx_sys, max_idx_sys, self.geom_sys, self.translations_sys_rest, self.rotations_sys_rest)
        return self.reference_dia, self.reference_sys

    def _make_contour(self,
                      idx: int,
                      elliptic_ratio: float,
                      area: float,
                      aortic_thickness: float,
                      pulmonary_thickness: float,
                      z: float,
                      translation: Tuple[float,float],
                      rotation: float) -> Contour:
        """
        Helper to construct a single Contour from the given parameters,
        apply translation & rotation, then recompute centroid.
        """
        c = Contour()
        c.idx = idx

        # create base ellipse
        x, y = self.create_ellipse(area, elliptic_ratio)
        c.points_x = x.tolist()
        c.points_y = y.tolist()
        c.points_z = float(z)
        c.aortic_thickness = aortic_thickness
        c.pulmonary_thickness = pulmonary_thickness

        # initial centroid
        c.calculate_centroid()

        # transform
        self._translate_single(c, translation)
        self._rotate_single(c, rotation)

        # updated centroid
        c.calculate_centroid()

        return c

    def calculate_z_coords(self):
        """
        Returns two NumPy arrays (z_coords_dia, z_coords_sys) giving the cumulative
        pullback distance (in mm) at each contour index. The first entry in each array is 0.0.
        """
        def _compute_z(idxs):
            frame_diffs = np.diff(idxs, prepend=idxs[0])
            dist_mm = frame_diffs / self.frame_rate * self.pullback_speed
            return np.cumsum(dist_mm)

        return _compute_z(self.idx_dia_rest_sorted), _compute_z(self.idx_sys_rest_sorted)

    def create_ellipse(self, area, elliptic_ratio, num_points=501) -> Tuple[np.ndarray, np.ndarray]:
        a = np.sqrt(area / (np.pi * elliptic_ratio))
        b = elliptic_ratio * a
        t = np.linspace(0, 2 * np.pi, num_points)
        return a * np.cos(t) + 4.5, b * np.sin(t) + 4.5

    def _translate_single(self, contour: Contour, translation: Tuple[float, float]):
        dx, dy = translation
        contour.points_x = (np.array(contour.points_x) + dx).tolist()
        contour.points_y = (np.array(contour.points_y) + dy).tolist()
        if contour.centroid is not None:
            cx, cy, cz = contour.centroid
            contour.centroid = (cx + dx, cy + dy, cz)

    def _rotate_single(self, contour: Contour, degree: float):
        θ = np.deg2rad(degree)
        c, s = np.cos(θ), np.sin(θ)
        if contour.centroid is None:
            contour.calculate_centroid()
        cx, cy, _ = contour.centroid

        x = np.array(contour.points_x) - cx
        y = np.array(contour.points_y) - cy
        x_rot = x * c - y * s
        y_rot = x * s + y * c

        contour.points_x = (x_rot + cx).tolist()
        contour.points_y = (y_rot + cy).tolist()

    def _plot_contour(self):
        if not self.geom_dia:
            return
        contour = self.geom_dia[len(self.geom_dia) - 1]
        plt.scatter(contour.points_x, contour.points_y)
        plt.scatter(self.reference_dia[1], self.reference_dia[2])
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.xlim(0, 9)
        plt.ylim(0, 9)
        plt.title(f'Contour at Z={contour.points_z:.2f} mm')
        plt.gca().set_aspect('equal', adjustable='box')
        plt.show()
