import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from geometry import GeometryAssembler


def main():
    ref_dia_rest, geom_dia_rest, ref_sys_rest, geom_sys_rest = GeometryAssembler(
        "rest"
    )()
    ref_dia_stress, geom_dia_stress, ref_sys_stress, geom_sys_stress = (
        GeometryAssembler("stress")()
    )

if __name__ == "__main__":
    main()
