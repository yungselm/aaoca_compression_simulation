import numpy as np


FRAME_RATE = 30
PULLBACK_SPEED = 1.0
START_DIA = 22.0
START_SYS = 27.0
NUM_POINTS_CONTOUR = 501
NUM_POINTS_CATHETER = 20

IDX_DIA_REST_SORTED = [126, 155, 183, 211, 238, 267, 296, 324, 353, 511, 377, 463, 489, 407, 436, 542, 596, 569, 622, 649, 775, 748, 678]
IDX_SYS_REST_SORTED = [111, 142, 170, 197, 224, 252, 280, 309, 339, 366, 393, 474, 501, 447, 421, 526, 581, 556, 635, 608, 661, 748, 678]

ELLIP_DIA_REST = [1.0, 1.1, 1.2, 1.2, 1.2, 1.2, 1.2, 1.0, 1.2, 1.3, 1.2, 1.3, 1.3, 1.4, 1.4, 1.4, 1.4, 1.5, 1.6, 1.8, 2.1, 1.9, 2.2]
ELLIP_SYS_REST = [1.2, 1.1, 1.0, 1.1, 1.1, 1.0, 1.1, 1.0, 1.0, 1.1, 1.4, 1.3, 1.4, 1.4, 1.8, 1.7, 3.4, 3.1, 2.8, 2.7, 3.8, 1.9, 2.2]

AREA_DIA_REST = [13.4, 13.1, 13.5, 14.0, 12.4, 12.0, 11.0, 13.2, 11.8, 8.1, 10.3, 8.9, 8.4, 8.8, 8.8, 7.3, 7.9, 7.4, 8.3, 9.1, 8.3, 8.8, 8.1]
AREA_SYS_REST = [12.6, 10.6, 10.1, 10.8, 10.2, 9.9, 9.6, 9.1, 9.5, 9.6, 7.5, 7.4, 7.4, 7.3, 6.8, 6.9, 6.7, 6.5, 6.4, 6.7, 5.5, 8.8, 8.1]

AORTIC_DIA_REST = [None, None, None, None, None, None, None, None, None, 1.6, 1.6, 1.5, 1.6, 1.5, 1.6, 1.5, 1.5, 1.5, 1.6, 1.5, 1.4, 1.5, 1.2]
AORTIC_SYS_REST = [None, None, None, None, None, None, None, 1.7, 1.6, 1.5, 1.5, 1.4, 1.4, 1.5, 1.5, 1.8, 1.3, 1.3, 1.3, 1.3, 1.3, 1.5, 1.2]

PULMONARY_DIA_REST = [None, None, None, None, None, None, None, None, 2.2, 1.9, 2.0, 1.9, 1.9, 1.8, 1.9, 1.9, 1.8, 1.8, 1.7, 1.9, 1.7, 1.9, 1.6]
PULMONARY_SYS_REST = [None, None, None, None, None, None, None, None, None, None, None, 2.0, 2.0, 1.9, 1.8, 1.6, 1.8, 1.8, 1.7, 1.9, 1.6, 1.5, 1.2]

# make a 5×5 grid of values between -2.0 and +2.0
_vals = np.linspace(-1.0, 1.0, 5)
_all_pairs = [(float(x), float(y)) for x in _vals for y in _vals]

# pick exactly as many as you need (23)
TRANSLATION_DIA_REST = _all_pairs[:len(IDX_DIA_REST_SORTED)]
TRANSLATION_SYS_REST = _all_pairs[:len(IDX_SYS_REST_SORTED)]
TRANSLATION_DIA_REST[-1] = (0, 0)

ROTATION_DIA_REST = [260, 270, 280, 290, 300, 310, 320, 330, 340, 350, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120]
ROTATION_SYS_REST = [230, 240, 250, 260, 270, 280, 290, 300, 310, 320, 330, 340, 350, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90]