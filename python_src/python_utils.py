import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Configuration
filename = "test_cl.txt"  # Update with your filename
# filename = r"C:\WorkingData\Documents\2_Coding\Rust\aaoca_compression_simulation\PDJ1W2CV_csv_files\Anonymous_Female_Centerline.txt"
# obj_filename = r"C:\WorkingData\Documents\2_Coding\Rust\aaoca_compression_simulation\PDJ1W2CV_csv_files\boeck.obj"
target_distance = 0.403333333
# target_distance = 0.91711111
obj_filename = "test.obj"  # Update with your OBJ filename

# Read and parse centerline data
points = []
with open(filename, 'r') as f:
    lines = f.readlines()

start_parsing = False
for line in lines:
    stripped = line.strip()
    
    if stripped.startswith('Px') and 'Py' in stripped and 'Pz' in stripped:
        start_parsing = True
        continue
    
    if start_parsing:
        if not stripped:
            continue
        
        parts = line.split()
        if len(parts) >= 3:
            try:
                px, py, pz = map(float, parts[:3])
                points.append([px, py, pz])
            except ValueError:
                continue

if not points:
    raise ValueError("No valid points found in the centerline file")

points = np.array(points)

# reverse the points
points = points[::-1]

# remove last 15 points
points = points[:-4]

# reverse the points
points = points[::-1]

# Read and reduce OBJ file points
obj_points = []
with open(obj_filename, 'r') as f:
    for line in f:
        if line.startswith('v '):  # Only process vertex lines
            parts = line.split()
            if len(parts) >= 4:
                try:
                    x = float(parts[1]) * 1000
                    y = float(parts[2]) * 1000
                    z = float(parts[3]) * 1000
                    obj_points.append([x, y, z])
                except ValueError:
                    continue

if not obj_points:
    raise ValueError("No valid vertices found in the OBJ file")

# Reduce OBJ points (take every 100th point)
# reduced_obj = np.array(obj_points[::2])
reduced_obj = np.array(obj_points)

# Calculate cumulative distances for centerline
diffs = np.diff(points, axis=0)
segment_distances = np.linalg.norm(diffs, axis=1)
cumulative_distances = np.zeros(len(points))
cumulative_distances[1:] = np.cumsum(segment_distances)

# Generate target sampling points
total_length = cumulative_distances[-1]
s_new = np.arange(0, total_length + target_distance, target_distance)
s_new = s_new[s_new <= total_length]

# Linear interpolation function
def interpolate_points(s):
    i = np.searchsorted(cumulative_distances, s, side='right') - 1
    i = max(0, min(i, len(points) - 2))
    t = (s - cumulative_distances[i]) / (cumulative_distances[i+1] - cumulative_distances[i] + 1e-9)
    return points[i] + t * (points[i+1] - points[i])

# Resample centerline points
resampled = np.array([interpolate_points(s) for s in s_new])

# Save resampled points
np.savetxt('resampled_centerline.txt', resampled, fmt='%.6f', delimiter=' ')

# Create 3D plot
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')

# Plot original centerline
ax.plot(points[:, 0], points[:, 1], points[:, 2], 
        'b-', linewidth=2, label='Original Centerline')

# Plot resampled centerline
ax.plot(resampled[:, 0], resampled[:, 1], resampled[:, 2], 
        'ro', markersize=4, alpha=0.7, label=f'Resampled ({target_distance} units)')

# Plot reduced OBJ points
ax.plot(reduced_obj[:, 0], reduced_obj[:, 1], reduced_obj[:, 2], 
        'g.', markersize=1, alpha=0.3, label='Reduced OBJ Points (every 100th)')

# Plot three points with different colors
ax.scatter(9.5840, -201.3643, 1750.4165, color='cyan', s=50, label='Point 1')
ax.scatter(11.9164, -202.1920, 1754.6775, color='magenta', s=50, label='Point 2')
ax.scatter(15.5806, -202.1920, 1750.1251, color='yellow', s=50, label='Point 3')

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.legend()
plt.title('Centerline Resampling with OBJ Comparison')
plt.show()

# path_reference = "input/rest_csv_files"

# reference_contour = pd.read_csv(os.path.join(path_reference, "diastolic_contours.csv"), sep="\t", names=["index", "x", "y", "z"])
# reference_point = pd.read_csv(os.path.join(path_reference, "diastolic_reference_points.csv"), sep="\t", names=["index", "x", "y", "z"])

# #keep only every 100th point in reference_contour
# reference_contour = reference_contour.iloc[::10]

# def find_farthest_points(contour):
#     max_dist = 0
#     farthest_pair = None

#     contour = contour.reset_index(drop=True)  # Reset index to ensure integer indexing

#     for i in range(0, len(contour)):
#         for j in range(i + 1, len(contour)):
#             dx = contour.loc[i, 'x'] - contour.loc[j, 'x']
#             dy = contour.loc[i, 'y'] - contour.loc[j, 'y']
#             dist = np.sqrt(dx * dx + dy * dy)
#             if dist > max_dist:
#                 max_dist = dist
#                 farthest_pair = (contour.iloc[i], contour.iloc[j])
    
#     return farthest_pair, max_dist

# def compute_centroid(contour):
#     center_x = contour['x'].mean()
#     center_y = contour['y'].mean()
#     return center_x, center_y

# def sort_contour_points(contour):
#     center_x, center_y = compute_centroid(contour)

#     # Compute angles in radians around centroid
#     contour = contour.copy()
#     contour['angle'] = np.arctan2(contour['y'] - center_y, contour['x'] - center_x)

#     # Sort by angle (counterclockwise)
#     contour = contour.sort_values(by='angle').reset_index(drop=True)

#     # Find the index of the point with the highest y-value
#     start_idx = contour['y'].idxmax()

#     # Rotate the DataFrame so the highest y-value point comes first
#     contour = pd.concat([contour.iloc[start_idx:], contour.iloc[:start_idx]]).reset_index(drop=True)

#     # Drop the helper column
#     contour = contour.drop(columns=['angle'])

#     return contour

# def rotate_point(x, y, angle, center):
#     cx, cy = center
#     # Translate point back to origin:
#     x -= cx
#     y -= cy
#     # Rotate point
#     x_new = x * np.cos(angle) - y * np.sin(angle)
#     y_new = x * np.sin(angle) + y * np.cos(angle)
#     # Translate point back:
#     x_rot = x_new + cx
#     y_rot = y_new + cy
#     return x_rot, y_rot

# def rotate_contour(contour, angle, center):
#     contour = contour.copy()  # avoid mutating original DataFrame
#     rotated = contour.apply(lambda row: rotate_point(row['x'], row['y'], angle, center), axis=1)
#     contour[['x', 'y']] = list(rotated)
#     return contour

# reference_contour = reference_contour[reference_contour.iloc[:, 0] == 775]

# farthest_pair, dist = find_farthest_points(reference_contour)

# print(farthest_pair)

# print(reference_contour)
# print(reference_point)

# plt.scatter(reference_contour['x'], reference_contour['y'], label='Contour Points')
# plt.scatter(reference_point['x'], reference_point['y'], label='Reference Points')
# plt.scatter(farthest_pair[0]['x'], farthest_pair[0]['y'], color='red', label='Farthest Pair')
# plt.scatter(farthest_pair[1]['x'], farthest_pair[1]['y'], color='red')

# # Annotate each point with its corresponding row index from the DataFrame
# for i, row in reference_contour.iterrows():
#     plt.text(row['x'], row['y'], str(row.name), fontsize=8, color='blue')

# for i, row in reference_point.iterrows():
#     plt.text(row['x'], row['y'], str(row.name), fontsize=8, color='green')

# plt.xlim(0, 10)
# plt.ylim(0, 10)
# plt.legend()
# plt.show()

# reference_contour = sort_contour_points(reference_contour)

# plt.scatter(reference_contour['x'], reference_contour['y'], label='Contour Points')
# plt.scatter(reference_point['x'], reference_point['y'], label='Reference Points')
# plt.scatter(farthest_pair[0]['x'], farthest_pair[0]['y'], color='red', label='Farthest Pair')
# plt.scatter(farthest_pair[1]['x'], farthest_pair[1]['y'], color='red')

# # Annotate each point with its corresponding row index from the DataFrame
# for i, row in reference_contour.iterrows():
#     plt.text(row['x'], row['y'], str(row.name), fontsize=8, color='blue')

# for i, row in reference_point.iterrows():
#     plt.text(row['x'], row['y'], str(row.name), fontsize=8, color='green')

# plt.xlim(0, 10)
# plt.ylim(0, 10)
# plt.legend()
# plt.show()

# reference_contour = rotate_contour(reference_contour, np.radians(30), compute_centroid(reference_contour))
# reference_point = reference_point.copy()
# rotated_ref_points = reference_point.apply(
#     lambda row: rotate_point(row['x'], row['y'], np.radians(30), (5.0, 5.0)), axis=1
# )
# reference_point[['x', 'y']] = list(rotated_ref_points)

# plt.scatter(reference_contour['x'], reference_contour['y'], label='Contour Points')
# plt.scatter(reference_point['x'], reference_point['y'], label='Reference Points')
# plt.scatter(farthest_pair[0]['x'], farthest_pair[0]['y'], color='red', label='Farthest Pair')
# plt.scatter(farthest_pair[1]['x'], farthest_pair[1]['y'], color='red')

# # Annotate each point with its corresponding row index from the DataFrame
# for i, row in reference_contour.iterrows():
#     plt.text(row['x'], row['y'], str(row.name), fontsize=8, color='blue')

# for i, row in reference_point.iterrows():
#     plt.text(row['x'], row['y'], str(row.name), fontsize=8, color='green')

# plt.xlim(0, 10)
# plt.ylim(0, 10)
# plt.legend()
# plt.show()