import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Configuration
filename = "test_cl.txt"  # Update with your filename
target_distance = 0.40439
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

# remove last 15 points
points = points[:-12]

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

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.legend()
plt.title('Centerline Resampling with OBJ Comparison')
plt.show()