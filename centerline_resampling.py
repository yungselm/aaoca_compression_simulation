import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Configuration
filename = "test_cl.txt"  # Update with your filename
target_distance = 0.40439

# Read and parse data
points = []
with open(filename, 'r') as f:
    lines = f.readlines()

start_parsing = False
for line in lines:
    stripped = line.strip()
    
    # Find the header line containing Px, Py, Pz
    if stripped.startswith('Px') and 'Py' in stripped and 'Pz' in stripped:
        start_parsing = True
        continue
    
    if start_parsing:
        if not stripped:  # Skip empty lines
            continue
        
        parts = line.split()
        if len(parts) >= 3:
            try:
                # Extract coordinates (first three columns)
                px, py, pz = map(float, parts[:3])
                points.append([px, py, pz])
            except ValueError:
                continue  # Skip lines with invalid data

if not points:
    raise ValueError("No valid points found in the file")

points = np.array(points)

# Calculate cumulative distances
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
    i = max(0, min(i, len(points) - 2))  # Clamp index
    
    t = (s - cumulative_distances[i]) / (cumulative_distances[i+1] - cumulative_distances[i] + 1e-9)
    return points[i] + t * (points[i+1] - points[i])

# Resample points
resampled = np.array([interpolate_points(s) for s in s_new])

# save as txt file
np.savetxt('resampled_centerline.txt', resampled, fmt='%.6f', delimiter=' ')

# Create 3D plot
fig = plt.figure(figsize=(10, 6))
ax = fig.add_subplot(111, projection='3d')

# Plot original points
ax.plot(points[:, 0], points[:, 1], points[:, 2], 
        'b-', linewidth=2, label='Original Centerline')

# Plot resampled points
ax.plot(resampled[:, 0], resampled[:, 1], resampled[:, 2], 
        'ro', markersize=4, alpha=0.7, label=f'Resampled ({target_distance} units)')

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.legend()
plt.title('Centerline Resampling')
plt.show()