import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Load data
# centerline = pd.read_csv("output/centerline.csv")
geometry_original = pd.read_csv("output/debugging/original_geometry_rest_dia.csv")
geometry_zadjusted = pd.read_csv("output/debugging/zadjusted_geometry_rest_dia.csv")
geometry_aligned = pd.read_csv("output/debugging/aligned_geometry_rest_dia.csv")
geometry_smoothed = pd.read_csv("output/debugging/smoothed_geometry_rest_dia.csv")
geometry_reloaded = pd.read_csv("output/debugging/reloaded_geometry_rest_dia.csv")
data_original = pd.read_csv("input/rest_csv_files/systolic_contours.csv", sep="\t", names=["contour_id", "x", "y", "z"])
data_original_dia = pd.read_csv("input/rest_csv_files/diastolic_contours.csv", sep="\t", names=["contour_id", "x", "y", "z"])
reference_point = pd.read_csv("input/rest_csv_files/systolic_reference_points.csv", sep="\t", names=["contour_id", "x", "y", "z"])
reference_point_dia = pd.read_csv("input/rest_csv_files/diastolic_reference_points.csv", sep="\t", names=["contour_id", "x", "y", "z"])
before_rotation = pd.read_csv("output/debugging/before_rotation_geometry.csv")
after_rotation = pd.read_csv("output/debugging/after_rotation_geometry.csv")

before_rotation = before_rotation[before_rotation['contour_id'] == 25]
after_rotation = after_rotation[after_rotation['contour_id'] == 25]

data_dia = data_original_dia[data_original_dia['contour_id']==803]
data_sys = data_original[data_original['contour_id']==661]
data_dia_reloaded = geometry_reloaded[geometry_reloaded['contour_id']==0]

plt.plot(data_sys['x'], data_sys['y'], 'o', label='diastole', alpha=0.6)
plt.plot(data_dia['x'], data_dia['y'], 'o', label='systole', alpha=0.6)
plt.plot(data_dia_reloaded['x'], data_dia_reloaded['y'], 'o', label='reloaded', alpha=0.6)
plt.scatter(reference_point['x'], reference_point['y'], c='red', marker='x', s=100, linewidths=1.5)
plt.scatter(reference_point_dia['x'], reference_point_dia['y'], c='blue', marker='x', s=100, linewidths=1.5)
plt.xlim(0, 9)
plt.ylim(0, 9)
plt.show()

# # Sort the dataframes by a specific column in descending order
# geometry_original = geometry_original.sort_values(by="contour_id", ascending=False)
# geometry_zadjusted = geometry_zadjusted.sort_values(by="contour_id", ascending=False)
# geometry_aligned = geometry_aligned.sort_values(by="contour_id", ascending=False)


# def display_all_contours_in_geometry_reloaded(geometry_reloaded):
#     # Get sorted list of unique contour IDs from geometry_reloaded
#     contour_ids = sorted(geometry_reloaded['contour_id'].unique(), reverse=True)
    
#     # Set up consistent styling
#     color = 'purple'
#     label = 'Reloaded'
    
#     # Create a single plot for all contours
#     fig = plt.figure(figsize=(12, 10))
#     ax = fig.add_subplot(111, projection='3d')
    
#     for contour_id in contour_ids:
#         # Get subset for current contour ID
#         subset = geometry_reloaded[geometry_reloaded['contour_id'] == contour_id]
        
#         if subset.empty:
#             continue  # Skip if no data for this contour
        
#         # Split into aortic and non-aortic points
#         non_aortic = subset[subset['aortic'] == False]
#         aortic = subset[subset['aortic'] == True]
        
#         # Plot non-aortic points
#         ax.scatter(non_aortic['x'], non_aortic['y'], non_aortic['z'],
#                    c=color, marker='o', label=label if contour_id == contour_ids[0] else "", alpha=0.6)
        
#         # Plot aortic points with annotations
#         if not aortic.empty:
#             ax.scatter(aortic['x'], aortic['y'], aortic['z'],
#                        c=color, marker='x', s=100, linewidths=1.5)
#             for _, row in aortic.iterrows():
#                 ax.text(row['x'], row['y'], row['z'], 
#                         str(row['point_index']), color=color, fontsize=8)
    
#     # Configure plot settings
#     ax.set_xlabel('X')
#     ax.set_xlim(0, 9)
#     ax.set_ylabel('Y')
#     ax.set_ylim(0, 9)
#     ax.set_zlabel('Z')
#     ax.set_zlim(20, 30)
#     ax.legend()
#     plt.title('3D View: All Contours in Geometry Reloaded')
#     plt.show()

# # Call the function to display all contours in geometry_reloaded
# display_all_contours_in_geometry_reloaded(geometry_reloaded)

# # # small z ofsett for geometry_reloaded to visualize better
# geometry_reloaded['z'] += 0.1

# # Get sorted list of unique contour IDs from original geometry
# contour_ids = sorted(geometry_original['contour_id'].unique(), reverse=True)

# # Set up consistent styling
# colors = ['blue', 'green', 'red', 'orange', 'purple']
# labels = ['Original', 'Z-adjusted', 'Aligned', 'Smoothed', 'Reloaded']

# # Create a separate plot for each contour ID
# for contour_id in contour_ids:
#     # Get subsets for current contour ID
#     original_subset = geometry_original[geometry_original['contour_id'] == contour_id]
#     zadjusted_subset = geometry_zadjusted[geometry_zadjusted['contour_id'] == contour_id]
#     aligned_subset = geometry_aligned[geometry_aligned['contour_id'] == contour_id]
#     smoothed_subset = geometry_smoothed[geometry_smoothed['contour_id'] == contour_id]
#     reloaded_subset = geometry_reloaded[geometry_reloaded['contour_id'] == contour_id]
    
#     # Create figure
#     fig = plt.figure(figsize=(12, 10))
#     ax = fig.add_subplot(111, projection='3d')
    
#     # Plot each geometry's data
#     for subset, color, label in zip([original_subset, zadjusted_subset, aligned_subset, smoothed_subset, reloaded_subset], colors, labels):
#         if subset.empty:
#             continue  # Skip if no data for this contour
        
#         # Split into aortic and non-aortic points
#         non_aortic = subset[subset['aortic'] == False]
#         aortic = subset[subset['aortic'] == True]
        
#         # Plot non-aortic points
#         ax.scatter(non_aortic['x'], non_aortic['y'], non_aortic['z'],
#                    c=color, marker='o', label=label, alpha=0.6)
#         for _, row in non_aortic.iterrows():
#             ax.text(row['x'], row['y'], row['z'], 
#                     str(row['point_index']), color=color, fontsize=8)
        
#         # Plot aortic points with annotations
#         if not aortic.empty:
#             ax.scatter(aortic['x'], aortic['y'], aortic['z'],
#                        c=color, marker='x', s=100, linewidths=1.5)
#             for _, row in aortic.iterrows():
#                 ax.text(row['x'], row['y'], row['z'], 
#                         str(row['point_index']), color=color, fontsize=8)
    
#     # Configure plot settings
#     ax.set_xlabel('X')
#     ax.set_ylabel('Y')
#     ax.set_zlabel('Z')
#     ax.legend()
#     plt.title(f'3D View: Geometries Comparison (Contour ID {contour_id})')
#     plt.show()

# Plot non-aortic points
plt.scatter(before_rotation[before_rotation['aortic'] == False]['x'], 
            before_rotation[before_rotation['aortic'] == False]['y'], 
            c='blue', marker='o', label='Before Rotation (Non-Aortic)', alpha=0.6)
for _, row in before_rotation[before_rotation['aortic'] == False].iterrows():
    plt.text(row['x'], row['y'], str(row['point_index']), color='blue', fontsize=8)

plt.scatter(after_rotation[after_rotation['aortic'] == False]['x'], 
            after_rotation[after_rotation['aortic'] == False]['y'], 
            c='red', marker='o', label='After Rotation (Non-Aortic)', alpha=0.6)
for _, row in after_rotation[after_rotation['aortic'] == False].iterrows():
    plt.text(row['x'], row['y'], str(row['point_index']), color='red', fontsize=8)

# Plot aortic points
plt.scatter(before_rotation[before_rotation['aortic'] == True]['x'], 
            before_rotation[before_rotation['aortic'] == True]['y'], 
            c='blue', marker='x', s=100, linewidths=1.5, label='Before Rotation (Aortic)')
for _, row in before_rotation[before_rotation['aortic'] == True].iterrows():
    plt.text(row['x'], row['y'], str(row['point_index']), color='blue', fontsize=8)

plt.scatter(after_rotation[after_rotation['aortic'] == True]['x'], 
            after_rotation[after_rotation['aortic'] == True]['y'], 
            c='red', marker='x', s=100, linewidths=1.5, label='After Rotation (Aortic)')
for _, row in after_rotation[after_rotation['aortic'] == True].iterrows():
    plt.text(row['x'], row['y'], str(row['point_index']), color='red', fontsize=8)

# Configure plot settings
plt.xlim(0, 9)
plt.ylim(0, 9)
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.title('Before and After Rotation Comparison')
plt.legend()
plt.show()