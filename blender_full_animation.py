import bpy
import os
from mathutils import Vector, Matrix, Quaternion
import math

# Configuration - UPDATE THIS PATH
obj_directory = "C:/WorkingData/Documents/2_Coding/Rust/aaoca_compression_simulation/output/rest"
# Object animation (only one mesh visible at a time) runs for frames 1 to 62.
object_end_frame = 62
# Increase z-axis orbit duration to twice as long as before.
z_orbit_frames = 124
# Let the y-axis (vertical tilt) phase last 124 frames as well.
y_orbit_frames = 124
total_frames = object_end_frame + z_orbit_frames + y_orbit_frames
frame_rate = 64

def clean_scene():
    bpy.ops.object.select_all(action='SELECT')
    bpy.ops.object.delete()

def import_sequence():
    imported_groups = []
    for i in range(32):
        file_path = os.path.join(obj_directory, f"mesh_{i:03}_rest.obj")
        bpy.ops.import_scene.obj(
            filepath=file_path,
            axis_forward='Y',
            axis_up='Z',
            use_split_objects=True,
            use_split_groups=True,
            use_image_search=True
        )
        imported_group = sorted(bpy.context.selected_objects.copy(), key=lambda obj: obj.name)
        imported_groups.append(imported_group)
        bpy.ops.object.select_all(action='DESELECT')
    return imported_groups

def import_catheter_sequence():
    imported_catheter_groups = []
    for i in range(32):
        file_path = os.path.join(obj_directory, f"catheter_{i:03}_rest.obj")
        bpy.ops.import_scene.obj(
            filepath=file_path,
            axis_forward='Y',
            axis_up='Z',
            use_split_objects=True,
            use_split_groups=True,
            use_image_search=True
        )
        imported_group = sorted(bpy.context.selected_objects.copy(), key=lambda obj: obj.name)
        imported_catheter_groups.append(imported_group)
        bpy.ops.object.select_all(action='DESELECT')
    return imported_catheter_groups

def setup_animation(groups, total_frames, frame_rate, object_end_frame):
    scene = bpy.context.scene
    scene.frame_start = 1
    scene.frame_end = total_frames
    scene.render.fps = frame_rate

    current_frame = 1

    while current_frame < total_frames:
        # Hide all objects at the start of the loop.
        for group in groups:
            for obj in group:
                obj.hide_viewport = True
                obj.hide_render = True
                obj.keyframe_insert(data_path="hide_viewport", frame=current_frame)
                obj.keyframe_insert(data_path="hide_render", frame=current_frame)

        # Animate visibility: each group gets visible at its designated frame.
        for obj_idx, group in enumerate(groups):
            start_frame = current_frame + obj_idx + 1
            end_frame = current_frame + object_end_frame + 1 - obj_idx

            for obj in group:
                for frame in [start_frame, end_frame]:
                    if frame < total_frames:
                        obj.hide_viewport = False
                        obj.hide_render = False
                        obj.keyframe_insert(data_path="hide_viewport", frame=frame)
                        obj.keyframe_insert(data_path="hide_render", frame=frame)

                        if frame + 1 < total_frames:
                            obj.hide_viewport = True
                            obj.hide_render = True
                            obj.keyframe_insert(data_path="hide_viewport", frame=frame + 1)
                            obj.keyframe_insert(data_path="hide_render", frame=frame + 1)

        # Set constant interpolation for all keyframes.
        for group in groups:
            for obj in group:
                if obj.animation_data and obj.animation_data.action:
                    for fcurve in obj.animation_data.action.fcurves:
                        if 'hide_viewport' in fcurve.data_path or 'hide_render' in fcurve.data_path:
                            for keyframe in fcurve.keyframe_points:
                                keyframe.interpolation = 'CONSTANT'

        # Move to the next cycle
        current_frame += object_end_frame

def setup_render_settings():
    bpy.context.scene.render.engine = 'BLENDER_EEVEE'
    bpy.context.scene.render.image_settings.file_format = 'FFMPEG'
    bpy.context.scene.render.ffmpeg.format = 'MPEG4'
    bpy.context.scene.render.ffmpeg.codec = 'H264'
    bpy.context.scene.render.filepath = os.path.join(obj_directory, "animation.mp4")
    bpy.context.scene.render.ffmpeg.constant_rate_factor = 'PERC_LOSSLESS'
    bpy.context.scene.render.resolution_percentage = 100
    bpy.context.scene.render.resolution_x = 1920
    bpy.context.scene.render.resolution_y = 1080

def setup_camera():
    min_coord = Vector((float('inf'),) * 3)
    max_coord = -min_coord.copy()
    for obj in bpy.data.objects:
        if obj.type == 'MESH':
            for corner in obj.bound_box:
                world_corner = obj.matrix_world @ Vector(corner)
                min_coord.x = min(min_coord.x, world_corner.x)
                min_coord.y = min(min_coord.y, world_corner.y)
                min_coord.z = min(min_coord.z, world_corner.z)
                max_coord.x = max(max_coord.x, world_corner.x)
                max_coord.y = max(max_coord.y, world_corner.y)
                max_coord.z = max(max_coord.z, world_corner.z)
    scene_center = (min_coord + max_coord) / 2
    scene_dimensions = max_coord - min_coord
    max_dimension = max(scene_dimensions)
    
    bpy.ops.object.empty_add(type='PLAIN_AXES', location=scene_center)
    camera_target = bpy.context.object
    
    bpy.ops.object.camera_add()
    camera = bpy.context.object
    camera_distance = max_dimension * 2.5
    camera.location = scene_center + Vector((0, -camera_distance, camera_distance * 0.6))
    
    constraint = camera.constraints.new('TRACK_TO')
    constraint.target = camera_target
    constraint.track_axis = 'TRACK_NEGATIVE_Z'
    constraint.up_axis = 'UP_Y'
    
    bpy.context.scene.camera = camera
    camera.data.lens = 35
    camera.data.clip_end = camera_distance * 10
    
    return camera, camera_target

# def animate_camera(camera, target, start_frame, z_orbit_frames, y_orbit_frames):
#     """
#     Animates the camera in two phases while keeping it locked onto the vessel (target):
#     Phase 1: Horizontal orbit (around global Z-axis) lasting z_orbit_frames.
#     Phase 2: Vertical tilt (rotate about the camera's right vector) lasting y_orbit_frames.
#     """
#     scene_center = target.location.copy()
#     initial_offset = camera.location - scene_center

#     # Phase 1: Rotate around Z-axis (horizontal orbit).
#     z_start = start_frame
#     z_end = start_frame + z_orbit_frames - 1
#     for f in range(z_start, z_end + 1):
#         t = (f - z_start) / (z_orbit_frames - 1)
#         angle = 2 * math.pi * t  # Full 360° orbit.
#         rot_z = Matrix.Rotation(angle, 4, 'Z')
#         new_offset = rot_z @ initial_offset
#         camera.location = scene_center + new_offset
#         camera.keyframe_insert(data_path="location", frame=f)

#     # Phase 2: Rotate around the camera's right vector to tilt vertically.
#     # Use the offset at the end of Phase 1 as base.
#     phase1_end_offset = camera.location - scene_center
#     # Compute the camera's right vector from its current orientation.
#     right_vector = camera.matrix_world.to_quaternion() @ Vector((1, 0, 0))
#     y_start = z_end + 1
#     y_end = y_start + y_orbit_frames - 1
#     for f in range(y_start, y_end + 1):
#         t = (f - y_start) / (y_orbit_frames - 1)
#         angle = 2 * math.pi * t  # Full 360° tilt.
#         rot = Matrix.Rotation(angle, 4, right_vector)
#         new_offset = rot @ phase1_end_offset
#         camera.location = scene_center + new_offset
#         camera.keyframe_insert(data_path="location", frame=f)

def animate_camera(camera, target, start_frame, z_orbit_frames, y_orbit_frames):
    """
    Animates the camera in two phases while keeping it locked onto the vessel (target):
    Phase 1: Horizontal orbit (around global Z-axis) lasting z_orbit_frames.
    Phase 2: Vertical tilt (rotate about the camera's right vector) lasting y_orbit_frames.
    """
    scene_center = target.location.copy()
    initial_offset = camera.location - scene_center

    # Phase 1: Rotate around Z-axis (horizontal orbit).
    z_start = start_frame
    z_end = start_frame + z_orbit_frames - 1
    for f in range(z_start, z_end + 1):
        t = (f - z_start) / (z_orbit_frames - 1)
        angle = 2 * math.pi * t  # Full 360° orbit.
        rot_z = Matrix.Rotation(angle, 4, 'Z')
        new_offset = rot_z @ initial_offset
        camera.location = scene_center + new_offset
        camera.keyframe_insert(data_path="location", frame=f)

    # Phase 2: Rotate around the camera's right vector to tilt vertically.
    # Use the offset at the end of Phase 1 as base.
    phase1_end_offset = camera.location - scene_center
    # Compute the camera's right vector from its current orientation.
    right_vector = camera.matrix_world.to_quaternion() @ Vector((1, 0, 0))
    y_start = z_end + 1
    y_end = y_start + y_orbit_frames - 1
    
    # We need to track the camera's current quaternion orientation and apply rotation.
    current_rotation = camera.rotation_euler.copy()
    for f in range(y_start, y_end + 1):
        t = (f - y_start) / (y_orbit_frames - 1)
        angle = 2 * math.pi * t  # Full 360° tilt.
        
        # Instead of directly rotating around the right vector, use quaternion-based rotation
        quaternion_rotation = Quaternion(right_vector, angle)  # Create a quaternion for the tilt
        rotated_offset = quaternion_rotation @ phase1_end_offset
        
        camera.location = scene_center + rotated_offset
        camera.rotation_euler = current_rotation  # Ensure the rotation matches the camera's orientation
        camera.keyframe_insert(data_path="location", frame=f)

        # Keep the rotation smooth and consistent for vertical tilt without flipping.
        current_rotation = camera.rotation_euler.copy()

if __name__ == "__main__":
    clean_scene()
    groups = import_sequence()
    catheter_groups = import_catheter_sequence()
    
    setup_animation(groups, total_frames, frame_rate, object_end_frame)
    setup_animation(catheter_groups, total_frames, frame_rate, object_end_frame)
    
    setup_render_settings()
    
    camera, target = setup_camera()
    
    # Animate the camera:
    # Original object animation: frames 1-62.
    # Then Phase 1 (horizontal orbit) from frame 63 to 62+124 = 186.
    # Then Phase 2 (vertical tilt) from frame 187 to 186+124 = 310.
    animate_camera(camera, target, object_end_frame + 1, z_orbit_frames, y_orbit_frames)
    
    bpy.ops.object.light_add(type='SUN', location=(0, -2, 5), rotation=(0.785, 0, 0))
    bpy.ops.object.light_add(type='SUN', location=(0, 2, 5), rotation=(-0.785, 0, 0))
    
    print("Rendering animation...")
    bpy.ops.render.render(animation=True)
    print(f"Animation saved to: {bpy.context.scene.render.filepath}")
