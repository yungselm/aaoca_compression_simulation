import bpy
import os
from mathutils import Vector

# Configuration - UPDATE THIS PATH
obj_directory = "/mnt/c/WorkingData/Documents/2_Coding/Rust/aaoca_compression_simulation/output/rest"
frame_rate = 64
total_frames = 62

def clean_scene():
    # Clear existing objects
    bpy.ops.object.select_all(action='SELECT')
    bpy.ops.object.delete()

def import_sequence():
    imported_groups = []
    for i in range(32):
        file_path = os.path.join(obj_directory, f"mesh_{i:03}_rest.obj")
        
        # Use legacy OBJ importer for Blender <3.3
        bpy.ops.import_scene.obj(
            filepath=file_path,
            axis_forward='Y',  # Different parameter name in older versions
            axis_up='Z',
            use_split_objects=True,
            use_split_groups=True,
            use_image_search=True
        )
        
        imported_group = sorted(bpy.context.selected_objects.copy(), key=lambda obj: obj.name)
        imported_groups.append(imported_group)
        bpy.ops.object.select_all(action='DESELECT')
    return imported_groups

def setup_animation(groups):
    # Set animation parameters
    scene = bpy.context.scene
    scene.frame_start = 1
    scene.frame_end = total_frames
    scene.render.fps = frame_rate

    # First hide all objects at frame 0
    for group in groups:
        for obj in group:
            obj.hide_viewport = True
            obj.hide_render = True
            obj.keyframe_insert(data_path="hide_viewport", frame=0)
            obj.keyframe_insert(data_path="hide_render", frame=0)

    # Create visibility animation
    for obj_idx, group in enumerate(groups):
        # Calculate visible frames
        if obj_idx == 0:
            frames = [1]
        elif obj_idx == 31:
            frames = [32]
        else:
            frames = [obj_idx + 1, 63 - obj_idx]

        for obj in group:
            for frame in frames:
                # Set visible at current frame
                obj.hide_viewport = False
                obj.hide_render = False
                obj.keyframe_insert(data_path="hide_viewport", frame=frame)
                obj.keyframe_insert(data_path="hide_render", frame=frame)
                
                # Set hidden at next frame
                if frame < total_frames:
                    obj.hide_viewport = True
                    obj.hide_render = True
                    obj.keyframe_insert(data_path="hide_viewport", frame=frame+1)
                    obj.keyframe_insert(data_path="hide_render", frame=frame+1)

    # Set constant interpolation for all keyframes
    for group in groups:
        for obj in group:
            if obj.animation_data and obj.animation_data.action:
                for fcurve in obj.animation_data.action.fcurves:
                    if 'hide_viewport' in fcurve.data_path or 'hide_render' in fcurve.data_path:
                        for keyframe in fcurve.keyframe_points:
                            keyframe.interpolation = 'CONSTANT'

def setup_render_settings():
    # Set render engine
    bpy.context.scene.render.engine = 'BLENDER_EEVEE'
    
    # Output settings
    bpy.context.scene.render.image_settings.file_format = 'FFMPEG'
    bpy.context.scene.render.ffmpeg.format = 'MPEG4'
    bpy.context.scene.render.ffmpeg.codec = 'H264'
    bpy.context.scene.render.filepath = "/mnt/c/WorkingData/Documents/2_Coding/Rust/aaoca_compression_simulation/output/animation.mp4"
    
    # Quality settings
    bpy.context.scene.render.ffmpeg.constant_rate_factor = 'PERC_LOSSLESS'
    bpy.context.scene.render.resolution_percentage = 100
    bpy.context.scene.render.resolution_x = 1920
    bpy.context.scene.render.resolution_y = 1080

def setup_camera():
    # Calculate bounding box of all objects
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

    # Calculate scene center and dimensions
    scene_center = (min_coord + max_coord) / 2
    scene_dimensions = max_coord - min_coord
    max_dimension = max(scene_dimensions)
    
    # Create camera target
    bpy.ops.object.empty_add(type='PLAIN_AXES', location=scene_center)
    camera_target = bpy.context.object
    
    # Create and position camera
    bpy.ops.object.camera_add()
    camera = bpy.context.object
    camera_distance = max_dimension * 2.5  # Adjusted multiplier
    camera.location = scene_center + Vector((0, -camera_distance, camera_distance * 0.6))
    
    # Add tracking constraint
    constraint = camera.constraints.new('TRACK_TO')
    constraint.target = camera_target
    constraint.track_axis = 'TRACK_NEGATIVE_Z'
    constraint.up_axis = 'UP_Y'
    
    # Set camera settings
    bpy.context.scene.camera = camera
    camera.data.lens = 35  # 35mm focal length
    camera.data.clip_end = camera_distance * 10  # Ensure far clip is enough
    
    return camera, camera_target

if __name__ == "__main__":
    clean_scene()
    groups = import_sequence()
    setup_animation(groups)
    setup_render_settings()
    
    # Setup camera and lighting
    camera, target = setup_camera()
    
    # Improved lighting
    bpy.ops.object.light_add(type='SUN', location=(0, -2, 5), rotation=(0.785, 0, 0))
    bpy.ops.object.light_add(type='SUN', location=(0, 2, 5), rotation=(-0.785, 0, 0))
    
    print("Rendering animation...")
    bpy.ops.render.render(animation=True)
    print(f"Animation saved to: {bpy.context.scene.render.filepath}")
