import bpy
import os
from mathutils import Vector, Matrix, Quaternion
import math

# Configuration - UPDATE THIS PATH
obj_directory = "C:/WorkingData/Documents/2_Coding/Rust/aaoca_compression_simulation/output/rest"
# Object animation (only one mesh visible at a time) runs for frames 1 to 62.
object_end_frame = 62
# Increase z-axis orbit duration to twice as long as before.
z_orbit_frames = 248
# Let the y-axis (vertical tilt) phase last 124 frames as well.
y_orbit_frames = 248
total_frames = object_end_frame + z_orbit_frames + y_orbit_frames
frame_rate = 64

def clean_scene():
    bpy.ops.object.select_all(action='SELECT')
    bpy.ops.object.delete()

def flip_normals():
    """Flip normals for all selected objects and print the direction of one normal before and after flipping."""
    for obj in bpy.context.selected_objects:
        bpy.context.view_layer.objects.active = obj
        bpy.ops.object.mode_set(mode='EDIT')
        
        # Get the direction of one normal before flipping
        bpy.ops.mesh.select_all(action='SELECT')
        bpy.ops.mesh.normals_make_consistent(inside=False)
        bpy.ops.object.mode_set(mode='OBJECT')
        normal_before = obj.data.polygons[0].normal.copy()
        print(f"Normal before flipping: {normal_before}")
        
        # Flip normals
        bpy.ops.object.mode_set(mode='EDIT')
        bpy.ops.mesh.flip_normals()
        bpy.ops.object.mode_set(mode='OBJECT')
        
        # Get the direction of one normal after flipping
        normal_after = obj.data.polygons[0].normal.copy()
        print(f"Normal after flipping: {normal_after}")

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
        # Flip normals after import
        flip_normals()
        imported_group = sorted(bpy.context.selected_objects.copy(), key=lambda obj: obj.name)
        imported_groups.append(imported_group)
        bpy.ops.object.select_all(action='DESELECT')

    print(f"Imported {len(imported_groups)} groups")
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
        # Flip normals after import
        flip_normals()
        imported_group = sorted(bpy.context.selected_objects.copy(), key=lambda obj: obj.name)
        imported_catheter_groups.append(imported_group)
        bpy.ops.object.select_all(action='DESELECT')
    return imported_catheter_groups

def setup_animation(groups, catheter_groups, total_frames, frame_rate, object_end_frame):
    scene = bpy.context.scene
    scene.frame_start = 0
    scene.frame_end = total_frames
    scene.render.fps = frame_rate

    num_objects = len(groups)  # Number of mesh/catheter pairs
    
    for frame in range(total_frames):
        current_index = frame % (2 * num_objects - 2)
        if current_index >= num_objects:
            current_index = 2 * num_objects - 2 - current_index
        
        # Hide all objects
        for mesh_group, catheter_group in zip(groups, catheter_groups):
            for obj in mesh_group + catheter_group:
                obj.hide_viewport = True
                obj.hide_render = True
                obj.keyframe_insert(data_path="hide_viewport", frame=frame)
                obj.keyframe_insert(data_path="hide_render", frame=frame)
        
        # Show only the current object pair
        for obj in groups[current_index] + catheter_groups[current_index]:
            obj.hide_viewport = False
            obj.hide_render = False
            obj.keyframe_insert(data_path="hide_viewport", frame=frame)
            obj.keyframe_insert(data_path="hide_render", frame=frame)
    
    # Set constant interpolation for visibility keyframes
    for group in groups + catheter_groups:
        for obj in group:
            if obj.animation_data and obj.animation_data.action:
                for fcurve in obj.animation_data.action.fcurves:
                    if 'hide_viewport' in fcurve.data_path or 'hide_render' in fcurve.data_path:
                        for keyframe in fcurve.keyframe_points:
                            keyframe.interpolation = 'CONSTANT'

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

def animate_camera(camera, target, start_frame, z_orbit_frames, y_orbit_frames):
    """
    Animates the camera in two phases while keeping it locked onto the vessel (target):

    Phase 1: Horizontal orbit (around global Z-axis) lasting z_orbit_frames.
    Phase 2: Vertical tilt with custom segments:
      - Segment 1: Tilt 89° in the current direction.
      - Segment 2: Tilt 179° in the opposite direction.
      - Segment 3: Hold this final position for 64 frames.
    """
    scene_center = target.location.copy()
    initial_offset = camera.location - scene_center

    # --- Phase 1: Horizontal Orbit ---
    z_start = start_frame
    z_end = start_frame + z_orbit_frames - 1
    for f in range(z_start, z_end + 1):
        t = (f - z_start) / (z_orbit_frames - 1)
        angle = 2 * math.pi * t  # Full 360° orbit.
        rot_z = Matrix.Rotation(angle, 4, 'Z')
        new_offset = rot_z @ initial_offset
        camera.location = scene_center + new_offset
        camera.keyframe_insert(data_path="location", frame=f)

    # --- Phase 2: Vertical Tilt ---
    # Determine how many frames to allocate to tilt motion versus hold.
    hold_frames = 64
    tilt_motion_frames = y_orbit_frames - hold_frames  # Total frames available for tilting.
    
    # Divide the tilt motion in two segments: 89° and 179°.
    # (Using the ratio of the angles: 89 : 179)
    segment1_frames = round(tilt_motion_frames * (89 / (89 + 179)))
    segment2_frames = tilt_motion_frames - segment1_frames

    # Calculate the right vector (axis of tilt rotation).
    right_vector = camera.matrix_world.to_quaternion() @ Vector((1, 0, 0))
    
    # Use the camera's position at the end of Phase 1 as the starting offset.
    current_offset = camera.location - scene_center

    # --- Segment 1: Tilt from 0° to +100° ---
    angle1_final = math.radians(105)
    seg1_start = z_end + 1
    for f in range(seg1_start, seg1_start + segment1_frames):
        t = (f - seg1_start) / (segment1_frames - 1) if segment1_frames > 1 else 1
        angle = t * angle1_final
        quat = Quaternion(right_vector, angle)
        new_offset = quat @ current_offset
        camera.location = scene_center + new_offset
        camera.keyframe_insert(data_path="location", frame=f)
    
    # Update the offset after segment 1.
    current_offset = camera.location - scene_center

    # --- Segment 2: Tilt from +89° to (89° - 179°) = -90° ---
    angle2_final = math.radians(-160)
    seg2_start = seg1_start + segment1_frames
    for f in range(seg2_start, seg2_start + segment2_frames):
        t = (f - seg2_start) / (segment2_frames - 1) if segment2_frames > 1 else 1
        angle = t * angle2_final
        quat = Quaternion(right_vector, angle)
        new_offset = quat @ current_offset
        camera.location = scene_center + new_offset
        camera.keyframe_insert(data_path="location", frame=f)

    # --- Segment 3: Hold the final position for hold_frames ---
    hold_start = seg2_start + segment2_frames
    final_offset = camera.location - scene_center
    for f in range(hold_start, hold_start + hold_frames):
        camera.location = scene_center + final_offset
        camera.keyframe_insert(data_path="location", frame=f)

def scatter_lights_around_object(center, radius, num_lights, target):
    """
    Scatter lights evenly around the object in a spherical pattern.
    This version uses a golden spiral distribution for more even spacing.
    """
    lights = []
    for i in range(num_lights):
        theta = math.acos(1 - 2 * (i + 0.5) / num_lights)
        phi = math.pi * (1 + 5**0.5) * (i + 0.5)
        x = center.x + radius * math.sin(theta) * math.cos(phi)
        y = center.y + radius * math.sin(theta) * math.sin(phi)
        z = center.z + radius * math.cos(theta)
        bpy.ops.object.light_add(type='SUN', location=(x, y, z))
        light = bpy.context.object
        # Add a tracking constraint so the light always points to the target
        constraint = light.constraints.new('TRACK_TO')
        constraint.target = target
        constraint.track_axis = 'TRACK_NEGATIVE_Z'
        constraint.up_axis = 'UP_Y'
        lights.append(light)
    return lights

def get_object_geometric_center(obj):
    # Transform each bound_box corner from local to global coordinates.
    coords = [obj.matrix_world @ Vector(corner) for corner in obj.bound_box]
    # Compute min and max points.
    min_coord = Vector((min(v[i] for v in coords) for i in range(3)))
    max_coord = Vector((max(v[i] for v in coords) for i in range(3)))
    return (min_coord + max_coord) / 2

if __name__ == "__main__":
    clean_scene()
    groups = import_sequence()
    catheter_groups = import_catheter_sequence()
    
    setup_animation(groups, catheter_groups, total_frames, frame_rate, object_end_frame)
    
    setup_render_settings()
    
    camera, target = setup_camera()
    
    animate_camera(camera, target, object_end_frame + 1, z_orbit_frames, y_orbit_frames)
    
    # Find the mesh_000_rest object to determine the lighting center.
    mesh_obj = None
    for obj in bpy.data.objects:
        if obj.name.startswith("mesh_000_rest"):
            mesh_obj = obj
            break
    if mesh_obj:
        mesh_center = get_object_geometric_center(mesh_obj)
    else:
        mesh_center = target.location

    
    # Scatter lights around the mesh in a spherical pattern.
    light_radius = 10.0  # Adjust the distance as needed.
    num_lights = 12      # Adjust the number of lights as needed.
    lights = scatter_lights_around_object(mesh_center, light_radius, num_lights, target)
    
    # Configure light properties.
    for light in lights:
        light.data.energy = 2.0  # Adjust light intensity.
        light.data.color = (1.0, 1.0, 1.0)  # White light
    
    print("Rendering animation...")
    bpy.ops.render.render(animation=True)
    print(f"Animation saved to: {bpy.context.scene.render.filepath}")
