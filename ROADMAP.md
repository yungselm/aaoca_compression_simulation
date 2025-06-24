# Code fixes to do:
- [] Wall calculate the thickness based on distance from centroid, but should take a 90° angle offset or something similar for every point (currently always + fixed distance to distance centroid point creates uneven wall, depending on angle)

# Implementation mesh alignment:
- [] Align also the wall mesh
- [] Extrude the ostium based on a function where ostium slice is repeated and gets bigger

# Implementation CCTA mesh manipulation
- [] Shrink or expand CCTA mesh to best match wallthickness
- [] Shrink or expand coronary mesh to best match distal segments
- Potentially: [] Recreate CCTA mesh by first creating a centerline then keeping only the points corresponding to this position and fit a function through all points
- [] Remove points CCTA mesh that are in range of IVUS mesh
- [] Fully automate mesh creation

# Error handling
- [] remove every unspecified Error handling with '?'

# Unit Tests
## Module: IO
- [] Input
- [] Load Geometry
- [] Mod
- [] Output
## Module: Processing
- [] Comparison
- [] Contours
- [] Geometries
- [] Mod
- [] Process case
- [] Walls
## Module: Texture
- [] Mod
- [] Texture
## Module: Utils
- [] Utils
## Module: Mesh to Centerline
- [] Operations
- [] Preprocessing
- [] Mod