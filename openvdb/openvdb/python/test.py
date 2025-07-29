import numpy as np
import openvdb as vdb  # or your module name
import trimesh

SIMPLE_EXAMPLE = False
MESH_PATH = "rubbertoy.obj"
OUTPUT_PATH = "output.vdb"

if SIMPLE_EXAMPLE:
    # Example mesh: a single triangle
    points = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0]
    ], dtype=np.float32)  # must be float32
    triangles = np.array([
        [0, 1, 2]
    ], dtype=np.uint32)   # must be uint32
    # No quads in this example
    quads = None
else:
    mesh = trimesh.load(MESH_PATH)
    points = np.array(mesh.vertices, dtype=np.float32)
    triangles = np.array(mesh.faces, dtype=np.uint32)
    quads = None
    
voxel_size = 0.01
transform = vdb.createLinearTransform(voxel_size)
ex_band_width = 3.0
in_band_width = 10000.0

grid_fill_interior = vdb.FloatGrid.createLevelSetFromPolygons(
    points,
    triangles,
    quads,
    transform,
    ex_band_width,
    in_band_width,
)
grid_fill_interior.name = "grid_fill_interior"

grid_narrow_band = vdb.FloatGrid.createLevelSetFromPolygons(
    points,
    triangles,
    quads,
    transform,
    ex_band_width,
)
grid_narrow_band.name = "grid_narrow_band"

vdb.write(OUTPUT_PATH, [grid_fill_interior, grid_narrow_band])
