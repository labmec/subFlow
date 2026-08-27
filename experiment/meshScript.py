# ------------------------------------------------------------------------------
#
#  Project:
#
# ------------------------------------------------------------------------------

import gmsh
import sys
import os

# Using OpenCASCADE as the built-in CAD kernel
partitions = 1 # 0 -> tri, 1 -> quad
mesh_type = 0 # 0 -> unstructured, 1 -> structured
mesh = "tri3" if (partitions == 0) else "quad4"
path = "/home/marina/programming/SubFlow-Research/subFlow/experiment"
model_name = "darcy"

# Initial parameters and mesh settings
lc = 0.5
L = 10.
h = 10.
shift = 0. # For AxisSymmetric
an_pts = [[2.,7.],[4.,7.],[6.,7.]]
obstruct = 2. # Assuming it starts from the top
num_divL = 20 # Number of divisions along the length of the beam
num_divH = 20

# Gmsh must be initialized
gmsh.initialize()


# Add a new model
gmsh.model.add(model_name)


# Element entities
# Points
gmsh.model.occ.addPoint(0.+shift, 0, 0, lc, 1)
gmsh.model.occ.addPoint(L+shift, 0, 0, lc, 2)
gmsh.model.occ.addPoint(L+shift, h, 0, lc, 3)
gmsh.model.occ.addPoint(0.+shift, h, 0, lc, 4)
gmsh.model.occ.addPoint(0.+shift, h-obstruct, 0, lc, 5)
# Measurements points
gmsh.model.occ.addPoint(an_pts[0][0]+shift, an_pts[0][1], 0, lc, 6)
gmsh.model.occ.addPoint(an_pts[1][0]+shift, an_pts[1][1], 0, lc, 7)
gmsh.model.occ.addPoint(an_pts[2][0]+shift, an_pts[2][1], 0, lc, 8)

# Curves
gmsh.model.occ.addLine(1, 2, 1)
gmsh.model.occ.addLine(2, 3, 2)
gmsh.model.occ.addLine(3, 4, 3)
gmsh.model.occ.addLine(4, 5, 4)
gmsh.model.occ.addLine(5, 1, 5)

# Curve loop defined by an ordered list of connected curves
gmsh.model.occ.addCurveLoop([1,2,3,4,5], 1)

# Surface, defined by a curve loop
gmsh.model.occ.addPlaneSurface([1], 1)

# CAD entities must be synchronized with the Gmsh model
gmsh.model.occ.synchronize()


# Physical Groups
gmsh.model.addPhysicalGroup(0, [6,7,8], name="measure")  
gmsh.model.addPhysicalGroup(1, [1,3], name="noflux")    
gmsh.model.addPhysicalGroup(1, [5], name="p_in")  
gmsh.model.addPhysicalGroup(1, [2], name="p_out")
gmsh.model.addPhysicalGroup(1, [4], name="obstruct")
gmsh.model.addPhysicalGroup(2, [1], name="domain")


if mesh_type:
    # Structured mesh using transfinite interpolation
    num_divh = int(num_divL * (h / L)/4)
    n_pointsL = num_divL + 1
    n_pointsh = num_divh + 1
    hel = L / num_divL  # Element size along the length of the beam
    gmsh.model.mesh.setTransfiniteCurve(1, n_pointsL)
    gmsh.model.mesh.setTransfiniteCurve(3, n_pointsL)
    gmsh.model.mesh.setTransfiniteCurve(2, n_pointsh)
    gmsh.model.mesh.setTransfiniteCurve(4, n_pointsh)

    gmsh.model.mesh.setTransfiniteSurface(1, "Left", [1, 2, 3, 4])

gmsh.model.occ.synchronize()


# Create quadrangles instead of triangles
if partitions:
    gmsh.model.mesh.setRecombine(2, 1)


# Apply an elliptic smoother to the grid to have a more regular mesh:
# gmsh.option.setNumber("Mesh.Smoothing", 100)


# Generate the mesh. (Dim 0 for points, 1 for lines, 2 for surfaces and 3 for volumes)
gmsh.model.mesh.generate(2)


# Save it to disk
gmsh.write(f"{path}/darcy_{mesh}.msh")


gmsh.option.setNumber("Mesh.SaveAll", 1)


# Launch the GUI to see the results:
# if '-nopopup' not in sys.argv:
#     gmsh.fltk.run()

# This should be called when you are done using the Gmsh Python API:
gmsh.finalize()