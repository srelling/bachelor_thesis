import matplotlib.pyplot as plt
from netgen.occ import *
from firedrake import *

def plot_mesh(mesh):
    fig, axes = plt.subplots()
    triplot(mesh, axes=axes)
    plt.gca().set_aspect('equal', adjustable='box')
    axes.legend()
    axes.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.show()

# Define the dimensions
rect_width = 10
rect_height = 4
hole_radius = 1
max_mesh_size = 1

# Define points for the rectangle boundary
rect_p1 = Pnt(0, -rect_height / 2, 0)
rect_p2 = Pnt(rect_width, -rect_height / 2, 0)
rect_p3 = Pnt(rect_width, rect_height / 2, 0)
rect_p4 = Pnt(0, rect_height / 2, 0)

# Create the rectangle boundary with segments
rect_seg1 = Segment(rect_p1, rect_p2)
rect_seg2 = Segment(rect_p2, rect_p3)
rect_seg3 = Segment(rect_p3, rect_p4)
rect_seg4 = Segment(rect_p4, rect_p1)
rect_wire = Wire([rect_seg1, rect_seg2, rect_seg3, rect_seg4])

# Define the center point of the circular hole
hole_center = Pnt(2.5, 0)  # Center point in 2D

# Create the circular hole
circle = Circle(hole_center, hole_radius)
circle_edge = circle.Edge()
circle_wire = Wire([circle_edge])

# Create the face by subtracting the circular hole from the rectangle
rect_face = Face(rect_wire)
hole_face = Face(circle_wire)
final_face = rect_face - hole_face

# Generate the mesh
geometry = OCCGeometry(final_face, dim=2)
ngmsh = geometry.GenerateMesh(maxh=max_mesh_size)
for i in range(5):
    ngmsh.Refine()
msh = Mesh(ngmsh)

plot_mesh(msh)