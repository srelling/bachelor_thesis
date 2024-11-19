from firedrake import *
import netgen
from netgen.occ import *
from netgen.geom2d import SplineGeometry
import matplotlib.pyplot as plt
import ufl
import os

def plot_mesh(mesh):
    fig, axes = plt.subplots()
    triplot(mesh, axes=axes)
    plt.gca().set_aspect('equal', adjustable='box')
    axes.legend()
    plt.show()

def generate_ng_square_circle1_exact( h_max ):
    geo = SplineGeometry()
    circle_pnts = [(0, 0), (1, 0), (1, 1),
            (0, 1), (-1, 1), (-1, 0),
            (-1, -1), (0, -1), (1, -1)]
    square_pnts = [(2, 2), (-2, 2), (-2, -2), (2, -2)]
    p1, p2, p3, p4, p5, p6, p7, p8, p9 = [geo.AppendPoint(*pnt) for pnt in circle_pnts]
    p10, p11, p12, p13 = [geo.AppendPoint(*pnt) for pnt in square_pnts]
    circle = [[["spline3", p2, p3, p4], "curve"],
              [["spline3", p4, p5, p6], "curve"],
              [["spline3", p6, p7, p8], "curve"],
              [["spline3", p8, p9, p2], "curve"]]
    square = [[["line", p10, p11], "line"],
              [["line", p11, p12], "line"],
              [["line", p12, p13], "line"],
              [["line", p13, p10], "line"]]
    [geo.Append(c, bc=bc, leftdomain=0, rightdomain=1) for c, bc in circle]
    [geo.Append(c, bc=bc, leftdomain=1, rightdomain=0) for c, bc in square]
    ngmsh = geo.GenerateMesh(maxh=h_max)
    for i in range(5):
        ngmsh.Refine()
    return ngmsh

from firedrake import *
import netgen
from netgen.geom2d import SplineGeometry

def generate_ng_square_circle1( h_max ):
    geo = SplineGeometry()
    circle_pnts = [(0, 0), (1, 0), (1, 1),
            (0, 1), (-1, 1), (-1, 0),
            (-1, -1), (0, -1), (1, -1)]
    square_pnts = [(2, 2), (-2, 2), (-2, -2), (2, -2)]
    p1, p2, p3, p4, p5, p6, p7, p8, p9 = [geo.AppendPoint(*pnt) for pnt in circle_pnts]
    p10, p11, p12, p13 = [geo.AppendPoint(*pnt) for pnt in square_pnts]
    circle = [[["spline3", p2, p3, p4], "curve"],
              [["spline3", p4, p5, p6], "curve"],
              [["spline3", p6, p7, p8], "curve"],
              [["spline3", p8, p9, p2], "curve"]]
    square = [[["line", p10, p11], "line"],
              [["line", p11, p12], "line"],
              [["line", p12, p13], "line"],
              [["line", p13, p10], "line"]]
    [geo.Append(c, bc=bc, leftdomain=0, rightdomain=1) for c, bc in circle]
    [geo.Append(c, bc=bc, leftdomain=1, rightdomain=0) for c, bc in square]
    ngmsh = geo.GenerateMesh(maxh=h_max)
    return ngmsh

def generate_ng_rect_circ1( h_max ):
    # Define the dimensions
    rect_width = 10
    rect_height = 4
    hole_radius = 1

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
    ngmsh = geometry.GenerateMesh(maxh=h_max)
    return ngmsh


def load_or_generate_mesh( mesh_name, h_max=0.4 ):
    vol_file = f"meshes/{mesh_name}_{str(h_max).replace('.', '')}.vol"
    if os.path.exists(vol_file):
        ngmsh = netgen.meshing.Mesh()
        ngmsh.Load(vol_file)
        print(f"loaded mesh from {vol_file}")
    else:
        mesh_function = globals().get(f"generate_{mesh_name}")
        if mesh_function is None:
            raise ValueError(f"no mesh or generator with the name {mesh_name} found.")
        ngmsh = mesh_function(h_max)
        ngmsh.Save(vol_file)
        print(f"saved mesh to {vol_file}")
    mesh = Mesh(ngmsh)
    return mesh

#def load_or_generate_mesh( mesh_name, h_max=0.4 ):
#    mesh_function = globals().get(f"generate_{mesh_name}")
#    if mesh_function is None:
#        raise ValueError(f"no mesh or generator with the name {mesh_name} found.")
#    ngmsh = mesh_function(h_max)
#    mesh = Mesh(ngmsh)
#    return mesh
    

if __name__ == "__main__":
    msh = load_or_generate_mesh( "ng_square_circle1_exact", h_max=0.4 )
    
    