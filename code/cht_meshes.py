from firedrake import *
import netgen
from netgen.geom2d import SplineGeometry
import matplotlib.pyplot as plt
import ufl


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
    msh = Mesh(ngmsh)
    return msh

def generate_mesh( mesh_name, h_max=0.4 ):
    mesh_function = globals().get(f"generate_{mesh_name}")
    if mesh_function is None:
        raise ValueError(f"no mesh with the name {mesh_name} found.")
    return mesh_function( h_max )

if __name__ == "__main__":
    mesh = generate_ng_square_circle1( 0.4 )
    fig, axes = plt.subplots()
    triplot(mesh, axes=axes)
    plt.gca().set_aspect('equal', adjustable='box')
    axes.legend()
    plt.show()
    
    