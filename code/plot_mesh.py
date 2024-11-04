import matplotlib.pyplot as plt
from firedrake import *
import netgen
from netgen.geom2d import SplineGeometry

def plot_mesh(mesh):
    fig, axes = plt.subplots()
    triplot(mesh, axes=axes)
    plt.gca().set_aspect('equal', adjustable='box')
    axes.legend()
    axes.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.show()
    

mesh_name = "ng_square_circle1"
h_max = 0.01
vol_file = f"meshes/{mesh_name}_{str(h_max).replace('.', '')}.vol"
ngmsh = netgen.meshing.Mesh()
ngmsh.Load(vol_file)
labels = [i+1 for i, name in enumerate(ngmsh.GetRegionNames(codim=1)) if name in ["line","curve"]]
print(labels)
msh = Mesh(ngmsh)
plot_mesh(msh)