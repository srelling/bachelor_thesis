import os
import datetime
import shutil
import matplotlib.pyplot as plt
from firedrake import *

from settings import *

daytime_string = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

def create_output_directory():
    try:
        os.makedirs("output/output" + daytime_string + "/plots")
        os.makedirs("output/output" + daytime_string + "/vtk")
        shutil.copy("cht_parameters.py", "output/output" + daytime_string + "/cht_parameters.py")
    except FileExistsError:
        pass
    
def plot_mesh(mesh, i):
    fig, axes = plt.subplots()
    triplot(mesh, axes=axes)
    plt.gca().set_aspect('equal', adjustable='box')
    axes.legend()
    axes.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.savefig("output/output" + daytime_string + f"/plots/mesh_{i}.pdf")       if settings["save_plot"]["mesh"] else None
    
def plot_velocity_and_pressure(u, p, i):
    fig, axes = plt.subplots(nrows=2, sharex=True, sharey=True)
    streamlines = streamplot(u, resolution=0.1,  seed=0, axes=axes[0])
    fig.colorbar(streamlines, ax=axes[0], fraction=0.046)
    axes[0].set_aspect("equal")
    axes[0].set_title("Velocity")
    contours = tricontourf(p, 30, axes=axes[1])
    fig.colorbar(contours, ax=axes[1], fraction=0.046)
    axes[1].set_aspect("equal")
    axes[1].set_title("Pressure")
    plt.gca().set_aspect('equal', adjustable='box')
    plt.savefig("output/output" + daytime_string + f"/plots/velocity_pressure_{i}.pdf")       if settings["save_plot"]["velocity_pressure"] else None
    
def plot_temperature(u, i):
    colorplot = tripcolor(u)
    plt.xlabel(r'$x$')
    plt.ylabel(r'$y$')
    plt.colorbar(colorplot)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.savefig("output/output" + daytime_string + f"/plots/temperature_{i}.pdf")       if settings["save_plot"]["temperature"] else None
    
def plot_solution_difference(mesh, bcs, u_sol, i):
    if "U_EXACT" not in bcs or bcs["U_EXACT"] is None:
        return
    u_exact_expr = bcs["U_EXACT"]
    V = u_sol.function_space() 
    u_exact = Function(V).interpolate(u_exact_expr)
    
    u_diff = Function(V)
    u_diff.assign(u_sol - u_exact)
    

    plt.figure()
    colorplot = tripcolor(u_diff)
    plt.xlabel(r'$x$')
    plt.ylabel(r'$y$')
    plt.title("Difference between computed and exact solution")
    plt.colorbar(colorplot)
    plt.savefig("output/output" + daytime_string + f"/plots/solution_difference_{i}.pdf") if settings["save_plot"]["solution_difference"] else None

def save_solution(u, p, T, i):
    if not any(settings["save_vtk"].values()):
        return
    VTKFile("output/output" + daytime_string + "/vtk/velocity.pvd", mode="a", adaptive=True).write(u)         if settings["save_vtk"]["velocity"] else None
    VTKFile("output/output" + daytime_string + "/vtk/pressure.pvd", mode="a", adaptive=True).write(p)         if settings["save_vtk"]["pressure"] else None
    VTKFile("output/output" + daytime_string + "/vtk/temperature.pvd", mode="a", adaptive=True).write(T)      if settings["save_vtk"]["temperature"] else None
    
    print("Solution saved to output/output" + daytime_string + "/vtk")
    