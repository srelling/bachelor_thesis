import os
import datetime
import shutil
import matplotlib.pyplot as plt
from firedrake import *

from settings import *

daytime_string = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

def create_output_directory(method):
    try:
        os.makedirs("output/" + daytime_string + "/" + method + "/plots")
        os.makedirs("output/" + daytime_string + "/" + method + "/vtk")
        shutil.copy("cht_parameters.py", "output/" + daytime_string + "/cht_parameters.py")
    except FileExistsError:
        pass
    output_dir = "output/" + daytime_string + "/" + method
    return output_dir
    
def plot_mesh(mesh, dir, i):
    fig, axes = plt.subplots()
    triplot(mesh, axes=axes)
    plt.gca().set_aspect('equal', adjustable='box')
    axes.legend()
    axes.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.savefig(f"{dir}/plots/mesh_{i}.pdf")
    
def plot_velocity_and_pressure(u, p, dir, i):
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
    plt.savefig(f"{dir}/plots/velocity_pressure_{i}.pdf")       
    
def plot_temperature(u, dir, i):
    colorplot = tripcolor(u)
    plt.xlabel(r'$x$')
    plt.ylabel(r'$y$')
    plt.colorbar(colorplot)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.savefig(f"{dir}/plots/temperature_{i}.pdf")
    
def plot_solution_difference(mesh, bcs, u_sol, dir, i):
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
    plt.savefig(f"{dir}/plots/solution_difference_{i}.pdf")

def plot_error_convergence(n_dofs_list, errors_list, names):
    plt.close("all")
    for n_dofs, errors, name in zip(n_dofs_list, errors_list, names):
        plt.loglog(n_dofs, errors, "-o", label=name)
    plt.xlabel("Number of DoFs")
    plt.ylabel("Error")
    plt.legend()
    plt.grid()
    plt.savefig("output/" + daytime_string + "/convergence.pdf")

def save_solution(u, p, T, dir, i):
    if not any(settings["save_vtk"].values()):
        return
    VTKFile(dir + "/vtk/velocity.pvd", mode="a", adaptive=True).write(u)         if settings["save_vtk"]["velocity"] else None
    VTKFile(dir + "/vtk/pressure.pvd", mode="a", adaptive=True).write(p)         if settings["save_vtk"]["pressure"] else None
    VTKFile(dir + "/vtk/temperature.pvd", mode="a", adaptive=True).write(T)      if settings["save_vtk"]["temperature"] else None
    
    print("Solution saved to " + dir + "/vtk")
    