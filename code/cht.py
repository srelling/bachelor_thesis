from firedrake import *
from firedrake.__future__ import interpolate
import matplotlib.pyplot as plt
from matplotlib import rc
import ufl

import datetime
import shutil
import os

from parameters import *
from settings import *

daytime_string = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

def create_output_directory():
    try:
        os.makedirs("output/output" + daytime_string + "/plots")
        os.makedirs("output/output" + daytime_string + "/vtk")
        shutil.copy("parameters.py", "output/output" + daytime_string + "/parameters.py")
    except FileExistsError:
        pass
    

def plot_mesh(mesh):
    fig, axes = plt.subplots()
    triplot(mesh, axes=axes)
    plt.gca().set_aspect('equal', adjustable='box')
    axes.legend()
    plt.savefig("output/output" + daytime_string + "/plots/mesh.pdf")       if settings["save_plot"]["mesh"] else None
    
def plot_velocity_and_pressure(u, p):
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
    plt.savefig("output/output" + daytime_string + "/plots/velocity_pressure.pdf")       if settings["save_plot"]["velocity_pressure"] else None
    
def plot_temperature(u):
    colorplot = tripcolor(u)
    plt.xlabel(r'$x$')
    plt.ylabel(r'$y$')
    plt.colorbar(colorplot)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.savefig("output/output" + daytime_string + "/plots/temperature.pdf")       if settings["save_plot"]["temperature"] else None
    

def plot_solution_difference(mesh, u_sol):
    # Define the exact solution: u_exact(x, y) = 0.25 * x + 0.25 * y^2
    x, y = SpatialCoordinate(mesh)
    u_exact_expr = 0.25 * x + 0.25 * y**2 + 1
    
    # Interpolate the exact solution into the same function space as u_sol
    V = u_sol.function_space()  # Use the same function space as the computed solution
    u_exact = Function(V).interpolate(u_exact_expr)
    
    # Compute the difference: u_diff = u_sol - u_exact
    u_diff = Function(V)
    u_diff.assign(u_sol - u_exact)
    
    # Plot the difference
    plt.figure()
    colorplot = tripcolor(u_diff)
    plt.xlabel(r'$x$')
    plt.ylabel(r'$y$')
    plt.title("Difference between computed and exact solution")
    plt.colorbar(colorplot)
    plt.savefig("output/output" + daytime_string + "/plots/solution_difference.pdf") if settings["save_plot"]["solution_difference"] else None


def save_solution(u, p, T):
    if not any(settings["save_vtk"].values()):
        return
    u_file = VTKFile("output/output" + daytime_string + "/vtk/velocity.pvd")         if settings["save_vtk"]["velocity"] else None
    p_file = VTKFile("output/output" + daytime_string + "/vtk/pressure.pvd")         if settings["save_vtk"]["pressure"] else None
    T_file = VTKFile("output/output" + daytime_string + "/vtk/temperature.pvd")      if settings["save_vtk"]["temperature"] else None
    u_file.write(u)                                                                  if settings["save_vtk"]["velocity"] else None
    p_file.write(p)                                                                  if settings["save_vtk"]["pressure"] else None
    T_file.write(T)                                                                  if settings["save_vtk"]["temperature"] else None
    print("Solution saved to output/output" + daytime_string + "/")
    

# ======================NAVIER-STOKES======================
def solveNS( msh ):
    V = VectorFunctionSpace( msh, "CG", 2 )
    W = FunctionSpace( msh, "CG", 1 )
    Z = V * W
    up = Function( Z )
    u, p = split( up )
    v, q = TestFunctions( Z )
    Re = Constant( Reynolds_number )
    F = 1.0 / Re * inner( grad( u ), grad( v ) ) * dx \
                              - p * div( v ) * dx \
                              + div( u ) * q * dx
    x, y = SpatialCoordinate( msh )
    u_inflow = as_vector( [( -(y**2) + 1 ) * Inflow_velocity, 0] ) 
    bcs = [
        DirichletBC( Z.sub( 0 ), assemble( interpolate( u_inflow, V ) ), 12 ),
        DirichletBC( Z.sub( 0 ), Constant( (0, 0) ), (11) ) # no Cylinder
        #DirichletBC( Z.sub( 0 ), Constant( (0, 0) ), (11, 14) ) # Cylinder
    ]
    solve( F == 0, up, bcs=bcs, 
            solver_parameters={
                "snes_monitor": None,
            },
        )
    u_init, p_init = up.subfunctions
    return u_init, p_init


# ====================ADVECTION-DIFFUSION====================
def solveAD( msh , u_init ):
    V = FunctionSpace( msh , "CG" , 2 )
    x, y = SpatialCoordinate( msh )
    g_left = 0.25 * y**2 + 0.75
    g_horEdges = 0.25 * x + 1.25
    g_cylinder = Constant( Cylinder_temperature )
    bc_left = DirichletBC( V , g_left , 12 )
    bc_horEdges = DirichletBC( V , g_horEdges , 11 )
    #bc_cylinder = DirichletBC( V , g_cylinder , 14 )
    dirichlet_condition = [ bc_left , bc_horEdges ]
    f = (-0.0025 * y**2 + 0.0025) - 0.5 * Diffusion_coefficient
    u = TrialFunction( V )
    v = TestFunction( V )
    k = Constant( Diffusion_coefficient )  # Diffusion
    b = u_init
    a = k * inner( grad( v ) , grad( u ) ) * dx + dot( b , grad( u ) ) * v * dx
    L = f * v * dx
    h_k = sqrt( 2 ) * CellVolume( msh ) / CellDiameter( msh )
    b_norm = norm(b)
    print( "b_norm = ", b_norm )
    Pe_k = m_k * b_norm * h_k / ( 2.0 * k )
    one = Constant( 1.0 )
    eps_k = conditional( gt( Pe_k , one ) , one , Pe_k )
    tau_k = h_k / ( 2.0 * b_norm ) * eps_k
    #a += inner( ( dot( b , grad( u ) ) - k * div( grad( u ) ) ) , tau_k * dot( b , grad( v ) ) ) * dx
    #L += f * tau_k * dot( b , grad( v ) ) * dx
    u_sol = Function( V )
    problem = LinearVariationalProblem( a , L , u_sol , dirichlet_condition )
    solver = LinearVariationalSolver( problem )
    solver.solve()
    return u_sol

def solveCHT( msh ):
    create_output_directory()
    plot_mesh( msh )
    print( "Solving Navier-Stokes equations..." )
    u_init, p_init = solveNS( msh )
    print( "Navier-Stokes equations solved." )
    print( "Solving Advection-Diffusion equations..." )
    plot_velocity_and_pressure( u_init, p_init )
    print( "Advection-Diffusion equations solved." )
    u_sol = solveAD( msh , u_init )
    plot_temperature( u_sol )
    plot_solution_difference( msh, u_sol )
    save_solution( u_init, p_init, u_sol )
    plt.show()                                                  if settings["show_plots"] else None
    return u_sol

# ======================MAIN======================
if __name__ == "__main__":
    msh = Mesh( '../meshes/square.msh' )
    u_sol = solveCHT( msh )