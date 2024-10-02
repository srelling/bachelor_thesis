from firedrake import *
from firedrake.__future__ import interpolate
import matplotlib.pyplot as plt
from matplotlib import rc
import datetime
import ufl

# Navier-Stokes parameters
Reynolds_number = 10
Inflow_velocity = 10**(-2)

# Advection-Diffusion parameters
Inflow_temperature = 0
Cylinder_temperature = 1
Source_temperature = 0.01
HorEdges_temperature = 0
Diffusion_coefficient = 1e-8
m_k = 1.0 / 12.

def plot_mesh(mesh):
    fig, axes = plt.subplots()
    triplot(mesh, axes=axes)
    plt.gca().set_aspect('equal', adjustable='box')
    axes.legend()
    plt.draw()
    
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
    plt.draw()
    
def plot_temperature(u):
    colorplot = tripcolor(u)
    plt.xlabel(r'$x$')
    plt.ylabel(r'$y$')
    plt.colorbar(colorplot)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.show()
    
def save_solution(u, p, T):
    daytime_string = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M")
    u_file = VTKFile("output/output" + daytime_string + "/velocity.pvd")
    p_file = VTKFile("output/output" + daytime_string + "/pressure.pvd")
    T_file = VTKFile("output/output" + daytime_string + "/temperature.pvd")
    u_file.write(u)
    p_file.write(p)
    T_file.write(T)
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
    u_inflow = as_vector( [((y+2)*(4-(y+2))/2.0)*Inflow_velocity, 0] ) 
    bcs = [
        DirichletBC( Z.sub( 0 ), assemble( interpolate( u_inflow, V ) ), 12 ),
        DirichletBC( Z.sub( 0 ), Constant( (0, 0) ), (11, 14) )
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
    g_left = Constant( Inflow_temperature )
    g_horEdges = Constant( HorEdges_temperature )
    g_cylinder = Constant( Cylinder_temperature )
    bc_left = DirichletBC( V , g_left , 12 )
    bc_horEdges = DirichletBC( V , g_horEdges , 11 )
    bc_cylinder = DirichletBC( V , g_cylinder , 14 )
    dirichlet_condition = [ bc_left , \
                            bc_cylinder, bc_horEdges ]
    f = Constant( Source_temperature )
    u = TrialFunction( V )
    v = TestFunction( V )
    k = Constant( Diffusion_coefficient )  # Diffusion
    b = u_init
    a = k * inner( grad( v ) , grad( u ) ) * dx + dot( b , grad( u ) ) * v * dx
    L = f * v * dx
    h_k = sqrt( 2 ) * CellVolume( msh ) / CellDiameter( msh )
    b_norm = sqrt( dot( b , b ) )
    Pe_k = m_k * b_norm * h_k / ( 2.0 * k )
    one = Constant( 1.0 )
    eps_k = conditional( gt( Pe_k , one ) , one , Pe_k )
    tau_k = h_k / ( 2.0 * b_norm ) * eps_k
    a += inner( ( dot( b , grad( u ) ) - k * div( grad( u ) ) ) , tau_k * dot( b , grad( v ) ) ) * dx
    L += f * tau_k * dot( b , grad( v ) ) * dx
    u_sol = Function( V )
    problem = LinearVariationalProblem( a , L , u_sol , dirichlet_condition )
    solver = LinearVariationalSolver( problem )
    solver.solve()
    return u_sol

def solveCHT( msh, plot=True, save=True ):
    plot_mesh( msh ) if plot else None
    print( "Solving Navier-Stokes equations..." )
    u_init, p_init = solveNS( msh )
    print( "Navier-Stokes equations solved." )
    print( "Solving Advection-Diffusion equations..." )
    plot_velocity_and_pressure( u_init, p_init ) if plot else None
    print( "Advection-Diffusion equations solved." )
    u_sol = solveAD( msh , u_init )
    plot_temperature( u_sol ) if plot else None
    save_solution( u_init, p_init, u_sol ) if save else None
    return u_sol

# ======================MAIN======================
if __name__ == "__main__":
    msh = Mesh('../meshes/immersed_domain.msh')
    u_sol = solveCHT(msh, plot=True, save=False)