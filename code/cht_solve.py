from firedrake import *
from firedrake.petsc import PETSc
from firedrake.__future__ import interpolate
import matplotlib.pyplot as plt
from matplotlib import rc
import ufl

import datetime
import shutil
import os

from cht_boundary_conditions import *
from cht_meshes import *
from cht_parameters import *
from cht_plot import *
from settings import *

# ======================NAVIER-STOKES======================
def solveNS( msh, bcs ):
    V = VectorFunctionSpace( msh, "CG", 2 )
    W = FunctionSpace( msh, "CG", 1 )
    Z = V * W
    up = Function( Z )
    u, p = split( up )
    v, q = TestFunctions( Z )
    Re = Constant( Reynolds_number )
    F = 1.0 / Re * inner( grad( u ), grad( v ) ) * dx \
                              + dot( grad( u ) * u, v ) * dx \
                              - p * div( v ) * dx \
                              + div( u ) * q * dx
    x, y = SpatialCoordinate( msh )
    ns_bcs = []
    for name, (func, tags) in bcs["NS_BCS"]:
            bc = DirichletBC( Z.sub( 0 ) , func , tags )
            ns_bcs.append( bc )

    #nullspace = MixedVectorSpaceBasis(
    #    Z, [Z.sub(0), VectorSpaceBasis(constant=True)])
    #appctx = {"Re": Re, "velocity_space": 0}
    #parameters = {"mat_type": "matfree",
    #              "snes_monitor": None,
    #              "ksp_type": "fgmres",
    #              "ksp_gmres_modifiedgramschmidt": None,
    #              "ksp_monitor_true_residual": None,
    #              "pc_type": "fieldsplit",
    #              "pc_fieldsplit_type": "schur",
    #              "pc_fieldsplit_schur_fact_type": "lower",
    #              "fieldsplit_0_ksp_type": "preonly",
    #              "fieldsplit_0_pc_type": "python",
    #              "fieldsplit_0_pc_python_type": "firedrake.AssembledPC",
    #              "fieldsplit_0_assembled_pc_type": "lu",
    #              "fieldsplit_1_ksp_type": "gmres",
    #              "fieldsplit_1_ksp_rtol": 1e-4,
    #              "fieldsplit_1_pc_type": "python",
    #              "fieldsplit_1_pc_python_type": "firedrake.PCDPC",
    #              "fieldsplit_1_pcd_Mp_ksp_type": "preonly",
    #              "fieldsplit_1_pcd_Mp_pc_type": "lu",
    #              "fieldsplit_1_pcd_Kp_ksp_type": "preonly",
    #              "fieldsplit_1_pcd_Kp_pc_type": "lu",
    #              "fieldsplit_1_pcd_Fp_mat_type": "matfree"}
    #
    #up.assign(0)
    #solve(F == 0, up, bcs=ns_bcs, nullspace=nullspace, solver_parameters=parameters,
    #        appctx=appctx)

        
    solve( F == 0, up, bcs=ns_bcs, 
            solver_parameters={
                "snes_monitor": None,
            },
        )
    u_init, p_init = up.subfunctions
    return u_init, p_init

# ====================ADVECTION-DIFFUSION====================
def solveAD( msh, bcs, u_init ):
    V = FunctionSpace( msh , "CG" , 2 )
    x, y = SpatialCoordinate( msh )
    boundary_conditions = []
    for name, (func, tags) in bcs["AD_BCS"]:
            bc = DirichletBC( V , func , tags )
            boundary_conditions.append( bc )
    f = bcs["AD_SOURCE"]
    u = TrialFunction( V )
    v = TestFunction( V )
    k = Constant( Diffusion_coefficient )  # Diffusion
    b = u_init
    a = k * inner( grad( v ) , grad( u ) ) * dx + dot( b , grad( u ) ) * v * dx
    L = f * v * dx
    h_k = sqrt( 2 ) * CellVolume( msh ) / CellDiameter( msh )
    b_norm = norm(b)
    Pe_k = m_k * b_norm * h_k / ( 2.0 * k )
    one = Constant( 1.0 )
    eps_k = conditional( gt( Pe_k , one ) , one , Pe_k )
    tau_k = h_k / ( 2.0 * b_norm ) * eps_k
    a += inner( ( dot( b , grad( u ) ) - k * div( grad( u ) ) ) , tau_k * dot( b , grad( v ) ) ) * dx
    L += f * tau_k * dot( b , grad( v ) ) * dx
    u_sol = Function( V  , name="u_sol" )
    problem = LinearVariationalProblem( a , L , u_sol , boundary_conditions )
    solver = LinearVariationalSolver( problem )
    solver.solve()
    return u_sol

# ====================CHT====================
def solveCHT( msh, bcs, dir, iteration=0 ):
    plot_mesh( msh, dir, iteration )                                 if settings["save_plot"]["mesh"] else None
    print( "Solving Navier-Stokes equations..." )                    
    u_init, p_init = solveNS( msh, bcs )
    print( "Navier-Stokes equations solved." )
    print( "Solving Advection-Diffusion equations..." )
    plot_velocity_and_pressure( u_init, p_init, dir, iteration )     if settings["save_plot"]["velocity_pressure"] else None
    print( "Advection-Diffusion equations solved." )
    u_sol = solveAD( msh , bcs, u_init )
    plot_temperature( u_sol, dir, iteration )                        if settings["save_plot"]["temperature"] else None
    plot_solution_difference( msh, bcs, u_sol, dir, iteration )      if settings["save_plot"]["solution_difference"] else None
    save_solution( u_init, p_init, u_sol, dir, iteration )
    return u_sol

# ======================MAIN======================
if __name__ == "__main__":
    msh = load_or_generate_mesh( "ng_square_circle1_exact", h_max=0.4 )
    msh.name = "ng_square_circle1_exact"
    bcs = get_bcs( msh, "ng_square_circle1_bc" )
    create_output_directory()
    u_sol = solveCHT( msh, bcs )
    with CheckpointFile("solutions/ng_square_circle1_001/temperature.h5", 'w') as afile:
        afile.save_mesh(msh)
        afile.save_function(u_sol)
    
    