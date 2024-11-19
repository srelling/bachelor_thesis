from firedrake import *
from firedrake.petsc import PETSc
from firedrake.__future__ import interpolate
import matplotlib.pyplot as plt
from matplotlib import rc
import ufl
import numpy as np

import datetime
import shutil
import os

from cht_solve import *
from cht_mark import *
from cht_meshes import *
from cht_boundary_conditions import *
from parameters import *
from cht_plot import *
from settings import *

import netgen
from netgen.geom2d import SplineGeometry


if __name__ == "__main__":
    j = 0
    
    # Exact Solution
    with CheckpointFile("solutions/ng_rect_circ1_exact/temperature.h5", 'r') as afile:
        msh_sol = afile.load_mesh("ng_rect_circ1_exact")
        u_sol = afile.load_function(msh_sol, "u_sol")
    
    V_exact = FunctionSpace(msh_sol, "CG", 2)
    coords = np.array(msh_sol.coordinates.dat.data)
    u_sol_at_nodes = np.array([u_sol.at(coord) for coord in coords])

    
    errors = []
    n_dofs = []
    
    dir = create_output_directory(j)
    msh = load_or_generate_mesh( mesh_name, h_max )

    for i in range(refinement_iterations):
        print( f"level {i}" )
        bcs = get_bcs( msh, bc_name )
        uh = solveCHT( msh, bcs, dir, i )
        uh_dest = Function( V_exact ).interpolate( uh )
        l2_error = errornorm( u_sol, uh_dest, norm_type="L2" )
        errors.append( l2_error )
        V = FunctionSpace( msh, "CG", 2 )
        n_dofs.append( V.dim() )
        mark = Mark( msh, uh )
        msh = msh.refine_marked_elements( mark )
        
    save_convergence_data( n_dofs, errors, dir )                            if settings["save_data"]["convergence"] else None
    plot_error_convergence( [n_dofs], [errors] )                            if settings["save_plot"]["convergence"] else None

        
        
        
    


