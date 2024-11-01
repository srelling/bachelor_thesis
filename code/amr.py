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
from cht_meshes import *
from cht_boundary_conditions import *
from cht_parameters import *
from cht_plot import *
from settings import *

import netgen
from netgen.geom2d import SplineGeometry
from netgen.geom2d import CSG2d, Circle, Rectangle

tolerance = 1e-16


def Mark(msh, uh):
    W = FunctionSpace(msh, "DG", 0)
    V = FunctionSpace(msh, "CG", 2)
    w = TestFunction(W)
    R_T = uh + div(grad(uh))
    n = FacetNormal(V.mesh())
    h = CellDiameter(msh)
    R_dT = dot(grad(uh), n)
    # Assembling the error indicator.
    eta = assemble(h**2*R_T**2*w*dx +
        (h("+")+h("-"))*(R_dT("+")-R_dT("-"))**2*(w("+")+w("-"))*dS)
    frac = .95
    delfrac = .05
    part = .2
    mark = Function(W)
    # Filling in the marked element vector using eta.
    with mark.dat.vec as markedVec:
        with eta.dat.vec as etaVec:
            sum_eta = etaVec.sum()
            if sum_eta < tolerance:
                return markedVec
            eta_max = etaVec.max()[1]
            sct, etaVec0 = PETSc.Scatter.toZero(etaVec)
            markedVec0 = etaVec0.duplicate()
            sct(etaVec, etaVec0)
            if etaVec.getComm().getRank() == 0:
                eta = etaVec0.getArray()
                marked = np.zeros(eta.size, dtype='bool')
                sum_marked_eta = 0.
                #Marking strategy
                while sum_marked_eta < part*sum_eta:
                    new_marked = (~marked) & (eta > frac*eta_max)
                    sum_marked_eta += sum(eta[new_marked])
                    marked += new_marked
                    frac -= delfrac
                markedVec0.getArray()[:] = 1.0*marked[:]
            sct(markedVec0, markedVec, mode=PETSc.Scatter.Mode.REVERSE)
    return mark  
    
    
if __name__ == "__main__":
    msh = load_or_generate_mesh( "ng_square_circle1", h_max=0.4 )    
    bcs = get_bcs( msh, "ng_square_circle1_bc" )
    create_output_directory()
    
    with CheckpointFile("solutions/ng_square_circle1_001/temperature.h5", 'r') as afile:
        msh_sol = afile.load_mesh("ng_square_circle1")
        u_sol = afile.load_function(msh_sol, "u_sol")
    
    V_exact = FunctionSpace(msh_sol, "CG", 2)
    coords = np.array(msh_sol.coordinates.dat.data)
    u_sol_at_nodes = np.array([u_sol.at(coord) for coord in coords])

    
    errors = []
    n_dofs = []
    
    for i in range(5):
        print(f"level {i}")
        uh = solveCHT(msh, bcs, i)
        uh_at_nodes = np.array([uh.at(coord) for coord in coords])
        l2_error = np.linalg.norm(u_sol_at_nodes - uh_at_nodes)
        errors.append(l2_error)
        V = FunctionSpace(msh, "CG", 2)
        n_dofs.append(V.dim())
        mark = Mark(msh, uh)
        msh = msh.refine_marked_elements(mark)
    
    print(errors)
    print(n_dofs)
    plt.clf()
    plt.loglog(n_dofs, errors, "-o")
    plt.xlabel("Number of DoFs")
    plt.ylabel("Error")
    plt.grid()
    plt.show()
        

        
        
        
    


