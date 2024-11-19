from firedrake import *
from firedrake.petsc import PETSc
import netgen
from parameters import *

def residual_based1( msh, uh ):
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
    return eta

def residual_based2( msh, uh ):
    W = FunctionSpace(msh, "DG", 0)
    V = FunctionSpace(msh, "CG", 2)
    w = TestFunction(W)
    R_T = uh + div(grad(uh))
    n = FacetNormal(V.mesh())
    h = CellDiameter(msh)
    R_dT = dot(grad(uh), n)
    # Assembling the error indicator.
    eta = assemble((h("+")+h("-"))*(R_dT("+")-R_dT("-"))**2*(w("+")+w("-"))*dS)
    return eta

def get_error_indicator( msh, uh ):
    error_indicator_function = globals().get(error_estimator)
    if error_indicator_function is None:
        raise ValueError(f"no error indicator with the name {error_indicator_function} found.")
    return error_indicator_function(msh, uh)