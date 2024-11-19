from firedrake import *
from firedrake.petsc import PETSc
import netgen
from parameters import *
from cht_error_indicator import *

def Mark_naive(msh, uh):
    W = FunctionSpace(msh, "DG", 0)
    mark = Function(W)

    # Assembling the error indicator.
    eta = get_error_indicator(msh, uh)
    
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

def Mark_uniform( msh, uh ):
    W = FunctionSpace(msh, "DG", 0)
    mark = Function(W)
    x, y = SpatialCoordinate(msh)
    # Fill in all of mark with 1
    with mark.dat.vec as markVec:
        markVec.set(1.0)
    return mark


def Mark( msh, uh ):
    marking_function = globals().get("Mark_" + marking_strategy)
    if marking_function is None:
        raise ValueError(f"no marking strategy with the name {marking_function} found.")
    return marking_function(msh, uh)