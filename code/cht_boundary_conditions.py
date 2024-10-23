from firedrake import *
import ufl
from cht_parameters import *

def get_square_exact_bc(msh):
    x, y = SpatialCoordinate(msh)
    square_exact_bc = {
        "NS_BCS": [
            ("inflow", ([(-(y**2) + 1) * Inflow_velocity, 0], (12))),
            ("no_slip", ([0, 0], (11)))
        ],
        "AD_BCS": [
            ("left", (0.25 * y**2 + 0.75, (12))),
            ("hor", (0.25 * x + 1.25, (11))),
            ("right", (0.25 * y**2 + 1.25, (13)))
        ],
        "AD_SOURCE": (-0.0025 * y**2 + 0.0025) - 0.5 * Diffusion_coefficient,
        "U_EXACT": 0.25 * x + 0.25 * y**2 + 1
    }
    return square_exact_bc

def get_bcs( msh, bc_name ):
    bc_function = globals().get(f"get_{bc_name}")
    if bc_function is None:
        raise ValueError(f"no boundry conditions with the name {bc_name} found.")
    return bc_function(msh)
