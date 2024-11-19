# Navier-Stokes parameters
Reynolds_number = 100
Inflow_velocity = 10**(-2)

# Advection-Diffusion parameters
Diffusion_coefficient = 1e-4
m_k = 1.0 / 12.

# AMR parameters
method_name = "Uniform"
marking_strategy = "uniform"
error_estimator = "residual_based2"
tolerance = 1e-16
refinement_iterations = 5

# Initial Mesh parameters
mesh_name = "ng_rect_circ1"
h_max = 1

# Boundary conditions parameters
bc_name = "ng_rect_circ1_bc"
