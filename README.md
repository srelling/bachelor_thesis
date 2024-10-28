# bachelor_thesis
Adaptive Mesh Refinement for the Conjugate Heat Transfer Equations Using Firedrake
# TODOLIST
- [x] Implement SolveCHT function
- [x] More Elegant swapping between meshes with different number of boundries
- [x] Add inflow velocity and source function to input parameters as firedrake functions
- [x] Set up test problem/mesh for which i can compare to a analytical solution
- [x] Add a manufactured solution with a degree 3 polynomial
- [x] Plot the convergence of CHTSolve on uniform refined mesh of exact solution of degree 3
- [ ] Simple uniform refinement implementation
- [ ] Create First Adaptive refinement Error indicator
- [ ] Make different meshes to try out
- [ ] Clean up output of a refinement run
# Problems
- [ ] Reynolds number only influences the pressure field and npt the velocity field so it doesn't make a difference for the final temperature field, does this make sense?
- [ ] The magnitude of the velocity field only seems to make a tiny difference in the final temperature field
- [ ] maybe I should use neuman boundry conditions on the surfaces of the domain to simulate the exchange of heat on the domain? But then the exchange of heat would be fix and not dependent on the flow of the liqiud?
- [ ] I don't understand the fluid dynamics well enough.
# Bugs