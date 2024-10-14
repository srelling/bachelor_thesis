# bachelor_thesis
Adaptive Mesh Refinement for the Conjugate Heat Transfer Equations Using Firedrake
# TODOLIST
- [x] Implement SolveCHT function
- [ ] Simple uniform refinement implementation
- [ ] Make different meshes to try out
- [ ] Set up test problem/mesh for which i can compare to a analytical solution
- [ ] More Elegant swapping between meshes with different number of boundries
- [ ] Add inflow velocity and source function to input parameters as firedrake functions
# Problems
- [ ] Reynolds number only influences the pressure field and npt the velocity field so it doesn't make a difference for the final temperature field, does this make sense?
- [ ] The magnitude of the velocity field does seem to make only a tiny difference in the final temperature field
- [ ] maybe I should use neuman boundry conditions on the surfaces of the domain to simulate the exchange of heat on the domain? But then the exchange of heat would be fix and not dependent on the flow of the liqiud?
# Bugs