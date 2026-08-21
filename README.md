# LDT_for_CL
The code attached was used to produce the figures for the 'Large Deviation Theory for Complex Langevin' submission to Physical Review D
sde_solver.m solves an SDE using the Euler-Heun method, implementing an adaptive step size
Complex_Langevin_Kernelled utilises that solver for a Monte Carlo Complex Langevin Simulation
Expectation values of different observables may then be calculated. In particular, boundary_terms.m finds boundary terms of different observables
Lefschetz_thimble.m can be used to find Lefschetz thimbles for a given action
Quasipotential_final.m can be used to solve the Hamilton-Jacobi PDE for the Quasipotential using the ray tracing method
If instead a Quasipotential at a particular point is desired, Quasipotential_two_ray.m should be employed 
Finally, best_kernel.m implements a direct scan of different kernels, plotting the Quasipotential against the kernel phase 
