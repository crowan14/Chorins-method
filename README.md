# Chorins-method

Numerical solution to incompressible Navier-Stokes equations on square domain with velocity boundary conditions. Pressure Poisson problem is solved with the finite element method and the velocity field is updated with the finite difference method. Chorin's method is the numerical technique to enforce the incompressibility condition on the flow, and is an operator splitting approach.

The finite element implementation of the pressure problem uses a MATLAB toolbox called calfem. Calfem can be downloaded at: https://github.com/CALFEM/calfem-matlab/tree/master/fem. 
