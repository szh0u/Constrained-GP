# Constrained GP model
Function estimation under shape constraints using Gaussian Process model.
# Model & Constraints:

Observations on f: {x_i, y_i}, i =1,...,n.

Model: y_i = f(x_i) + u_i,  u_i i.i.d. N(0, sig), sig unknown.

Constraints: 

             (1) f(0) = 1;

             (2) f'(x) < 0;
             
             (3)f''(x) > 0.
             
Priors: 
    f|x ~ GP(0, tau K(x,x', nu, l)), K matern kernel with length-scale parameter l and smoothness parameter nu.

   p(tau) = 1/tau.
   
   p(sig^2) = 1/sig^2.

# Simulation examples 

function: f(x) = 

# R functions of constrained GP model

cGP: GP with 

