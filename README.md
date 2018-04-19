# Constrained GP model
Function estimation under shape constraints using Gaussian Process model.
# Model & Constrstraints:
Model: y = f(x) + u,  u ~ N(0, sig), sig unknown.

Constraints: f(0) = 1;
             f'(x) < 0;
             f''(x) > 0.
             
Priors: 
    f|x ~ GP(0, tau K(x,x', nu, l)), K matern kernel with length-scale parameter l and smoothness parameter nu.

   p(tau) = 1/tau.
   
   p(sig^2) = 1/sig^2.

# Gassian 
