# Constrained GP model
Function estimation under shape constraints using Gaussian Process model.
# Model & Constraints:

Observations: 
    
    {x_i, y_i}, i = 1,...,n.

Model: 

     y_i = f(x_i) + u_i,  u_i i.i.d. N(0, sig), where sig is set to be unknown.

Constraints: 

             (1) f(0) = 1;

             (2) f'(x) < 0;
             
             (3) f''(x) > 0.
             
Priors: 

    f|x ~ GP(0, tau*K(x,x', nu, l)), K matern kernel with length-scale parameter l and smoothness parameter nu.

    p(tau) = 1/tau.
   
    p(sig^2) = 1/sig^2.
   
   
# R functions of constrained GP model used in the simulation example 

     cGP:  GP model with constraints (1), (2), (3) ;
     
     c0GP: GP model with constraints (1) ;
     
     c1GP: GP model with constraints (2), (3) ;
     
     uGP:  GP model with no constraint.
     

# Simulation example

function: 
         
        f(x) = 1/(1 + 0.1176*x/2)^2. 

setting: 

     n = 200,  x in [0, 10], sig = 0.0002, l = 20, nu = 2.5 
        


