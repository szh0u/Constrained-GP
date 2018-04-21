# Constrained GP model
Function estimation under shape constraints using Gaussian Process model. 

# Model & Constraints

Model: 

     y_i = f(x_i) + u_i,  u_i i.i.d. N(0, sigma), where sigma is set to be unknown.

Constraints: 

             (1) f(0) = 1;

             (2) f'(x) < 0;
             
             (3) f''(x) > 0.
             
Priors: 

    f|x ~ GP(0, tau*K(x,x', nu, l)), K matern kernel with length-scale parameter l and smoothness parameter nu.

    p(tau) = 1/tau.
   
    p(sigma^2) = 1/sigma^2.
   
   
# R functions used in the simulation example 

     maternCov: generate the matern covariance matrix with arbitrary smoothness parameter nu and length-scale parameter l ; 

     cGP:  GP model with constraints (1), (2), (3) ;
     
     c0GP: GP model with constraints (1) ;
     
     c1GP: GP model with constraints (2), (3) ;
     
     uGP:  GP model with no constraint.
     

# Output 

!(Constrained-GP/example.png)




