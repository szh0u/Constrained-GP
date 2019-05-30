## Constrained GP model
Function estimation under shape constraints using Gaussian Process (Physical Review C: https://journals.aps.org/prc/abstract/10.1103/PhysRevC.99.055202).

## Model & Constraints

Model: 

     y_i = f_1(x_i) + e_i,  e_i iid N(0, sigma), where sigma is set to be unknown.
     
     f_1(x) = f(0) + f'(0)*x + \sum_{j=0}^n f''(u_j)*h_j(x), with pre-defined knots {u_j} and basis functions {h_j}.
     
     Denote parameters F = [f(0), f'(0), f''(u_0), ..., f''(u_n)] to be updated.

Constraints: 

    (1) f(0) = 1;

    (2) f'(x) < 0;
             
    (3) f''(x) > 0.
             
Priors: 

    F|x ~ GP(0, tau*K(x,x', nu, l)), K matern kernel with length-scale parameter l and smoothness parameter nu.

    p(tau) = 1/tau.
   
    p(sigma^2) = 1/sigma^2.
   
   
## R functions used in the simulation example 

     maternCov: generate the matern covariance matrix with arbitrary smoothness parameter nu and length-scale parameter l ; 

     cGP:  GP model with constraints (1), (2), (3) ;
     
     c0GP: GP model with constraints (1) ;
     
     c1GP: GP model with constraints (2), (3) ;
     
     uGP:  GP model with no constraint.
     
     
     

## Output

![szh0u\Constrained-GP](simu_plot.png)

red dashed line -- true function; 
blue solid line -- cGP (darker shade 95% CI);
green solid line -- c0GP (light shade 95% CI);
purple solid line -- uGP (purple dashed lines 95% CI).



