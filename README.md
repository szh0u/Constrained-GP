# Constrained GP model
Function estimation under shape constraints using Gaussian Process model.
# Model & Constraints:

Observations: 
    
    {x_i, y_i}, i = 1,...,n.

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
   
   
# R functions of constrained GP model used in the simulation example 

     cGP:  GP model with constraints (1), (2), (3) ;
     
     c0GP: GP model with constraints (1) ;
     
     c1GP: GP model with constraints (2), (3) ;
     
     uGP:  GP model with no constraint.
     

# Simulation example

  
     ## simulate data ##
     a1 = -0.1176 
     D = function(x) 1/(1-a1*x/2)^2 
     n = 250
     x = runif(n, min = 0, max = 10)
     sigma = 0.002
     err = sigma*rnorm(n,0,1)
     y = D(x) + err
     N = floor(n/4)+1   # number of knots  
     Niter = 500
     ## basis matrix ##
     C = max(x)
     u = seq(0,C, length.out = N)
     dN = C/(N-1)
     Phi_x = matrix(nrow = n, ncol = N)
     for(i in 1:n){
     Phi_x[i,1] = psi_1(x[i])
     Phi_x[i,N] = psi_N(x[i])
     for(j in 2:(N-1)){Phi_x[i,j] = psi(x[i],j)}}
     Phi = cbind(as.matrix(rep(1,n),nrow=n), x, Phi_x)

    ## transform coeff matrix ##
    max_phi = c(dN/2, rep(dN, N-2), dN/2)
    trans_mat = matrix(nrow = N+1, ncol = N+1)
    trans_mat[1,] = -c(1, max_phi)
    trans_mat[-1,] = cbind(as.matrix(rep(0,N),ncol=1),diag(N)) 
    trans_mat1 = as.matrix(bdiag(1, trans_mat)) # transform matrix for c1GP

    # cGP #
    r0 = constrGP3(x, y,  nu = nu,  l,
               niter = Niter, u = u, Phi = Phi, trans_mat = trans_mat)

    # c0GP #
    r0_n = constrGP_n0(x, y,  nu = nu , l,
                   niter = Niter, u = u, Phi = Phi, trans_mat = trans_mat)

    # c1GP #
    r_n1 = constrGP_n1(x, y,  nu = nu , l,
                   niter = Niter, u = u, Phi = Phi, trans_mat = trans_mat1)

    # uGP #
    ru_n = constrGP_n(x, y,  nu = nu , l,
                  niter = Niter, u = u, Phi = Phi, trans_mat = trans_mat)



        


