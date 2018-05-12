#########################################################
#
# GP model incorporating constraints (1), (2), (3).
#
# Input: x:         Covariate, a n*1 vector. 
#        y:         Response, a n*1 vector.  
#        nu:        Smoothness parameter in matern kernel. 
#        l:         Length-scale parameter in matern kernel.
#        niter:     Iteration numvers of MC. 
#        u:         Knots of basis functions.
#        Phi:       Basis function evaluation, a n*(N+2) matrix.
#        trans_mat: tansformation matrix used in the constrained GP methods.
#
# Output: f:        mcmc samples of parameters F.  
#         tau:      mcmc samples of signal-to-noise level in matern kernel. 
#         sigma:    mcmc samples of noise level.
#
##########################################################


cGP = function(x, y, nu, l, niter, u, Phi, trans_mat){

  n=length(x)
  
  mcov1 = maternCov(u, l, nu = nu)
  mcov0 = mcov1[-1,]
  mcov = mcov0[,-1]
  D = eigen(mcov)$values
  U = eigen(mcov)$vectors
  for(i in 1:length(D)){ if(D[i] < 0) D[i] = -D[i]}
  mcov = U%*%diag(D)%*%t(U)
  mcov_inv = chol2inv(chol(mcov))
  
  # set initial values of [f, tau, sigma] # 
  tau = c()
  tau[1] = 0.01
  sigma = c()
  sigma[1] = 0.01
  f = matrix(nrow = N+1, ncol = niter)
  f[,1] = t(chol(mcov))%*%rnorm(N+1, 0, 1)*sqrt(tau[1])
  
  L = rep(0, N+1) 
  U = rep(Inf, N+1)
  re = list()
  
  # mcmc # 
  for(i in 2:niter){
    
    if (i%%10 == 0) print(i)
    
    # update f #
    Sig = tau[i-1]*t(Phi) %*% Phi + sigma[i-1]*mcov_inv
    Siginv = tau[i-1]*sigma[i-1]*chol2inv(chol(Sig))
    postmean = Siginv%*%t(Phi)%*%y/sigma[i-1]

    tran_Sig = trans_mat%*%Siginv%*%t(trans_mat) 
    f_trans = trans_mat %*% postmean + mvrandn((L-trans_mat%*%postmean), U, tran_Sig, 1)
    f_cond = solve(trans_mat, f_trans)
    f[,i] = as.matrix(f_cond, ncol = 1)
    
    # update sigma #
    sigma[i] = rinvgamma(1, shape = n/2, rate = t(y-Phi%*%f[,i])%*%(y-Phi%*%f[,i])/2)
    
    # update tau #
    tau[i] = rinvgamma(1,shape = (N+1)/2, rate = t(f[,i])%*%mcov_inv%*%f[,i]/2)
  }
  
  re = list(f, tau, sigma)
  
  return(re)
  
}



