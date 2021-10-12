#------------------------------------------------------------------------------#
# File:        kalman.standard.errors.R
#
# Description: This file computes confidence intervals and corresponding 
#              standard errors for the estimates of the states using
#              Hamilton's (1986) Monte Carlo procedure that accounts for
#              both filter and parameter uncertainty. See footnote 7 in HLW.
#------------------------------------------------------------------------------#
kalman.standard.errors <- function(T, states, theta, y.data, x.data, stage,
                                   lambda.g, lambda.z, xi.00, P.00, niter = 5000,
                                   a3.constraint=NA, b2.constraint=NA) {

    print('Computing Standard Errors')
    
    # Set a3.constraint to -0.0025 if a constraint is not specified in stage 3
    if (is.na(a3.constraint)) {
        a3.constraint <- -0.0025
    }
    # Set b2.constraint to 0.025 if a constraint is not specified in stage 3
    if (is.na(b2.constraint)) {
        b2.constraint <- 0.025
    }

    print("Standard Error Procedure: a3.constraint")
    print(a3.constraint)
    
    print("Standard Error Procedure: b2.constraint")
    print(b2.constraint)

    n.params <- length(theta)
    n.state.vars <- length(xi.00)

    # Return vector of log likelihood values at each time t 
    log.likelihood.estimated.vector <- log.likelihood.wrapper(theta, y.data, x.data, stage = 3,lambda.g=lambda.g, lambda.z=lambda.z, xi.00=xi.00, P.00=P.00)$ll.vec
    stin <- states$smoothed$xi.tT[1,] # First smoothed state vector
    pp1  <- states$filtered$P.ttm1[1:n.state.vars,] # First covariance matrix
    eigenstuff.pp1   <- eigen(pp1)
    eigenvectors.pp1 <- eigenstuff.pp1$vectors # Eigenvectors of first covariance matrix
    # Eigenvectors without a positive first entry are multiplied by -1 to ensure
    # consistency across different versions of R, which choose the sign differently
    for (l in 1:n.state.vars) {
        if (eigenvectors.pp1[1,l] < 0 ) { eigenvectors.pp1[,l] <- -eigenvectors.pp1[,l] }
    } 
    eigenvalues.pp1  <- eigenstuff.pp1$value   # Eigenvalues of first covariance matrix
    dg   <- diag(x = eigenvalues.pp1) 
    hh2  <- eigenvectors.pp1 %*% sqrt(dg)    

    # Compute information matrix from difference in gradients of the likelihood function
    # from varying theta (parameter vector) values
    likelihood.gradient <- matrix(NA,T,n.params)
    for (i in 1:n.params){
        delta   <- max(theta[i]*1e-6, 1e-6)
        d.theta <- theta
        d.theta[i] <- theta[i] + delta
        likelihood.gradient[,i] <-  (log.likelihood.wrapper(d.theta, y.data, x.data, stage = 3,lambda.g=lambda.g, lambda.z=lambda.z, xi.00=xi.00, P.00=P.00)$ll.vec - 
                                         log.likelihood.estimated.vector)/delta
    }
    info <- solve(t(likelihood.gradient) %*% likelihood.gradient) # Information matrix
    bse <- sqrt(diag(info))
    t.stats <- abs(theta) / bse

    # Smoothed estimates
    g      <- 4 * states$smoothed$xi.tT[,4]
    ypot   <- states$smoothed$xi.tT[,1]
    z      <- states$smoothed$xi.tT[,6]
    rstar  <- g + z

    # cum1 cumulates terms for parameter uncertainty;
    # cum2 cumulates terms for filter uncertainty
    cum1 <- matrix(0,T,3)
    cum2 <- matrix(0,T,3)
    eigenstuff.info   <- eigen(info)
    eigenvectors.info <- eigenstuff.info$vectors # Eigenvectors of information matrix
    # Eigenvectors without a positive first entry are multiplied by -1 to ensure
    # consistency across different versions of R, which choose the sign differently
    for (l in 1:n.params) {
        if (eigenvectors.info[1,l] < 0 ) { eigenvectors.info[,l] <- -eigenvectors.info[,l] }
    }
    eigenvalues.info  <- eigenstuff.info$value # Eigenvalues of information matrix
    dg <- diag(x = eigenvalues.info)
    hh <- eigenvectors.info %*% sqrt(dg)

    set.seed(50)

    # Store the number of draws excluded for violating constraints
    good.draws                 <- 0
    excluded.draw.counter      <- 0
    excluded.draw.counter.a3   <- 0
    excluded.draw.counter.b2   <- 0
    excluded.draw.counter.a1a2 <- 0

    # See HLW footnote 7 for description of procedure
    # niter is the number of iterations; we discard draws that violate constraints
    while (good.draws < niter) {
      theta.i <- (hh %*% rnorm(n.params) + theta)[,1] 
      if ( (theta.i[3] <= a3.constraint) & (theta.i[5] >= b2.constraint) & (theta.i[1] + theta.i[2] < 1) ) {
          xi.00.i  <- c(t(hh2 %*% rnorm(n.state.vars) + stin)) 
          states.i <- kalman.states.wrapper(theta.i, y.data, x.data, stage, lambda.g, lambda.z, xi.00.i, pp1)

          g.i    <- 4 * states.i$smoothed$xi.tT[,4]
          ypot.i <- states.i$smoothed$xi.tT[,1]
          z.i    <- states.i$smoothed$xi.tT[,6]
          r.i    <- g.i + z.i          
                    
          cum1[,1] <- cum1[,1]+(ypot.i-ypot)^2
          cum1[,2] <- cum1[,2]+(r.i-rstar)^2
          cum1[,3] <- cum1[,3]+(g.i-g)^2

          P.ttm1.i   <- states.i$smoothed$P.tT
          P.ttm1.i.f <- states.i$filtered$P.tt
          for (j in 1:(T-1)){
              cum2[j,1]  <- cum2[j,1] + P.ttm1.i[(j * n.state.vars +1),1]
              cum2[j,2]  <- cum2[j,2] + 16 * P.ttm1.i[(j*n.state.vars+4),4] + P.ttm1.i[(j*n.state.vars+6),6]
              cum2[j,3]  <- cum2[j,3] + P.ttm1.i[(j*n.state.vars+4),4]
          }
          cum2[T,1] <- cum2[T,1] + P.ttm1.i.f[((T-1)*n.state.vars+1),1]
          cum2[T,2] <- cum2[T,2] + (16 * P.ttm1.i.f[((T-1)*n.state.vars+4),4] + P.ttm1.i.f[((T-1)*n.state.vars+6),6])
          cum2[T,3] <- cum2[T,3] + P.ttm1.i.f[((T-1)*n.state.vars+4),4]
          good.draws <- good.draws + 1
      } else {
          excluded.draw.counter <- excluded.draw.counter + 1
          if (theta.i[3] > a3.constraint) {
              excluded.draw.counter.a3 <- excluded.draw.counter.a3 + 1
          }
          if (theta.i[5] < b2.constraint) {
              excluded.draw.counter.b2 <- excluded.draw.counter.b2 + 1
          }
          if ((theta.i[1] + theta.i[2]) >= 1) {
              excluded.draw.counter.a1a2 <- excluded.draw.counter.a1a2 + 1
          }
      }
    }
    cum1 <- cum1/niter # Measure of parameter uncertainty
    cum2 <- cum2/niter # Measure of filter uncertainty
    cum2[,3] <- 16*cum2[,3] # Variance for growth at an annualized rate

    # Standard errors for estimates of the states
    # Order: y*, r*, g
    se <- sqrt(cum1 + cum2)

    rm(.Random.seed)
    return(list("se.mean"=colMeans(se),
                "se"=se,"t.stats"=t.stats,"bse"=bse,
                "number.excluded"=excluded.draw.counter,
                "number.excluded.a3"=excluded.draw.counter.a3,
                "number.excluded.b2"=excluded.draw.counter.b2,
                "number.excluded.a1a2"=excluded.draw.counter.a1a2))
}    
