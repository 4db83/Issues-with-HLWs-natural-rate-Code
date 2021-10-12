#------------------------------------------------------------------------------#
# File:        rstar.stage3.R
#
# Description: This file runs the model in the third stage of the HLW estimation.
#------------------------------------------------------------------------------#
rstar.stage3 <- function(log.output,
                         inflation,
                         real.interest.rate,
                         nominal.interest.rate,
                         lambda.g,
                         lambda.z,
                         a3.constraint=NA,
                         b2.constraint=NA,
                         run.se = TRUE) {

  stage <- 3
  
  # Data must start 4 quarters before the estimation period
  T <- length(log.output) - 4

  # Original output gap estimate
  x.og <- cbind(rep(1,T+4), 1:(T+4))
  y.og <- log.output
  output.gap <- (y.og - x.og %*% solve(t(x.og) %*% x.og, t(x.og) %*% y.og)) * 100

  # Initialization of state vector for Kalman filter using HP trend of log output
  log.output.hp.trend <- hpfilter(log.output,freq=36000,type="lambda",drift=FALSE)$trend
  g.pot <- log.output.hp.trend[(g.pot.start.index):length(log.output.hp.trend)]
  g.pot.diff <- diff(g.pot)
# ---------------------------------------------------------------------------------------------
  # INITIAL CONDITION FOR STATE VECTOR
  xi.00 <- c(100*g.pot[3:1],100*g.pot.diff[2:1],0,0) # their default starting values
#  xi.00 <- c(100*g.pot[3:1],100*g.pot.diff[2:1],2,2)
#  xi.00 <- c(100*g.pot[3:1],0,0,2,2)
# ---------------------------------------------------------------------------------------------
  # IS curve
  y.is <- output.gap[5:(T+4)]
  x.is <- cbind(output.gap[4:(T+3)], output.gap[3:(T+2)],
                (real.interest.rate[4:(T+3)] + real.interest.rate[3:(T+2)])/2,
                rep(1,T))
  b.is <- solve(t(x.is) %*% x.is, t(x.is) %*% y.is)
  r.is <- as.vector(y.is - x.is %*% b.is)
  s.is <- sqrt(sum(r.is^2) / (length(r.is)-(dim(x.is)[2])))

  # Phillips curve
  y.ph <- inflation[5:(T+4)]
  x.ph <- cbind(inflation[4:(T+3)],
                (inflation[3:(T+2)]+inflation[2:(T+1)]+inflation[1:T])/3,
                output.gap[4:(T+3)])
  b.ph <- solve(t(x.ph) %*% x.ph, t(x.ph) %*% y.ph)
  r.ph <- y.ph - x.ph %*% b.ph
  s.ph <- sqrt(sum(r.ph^2) / (length(r.ph)-(dim(x.ph)[2])))
  
  y.data <- cbind(100 * log.output[5:(T+4)],
                  inflation[5:(T+4)])
  x.data <- cbind(100 * log.output[4:(T+3)],
                  100 * log.output[3:(T+2)],
                  real.interest.rate[4:(T+3)],
                  real.interest.rate[3:(T+2)],                  
                  inflation[4:(T+3)],
                  (inflation[3:(T+2)]+inflation[2:(T+1)]+inflation[1:T])/3)

  # Starting values for the parameter vector  
  initial.parameters <- c(b.is[1:3], b.ph[1], b.ph[3], s.is, s.ph, 0.7) 

  # Set an upper and lower bound on the parameter vectors:
  # The vector is unbounded unless values are otherwise specified  
  theta.lb <- c(rep(-Inf,length(initial.parameters)))
  theta.ub <- c(rep( Inf,length(initial.parameters)))

  # Set a lower bound for the Phillips curve slope (b_2) of b2.constraint, if not NA
  # In HLW, b2.constraint = 0.025  
  if (!is.na(b2.constraint)) {
#      print(paste0("Setting a lower bound on b_2 of ",as.character(b2.constraint)))
      if (initial.parameters[5] <  b2.constraint) {
          initial.parameters[5] <- b2.constraint
      }
      theta.lb[5] <- b2.constraint
  }

  # Set an upper bound for the IS curve slope (a_3) of a3.constraint, if not NA
  # In HLW, a3.constraint = -0.0025  
  if (!is.na(a3.constraint)) {
#      print(paste0("Setting an upper bound on a_3 of ",as.character(a3.constraint)))
      if (initial.parameters[3] >  a3.constraint) {
          initial.parameters[3] <- a3.constraint
      }
      theta.ub[3] <- a3.constraint      
  }

  # Set the initial covariance matrix (see footnote 6)   
  P.00 <- calculate.covariance(initial.parameters, theta.lb, theta.ub, y.data, x.data, stage, lambda.g, lambda.z, xi.00)
  
  # Get parameter estimates via maximum likelihood
  f <- function(theta) {return(-log.likelihood.wrapper(theta, y.data, x.data, stage, lambda.g, lambda.z, xi.00, P.00)$ll.cum)}
  nloptr.out <- nloptr(initial.parameters, f, eval_grad_f=function(x) {gradient(f, x)},
                       lb=theta.lb,ub=theta.ub,
                       opts=list("algorithm"="NLOPT_LD_LBFGS","xtol_rel"=1.0e-8))
#                       opts=list("algorithm"="NLOPT_GN_CRS2_LM","xtol_rel"=1.0e-8,"print_level" = 1))
  theta <- nloptr.out$solution
  
  log.likelihood <- log.likelihood.wrapper(theta, y.data, x.data, stage, lambda.g, lambda.z, xi.00, P.00)$ll.cum

#browser()

  # Get state vectors (xi.tt, xi.ttm1, xi.tT, P.tt, P.ttm1, P.tT) via Kalman filter
  states <- kalman.states.wrapper(theta, y.data, x.data, stage, lambda.g, lambda.z, xi.00, P.00)

  # If run.se = TRUE, compute standard errors for estimates of the states (see footnote 7) and report run time
  if (run.se) {
      ptm <- proc.time()
      se <- kalman.standard.errors(T, states, theta, y.data, x.data, stage, lambda.g, lambda.z, xi.00, P.00, niter, a3.constraint, b2.constraint)
      print("Standard error procedure run time")
      print(proc.time() - ptm)
  }
  
  ########## I have added this ##############
  log.likelihood.initialvals <- log.likelihood.wrapper(initial.parameters, y.data, x.data, stage, lambda.g, lambda.z, xi.00, P.00)$ll.cum
  
  # One-sided (filtered) estimates
  trend.filtered      <- states$filtered$xi.tt[,4] * 4
  z.filtered          <- states$filtered$xi.tt[,6]
  rstar.filtered      <- trend.filtered + z.filtered
  potential.filtered  <- states$filtered$xi.tt[,1]/100
  output.gap.filtered <- y.data[,1] - (potential.filtered * 100)

	########### I have added this ################
  # Two-sided (smoothed) estimates (I HAVE FIXED THESE AS THEY WERE ONLY RETURNING FILTERED SERIES)
  trend.smoothed      <- states$smoothed$xi.tT[,4] * 4
  z.smoothed          <- states$smoothed$xi.tT[,6]
  rstar.smoothed      <- trend.smoothed + z.smoothed
  potential.smoothed  <- states$smoothed$xi.tT[,1]/100
  output.gap.smoothed <- y.data[,1] - (potential.smoothed * 100)
  
#  browser()
   
  # Save variables to return
  return.list <- list()
  return.list$rstar.filtered      <- rstar.filtered
  return.list$trend.filtered      <- trend.filtered
  return.list$z.filtered          <- z.filtered
  return.list$potential.filtered  <- potential.filtered
  return.list$output.gap.filtered <- output.gap.filtered
  return.list$rstar.smoothed      <- rstar.smoothed
  return.list$trend.smoothed      <- trend.smoothed
  return.list$z.smoothed          <- z.smoothed
  return.list$potential.smoothed  <- potential.smoothed
  return.list$output.gap.smoothed <- output.gap.smoothed
  return.list$theta               <- theta
  return.list$log.likelihood      <- log.likelihood
  return.list$states              <- states
  return.list$xi.00               <- xi.00
  return.list$P.00                <- P.00
  return.list$y.data              <- y.data
  ########### I have added this #############################################
  return.list$initial.parameters          <- initial.parameters
  return.list$log.likelihood.initialvals  <- log.likelihood.initialvals   
  if (run.se) { return.list$se    <- se }
  return(return.list)
}
