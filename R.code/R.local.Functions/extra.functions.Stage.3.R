#------------------------------------------------------------------------------#
# File: extra.functions.Stage.3.R
#
# Description: This file adds various Stage 3 model versions to the HLW estimation.


#***************************************************************************************************************  
#                            Stage 3 with \sigma_g and \sigma_z estimated by MLE
#***************************************************************************************************************
fit.stage3.gz <- function( log.output, inflation, real.interest.rate,
                          lambda.g, lambda.z, a3.constraint=NA, b2.constraint=NA,
                          P.00, initial.parameters ) {
  stage <- 3
  # Data must start 4 quarters before the estimation period  
  T <- length(log.output) - 4
 
  # Initialization of state vector for Kalman filter using HP trend of log output
  log.output.hp.trend <- hpfilter(log.output,freq=36000,type="lambda",drift=FALSE)$trend
  g.pot <- log.output.hp.trend[(g.pot.start.index):length(log.output.hp.trend)]
  g.pot.diff <- diff(g.pot)
  xi.00 <- c(100*g.pot[3:1],100*g.pot.diff[2:1],0,0) # their default starting values

	#	MAKE Y AND X DATA MATRICES FROM INPUT DATA ------------------------------------------------
  y.data <- cbind(100 * log.output[5:(T+4)],
                  inflation[5:(T+4)])
  x.data <- cbind(100 * log.output[4:(T+3)],
                  100 * log.output[3:(T+2)],
                  real.interest.rate[4:(T+3)],
                  real.interest.rate[3:(T+2)],                  
                  inflation[4:(T+3)],
                  (inflation[3:(T+2)]+inflation[2:(T+1)]+inflation[1:T])/3)

  # Starting values for the parameter vector: THESE ARE PASSED INTO THE FUNCTION NOW
  # initial.parameters <- c(b.is[1:3], b.ph[1], b.ph[3], s.is, s.ph, 0.7) 

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

  # Set the initial covariance matrix (see footnote 6) : THESE ARE PASSED INTO THE FUNCTION NOW  
  # P.00 <- calculate.covariance(initial.parameters, theta.lb, theta.ub, y.data, x.data, stage, lambda.g, lambda.z, xi.00)
  
  # Get parameter estimates via maximum likelihood
  f <- function(theta) {return(-LL.S3.gz(theta, y.data, x.data, stage, lambda.g, lambda.z, xi.00, P.00)$ll.cum)}

  nloptr.out <- nloptr(initial.parameters, f, eval_grad_f=function(x) {gradient(f, x)},
                       lb=theta.lb,ub=theta.ub,
                       opts=list("algorithm"="NLOPT_LD_LBFGS","xtol_rel"=1.0e-8))

  theta <- nloptr.out$solution
  
  log.likelihood <- LL.S3.gz(theta, y.data, x.data, stage, lambda.g, lambda.z, xi.00, P.00)$ll.cum

  # Get state vectors (xi.tt, xi.ttm1, xi.tT, P.tt, P.ttm1, P.tT) via Kalman filter
  states <- KSW.S3.gz(theta, y.data, x.data, stage, lambda.g, lambda.z, xi.00, P.00)

  ########## I have added this ##############
  log.likelihood.initialvals <- LL.S3.gz(initial.parameters, y.data, x.data, stage, lambda.g, lambda.z, xi.00, P.00)$ll.cum
  
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
  return(return.list)
}

#--------------------------------------------------------------------------------------------------------------#
# Stage3.gz model: log.likelihood.wrapper 
#--------------------------------------------------------------------------------------------------------------#
LL.S3.gz <- function(parameters, y.data, x.data, stage=NA, lambda.g=NA, lambda.z=NA, xi.00=NA, P.00=NA){	
  # only this line below changes for the 3 different Stage 3 models with estimated g,z	
  out <- UP.S3.gz(parameters, y.data, x.data, lambda.g, lambda.z, xi.00, P.00)
  # this remains the same
  return(kalman.log.likelihood( out$xi.00, out$P.00, out$F, out$Q, out$A, 
                                out$H, out$R, out$y.data, out$x.data) )
}

#--------------------------------------------------------------------------------------------------------------#
# Stage3.gz model: Unpack parameters as need
#--------------------------------------------------------------------------------------------------------------#
UP.S3.gz <- function(parameters, y.data, x.data, lambda.g, lambda.z, xi.00, P.00) {
# unpack.parameters.stage3 <- function(parameters, y.data, x.data, lambda.g, lambda.z, xi.00, P.00) {
  A         <- matrix(0, 2, 6)   
  A[1, 1:2] <- parameters[1:2]   # a_y,1, a_y,2
  A[1, 3:4] <- parameters[3]/2   # a_r/2
  A[2, 1]   <- parameters[5]     # b_y
  A[2, 5]   <- parameters[4]     # b_pi
  A[2, 6]   <- 1 - parameters[4] # 1 - b_pi
  A         <- t(A)

  H         <- matrix(0, 2, 7)
  H[1, 1]   <- 1
  H[1, 2:3] <- -parameters[1:2]       # a_y,1, a_y,2
  H[1, 4:5] <- -parameters[3] * 2     # -a_r/2 (annualized)
  H[1, 6:7] <- -parameters[3]/2       # -a_r/2
  H[2, 2]   <- -parameters[5]         # -b_y
  H         <- t(H)

  R         <- diag(c(parameters[6]^2, parameters[7]^2)) # sigma_y~, sigma_pi

  Q         <- matrix(0, 7, 7)
  Q[1, 1]   <- parameters[8]^2 + parameters[9]^2          # sigma2_y* + sigma2_g
  Q[1, 4]   <- Q[4, 1] <- Q[4, 4]<- (parameters[9])^2 		# sigma2_g
  Q[6, 6]   <- (parameters[10])^2                         # sigma2_z

  F <- matrix(0, 7, 7)
  F[1, 1] <- F[1, 4] <- F[2, 1] <- F[3, 2] <- F[4,4] <- F[5,4]<- F[6,6] <- F[7,6] <- 1
 
  return(list("xi.00"=xi.00, "P.00"=P.00, "F"=F, "Q"=Q, "A"=A, "H"=H, "R"=R, "x.data"=x.data, "y.data"=y.data))
}

#--------------------------------------------------------------------------------------------------------------#
# Stage3.gz model: kalman.states.wrapper.R
#--------------------------------------------------------------------------------------------------------------#
KSW.S3.gz <- function(parameters, y.data, x.data, stage=NA, lambda.g=NA, lambda.z=NA, xi.00=NA, P.00=NA){
  # only this line below changes for the 3 different Stage 3 models with estimated g,z
  out <- UP.S3.gz(parameters, y.data, x.data, lambda.g, lambda.z, xi.00, P.00)
  # this remains the same
  states <- kalman.states(out$xi.00, out$P.00, out$F, out$Q, out$A, out$H, out$R, out$y.data, out$x.data)
  # browser()
  return(states)
}

#***************************************************************************************************************  
#                            Stage 3 with \sigma_g estimated by MLE, only lambda.z fixed at some exogenous value
#***************************************************************************************************************
fit.stage3.g <- function( log.output, inflation, real.interest.rate,
                          lambda.g, lambda.z, a3.constraint=NA, b2.constraint=NA,
                          P.00, initial.parameters ) {
  stage <- 3
  # Data must start 4 quarters before the estimation period  
  T <- length(log.output) - 4
 
  # Initialization of state vector for Kalman filter using HP trend of log output
  log.output.hp.trend <- hpfilter(log.output,freq=36000,type="lambda",drift=FALSE)$trend
  g.pot <- log.output.hp.trend[(g.pot.start.index):length(log.output.hp.trend)]
  g.pot.diff <- diff(g.pot)
  xi.00 <- c(100*g.pot[3:1],100*g.pot.diff[2:1],0,0) # their default starting values

	#	MAKE Y AND X DATA MATRICES FROM INPUT DATA ------------------------------------------------
  y.data <- cbind(100 * log.output[5:(T+4)],
                  inflation[5:(T+4)])
  x.data <- cbind(100 * log.output[4:(T+3)],
                  100 * log.output[3:(T+2)],
                  real.interest.rate[4:(T+3)],
                  real.interest.rate[3:(T+2)],                  
                  inflation[4:(T+3)],
                  (inflation[3:(T+2)]+inflation[2:(T+1)]+inflation[1:T])/3)

  # Starting values for the parameter vector: THESE ARE PASSED INTO THE FUNCTION NOW
  # initial.parameters <- c(b.is[1:3], b.ph[1], b.ph[3], s.is, s.ph, 0.7) 

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

  # Set the initial covariance matrix (see footnote 6) : THESE ARE PASSED INTO THE FUNCTION NOW  
  # P.00 <- calculate.covariance(initial.parameters, theta.lb, theta.ub, y.data, x.data, stage, lambda.g, lambda.z, xi.00)
  
  # Get parameter estimates via maximum likelihood
  f <- function(theta) {return(-LL.S3.g(theta, y.data, x.data, stage, lambda.g, lambda.z, xi.00, P.00)$ll.cum)}

  nloptr.out <- nloptr(initial.parameters, f, eval_grad_f=function(x) {gradient(f, x)},
                       lb=theta.lb,ub=theta.ub,
                       opts=list("algorithm"="NLOPT_LD_LBFGS", "xtol_rel"=1.0e-8) )
  
                       # opts=list( "algorithm"="NLOPT_LD_LBFGS",
                       #            "print_level" = 3,
                       #            "xtol_rel"=1.0e-8))

  theta <- nloptr.out$solution
  
  log.likelihood <- LL.S3.g(theta, y.data, x.data, stage, lambda.g, lambda.z, xi.00, P.00)$ll.cum

  # Get state vectors (xi.tt, xi.ttm1, xi.tT, P.tt, P.ttm1, P.tT) via Kalman filter
  states <- KSW.S3.g(theta, y.data, x.data, stage, lambda.g, lambda.z, xi.00, P.00)

  ########## I have added this ##############
  log.likelihood.initialvals <- LL.S3.g(initial.parameters, y.data, x.data, stage, lambda.g, lambda.z, xi.00, P.00)$ll.cum
  
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
  return(return.list)
}

#--------------------------------------------------------------------------------------------------------------#
# Stage3.g model: log.likelihood.wrapper 
#--------------------------------------------------------------------------------------------------------------#
LL.S3.g <- function(parameters, y.data, x.data, stage=NA, lambda.g=NA, lambda.z=NA, xi.00=NA, P.00=NA){	
  # only this line below changes for the 3 different Stage 3 models with estimated g,z	
  out <- UP.S3.g(parameters, y.data, x.data, lambda.g, lambda.z, xi.00, P.00)
  # this remains the same
  return(kalman.log.likelihood( out$xi.00, out$P.00, out$F, out$Q, out$A, 
                                out$H, out$R, out$y.data, out$x.data) )
}

#--------------------------------------------------------------------------------------------------------------#
# Stage3.g model: Unpack parameters as need
#--------------------------------------------------------------------------------------------------------------#
UP.S3.g <- function(parameters, y.data, x.data, lambda.g, lambda.z, xi.00, P.00) {
# unpack.parameters.stage3 <- function(parameters, y.data, x.data, lambda.g, lambda.z, xi.00, P.00) {
  A         <- matrix(0, 2, 6)   
  A[1, 1:2] <- parameters[1:2]   # a_y,1, a_y,2
  A[1, 3:4] <- parameters[3]/2   # a_r/2
  A[2, 1]   <- parameters[5]     # b_y
  A[2, 5]   <- parameters[4]     # b_pi
  A[2, 6]   <- 1 - parameters[4] # 1 - b_pi
  A         <- t(A)

  H         <- matrix(0, 2, 7)
  H[1, 1]   <- 1
  H[1, 2:3] <- -parameters[1:2]       # a_y,1, a_y,2
  H[1, 4:5] <- -parameters[3] * 2     # -a_r/2 (annualized)
  H[1, 6:7] <- -parameters[3]/2       # -a_r/2
  H[2, 2]   <- -parameters[5]         # -b_y
  H         <- t(H)

  R         <- diag(c(parameters[6]^2, parameters[7]^2)) # sigma_y~, sigma_pi

  Q         <- matrix(0, 7, 7)
  Q[1, 1]   <- parameters[8]^2 + parameters[9]^2          # sigma2_y* + sigma2_g
  Q[1, 4]   <- Q[4, 1] <- Q[4, 4]<- (parameters[9])^2 		# sigma2_g
  Q[6, 6]   <- (lambda.z*parameters[6]/parameters[3])^2   # sigma_y~/a_r

  F <- matrix(0, 7, 7)
  F[1, 1] <- F[1, 4] <- F[2, 1] <- F[3, 2] <- F[4,4] <- F[5,4]<- F[6,6] <- F[7,6] <- 1
 
  return(list("xi.00"=xi.00, "P.00"=P.00, "F"=F, "Q"=Q, "A"=A, "H"=H, "R"=R, "x.data"=x.data, "y.data"=y.data))
}

#--------------------------------------------------------------------------------------------------------------#
# Stage3.g model: kalman.states.wrapper.R
#--------------------------------------------------------------------------------------------------------------#
KSW.S3.g <- function(parameters, y.data, x.data, stage=NA, lambda.g=NA, lambda.z=NA, xi.00=NA, P.00=NA){
  # only this line below changes for the 3 different Stage 3 models with estimated g,z
  out <- UP.S3.g(parameters, y.data, x.data, lambda.g, lambda.z, xi.00, P.00)
  # this remains the same
  states <- kalman.states(out$xi.00, out$P.00, out$F, out$Q, out$A, out$H, out$R, out$y.data, out$x.data)
  # browser()
  return(states)
}
