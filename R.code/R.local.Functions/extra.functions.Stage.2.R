#------------------------------------------------------------------------------#
# File: extra.functions.Stage.2.R
#
# Description: This file adds various Stage 2 model versions to the HLW estimation.

#*************************************************************************************************************  
#                            Correct Stage 2 M0 with \sigma_g estimated
#*************************************************************************************************************  

#--------------------------------------------------------------------------------------------------------------#
# Stage2.M0g model: fit the Correct Stage 2 M0 with \sigma_g estimated
#--------------------------------------------------------------------------------------------------------------#
fit.stage2.M0g <- function( log.output, inflation, real.interest.rate,
                            lambda.g, a3.constraint=NA, b2.constraint=NA,
                            P.00, initial.parameters) {
  stage <- 2
  # Data must start 4 quarters before the estimation period  
  T <- length(log.output) - 4
  # pass these in through function input
  # initial.parameters 

  y.data <- cbind(100 * log.output[5:(T+4)],
                  inflation[5:(T+4)])
  x.data <- cbind(100 * log.output[4:(T+3)],
                  100 * log.output[3:(T+2)],
                  real.interest.rate[4:(T+3)],
                  real.interest.rate[3:(T+2)],                  
                  inflation[4:(T+3)],
                 (inflation[3:(T+2)]+inflation[2:(T+1)]+inflation[1:T])/3)
  # ,rep(1,T)) # remove the last column of ones from the x.data input vector

  # Initialization of state vector for Kalman filter using HP trend of log output
  log.output.hp.trend <- hpfilter(log.output,freq=36000,type="lambda",drift=FALSE)$trend
  g.pot <- log.output.hp.trend[(g.pot.start.index):length(log.output.hp.trend)]
  g.pot.diff <- diff(g.pot)

  # ---------------------------------------------------------------------------------------------
  # INITIAL CONDITION FOR STATE VECTOR
  # xi.00 <- c(100*g.pot[3:1],100*g.pot.diff[2:1],0,0) # their default starting values
  xi.00 <- c(100*g.pot[3:1],100*g.pot.diff[2:1])

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
    # print(paste0("Setting an upper bound on a_3 of ",as.character(a3.constraint)))
      if (initial.parameters[3] >  a3.constraint) {
          initial.parameters[3] <- a3.constraint
      }
      theta.ub[3] <- a3.constraint      
  }

  # Get parameter estimates via maximum likelihood this minimizes
  f <- function(theta) {return( -LL.S2.M0g(theta, y.data, x.data, xi.00, P.00)$ll.cum)}
  
  nloptr.out <- nloptr(initial.parameters, f, eval_grad_f=function(x) {gradient(f, x)},
                       lb=theta.lb,ub=theta.ub,
                       opts=list("algorithm"="NLOPT_LD_LBFGS", "xtol_rel"=1.0e-8) )

# opts=list("algorithm"="NLOPT_LD_LBFGS", "xtol_rel"=1.0e-8) )  
  
  theta <- nloptr.out$solution

  log.likelihood <- LL.S2.M0g(theta, y.data, x.data, xi.00, P.00)$ll.cum

  ##  # log-likelihood at initial values  ----
  log.likelihood.initialvals <- LL.S2.M0g(initial.parameters, y.data, x.data, xi.00, P.00)$ll.cum

# browser() 
  # Get state vectors (xi.tt, xi.ttm1, xi.tT, P.tt, P.ttm1, P.tT) via Kalman filter
  # states <- kalman.states.wrapper(theta, y.data, x.data, stage, lambda.g, NA, xi.00, P.00)
  states <- KSW.S2.M0g(theta, y.data, x.data, stage, lambda.g, NA, xi.00, P.00)

  # Two-sided (smoothed) estimates
  trend.smoothed      <- states$smoothed$xi.tt[,4] * 4
  potential.smoothed  <- c(states$smoothed$xi.tT[1, 3:2], states$smoothed$xi.tT[,1])
  output.gap.smoothed <- 100 * log.output[3:(T+4)] - potential.smoothed

  # Inputs for median.unbiased.estimator.stage2.R, these need to be adapted to the Stage2.M0 model
  y <- output.gap.smoothed[3:length(output.gap.smoothed)]

  x <- cbind( output.gap.smoothed[2:(length(output.gap.smoothed)-1)],
              output.gap.smoothed[1:(length(output.gap.smoothed)-2)],
              (x.data[,3]+x.data[,4])/2 - 4*( states$smoothed$xi.tT[,4]+states$smoothed$xi.tT[,5] )/2 )


  # One-sided (filtered) estimates  
  trend.filtered      <- states$filtered$xi.tt[,4] * 4
  potential.filtered  <- states$filtered$xi.tt[,1]/100
  output.gap.filtered <- y.data[,1] - (potential.filtered * 100)
  
  # Save variables to return -----
  return.list <- list()
  return.list$y                           <- y
  return.list$x                           <- x
  return.list$theta                       <- theta
  return.list$log.likelihood              <- log.likelihood
  return.list$states                      <- states
  return.list$xi.00                       <- xi.00
  return.list$P.00                        <- P.00
  return.list$trend.filtered              <- trend.filtered
  return.list$potential.filtered          <- potential.filtered
  return.list$output.gap.filtered         <- output.gap.filtered
  return.list$trend.smoothed              <- trend.smoothed
  return.list$potential.smoothed          <- potential.smoothed
  return.list$output.gap.smoothed         <- output.gap.smoothed
  ########### I have added this #############################################
  return.list$initial.parameters          <- initial.parameters
  return.list$log.likelihood.initialvals  <- log.likelihood.initialvals   
  return(return.list)
} 

#--------------------------------------------------------------------------------------------------------------#
# Stage2.M0g model: log.likelihood.wrapper.R
#--------------------------------------------------------------------------------------------------------------#
LL.S2.M0g <- function(parameters, y.data, x.data, 
                          xi.00=NA, P.00=NA){
  out <- UP.S2.M0g(parameters, y.data, x.data, xi.00, P.00)
  # this remains the same
  return(kalman.log.likelihood( out$xi.00, out$P.00, out$F, out$Q, out$A, 
                                out$H, out$R, out$y.data, out$x.data))
}

#--------------------------------------------------------------------------------------------------------------#
# Stage2.M0g model: Unpack parameters as need
#--------------------------------------------------------------------------------------------------------------#
UP.S2.M0g <- function(parameters, y.data, x.data, xi.00=NA, P.00=NA) {
# this is the stage 3 model
  A         <- matrix(0, 2, 6)   
  A[1, 1:2] <- parameters[1:2]   # a_y,1, a_y,2
  A[1, 3:4] <- parameters[3]/2   # a_r/2
  A[2, 1]   <- parameters[5]     # b_y
  A[2, 5]   <- parameters[4]     # b_pi
  A[2, 6]   <- 1 - parameters[4] # 1 - b_pi
  A         <- t(A)

  H         <- matrix(0, 2, 5)
  H[1, 1]   <- 1
  H[1, 2:3] <- -parameters[1:2]       # a_y,1, a_y,2
  H[1, 4:5] <- -parameters[3] * 2     # -a_r/2 (annualized)
  H[2, 2]   <- -parameters[5]         # -b_y
  H         <- t(H)

  R         <- diag(c(parameters[6]^2, parameters[7]^2)) # sigma_y~, sigma_pi

  Q         <- matrix(0, 5, 5)
  Q[1, 1]   <- (parameters[8]^2 + parameters[9]^2)    # sigma_y* + 
  Q[1, 4]   <- Q[4, 1] <- Q[4, 4]<- parameters[9]^2   # sigma_y*

  # browser()
  
  F <- matrix(0, 5, 5)
  F[1, 1] <- F[1, 4] <- F[2, 1] <- F[3, 2] <- F[4,4] <- F[5,4] <- 1

  # return values are as before  
  return(list("xi.00"=xi.00, "P.00"=P.00, "F"=F, "Q"=Q, "A"=A, "H"=H, "R"=R, "x.data"=x.data, "y.data"=y.data))
}

#--------------------------------------------------------------------------------------------------------------#
# Stage2.M0g model: kalman.states.wrapper.R
#--------------------------------------------------------------------------------------------------------------#
KSW.S2.M0g <- function(parameters, y.data, x.data, stage = NA,
                                  lambda.g=NA, lambda.z=NA, xi.00=NA, P.00=NA){
  out <- UP.S2.M0g(parameters, y.data, x.data, xi.00, P.00)
  # this remains the same
  states <- kalman.states(out$xi.00, out$P.00, out$F, out$Q, out$A, out$H, out$R, out$y.data, out$x.data)
  # browser()
  return(states)
}


#***************************************************************************************************************  
#                            Stage 2 with \sigma_g estimated
#***************************************************************************************************************


#--------------------------------------------------------------------------------------------------------------#
# Stage2.g model: fit the model Stage 2 with \sigma_g estimated
#--------------------------------------------------------------------------------------------------------------#
fit.stage2.g <- function( log.output, inflation, real.interest.rate,
                          lambda.g, a3.constraint=NA, b2.constraint=NA,
                          P.00, initial.parameters ) {

  stage <- 2
  # Data must start 4 quarters before the estimation period  
  T <- length(log.output) - 4
 
  # # Initialization of state vector for Kalman filter using HP trend of log output
  log.output.hp.trend <- hpfilter(log.output,freq=36000,type="lambda",drift=FALSE)$trend
  g.pot       <- log.output.hp.trend[(g.pot.start.index):length(log.output.hp.trend)]
  g.pot.diff  <- diff(g.pot)
  xi.00       <- c(100*g.pot[3:1],100*g.pot.diff[2])

  # make y and x data matrices from input data
  y.data <- cbind(100 * log.output[5:(T+4)],
                  inflation[5:(T+4)])
  x.data <- cbind(100 * log.output[4:(T+3)],
                  100 * log.output[3:(T+2)],
                  real.interest.rate[4:(T+3)],
                  real.interest.rate[3:(T+2)],                  
                  inflation[4:(T+3)],
                  (inflation[3:(T+2)]+inflation[2:(T+1)]+inflation[1:T])/3,
                  rep(1,T))

  # Starting values for the parameter vector
  # sd.g0 = sd(g.pot.diff)
  # initial.parameters <- c(b.is, -b.is[3], b.ph[1], b.ph[3], s.is, s.ph, 0.5, sd.g0)

  # browser()
  
  # Set an upper and lower bound on the parameter vectors:
  # The vector is unbounded unless values are otherwise specified
  theta.lb <- c(rep(-Inf,length(initial.parameters)))
  theta.ub <- c(rep( Inf,length(initial.parameters)))

  # Set a lower bound for the Phillips curve slope (b_2) of b2.constraint, if not NA
  # In HLW, b2.constraint = 0.025
  if (!is.na(b2.constraint)) {
      if (initial.parameters[7] <  b2.constraint) {
          initial.parameters[7] <- b2.constraint
      }
      theta.lb[7] <- b2.constraint
  }

  # Set an upper bound for the IS curve slope (a_3) of a3.constraint, if not NA
  # In HLW, a3.constraint = -0.0025
  if (!is.na(a3.constraint)) {
      if (initial.parameters[3] >  a3.constraint) {
          initial.parameters[3] <- a3.constraint
      }
      theta.ub[3] <- a3.constraint      
  }

  # Set the initial covariance matrix (see footnote 6)
  # P.00 <- calculate.covariance(initial.parameters, theta.lb, theta.ub, y.data, x.data, stage, lambda.g, NA, xi.00)

  # Get parameter estimates via maximum likelihood this minimizes
  f <- function(theta) {return( -LL.S2.g(theta, y.data, x.data, xi.00, P.00)$ll.cum)}
  
  nloptr.out <- nloptr(initial.parameters, f, eval_grad_f=function(x) {gradient(f, x)},
                       lb=theta.lb,ub=theta.ub,
                       opts=list("algorithm"="NLOPT_LD_LBFGS","xtol_rel"=1.0e-8))
  
  theta <- nloptr.out$solution

  # browser() 

  log.likelihood <- LL.S2.g(theta, y.data, x.data, xi.00, P.00)$ll.cum

  ##  # log-likelihood at initial values  ----
  log.likelihood.initialvals <- LL.S2.g(initial.parameters, y.data, x.data, xi.00, P.00)$ll.cum

  # Get state vectors (xi.tt, xi.ttm1, xi.tT, P.tt, P.ttm1, P.tT) via Kalman filter
  # states <- kalman.states.wrapper(theta, y.data, x.data, stage, lambda.g, NA, xi.00, P.00)
  states <- KSW.S2.g(theta, y.data, x.data, stage, lambda.g, NA, xi.00, P.00)

  # Two-sided (smoothed) estimates
  trend.smoothed      <- states$smoothed$xi.tt[,4] * 4
  potential.smoothed  <- c(states$smoothed$xi.tT[1, 3:2], states$smoothed$xi.tT[,1])
  output.gap.smoothed <- 100 * log.output[3:(T+4)] - potential.smoothed

  # Inputs for median.unbiased.estimator.stage2.R
  y <- output.gap.smoothed[3:length(output.gap.smoothed)]
  x <- cbind(output.gap.smoothed[2:(length(output.gap.smoothed)-1)],
             output.gap.smoothed[1:(length(output.gap.smoothed)-2)],
             (x.data[,3]+x.data[,4])/2,
             states$smoothed$xi.tT[,4],
             rep(1,T))

  # One-sided (filtered) estimates  
  trend.filtered      <- states$filtered$xi.tt[,4] * 4
  potential.filtered  <- states$filtered$xi.tt[,1]/100
  output.gap.filtered <- y.data[,1] - (potential.filtered * 100)
  
  # Save variables to return -----
  return.list <- list()
  return.list$y                           <- y
  return.list$x                           <- x
  return.list$theta                       <- theta
  return.list$log.likelihood              <- log.likelihood
  return.list$states                      <- states
  return.list$xi.00                       <- xi.00
  return.list$P.00                        <- P.00
  return.list$trend.filtered              <- trend.filtered
  return.list$potential.filtered          <- potential.filtered
  return.list$output.gap.filtered         <- output.gap.filtered
  return.list$trend.smoothed              <- trend.smoothed
  return.list$potential.smoothed          <- potential.smoothed
  return.list$output.gap.smoothed         <- output.gap.smoothed
  ########### I have added this #############################################
  return.list$initial.parameters          <- initial.parameters
  return.list$log.likelihood.initialvals  <- log.likelihood.initialvals   
  return(return.list)
} 

#--------------------------------------------------------------------------------------------------------------#
# Stage2.g model: log.likelihood.wrapper 
#--------------------------------------------------------------------------------------------------------------#
LL.S2.g <- function(parameters, y.data, x.data, 
                          xi.00=NA, P.00=NA){
  out <- UP.S2.g(parameters, y.data, x.data, xi.00, P.00)
  # this remains the same
  return(kalman.log.likelihood( out$xi.00, out$P.00, out$F, out$Q, out$A, 
                                out$H, out$R, out$y.data, out$x.data))
}

#--------------------------------------------------------------------------------------------------------------#
# Stage2.g model: Unpack parameters as need
#--------------------------------------------------------------------------------------------------------------#
UP.S2.g <- function(parameters, y.data, x.data, xi.00=NA, P.00=NA) {
  A         <- matrix(0, 2, 7)
  A[1, 1:2] <- parameters[1:2]   # a_y,1, a_y,2
  A[1, 3:4] <- parameters[3]/2   # a_r/2
  A[1, 7]   <- parameters[4]     # a_0
  A[2, 1]   <- parameters[7]     # b_y
  A[2, 5]   <- parameters[6]     # b_pi
  A[2, 6]   <- 1 - parameters[6] # 1 - b_pi
  A         <- t(A)
  
  H         <- matrix(0, 2, 4)
  H[1, 1  ] <- 1
  H[1, 2:3] <- -parameters[1:2] # -a_y,1, -a_y,2
  H[1, 4  ] <- parameters[5]    # a_g
  H[2, 2]   <- -parameters[7]   # -b_y
  H         <- t(H)

  R         <- diag(c(parameters[8]^2, parameters[9]^2)) # sigma_y~, sigma_pi
  Q         <- matrix(0, 4, 4)
  Q[1, 1]   <- parameters[10]^2              # sigma_y*^2
  # Q[4, 4]   <- (lambda.g * parameters[10])^2 # sigma_y*
  Q[4, 4]   <- parameters[11]^2              # sigma_g^2

  F <- matrix(0, 4, 4)
  F[1, 1] <- F[1, 4] <- F[2, 1] <- F[3, 2] <- F[4,4] <- 1
  
  return(list("xi.00"=xi.00, "P.00"=P.00, "F"=F, "Q"=Q, "A"=A, "H"=H, "R"=R, "x.data"=x.data, "y.data"=y.data))
}

#--------------------------------------------------------------------------------------------------------------#
# Stage2.g model: kalman.states.wrapper.R
#--------------------------------------------------------------------------------------------------------------#
KSW.S2.g <- function(parameters, y.data, x.data, stage = NA,
                                  lambda.g=NA, lambda.z=NA, xi.00=NA, P.00=NA){
  out <- UP.S2.g(parameters, y.data, x.data, xi.00, P.00)
  # this remains the same
  states <- kalman.states(out$xi.00, out$P.00, out$F, out$Q, out$A, out$H, out$R, out$y.data, out$x.data)
  # browser()
  return(states)
}