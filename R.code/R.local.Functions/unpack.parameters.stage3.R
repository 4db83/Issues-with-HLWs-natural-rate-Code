#------------------------------------------------------------------------------#
# File:        unpack.parameters.stage3.R
#
# Description: This file generates coefficient matrices for the stage 3
#              state-space model for the given parameter vector.
#
# Stage 3 parameter vector: [a_y,1, a_y,2, a_r, b_pi, b_y, sigma_y~, sigma_pi, sigma_y*]  
#------------------------------------------------------------------------------#
unpack.parameters.stage3 <- function(parameters, y.data, x.data, lambda.g, lambda.z, xi.00, P.00) {
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
  Q[1, 1]   <- (1+lambda.g^2)*parameters[8]^2                  # sigma_y*
  Q[1, 4]   <- Q[4, 1] <- Q[4, 4]<- (lambda.g*parameters[8])^2 # sigma_y*
  Q[6, 6]   <- (lambda.z*parameters[6]/parameters[3])^2        # sigma_y~/a_r

  F <- matrix(0, 7, 7)
  F[1, 1] <- F[1, 4] <- F[2, 1] <- F[3, 2] <- F[4,4] <- F[5,4]<- F[6,6] <- F[7,6] <- 1
 
  return(list("xi.00"=xi.00, "P.00"=P.00, "F"=F, "Q"=Q, "A"=A, "H"=H, "R"=R, "x.data"=x.data, "y.data"=y.data))
}



