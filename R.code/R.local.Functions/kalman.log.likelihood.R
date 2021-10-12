#------------------------------------------------------------------------------#
# File:        kalman.log.likelihood.R
#
# Description: This function takes as input the coefficient matrices of the
#              given state-space model and the conditional expectation and
#              covariance matrix of the initial state and returns a vector
#              with the log likelihood value from each period,
#              as well as the cumulative sum.
# --------------------------------------------------------------------------------------------------- #
# I Have replaced log(2 * atan(1) * 4) with log_2_pi = 1.837877066409345483560659472811235279722
# outside the TS loop for speed, and have precomputed some other values as well to speed up the 
# execution. Otherwise, numerical results are exactly the same as in the original kalman.log.likelihood.R
# file of HLW, which I have renamed to kalman.log.likelihood_ORIGINAL.R
# --------------------------------------------------------------------------------------------------- #
kalman.log.likelihood <- function(xi.tm1tm1, P.tm1tm1, F, Q, A, H, R, y, x) {
    T <- dim(y)[1]
    n <- dim(y)[2]
    ll.vec  <- matrix(NA,T,1)
    ll.cum  <- 0
    xi.tt   <- xi.tm1tm1
    P.tt    <- P.tm1tm1    
    # I have added this line to replace log(2 * atan(1) * 4)
    log_2_pi = 1.837877066409345483560659472811235279722
    # log_2_pi = log(2 * atan(1) * 4)

    # MAKE THE TRANSPOSES OUTSIDE THE LOOP
    t_F = t(F) 
    t_H = t(H)
    t_A = t(A)

# browser()

    for (t in 1:T){

        xi.ttm1 <- F %*% xi.tt
        P.ttm1  <- F %*% P.tt %*% t_F + Q
        
        prediction.error <- as.vector( (y[t,]) - ( t_A %*% (x[t,])) - (t_H %*% xi.ttm1) )
        # prediction.error <- ( y[t,] - as.vector(t(A) %*% as.vector(x[t,])) - as.vector(t_H %*% xi.ttm1))

        tHPttm1 <- t_H %*% P.ttm1
        HPHR    <- tHPttm1 %*% H + R

        # pre-compute once, this is expensive
        inv_HPHR_hat = solve(HPHR, prediction.error)
        # ( rcond(inv_HPHR_hat) )

        # ll.vec[t] <- drop( -(n/2)*log_2_pi -0.5*log( det(HPHR) ) -0.5*prediction.error %*% solve(HPHR, prediction.error) )
        ll.vec[t] <- drop( -(n/2)*log_2_pi -0.5*log( det(HPHR) ) -0.5*prediction.error %*% inv_HPHR_hat )
        ll.cum  <- ll.cum + ll.vec[t]

        xi.tt   <- xi.ttm1 + P.ttm1 %*% H %*% inv_HPHR_hat
        
        # xi.tt   <- xi.ttm1 + P.ttm1 %*% H %*% solve(HPHR, prediction.error)
        # print( cbind(solve(HPHR, t_H %*% P.ttm1) ) )

        P.tt    <-  P.ttm1 - P.ttm1 %*% H %*% solve(HPHR, tHPttm1)
    }
  
  # print( cbind(solve(HPHR, t_H %*% P.ttm1) ) )
    return(list("ll.vec"=ll.vec,"ll.cum"=ll.cum))

}
