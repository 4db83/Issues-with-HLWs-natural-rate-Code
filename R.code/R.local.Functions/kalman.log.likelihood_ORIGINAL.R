#------------------------------------------------------------------------------#
# File:        kalman.log.likelihood.R
#
# Description: This function takes as input the coefficient matrices of the
#              given state-space model and the conditional expectation and
#              covariance matrix of the initial state and returns a vector
#              with the log likelihood value from each period,
#              as well as the cumulative sum.
# I Have replaced log(2 * atan(1) * 4) with log_2_pi = 1.837877066409345483560659472811235279722
# outside the TS loop for speed.
#------------------------------------------------------------------------------#
kalman.log.likelihood_ORIGINAL <- function(xi.tm1tm1, P.tm1tm1, F, Q, A, H, R, y, x) {
    T <- dim(y)[1]
    n <- dim(y)[2]
    ll.vec <- matrix(NA,T,1)
    ll.cum <- 0
    xi.tt <- xi.tm1tm1
    P.tt  <- P.tm1tm1    
    # I have added this line to replace log(2 * atan(1) * 4)
    log_2_pi = 1.837877066409345483560659472811235279722
    # log_2_pi = log(2 * atan(1) * 4)

    for (t in 1:T){

        xi.ttm1 <- F %*% xi.tt
        P.ttm1  <- F %*% P.tt %*% t(F) + Q
        prediction.error <- (as.vector(y[t,]) - as.vector(t(A) %*% as.vector(x[t,])) - as.vector(t(H) %*% xi.ttm1))
        HPHR    <- t(H) %*% P.ttm1 %*% H + R
         # cat( det(HPHR), "\n")

        ll.vec[t] <- drop( -(n/2)*log_2_pi -0.5*log( det(HPHR) ) -0.5*prediction.error %*% solve(HPHR, prediction.error) )
        ll.cum  <- ll.cum + ll.vec[t]

        xi.tt   <- xi.ttm1 + P.ttm1 %*% H %*% solve(HPHR, prediction.error)
        P.tt    <-  P.ttm1 - P.ttm1 %*% H %*% solve(HPHR, t(H) %*% P.ttm1)
    }
    return(list("ll.vec"=ll.vec,"ll.cum"=ll.cum))
}
