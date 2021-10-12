#------------------------------------------------------------------------------#
# File:        median.unbiased.estimator.stage1.R
#
# Description: This file implements the median unbiased estimation of the
#              signal-to-noise ratio lambda_g following Stock and Watson (1998).
#------------------------------------------------------------------------------#
function.median.unbiased.estimator.stage1.correct.order <- function(series) {
    T <- length(series)
    y <- 400 * diff(series)

# initialize these anywhere
    stat <- rep(T-2*4)
    Fstat = NA
    bhats = NA
    
    for (i in 4:(T-5)) {
        xr <- cbind(rep(1, T-1), c(rep(0,i),rep(1,T-i-1)))
        xi <- solve(t(xr) %*% xr)
        b  <- solve(t(xr) %*% xr, t(xr) %*% y)
        s3 <- sum((y-xr%*%b)^2)/(T-2-1)
        stat[i+1-4] = b[2]/sqrt(s3*xi[2,2])
        # xr already includes a constant in the model (does not change the Fstat)
        # what is evident though from a plot of y is that there is a big drop in the variance of y
        olsreg = lm(y ~ 0 + xr)
        # F-stat = tstat^2
        Fstat[i+1-4] <- (summary(olsreg)$coefficients[2,3])^2 
#        bhats[i+1-4,2] = olsreg$coefficients
    }
##      	# run OLS regression of y(t) on an intercept as well as a vector of 0,1 with mean shifts.
##      	olsreg = lm(y ~ 0 + xr)
##      	Fstat[i+1-4] <- (summary(olsreg)$coefficients[2,3])^2 # F-stat = tstat^2
##      	
##      	# now do the Stock and Watson way
##      	# first half of sample
##      	y1 = y[xr[,2]==0]
##      	e1 = y1-mean(y1);
##      	# second half of sample
##      	y2 = y[xr[,2]==1]
##      	e2 = y2-mean(y2);
##    	
##      	SSB[i+1-4] = sum(e1^2) + sum(e2^2)
##    }
    #     ew <- 0
    # for (i in 1:length(stat)) {
    #     ew <- ew+exp(stat[i]^2/2)
    # }

    t0 = floor(T*.15)
    t1 = floor(T*.85)
    
    stat2 = stat[t0:t1]

    # this is much simpler    
    ew    <- log(mean(exp(0.5*stat^2)))
    ew2   <- log(mean(exp(0.5*stat2^2)))
    ewF   <- log(mean(exp(0.5*Fstat)))
    # ew  <- log(ew/length(stat))

    # this is much simpler
    mw    <- mean(stat^2)
    mw2   <- mean(stat2^2)
    mwF   <- mean(Fstat)
    # mw  <- sum(stat^2) / length(stat)

    qlr   <- max(stat^2)
    qlr2  <- max(stat2^2)
    qlrF  <- max(Fstat)
    
#    browser()

    # Values are from Table 3 in Stock and Watson (1998)
    # Test Statistic: Exponential Wald (EW)
    valew <- c(0.426, 0.476, 0.516, 0.661, 0.826, 1.111,
               1.419, 1.762, 2.355, 2.91,  3.413, 3.868, 4.925,
               5.684, 6.670, 7.690, 8.477, 9.191, 10.693, 12.024,
               13.089, 14.440, 16.191, 17.332, 18.699, 20.464,
               21.667, 23.851, 25.538, 26.762, 27.874)
    # Test Statistic: Mean Wald (MW)
    valmw <- c(0.689, 0.757, 0.806, 1.015, 1.234, 1.632,
               2.018, 2.390, 3.081, 3.699, 4.222, 4.776, 5.767,
               6.586, 7.703, 8.683, 9.467, 10.101, 11.639, 13.039,
               13.900, 15.214, 16.806, 18.330, 19.020, 20.562,
               21.837, 24.350, 26.248, 27.089, 27.758)
    # Test Statistic: QLR
    valql <- c(3.198, 3.416, 3.594, 4.106, 4.848, 5.689,
               6.682, 7.626, 9.16,  10.66, 11.841, 13.098, 15.451,
               17.094, 19.423, 21.682, 23.342, 24.920, 28.174, 30.736,
               33.313, 36.109, 39.673, 41.955, 45.056, 48.647, 50.983,
               55.514, 59.278, 61.311, 64.016)
    
    lame <- NA
    lamm <- NA
    lamq <- NA

    # Median-unbiased estimator of lambda_g for given values of the test
    # statistics are obtained using the procedure described in the 
    # footnote to Stock and Watson (1998) Table 3.
    if (ew <= valew[1]) {
        lame <- 0
    } else {
        for (i in 1:(length(valew)-1)) {
            if ((ew > valew[i]) & (ew <= valew[i+1])) {
                lame <- i-1+(ew-valew[i])/(valew[i+1]-valew[i])
            }
        }
    }

    if (mw <= valmw[1]) {
        lamm <- 0
    } else {
        for (i in 1:(length(valmw)-1)) {
            if ((mw > valmw[i]) & (mw <= valmw[i+1])) {
                lamm <- i-1+(mw-valmw[i])/(valmw[i+1]-valmw[i])
            }
        }
    }

    if (qlr <= valql[1]) {
        lamq <- 0
    } else {
        for (i in 1:(length(valql)-1)) {
            if ((qlr > valql[i]) & (qlr <= valql[i+1])) {
                lamq <- i-1+(qlr-valql[i])/(valql[i+1]-valql[i])
            }
        }
    }
    if (is.na(lame) | is.na(lamm) | is.na(lamq)) {
        print("At least one statistic has an NA value. Check to see if your EW, MW, and/or QLR value is outside of Table 3.")
    }
    
    # browser()
    
    stats <- c(ew, mw, qlr)
    lams  <- c(lame, lamm, lamq)
#    return(lame/(T-1))
    return(
      list(
#      "Statistics"  = rbind(ew,mw,qlr),
#      "Statistics2" = rbind(ew2,mw2,qlr2),
#      "Lambdas"     = rbind(lame/(T-1),lamm/(T-1),lamq/(T-1))
# I have changed the order so it fits the Table results in the paper
      "Statistics"  = rbind(mw,ew,qlr),
      "Lambdas"     = rbind(lamm/(T-1),lame/(T-1),lamq/(T-1))
      )
    )
}



