#------------------------------------------------------------------------------#
# File:        format.output.R
#
# Description: This generates a dataframe to be written to a CSV containing
#              one-sided estimates, parameter values, standard errors,
#              and other statistics of interest.
#------------------------------------------------------------------------------#
format.output <- function(country.estimation, one.sided.est.country, real.rate.country, start, end, run.se = TRUE) {
    output.country <- data.frame(matrix(NA,dim(one.sided.est.country)[1],22))
    
    output.country[,1]   <- seq(from = (as.Date(ti(shiftQuarter(start,-1),'quarterly'))+1), to = (as.Date(ti(shiftQuarter(end,-1),tif='quarterly'))+1), by = 'quarter')
    output.country[,2:5] <- one.sided.est.country

    output.country[1,7]    <- "Parameter Point Estimates"
    output.country[2,7:15] <- c("a_1","a_2","a_3","b_1","b_2","sigma_1","sigma_2","sigma_4","a_1 + a_2")
    output.country[3,7:14] <- country.estimation$out.stage3$theta
    output.country[3,15]   <- country.estimation$out.stage3$theta[1] + country.estimation$out.stage3$theta[2]
    # Include standard errors in output only if run.se switch is TRUE
    if (run.se) {
        output.country[4,7]    <- "T Statistics"
        output.country[5,7:14] <- country.estimation$out.stage3$se$t.stats
        
        output.country[8,7]    <- "Average Standard Errors"
        output.country[9,7:9]  <- c("y*","r*","g")
        output.country[10,7:9] <- country.estimation$out.stage3$se$se.mean
        
        output.country[12,7]   <- "Restrictions on MC draws: a_3 < -0.0025; b_2 > 0.025; a_1 + a_2 < 1"
        output.country[13,7]   <- "Draws excluded:"; output.country[13,9] <- country.estimation$out.stage3$se$number.excluded
        output.country[13,10]  <- "Total:"; output.country[13,11] <- niter
        output.country[14,7]   <- "Percent excluded:"; output.country[14,9] <- as.numeric(output.country[13,9]) / (as.numeric(output.country[13,9]) + as.numeric(output.country[13,11]))
        output.country[15,7]   <- "Draws excluded because a_3 > -0.0025:"; output.country[15,11] <- country.estimation$out.stage3$se$number.excluded.a3
        output.country[16,7]   <- "Draws excluded because b_2 <  0.025:"; output.country[16,11] <- country.estimation$out.stage3$se$number.excluded.b2
        output.country[17,7]   <- "Draws excluded because a_1 + a_2 < 1:"; output.country[17,11] <- country.estimation$out.stage3$se$number.excluded.a1a2
    }
    
    output.country[19,7] <- "Signal-to-noise Ratios"
    output.country[20,7] <- "lambda_g"; output.country[20,8] <- country.estimation$lambda.g
    output.country[21,7] <- "lambda_z"; output.country[21,8] <- country.estimation$lambda.z
    output.country[19,11] <- "Log Likelihood"; output.country[20,11] <- country.estimation$out.stage3$log.likelihood

    output.country[24,7] <- "State vector: [y_{t}* y_{t-1}* y_{t-2}* g_{t-1} g_{t-2} z_{t-1} z_{t-2}]"
    output.country[25,7] <- "Initial State Vector"
    output.country[26,7:13] <- country.estimation$out.stage3$xi.00
    output.country[28,7] <- "Initial Covariance Matrix"
    output.country[29:35,7:13] <- country.estimation$out.stage3$P.00

    if (run.se) {
        output.country[,17]   <- seq(from = (as.Date(ti(shiftQuarter(start,-1),'quarterly'))+1), to = (as.Date(ti(shiftQuarter(end,-1),tif='quarterly'))+1), by = 'quarter')
        output.country[,18:20] <- country.estimation$out.stage3$se$se
    }
    
    output.country[,22] <- real.rate.country[5:length(real.rate.country)] - country.estimation$out.stage3$rstar.filtered
  
    return(output.country)
}
