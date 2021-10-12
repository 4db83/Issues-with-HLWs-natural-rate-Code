#------------------------------------------------------------------------------
# File:        run.hlw.estimation.R
#
# Description: Runs the three stages of the HLW estimation using country-specific
#              inputs and returns output from each stage.
#------------------------------------------------------------------------------#
run.hlw.estimation <- function(log.output, inflation, real.interest.rate, nominal.interest.rate,
                               a3.constraint = NA, b2.constraint = NA, run.se = TRUE) {
    print("Estimating Stage 1")
    # Running the stage 1 model
    out.stage1 <- rstar.stage1(log.output,
                               inflation,
                               b2.constraint)
    print("Done Estimating Stage 1, --> MUE Stage 1")
    # Median unbiased estimate of lambda_g
    lambda.g <- median.unbiased.estimator.stage1(out.stage1$potential.smoothed)
    
    print("Done with MUE Stage 1. Lambda_g = ")
    print(lambda.g)

    print("running Stage 2")
    # Running the stage 2 model
    out.stage2 <- rstar.stage2(log.output,
                               inflation,
                               real.interest.rate,
                               lambda.g,
                               a3.constraint,
                               b2.constraint)
    print("Done Estimating Stage 2, --> MUE Stage 2")
    # Median unbiased estimate of lambda_z
    lambda.z <- median.unbiased.estimator.stage2(out.stage2$y, out.stage2$x)
#    if (lambda.z == 0) {
#       browser()
#      lambda.z = 1e-03
#    }
    print("Done with MUE Stage 2. lambda_z = ")
    print(lambda.z)

    print("running Stage 3")
    # Running the stage 3 model
    out.stage3 <- rstar.stage3(log.output,
                               inflation,
                               real.interest.rate,
                               nominal.interest.rate,
                               lambda.g,
                               lambda.z,
                               a3.constraint,
                               b2.constraint,
                               run.se)
                               
    print("Done Estimating Stage 3")                               
                               
##    # Running the stage 3 model with diffuse prior
##    out.stage3 <- rstar.stage3(log.output,
##                               inflation,
##                               real.interest.rate,
##                               nominal.interest.rate,
##                               lambda.g,
##                               lambda.z,
##                               a3.constraint,
##                               b2.constraint,
##                               run.se)                               

    return(list(out.stage1=out.stage1,
                out.stage2=out.stage2,
                out.stage3=out.stage3,
                lambda.g=lambda.g,lambda.z=lambda.z))

##    return(list(out.stage1=out.stage1,
##                out.stage2=out.stage2,
##                out.stage3=out.stage3,
##                out.stage3.diffuse=out.stage3.diffuse,
##                lambda.g=lambda.g,lambda.z=lambda.z))                
}
