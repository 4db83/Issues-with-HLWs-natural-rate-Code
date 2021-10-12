# ---------------------------------------------------------------------------------------------            
# File:        run.hlw.R modfied, to run over a longer period with extra output
# ---------------------------------------------------------------------------------------------            
# THIS THE SAME AS THE SCRIPT.RUN.HLW.R BUT ALSO DOES PLOTTING 
# DB: FIRST RUN SCRIPT SCRIPT.PREPARE.RSTAR.DATA.US.R WITH SAVE DATA STATEMENT SET TO 1.
# ORIGINAL SAMPLE *****************************************************************************                                            
# sample.start      <- c(1961,1)            ## "1960-01-01"  1960-01-01 1960:Q1   
# sample.end        <- c(2017,1)            ## "1960-04-01"  1960-04-01 1960:Q2   
# sample.end.string <- '2017Q1'             ## "1960-07-01"  1960-07-01 1960:Q3   
# "1960-10-01"  1960-10-01 1960:Q4   
# for(lib in pkgs) {remove.packages(lib) }  # to remove all pkgs 
###############################################################################################################

rm(list = ls())
# CLEAR THE CONSOLE
cat("\014")
# define some colors for plotting lager
grn1 = rgb(r=.2,g=.99,b=.5); grn2 = rgb(r=.2,g=.5,b=.1); blu1 = rgb(r=.2,g=.5,b=.99); red1 = rgb(r=.99,g=.4,b=.4);

# ---------------------------------------------------------------------------------------------            
# GET REQUIRED LIBRARIES 
# ---------------------------------------------------------------------------------------------            
# List of required packages
pkgs <- c('mFilter','tis','nloptr','tictoc','lubridate','zoo','MASS')
# get packages that are not installed
get.pkgs <- pkgs[!pkgs %in% installed.packages()]
if(!length(get.pkgs)==0){cat('Getting packages\n'); for(lib in get.pkgs) {install.packages(lib, dependencies = TRUE) } }
sapply( get.pkgs, require, character = TRUE )
# load packages
cat('  Loading packages ... \n')
for(ii in pkgs) library(ii, character.only = TRUE)

# ---------------------------------------------------------------------------------------------            
# LOAD HELPER FUNCTIONS STORE IN PATH.2.FUNCTIONS
# ---------------------------------------------------------------------------------------------            
path.2.functions  = c("./R.local.Functions/")
source_files      = list.files(path.2.functions, "*.R")
sapply( paste0(path.2.functions, source_files), source )

# ---------------------------------------------------------------------------------------------            
# START THE TIMER
# ---------------------------------------------------------------------------------------------            
tic("all Stage")
# increase max print option
options(max.print=10000)
options(scipen = 999)
options(digits = 14)

# ---------------------------------------------------------------------------------------------            
# SET TO ONE TO SAVE TO CSV FILE THE REQUIRED INPUTS
SAVE.RESULTS.2.FILE = 0
PLOT.OUTPUT         = 0
WRITE.OUTPUT.TXT    = 0

# write the results to txt. file 
if (WRITE.OUTPUT.TXT) { sink("full_estimation_output.txt") }
# if (WRITE.OUTPUT.TXT) { sink("short_estimation_output_beg_sample_1972.txt") }

# CHOOSE THE COUNTRY TO FIT +++++++ 
COUNTRIES	= c('US','EA','UK','CA')
# COUNTRIES = 'US'

for(country in COUNTRIES) {
  
ptm <- proc.time()

# estimation dates, date input etc. 
# ---------------------------------------------------------------------------------------------            
  # SET END DATES OF THE ESTIMATION SAMPLE (FORMAT IS C(YEAR,QUARTER))
  sample.end  <- c(2019,4)            ## "1960-04-01"  1960-04-01 1960:Q2
  
  # OUTPUT FILE DIRECTORY 
  if (!dir.exists('../data/')) {dir.create('../data/')}
  RESULTS.OUTPUT.DIR  = '../data/R.HLW.results/'
  # CHECK IF EXISTS, IF NOT MAKE, IT.
  if (!dir.exists(RESULTS.OUTPUT.DIR)){dir.create(RESULTS.OUTPUT.DIR)}
  
  # SPECIFY THE XLS DATA VINTAGE FILE TO READ IN ************************************************
  data.end.string = '2019Q4'
  FILE.INPUT.DIR  = '../data/R.data.for.estimation.2020.May.28/'
    
  # now load the data 
  if (country=='EA') { DATA.FILE.2.BE.USED = paste0(FILE.INPUT.DIR, 'EA.data.csv'); sample.start <- c(1972,1) }
  if (country=='US') { DATA.FILE.2.BE.USED = paste0(FILE.INPUT.DIR, 'US.data.csv'); sample.start <- c(1961,1) }
  if (country=='UK') { DATA.FILE.2.BE.USED = paste0(FILE.INPUT.DIR, 'UK.data.csv'); sample.start <- c(1961,1) }
  if (country=='CA') { DATA.FILE.2.BE.USED = paste0(FILE.INPUT.DIR, 'CA.data.csv'); sample.start <- c(1961,1) }

  ##  set sample start to 1972 for all countries to show influence of z(t)
  # sample.start <- c(1972,1)
  # ---------------------------------------------------------------------------------------------
  # now in string and date formats
  # ---------------------------------------------------------------------------------------------
  sample.start.string <- paste0( toString(sample.start[1]),"Q",toString(sample.start[2]) )
  sample.end.string   <- paste0( toString(sample.end[1]),"Q",toString(sample.end[2]) )            
  asDate.sample.end.string = as.Date(as.yearqtr(sample.end.string, format = "%YQ%q"))
  # THE ESTIMATION PROCESS USES DATA BEGINNING 4 QUARTERS PRIOR TO THE SAMPLE START
  pre.sample.start    <- shiftQuarter(sample.start,-4)
  # print(data.start)
  # MAKE SAMPLE START STRING TO BE USED TO TRIM THE SAMPLE BELOW
  pre.sample.start.string <- paste0( toString(pre.sample.start[1]),"Q",toString(pre.sample.start[2]) )
  asDate.pre.sample.start.string =  as.Date(as.yearqtr(pre.sample.start.string, format = "%YQ%q"))
  # as.yearqtr(sample.start.string, format = "%YQ%q")
  cat("--------------------------------------------------------------------------------------------\n")
  cat(paste0("   (pre-sample: ", pre.sample.start.string, ") Country: ", country ," Estimation Period from ", sample.start.string, " to ", sample.end.string))
  cat("\n--------------------------------------------------------------------------------------------\n")
  
  # FILE NAME OF OUTPUT FILE 
  Name.end = sample.start.string
  # Name.end      = '.Lambda.z'  
  # FILE.NAME.OUT = paste0(data.vintage.to.use, Name.end)
  # FILE.NAME.OUT = paste0(Name.end,'.to.',sample.end.string,'.dataVintage.',data.vintage.to.use)
  FILE.NAME.OUT = paste0(Name.end,'.to.',sample.end.string)
  # FILE.NAME.OUT = sample.end.string
  
  #**********************************************************************************************
  
  # ---------------------------------------------------------------------------------------------            
  # CHECK IF EXISTS, IF NOT MAKE, IT.
  DIR.NAME.OUT = paste0(RESULTS.OUTPUT.DIR, country, '/')
  if (!dir.exists(DIR.NAME.OUT)){dir.create(DIR.NAME.OUT)}
  # OUTPUT DIR NAMES:
    # OTHER OUTPUT NAMES
  Factors_one_sided         = paste0(DIR.NAME.OUT, 'Stage3.Filtered.',                        FILE.NAME.OUT, '.csv')
  Factors_two_sided         = paste0(DIR.NAME.OUT, 'Stage3.Smoothed.',                        FILE.NAME.OUT, '.csv')
  HLW_results_all           = paste0(DIR.NAME.OUT, 'Stage3.Full.Estimation.Results.',         FILE.NAME.OUT, '.csv')
  # OUTPUT FILE NAMES FOR DIFFUSE RESULTS:
  Factors_one_sided_diffuse = paste0(DIR.NAME.OUT, 'Stage3.Filtered.diffuse.',                FILE.NAME.OUT, '.csv')
  Factors_two_sided_diffuse = paste0(DIR.NAME.OUT, 'Stage3.Smoothed.diffuse.',                FILE.NAME.OUT, '.csv')
  HLW_results_all_diffuse   = paste0(DIR.NAME.OUT, 'Stage3.Full.Estimation.Results.diffuse.', FILE.NAME.OUT, '.csv')
  
  #**********************************************************************************************
  # READ IN DATA, RUN ESTIMATION, AND SAVE OUTPUT --
  #**********************************************************************************************
  data.Full <- read.table(  sep = ',', na.strings = ".", header=TRUE, stringsAsFactors=FALSE,
                              file = DATA.FILE.2.BE.USED )
  # head2tail(data.Full)
  
  #================================================================================================
  # sometimes the dates read by read.table are mixed around, leading to different formats for the 
  # dates between the different data files although they are the same format in the excel/csv file. 
  # this can happen when opening the csv file in excel and then saving it there and closing.
  #================================================================================================
  
  # CHANGE THE DATE FORMAT AS IN THE EXCEL FILE ----
  # dates.csv = as.Date(final.data.Full$Date,"%d.%m.%Y")   # use this with the OLD data 
  dates.csv = as.Date(data.Full$Date,"%Y-%m-%d")     # with the NEW data 
  # head2tail(dates.csv)
  
  # NOW TRIM/SELECT A SUBSAMPLE OF THE DATA testset[date>="2013-08-02" & date<="2013-11-01"]
  Subset.data = dates.csv >= as.Date(asDate.pre.sample.start.string) & 
                dates.csv <= as.Date(asDate.sample.end.string)

  # define the subset of data to be used
  final.data = data.Full[Subset.data,]
  # head2tail(final.data)
 
  # DEFINE VARIABLES
  log.output               <- final.data$gdp.log
  inflation                <- final.data$inflation
  inflation.expectations   <- final.data$inflation.expectations
  nominal.interest.rate    <- final.data$interest.rate
  real.interest.rate       <- nominal.interest.rate - inflation.expectations
  
    # Upper bound on a_3 parameter (slope of the IS curve)
  a3.constraint <- -0.0025
  # Lower bound on b_2 parameter (slope of the Phillips curve)
  b2.constraint <-  0.0250

  # Set start index for y
  g.pot.start.index <- 1 + ti(shiftQuarter(sample.start,-3),'quarterly')-ti(pre.sample.start,'quarterly')
  # Set column names for CSV output
  output.col.names  <- c("Date","rstar","g","z","output gap","","All results are output from the Stage 3 model.",rep("",8),"Standard Errors","Date","y*","r*","g","","rrgap")
  
  # set to false to not do the simulations.
  run.se = FALSE
  
  # ---------------------------------------------------------------------------------------------            
  # RUN HLW ESTIMATION FOR THE US MAIN CALLING FUNCTION ----
  # ---------------------------------------------------------------------------------------------            
  # MAIN RESULTS
  tic('Run all')
  hlw.full  <- run.hlw.estimation(log.output, inflation, real.interest.rate, nominal.interest.rate,
                  a3.constraint = a3.constraint, b2.constraint = b2.constraint, run.se = run.se)
  toc()
  
  # One-sided (filtered) estimates 
  one.side.est <- cbind(hlw.full$out.stage3$rstar.filtered,
                        hlw.full$out.stage3$trend.filtered,
                        hlw.full$out.stage3$z.filtered,
                        hlw.full$out.stage3$output.gap.filtered)
  
  # Two-sided (smoothed) estimates
  two.side.est <- cbind(hlw.full$out.stage3$rstar.smoothed,
                        hlw.full$out.stage3$trend.smoothed,
                        hlw.full$out.stage3$z.smoothed,
                        hlw.full$out.stage3$output.gap.smoothed)
   
  # ---------------------------------------------------------------------------------------------            
  # MAKE DATE SEQUENCE AND OUTPUT FORMAT
  # ---------------------------------------------------------------------------------------------            
  dates <- seq(from=( as.Date(ti(shiftQuarter(sample.start,-1),'quarterly'))+1), to = (
                      as.Date(ti(shiftQuarter(sample.end,-1),tif='quarterly'))+1), by = 'quarter')
  
  # MAKE FORMATED OUTPUT FOR PRINTING
  output.all <- format.output(hlw.full, one.side.est, real.interest.rate, sample.start, sample.end, run.se = run.se)
  
  # ---------------------------------------------------------------------------------------------            
  # Make all output required for MATLAB input 
  # ---------------------------------------------------------------------------------------------            
  # STAGE 1 ----
  Stage1.names = c( " a_y1            ",
                    " a_y2            ",
                    " b_pi            ",
                    " b_y             ",
                    " gama            ",
                    " sigma_y~        ",
                    " sigma_pi        ",
                    " sigma_y*        ",
                    " Log-Likelihood  ",                 
                    " Lambda.g        ")
  
  # combine Stage1 output
  Stage1.output = cbind (
    c(hlw.full$out.stage1$initial.parameters , hlw.full$out.stage1$log.likelihood.initialvals , 0           ),
    c(hlw.full$out.stage1$theta						   ,	hlw.full$out.stage1$log.likelihood  					, hlw.full$lambda.g)
  ) 
  # add row/colnames
  rownames(Stage1.output) = Stage1.names
  colnames(Stage1.output) = c(" init.Values", " HLW")
  cat(" -----------------------------------------------------------------------------------------------------\n")
  cat("                                         Stage 1 results")
  cat("\n -----------------------------------------------------------------------------------------------------\n")
  print( round(Stage1.output, digits = 10) )

  # make trend.smoothed and  get the full Lambda.g output and MUE.Stats, 
  trend.smoothed = hlw.full$out.stage1$potential.smoothed
  Dates = data.frame(dates);

  Stage1.lambdas  <- function.median.unbiased.estimator.stage1.correct.order(trend.smoothed)
  Stage1.MUE      <- cbind(Stage1.lambdas$Lambdas,Stage1.lambdas$Statistics)
  colnames(Stage1.MUE) = c("Lambda.g","MUE.Stat")

  # STAGE 2 ----------
  Stage2.names = c( " a_y1            ",  # 1   
                    " a_y2            ",  # 2
                    " a_r             ",  # 3
                    " a_0             ",  # 4
                    " a_g             ",  # 5
                    " b_pi            ",  # 6
                    " b_y             ",  # 7
                    " sigma_y~        ",  # 8
                    " sigma_pi        ",  # 9
                    " sigma_y*        ",  # 10
                    " sigma_g         ",  # 11
                    " Log-Likelihood  ",                 
                    " Lambda.g        ")
  
  #*************************************************************************************************************  
  ## 																	 Stage 2 model \sigma_g by MLE 
  #*************************************************************************************************************  
  options(scipen = 999)
  # tic("Stage 2 only")
  lambda.g = hlw.full$lambda.g
  # source('extra.functions.Stage.2.R')  
  
  # Stage(2) model with \sigma_g estimated by MLE
  stage2.g 	<- fit.stage2.g(  log.output, inflation, real.interest.rate,
                              NA, a3.constraint, b2.constraint, 
                              hlw.full$out.stage2$P.00, c(hlw.full$out.stage2$theta, 0.01) )

  # correct Stage(2) model with \sigma_g estimated by MLE  
  # use the Stage(2).g theta as initial values, remove a_0 and a_g estimates from stage2.g$theta
  stage2.M0g <- fit.stage2.M0g( log.output, inflation, real.interest.rate,
                              	NA, a3.constraint, b2.constraint,
                              	hlw.full$out.stage3$P.00[1:5,1:5], stage2.g$theta[-c(4,5)] ) 

  # combine Stage2 output for printing 
  #*************************************************************************************************************  
  # make inline functions to compute \sigma_g = \lambda_g*sigma_y* and \lambda_g = \sigma_g/\sigma_y* (use NaN insteast of NA for matlab compatibility in reading)
  sg2 = function(theta,Lambda_g){theta[10]*Lambda_g}
  Lg2 = function(theta){theta[11]/theta[10]}
  # compbine Stage 2 output from the different models
    Stage2.output = cbind (
    c(hlw.full$out.stage2$initial.parameters, NaN, hlw.full$out.stage2$log.likelihood.initialvals, NaN           ),
    c(hlw.full$out.stage2$theta	,	sg2(hlw.full$out.stage2$theta, hlw.full$lambda.g), hlw.full$out.stage2$log.likelihood , hlw.full$lambda.g),
    c(stage2.g$theta      , stage2.g$log.likelihood   , Lg2(stage2.g$theta)),
    c(stage2.M0g$theta[1:3], NaN,NaN, stage2.M0g$theta[4:9], stage2.M0g$log.likelihood   , stage2.M0g$theta[9]/stage2.M0g$theta[8])    
  )
  
  # add row/colnames
  rownames(Stage2.output) = Stage2.names
  colnames(Stage2.output) = c(" init.Values", " HLW-Baseline", "MLE.sigma_g", "MLE.M0.sigma_g")
  
  # print to screen with rounding 
  cat(" -----------------------------------------------------------------------------------------------------\n")
  cat("                                         Stage 2 results")
  cat("\n -----------------------------------------------------------------------------------------------------\n")
  print( round(Stage2.output, digits = 10) )
  
  # MAKE THE SMOOTHED STATE VECTOR, IF NEEDED FOR LATER ANALYSIS.
  Stage2.states.smoothed  = hlw.full$out.stage2$states$smoothed$xi.tT
  colnames(Stage2.states.smoothed) = c("y*(t)", "y*(t-1)", "y*(t-2)", "g(t-1)")
  rownames(Stage2.states.smoothed) = as.matrix(final.data[1])[5:dim(final.data)[1]]   # adding the dates

  #*************************************************************************************************************  
  #				TIME VARYING PHI VERSION: MAKE  STAGE 2 THETA VECTOR AND BASELINE GY AND XX SERIES AS IN HLW 
  #*************************************************************************************************************
  #                             function.MUE.S2 works correctly!!!
  #*************************************************************************************************************
  ### this reads in the yx data from Matlabs output as a check to the function.MUE.S2
  # check.yx = read.table(file = '../data/yx_M0g.csv', sep = ',', header=FALSE)
  # check.MUE.S2.TV.out = function.MUE.S2( as.matrix(check.yx[,1]), as.matrix(check.yx[,c(2,3,4)]) )
  # check.out = cbind(check.MUE.S2.TV.out$Lambdas, check.MUE.S2.TV.out$Statistics)
  # (round(check.out, digits = 8))
  #*************************************************************************************************************

  # Y and X inputs from their baseline Stage 2 model; --> input into the function.MUE.S2
  MUE.S2.GY 		= hlw.full$out.stage2$y;  
  MUE.S2.XX 		= hlw.full$out.stage2$x;
  # Y and X inputs from their baseline Stage 2 model (estimated sigma_g); --> input into the function.MUE.S2
  MUE.S2.g.GY 	= stage2.g$y;     
  MUE.S2.g.XX 	= stage2.g$x;
  # Y and X inputs from their CORRECT (M0) Stage 2 model (estimated sigma_g); --> input into the function.MUE.S2
  MUE.S2.M0g.GY = stage2.M0g$y;        
  MUE.S2.M0g.XX = stage2.M0g$x;
  
  ## Call to function.MUE.S2
  MUE.S2.TV.out      = function.MUE.S2(MUE.S2.GY, MUE.S2.XX)
  MUE.S2.TV.g.out    = function.MUE.S2(MUE.S2.g.GY, MUE.S2.g.XX)
  MUE.S2.TV.M0g.out  = function.MUE.S2(MUE.S2.M0g.GY, MUE.S2.M0g.XX)
  
  # MAKE OUTPUT FOR PRINTING (lambdas and F-stats)
  MUE.S2.Lambda_z.TV.out = cbind( MUE.S2.TV.out$Lambdas,
                                  MUE.S2.TV.g.out$Lambdas,
                                  MUE.S2.TV.M0g.out$Lambdas )
  colnames(MUE.S2.Lambda_z.TV.out) = c("Original","S2.MLE(g)","S2.M0.MLE(g)")  
  
  cat(" -----------------------------------------------------------------------------------------------------\n")
  cat("                     MUE Stage 2 Lambda_z Time varying Phi as in HLW")
  cat("\n -----------------------------------------------------------------------------------------------------\n")
  print( round(MUE.S2.Lambda_z.TV.out, digits = 10) )
  
  MUE.S2.Fstats.TV.out = cbind( MUE.S2.TV.out$Statistics,
                                MUE.S2.TV.g.out$Statistics,
                                MUE.S2.TV.M0g.out$Statistics )
  colnames(MUE.S2.Fstats.TV.out) = colnames(MUE.S2.Lambda_z.TV.out)
  
  cat(" -----------------------------------------------------------------------------------------------------\n")
  cat("                     MUE Stage 2 Structural Break F-Stats Time varying Phi as in HLW")
  cat("\n -----------------------------------------------------------------------------------------------------\n")
  print( round(MUE.S2.Fstats.TV.out, digits = 10) )
  cat(" -----------------------------------------------------------------------------------------------------\n")
  
  # -------------------------------------------------------------------------------------------------------------
  # 				CONSTANT PHI VERSION: MAKE  STAGE 2 THETA VECTOR AND BASELINE GY AND XX SERIES AS IN HLW 
  # -------------------------------------------------------------------------------------------------------------
  # Stage2.names = c( " a_y1            ", 1
  #                   " a_y2            ", 2 
  #                   " a_r             ", 3 
  #                   " a_0             ", 4
  #                   " a_g             ", 5
  # Constant phi -> GY is made ONCE outside the Structural break loop from fixed theta2 vector, need to re-order theta
  # NOTE:theta[c(1,2,3,5,4)] due to parameter ordering (must re-order), since x = [y*(t-1), y*(t-2), r(t-1) g(t-1) 1]
  # -------------------------------------------------------------------------------------------------------------
  
  # MAKE A VECTOR OF ONES AS THE CONSTANT/INTERCEPT PARAMETER IN THE CONSTANT PHI REGRESSIONS. 
  # This will be the same for all three scenarios
  II = as.matrix(hlw.full$out.stage2$x[,5])
  
  # construct the various GY inputs from the constant phi models
  MUE.S2.Cnst.GY  		= hlw.full$out.stage2$y - hlw.full$out.stage2$x %*% hlw.full$out.stage2$theta[c(1,2,3,5,4)]
  MUE.S2.Cnst.g.GY  	= stage2.g$y 						- stage2.g$x 						%*% stage2.g$theta[						c(1,2,3,5,4)]
  MUE.S2.Cnst.M0g.GY  = stage2.M0g$y 					- stage2.M0g$x 					%*% stage2.M0g$theta[					c(1,2,3)		]

 	# compute the MUE.2 output from the different constant Phi models  
  MUE.S2.Cnst.out 		= function.MUE.S2(MUE.S2.Cnst.GY,    	II)
  MUE.S2.Cnst.g.out 	= function.MUE.S2(MUE.S2.Cnst.g.GY,  	II)
  MUE.S2.Cnst.M0g.out = function.MUE.S2(MUE.S2.Cnst.M0g.GY, II)

  # MAKE OUTPUT FOR PRINTING (lambdas and F-stats)
  MUE.S2.Lambda_z.Cnst.out = cbind( MUE.S2.Cnst.out$Lambdas,
                                    MUE.S2.Cnst.g.out$Lambdas,
                                    MUE.S2.Cnst.M0g.out$Lambdas
                                  )

  colnames(MUE.S2.Lambda_z.Cnst.out) = colnames(MUE.S2.Lambda_z.TV.out)
 	cat("                             MUE Stage 2 Lambda_z Constant Phi")
 	cat("\n -----------------------------------------------------------------------------------------------------\n")
  print( round(MUE.S2.Lambda_z.Cnst.out, digits = 10) )
  
  MUE.SE.Fstats.Cnst.out <- cbind( MUE.S2.Cnst.out$Statistics, 
                                   MUE.S2.Cnst.g.out$Statistics,
                                   MUE.S2.Cnst.M0g.out$Statistics
                                  )
  
  colnames(MUE.SE.Fstats.Cnst.out) = colnames(MUE.S2.Lambda_z.TV.out)

  cat(" -----------------------------------------------------------------------------------------------------\n")
  cat("                             MUE Stage 2 Structural Break F-Stats Constant Phi")
  cat("\n -----------------------------------------------------------------------------------------------------\n")
	print( round(MUE.SE.Fstats.Cnst.out, digits = 10) )
	cat(" -----------------------------------------------------------------------------------------------------\n")

	#***************************************************************************************************************************
  # SOME PLOTS OF THE FSTATISTICS 
  #***************************************************************************************************************************
	# baseline TV Phi versions Stage 2 and M0g models----
	par(fig = c(0,1,.45,1))
  plot(	MUE.S2.TV.out$Fstat  ,   x = dates[4:(length(dates)-4)], col=red1, type="l", lwd=3, lty=1, ylim = c(0,14), xlab="", ylab= "", las=1, panel.first=grid() )
  mtext(side = 1, text = "HLW Misspecified Stage 2 Model", line = 2.2)
  lines(MUE.S2.Cnst.out$Fstat, x = dates[4:(length(dates)-4)], col=blu1, type="l", lwd=3, lty=1)
  legend( "topright",legend = c("Time varying", "Constant"), col=c(red1,blu1), lty=1:1, lwd=3, cex=1.0 )
  
  # baseline TV Phi versions Stage 2 and M0g models
  par(fig = c(0,1,.05,0.55), new = TRUE) 
  plot(	MUE.S2.TV.M0g.out$Fstat  ,   x = dates[4:(length(dates)-4)], col=red1, type="l", lwd=3, lty=1, ylim = c(0,07), xlab="", ylab= "", las=1, panel.first=grid())
  mtext(side = 1, text = "Correctly specified Stage 2 Model", line = 2.2)
  lines(MUE.S2.Cnst.M0g.out$Fstat, x = dates[4:(length(dates)-4)], col=blu1, type="l", lwd=3, lty=1,)
  legend( "topright",legend = c("Time varying", "Constant"), col=c(red1,blu1), lty=1:1, lwd=3, cex=1.0 )

  #***************************************************************************************************************************
  # 																					STAGE 3                                                                               ----
  #***************************************************************************************************************************
  Stage3.names = c( " a_y1              ",      # 1
                    " a_y2              ",      # 2
                    " a_r               ",      # 3
                    " b_pi              ",      # 4
                    " b_y               ",      # 5
                    " sigma_y~          ",      # 6
                    " sigma_pi          ",      # 7
                    " sigma_y*          ",      # 8
                    " sigma_g(implied)  ",      # 9
                    " sigma_z(implied)  ",      # 10                  
                    " Log-Likelihood    ",                 
                    " Lambda.g(implied) ",
                    " Lambda.z(implied) "  )
  
  # source('extra.functions.Stage.3.R')  
  
  # DEFINE WHICH LAMBDA.Z TO USE, FIRST USE THE HLW ESTIMATES
  lambda.z.HLW  = hlw.full$lambda.z
  lambda.z.M0g  = MUE.S2.Cnst.M0g.out$Lambdas[2]
  
	# --------------------------------------------------------------------------------------------------- #	
  # Stage 3 model with \sigma_g estimated by MLE and lambda.z from wrong HLW Stage 2 MUE
  # --------------------------------------------------------------------------------------------------- #
  stage3.g    <- fit.stage3.g(  log.output, inflation, real.interest.rate,
                                NA, lambda.z.HLW, a3.constraint, b2.constraint, 
                                hlw.full$out.stage3$P.00, c(hlw.full$out.stage3$theta, .03) )
	# --------------------------------------------------------------------------------------------------- #
  # Stage 3 model with \sigma_g estimated by MLE and lambda.z from correct M0 Stage 2 MUE (uses the same fit.stage3.g function)
  # --------------------------------------------------------------------------------------------------- #
  stage3.M0g  <- fit.stage3.g(  log.output, inflation, real.interest.rate,
                                NA, lambda.z.M0g, a3.constraint, b2.constraint, 
                                hlw.full$out.stage3$P.00, stage3.g$theta )
  # --------------------------------------------------------------------------------------------------- #
  # Stage 3 model with \sigma_g and \sigma_z estimated by MLE, no use of MUE in Stages 1 or 2.
  # --------------------------------------------------------------------------------------------------- #
  stage3.gz   <- fit.stage3.gz( log.output, inflation, real.interest.rate,
                                NA, NA, a3.constraint, b2.constraint, 
                                hlw.full$out.stage3$P.00, c(stage3.g$theta, .01) )
  
  # cbind(stage3.g$theta, stage3.M0g$theta, stage3.gz$theta)

  # COMBINE STAGE3 OUTPUT  ----
  Stage3.output = cbind (
    c(hlw.full$out.stage3$initial.parameters , NA, NA, hlw.full$out.stage3$log.likelihood.initialvals , NA ,NA),
    c(hlw.full$out.stage3$theta,	
      hlw.full$lambda.g*hlw.full$out.stage3$theta[8], 
      hlw.full$lambda.z*hlw.full$out.stage3$theta[6]/abs(hlw.full$out.stage3$theta[3]),
      hlw.full$out.stage3$log.likelihood, hlw.full$lambda.g , hlw.full$lambda.z),
    c(stage3.g$theta, lambda.z.HLW*stage3.g$theta[6]/abs(stage3.g$theta[3]),
      stage3.g$log.likelihood , stage3.g$theta[9]/stage3.g$theta[8], lambda.z.HLW),
    c(stage3.M0g$theta, lambda.z.HLW*stage3.M0g$theta[6]/abs(stage3.M0g$theta[3]),
      stage3.M0g$log.likelihood , stage3.M0g$theta[9]/stage3.M0g$theta[8], lambda.z.M0g),
    c(stage3.gz$theta, stage3.gz$log.likelihood , stage3.gz$theta[9]/stage3.gz$theta[8],
      abs( stage3.gz$theta[10]*stage3.gz$theta[3] / stage3.gz$theta[6] ) )
  )
  
  # add row/colnames
  rownames(Stage3.output) = Stage3.names
  colnames(Stage3.output) = c(" init.Values", "HLW-replicated", "S3_MLE.g", "S3_MLE.M0g", "S3_MLE.gz")
  # colnames(Stage3.output) = c(" init.Values", "HLW-replicated", "S3_MLE.g")
  cat("                                         Stage 3 results")
  cat("\n -----------------------------------------------------------------------------------------------------\n")
  print( round(Stage3.output, digits = 10) )
  cat(" -----------------------------------------------------------------------------------------------------\n")

  
  # WRITE Y AND X USED IN MUE STAGE 2 TO CSV FILE 
  Stage2yX4MUE = cbind(hlw.full$out.stage2$y, hlw.full$out.stage2$x)
  colnames(Stage2yX4MUE) = c( "y~(t)","y~(t-1)","y~(t-2)","[r(t-1)+r(t-2)]/2","g(t-1)","C")
  rownames(Stage2yX4MUE) = as.matrix(final.data[1])[5:dim(final.data)[1]]   # adding the dates
  
  # MAKE STAGE2 MUE OUTPUT FOR PRINTING TO CSV FILE
  Stage2.MUE <- cbind(MUE.S2.TV.out$Lambdas,MUE.S2.TV.out$Statistics)
  colnames(Stage2.MUE) = c("Lambda.z","MUE.Stat")

  # PRINT TO CSV FILE -----------------
  if (SAVE.RESULTS.2.FILE) {
    ## STAGE1.OUTPUT ## ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    write.table(Stage1.output, row.names=TRUE, col.names = TRUE, sep = ",",
                file = paste0(DIR.NAME.OUT, 'Stage1.theta.', FILE.NAME.OUT, '.csv') )
    # Stage1.P00
    write.table(hlw.full$out.stage1$P.00, row.names = FALSE, col.names = FALSE, sep = ",",
                file = paste0(DIR.NAME.OUT, 'Stage1.P00.', FILE.NAME.OUT, '.csv') )
    # Stage1.xi00
    write.table(hlw.full$out.stage1$xi.00, row.names = FALSE, col.names = FALSE, sep = ",",
                file = paste0(DIR.NAME.OUT, 'Stage1.xi00.', FILE.NAME.OUT, '.csv') )
    # Stage1.Trend.smoothed with dates
    write.table(cbind(Dates,trend.smoothed), row.names = FALSE, col.names = TRUE, sep = ",",
                file = paste0(DIR.NAME.OUT, 'Stage1.trend.smoothed.', FILE.NAME.OUT, '.csv') )
    # Stage1.Lambda.z and MUE.stats
    write.table(Stage1.MUE, row.names = TRUE, col.names = TRUE, sep = ",",
                file = paste0(DIR.NAME.OUT, 'Stage1.Lambda.g.', FILE.NAME.OUT, '.csv') )
    ## STAGE2.OUTPUT ## ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    write.table(Stage2.output, row.names=TRUE, col.names = TRUE, sep = ",",
                file = paste0(DIR.NAME.OUT, 'Stage2.theta.', FILE.NAME.OUT, '.csv') )
    # Stage2.P00
    write.table(hlw.full$out.stage2$P.00, row.names = FALSE, col.names = FALSE, sep = ",",
                file = paste0(DIR.NAME.OUT, 'Stage2.P00.', FILE.NAME.OUT, '.csv') )
    # Stage2.xi00
    write.table(hlw.full$out.stage2$xi.00, row.names = FALSE, col.names = FALSE, sep = ",",
                file = paste0(DIR.NAME.OUT, 'Stage2.xi00.', FILE.NAME.OUT, '.csv') )
    # Stage2.Stage2.states.smoothed
    write.table(cbind(Dates,Stage2.states.smoothed), row.names = FALSE, col.names = TRUE, sep = ",",
                file = paste0(DIR.NAME.OUT, 'Stage2.states.smoothed.', FILE.NAME.OUT, '.csv') )
    # Stage2.Lambda.z and MUE.stats
    write.table(Stage2.MUE, row.names = TRUE, col.names = TRUE, sep = ",",  
                file = paste0(DIR.NAME.OUT, 'Stage2.Lambda.z.', FILE.NAME.OUT, '.csv') )  
    # Stage2.yX.4.MUE
    write.table(cbind(Dates,Stage2yX4MUE), row.names = FALSE, col.names = TRUE, sep = ",", 
                file = paste0(DIR.NAME.OUT, 'Stage2.yX.4.MUE.', FILE.NAME.OUT, '.csv') )
    ## STAGE3.OUTPUT ## ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    write.table(Stage3.output, row.names=TRUE, col.names = TRUE, sep = ",",
                file = paste0(DIR.NAME.OUT, 'Stage3.theta.', FILE.NAME.OUT, '.csv') )
    # Stage3.P00
    write.table(hlw.full$out.stage3$P.00, row.names = FALSE, col.names = FALSE, sep = ",",
                file = paste0(DIR.NAME.OUT, 'Stage3.P00.', FILE.NAME.OUT, '.csv') )
    # Stage3.xi00
    write.table(hlw.full$out.stage3$xi.00, row.names = FALSE, col.names = FALSE, sep = ",",
                file = paste0(DIR.NAME.OUT, 'Stage3.xi00.', FILE.NAME.OUT, '.csv') )
                
    # *****************************************************************************************************
    # SAVE ALL STATE 3 ESTIMATES TO CSV THESE ARE THE BASELINE HLW ESTIMATES THAT ARE PRINTED
    # *****************************************************************************************************
    # one sided estimation filtered results
    write.table(one.side.est, Factors_one_sided, row.names = dates, quote = FALSE, sep = ',', na = ".", 
                col.names = c("rstar","g","z","output gap"))
    
    # two sided estimation smoothed results
    write.table(two.side.est, Factors_two_sided, row.names = dates, quote = FALSE, sep = ',', na = ".", 
                col.names = c("rstar","g","z","output gap"))
    
    # factors all results
    write.table(output.all, HLW_results_all, quote=FALSE, row.names=FALSE, sep = ',', na = '', 
                col.names = output.col.names)
    
    cat("Stage 3 Estimates written to .csv file \n")
  }
  
  ## ***************************************************************************************************
  ## Make smoothed plots ------
  ## ***************************************************************************************************
  if (PLOT.OUTPUT) {
  YLims = c(-2,6)
  # YLims = c(-2,2)
  
  ## par(mfrow=c(3,1))
  #op <- par(mfrow = c(3,1),
  #          oma = c(5,4,0,0) + 0.5,
  #          mar = c(0,0,-1.5,1) + 3.0,
  #          mgp = c(2, .8, 0))   
  #

  smoothed.rstar  = hlw.full$out.stage3$rstar.smoothed
  filtered.rstar  = hlw.full$out.stage3$rstar.filtered
  smoothed.g      = hlw.full$out.stage3$trend.smoothed
  filtered.g      = hlw.full$out.stage3$trend.filtered
  smoothed.z      = hlw.full$out.stage3$z.smoothed
  filtered.z      = hlw.full$out.stage3$z.filtered

  # now add the diffuse results.
  # diff.filtered = diffuse.prior.stage3trendd.filtered
  # diff.smoothed = diffuse.prior.stage3$trend.smoothed
  
  # NOW DO THE PLOTTING RSTAR
  pdf( paste0(country, "_", sample.start.string, "_plots.pdf") , width=8, height=5)
  ## smoothed
  # plot(	y = smoothed.rstar   ,  x = dates, col = blu1, type="l", lwd = 3, lty=1, xlab = "", ylab = "", xaxt = "n" , ylim = YLims, las=1) ## the last bit xaxt = "n"  kills the x-axis labels
  # lines(y = smoothed.g       ,  x = dates, col = red1, type="l", lwd = 3, lty=1)
  # lines(y = smoothed.z       ,  x = dates, col = grn2, type="l", lwd = 3, lty=1)
  ## filtered
  plot(	y = filtered.rstar   ,  x = dates, col = blu1, type="l", lwd = 3, lty=1, xlab = "", ylab = "", xaxt = "n" , ylim = YLims, las=1) ## the last bit xaxt = "n"  kills the x-axis labels
  lines(y = filtered.g       ,  x = dates, col = red1, type="l", lwd = 3, lty=2)
  lines(y = filtered.z       ,  x = dates, col = grn2, type="l", lwd = 3, lty=1)
  
  ###lines(y = smoothed.rstar2  ,  x = dates, col = blu1, type="l", lwd = 3, lty=2)
  ###lines(y = smoothed.g2      ,  x = dates, col = red1, type="l", lwd = 3, lty=2)
  ###lines(y = smoothed.z2      ,  x = dates, col = grn2, type="l", lwd = 3, lty=2)
  
  # date labels every 7 or so quarters
  labDates <- seq(head(dates,1), tail(dates, 1), by = "8 quarters")
  axis.Date(side = 1, dates, at = labDates, format = "%Y", las = 1)
  grid(29,NULL)
  abline(h = 0 , lwd = 1)
  # add legend
  legend("topright", legend = c("rstar","trend growth","other factor"),  
         col=c(blu1, red1, grn2), lty=1, lwd = 3)
  
  dev.off()
  
  }

  print( proc.time() - ptm )
  
# close the loop over the countries
}
  

# end sink -> writing to .txt file results
if (WRITE.OUTPUT.TXT) { sink() }

print("Done")    



  
  
  
    
  
## EOF ---------------------------