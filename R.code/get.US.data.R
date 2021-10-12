#------------------------------------------------------------------------------#
# File:        get.US.data.R
###############################################################################################################
# Description: This file gets and prepares the data used in HLW for the US.
# DB: source("modified.utilities.R") HAS MODFIED THEIR CALL TO GETFREDWGET AS IT DOES NOT WORK ON MY Machine
###############################################################################################################




rm(list = ls())
# CLEAR THE CONSOLE
cat("\014")
# load helper functions
path.2.functions = "./R.local.Functions/"
source(paste0(path.2.functions,"./modified.utilities.R"))

# LOAD REQUIRED LIBRARIES ----------------------------------------------------------------------------------------------
# List of required packages
pkgs <- c('zoo','tis')
# get packages that are not installed
get.pkgs <- pkgs[!pkgs %in% installed.packages()]
if(!length(get.pkgs)==0){cat('Getting packages\n'); for(lib in get.pkgs) {install.packages(lib, dependencies = TRUE) } }
sapply( get.pkgs, require, character = TRUE )
# LOAD PACKAGES  
cat('Loading packages\n')
for(ii in pkgs) library(ii, character.only = TRUE)

# SET TO ONE TO SAVE DATA TO CSV FILE
save.data.2.csv     <- 0
get.new.FRED.data   <- 0
download.data.vintage = '2020.Oct.5'


#------------------------------------------------------------------------------#
# Get Raw Data
#------------------------------------------------------------------------------#
# SET THE START AND END DATES OF THE DATA USED IN THE ESTIMATION
data.start      <- c(1947,1)
data.end        <- c(2020,2)
# Covert to strings
data.beg.string <- paste0( toString(data.start[1]),"Q",toString(data.start[2])  )
data.end.string <- paste0( toString(data.end[1])  ,"Q",toString(data.end[2])    )

# OUTPUT FILE DIRs and NAMEs
# main data store dir
if (!dir.exists('../data/')) {dir.create('../data/')}
# raw source dir in main data store dir
raw.source.dir  <- paste0('../data/R.source.data.',download.data.vintage,'/')
if (!dir.exists(raw.source.dir)) {dir.create(raw.source.dir)}
# data store dir for Matlab and R
store.data.dir  <- paste0('../data/R.data.for.estimation.',download.data.vintage,'/')
save.file.name  <- paste0(store.data.dir, 'US.data.csv')

#------------------------------------------------------------------------------#
#{ Import data using the function getFRED() in utilities.R
## ***************************************************************************************************************************************
## NOTE: THIS SCRIPT BELOW CREATES AN AVERAGED MONTHLY SERIES FOR THE QUARTERLY SERIES,
## NOT END OF MONTH VALUES. ALSO THE WAY THIS IS DONE IS SOMEWHAT DIFFERENT THAN NORMAL
## MONTHLY AVERAGE, WHICH IS CONTROLLED BY THE method = 'discrete' OR method = 'constant' OPTIONS.
## ***************************************************************************************************************************************
if (get.new.FRED.data) {
  if (!dir.exists(raw.source.dir)){dir.create(raw.source.dir)}
  gdp             <- getFREDnoWget('https://fred.stlouisfed.org/data/GDPC1.txt'       	  , paste0(raw.source.dir, '/GDPC1.txt'))
  pce.index       <- getFREDnoWget('https://fred.stlouisfed.org/data/PCEPILFE.txt'		    ,	paste0(raw.source.dir, '/PCEPILFE.txt'))
  ## Careful with UPDATING ny.discount as the new data seems to be of much shorter lengths. so just read in the one downloaded in 2017
  ny.discount     <- getFREDnoWget('https://fred.stlouisfed.org/data/INTDSRUSM193N.txt'   ,	paste0(raw.source.dir, '/INTDSRUSM193N.txt'))
  ## this is the 'new' series that I use which goes back to Nov 1914 to Jul 1969 (2012-08-20)
  new.ny.discount <- getFREDnoWget('https://fred.stlouisfed.org/data/M13009USM156NNBR.txt',	paste0(raw.source.dir, '/M13009USM156NNBR.txt'))
  fed.funds       <- getFREDnoWget('https://fred.stlouisfed.org/data/FEDFUNDS.txt'		    , paste0(raw.source.dir, '/FEDFUNDS.txt'))
  # (TOTAL) Personal Consumption Expenditures: Chain-type Price Index (PCECTPI)	Quarterly series INCLUDES VOLATLIE ITMES
  # this series is used in LW03 to splice data back to 1948.
  pce.lw03        <- getFREDnoWget('https://fred.stlouisfed.org/data/PCECTPI.txt'         ,	paste0(raw.source.dir, '/PCECTPI.txt'))
} else {
  # LOAD FROM FILES DOWNLOADED TO raw.source.dir
  raw.data.dir    <- raw.source.dir
  gdp             <- readFREDrawData(paste0(raw.data.dir, 'GDPC1.txt'))
  pce.index       <- readFREDrawData(paste0(raw.data.dir, 'PCEPILFE.txt'))
  ny.discount     <- readFREDrawData(paste0(raw.data.dir, 'INTDSRUSM193N.txt'))
  new.ny.discount <- readFREDrawData(paste0(raw.data.dir, 'M13009USM156NNBR.txt'))
  fed.funds       <- readFREDrawData(paste0(raw.data.dir, 'FEDFUNDS.txt'))
  # (TOTAL) Personal Consumption Expenditures: Chain-type Price Index (PCECTPI)	Quarterly series INCLUDES VOLATLIE ITMES # this series is used in LW03 to splice data back to 1948.
  pce.lw03        <- readFREDrawData(paste0(raw.data.dir, 'PCECTPI.txt'))
}

#------------------------------------------------------------------------------#
## Prepare Data ----
#------------------------------------------------------------------------------#
# Take log of real GDP
gdp.log         <- log(gdp)

# Create an annualized inflation series 
inflation.hlw   <- 400*log(pce.index/Lag(pce.index, k=1)) # HLW CORE PCE ONLY
inflation.pce   <- 400*log(pce.lw03/Lag(pce.lw03, k=1))   # LW03 use also PCE to go back in time
# merge into one series. 
inflation       <- mergeSeries(window(inflation.pce, end    = c(1959,1)) ,
                               window(inflation.hlw, start  = c(1959,2)) )
# HLW: Inflation expectations measure: 4-quarter moving average of past inflation
inflation.expectations  <- (    inflation       +
                            Lag(inflation, k=1) +
                            Lag(inflation, k=2) +
                            Lag(inflation, k=3))/4

## EXPRESS INTEREST RATE DATA ON A 365-DAY BASIS FROM LW -----
ny.discount.eff     <- 100*((1+ny.discount/36000)^365 -1) # this is not really needed
new.ny.discount.eff <- 100*((1+new.ny.discount/36000)^365 -1)
fed.funds.eff       <- 100*((1+fed.funds/36000)^365 -1)
# HLW use NYFed discount rate is used prior to 1965; thereafter, use the effective federal funds rate
interest.rate <- mergeSeries(window(new.ny.discount.eff,  end   = c(1964,4)) ,
                             window(fed.funds.eff,        start = c(1965,1)) )
# NOW MAKE THE REAL RATE SERIES
real.rate     <- interest.rate - inflation.expectations

#------------------------------------------------------------------------------#
# Output Data
#------------------------------------------------------------------------------#
# COMBINE THE DATA. THE CBIND FUNCTION ALREADY OPERATES ON A TIS OBJECT, SO JUST ADD UNION TRUE AND IT WILL COMBINE WITH NANS
US.data <- cbind(gdp.log,
                 real.rate,
                 interest.rate,
                 inflation,
                 inflation.expectations,
            union = TRUE)
head(US.data);tail(US.data)

# NOW TRIM TO DESIRED LENGTH: START.DATE TO END.DATE
US.data.out  <- window(US.data, start = data.start, end = data.end)
head(US.data.out); tail(US.data.out)

## PLOTTING -----
# tisPlot( US.data.out,
#          color = c("red", "blue", "green", "gray"),
#          lineType = c(1,2,4,4),labelLeftTicks = TRUE,
#          xTickFreq = "NonRe",
#          yearLabels = TRUE,
#          lineWidth = c(6,4,4,4),tck = .01,
#          xMinorTickFreq = "annual", nberShade = TRUE)
# abline(h=0,lwd = 2)
# # tisLegend(legend = colnames(series.to.plot), boxType ="y", xrel = 0, yrel = 0)
# legend( x="topleft",
#         legend = colnames(US.data.out),
#         col = c("red", "blue", "green", "gray"),

#         lty = c(1,2,4), lwd = 3)
# box(which = "plot", lty = "solid", col = 'black', lwd = 2)

###############################################################################################
# WRITE DATA WITH DATES. 
# this format below takes 1.Jan.1960 as date for 1960:Q1, as opposed to 31.3.1960. I prefer the latter as also in the US.data.out file
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
date.sequence <- seq(from=(as.Date(ti(shiftQuarter(data.start,-1),  'quarterly'))+1), to =
                          (as.Date(ti(shiftQuarter(data.end,-1),tif='quarterly'))+1), by = 'quarter')

# ADD DATES TO DATA.OUT TIS FRAME
US.data.withdates       <- data.frame(matrix(NA, dim(US.data.out)[1], dim(US.data.out)[2]+1))
US.data.withdates[,1]   <- date.sequence
US.data.withdates[,2:(dim(US.data.out)[2]+1)] <- US.data.out
# check the series
head(US.data.withdates$X1)
 
# # SET COLUMN NAMES FOR CSV OUTPUT
output.col.names <- c("Date",colnames(US.data.out))
colnames(US.data.withdates) <- output.col.names

print(head2tail(US.data.withdates))
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# NOW SAVE TO CSV FILE IF save.data.2.csv == 1
if (save.data.2.csv) {
  # check if store.data.dir exists, if not make it
  if (!dir.exists(store.data.dir)){dir.create(store.data.dir)}
	write.table(US.data.withdates, save.file.name, col.names = output.col.names,
	            quote=FALSE, row.names=FALSE, sep = ',', na = '')
  print( paste0("DONE! all data saved to ", store.data.dir) )
}
print("ALL DONE")

















# EOF