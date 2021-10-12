#------------------------------------------------------------------------------#
# File:        get.EA.data.R
###############################################################################################################
# Description: This file gets and prepares the data used in HLW for the EuroArea.
# DB: source("modified.utilities.R") HAS MODFIED THEIR CALL TO GETFREDWGET AS IT DOES NOT WORK ON MY Machine
# The ECB package is used to get the individual series,
###############################################################################################################



rm(list = ls())
# CLEAR THE CONSOLE
cat("\014")
# load helper functions
path.2.functions = "./R.local.Functions/"
source(paste0(path.2.functions,"./modified.utilities.R"))

# LOAD REQUIRED LIBRARIES ----------------------------------------------------------------------------------------------
# List of required packages
pkgs <- c('zoo','tis','seasonal','ecb','RCurl')
# get packages that are not installed
get.pkgs <- pkgs[!pkgs %in% installed.packages()]
if(!length(get.pkgs)==0){cat('Getting packages\n'); for(lib in get.pkgs) {install.packages(lib, dependencies = TRUE) } }
sapply( get.pkgs, require, character = TRUE )
# LOAD PACKAGES  
cat('Loading packages\n')
for(ii in pkgs) library(ii, character.only = TRUE)

# SET TO ONE TO SAVE DATA TO CSV FILE
save.data.2.csv     <- 0
get.new.AWMD.data   <- 0
download.data.vintage = '2020.Oct.5'


#------------------------------------------------------------------------------#
# Get Raw Data
#------------------------------------------------------------------------------#
# SET THE START AND END DATES OF THE DATA USED IN THE ESTIMATION
data.start 	    <- c(1971,1)
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
save.file.name  <- paste0(store.data.dir, 'EA.data.csv')

#------------------------------------------------------------------------------#
# Identify the date at which the core CPI series begins (Mnemonic: HEX)
core.cpi.start <- c(1987,4)

# IMPORT AREA WIDE MODEL DATA
# NOTE: HLW wrote: We seasonally adjust data rather than using HEXSA, HICPSA because of availability
awm.file <- paste0(raw.source.dir, 'awm19up18.csv')

# if (!file.exists(awm.file)){
if (get.new.AWMD.data) {
  print(" Downloading data from eabcn.org network --> DONE")
  download.file(url='https://eabcn.org/sites/default/files/awm19up18.csv', method='curl', destfile = awm.file )
}

#------------------------------------------------------------------------------#
# GET DATA FROM EABCN.ORG: must be connected to internet to download
#------------------------------------------------------------------------------#
# https://eabcn.org/sites/default/files/awm19up18.csv
awm.data <- read.table(awm.file, header=TRUE,sep=',')
# head(awm.data)
awm.start    <- c(as.numeric(substr(awm.data[1,'X'],1,4)),
                  as.numeric(substr(awm.data[1,'X'],6,7)))

gdp          <- tis(awm.data$YER,   start = awm.start, tif = 'quarterly')
cpi.nsa      <- tis(awm.data$HICP,  start = awm.start, tif = 'quarterly')
core.cpi.nsa <- tis(awm.data$HEX,   start = awm.start, tif = 'quarterly')
core.cpi.nsa <- window(core.cpi.nsa,start = core.cpi.start)
interest.q   <- tis(awm.data$STN,   start = awm.start, tif = 'quarterly')

# IMPORT ECB SDW SERIES FROM WAREHOUSE
### GDP ---------------------------------------------------------------------------
gdp.api       <- get_data("MNA.Q.Y.I8.W2.S1.S1.B.B1GQ._Z._Z._Z.EUR.LR.N")
gdp.start.ecb <- c(1995,1)
gdp.ecb       <- tis(gdp.api$obsvalue, start = gdp.start.ecb, tif = 'quarterly')
gdp.spliced   <- splice(gdp, gdp.ecb, gdp.start.ecb, "Quarterly")
# Take log of real GDP
gdp.log       <- log(gdp.spliced)

### CPI ---------------------------------------------------------------------------
core.cpi.api          <- get_data("ICP.M.U2.N.XE0000.4.INX")
core.cpi.start.ecb    <- c(1996,1)
core.cpi.splice.date  <- c(2018,1)
core.cpi.nsa.m.ecb    <- tis(core.cpi.api$obsvalue, start = core.cpi.start.ecb, tif = 'monthly')


# Aggregate ECB core CPI data to quarterly from monthly frequency
core.cpi.nsa.ecb      <- convert(core.cpi.nsa.m.ecb, tif = 'quarterly', observed = 'averaged')
# Splice ECB core CPI data with AWM core CPI data in 2015q4
core.cpi.nsa.spliced  <- splice(core.cpi.nsa, core.cpi.nsa.ecb, shiftQuarter(core.cpi.splice.date, -1), "Quarterly")
# Seasonally adjust CPI and core CPI data and re-format as tis series
cpi                   <- final(seas(as.ts(naWindow(cpi.nsa),freq=4)))
cpi                   <- as.tis(cpi,start=awm.start,tif='quarterly')
core.cpi              <- final(seas(as.ts(naWindow(core.cpi.nsa.spliced),freq=4)))
core.cpi              <- as.tis(core.cpi,start=core.cpi.start,tif='quarterly')

# Create annualized core inflation and inflation series using the price indices
core.inflation.q      <- 400*log(core.cpi/Lag(core.cpi, k=1))
inflation.q           <- 400*log(cpi/Lag(cpi, k=1))

# Final inflation series: CPI series is used prior to 1988q1;
# thereafter, use core CPI series
inflation             <- mergeSeries(window(inflation.q,      end   = core.cpi.start) ,
                                     window(core.inflation.q, start = shiftQuarter(core.cpi.start,1)) )

# compare.inflation = cbind(inflation.q, core.inflation.q, inflation, union = TRUE)
# # # check with a plot if needed
# tisPlot( compare.inflation,
#          color = c("red", "blue","green"),
#          lineType = c(1,2,2),labelLeftTicks = TRUE,
#          yearLabels = TRUE,
#          lineWidth = c(6,5,4),tck = .01,
#          xMinorTickFreq = "annual", nberShade = TRUE)

# Inflation expectations measure: 4-quarter moving average of past inflation
inflation.expectations <- (inflation + Lag(inflation, k=1) + Lag(inflation, k=2) + Lag(inflation, k=3))/4

### interest rate ------------------------------------------------------------------
interest.api          <- get_data("FM.Q.U2.EUR.RT.MM.EURIBOR3MD_.HSTA")
interest.start.ecb    <- c(1994,1)
interest.splice.date  <- c(2018,1)
interest.q.ecb        <- tis(interest.api$obsvalue, start = interest.start.ecb, tif = 'quarterly')
interest.spliced      <- mergeSeries( window(interest.q,     end    = shiftQuarter(interest.splice.date,-1) ),
                                      window(interest.q.ecb, start  = interest.splice.date))

# Express interest rate data on a 365-day basis
interest.rate <- 100*((1+interest.spliced/36000)^365 -1)
real.rate     <- interest.rate - inflation.expectations

#------------------------------------------------------------------------------#
# Output Data
#------------------------------------------------------------------------------#
# COMBINE THE DATA. THE CBIND FUNCTION ALREADY OPERATES ON A TIS OBJECT, SO JUST ADD UNION TRUE AND IT WILL COMBINE WITH NANS
EA.data <- cbind(gdp.log,
                 real.rate,
                 interest.rate,
                 inflation,
                 inflation.expectations,
            union = TRUE)
# head(EA.data); tail(EA.data)

# NOW TRIM TO DESIRED LENGTH: START.DATE TO END.DATE
EA.data.out  <- window( EA.data,  start = data.start, end = data.end)
# head(EA.data.out); tail(EA.data.out)

###############################################################################################
# WRITE DATA WITH DATES.
# this format below takes 1.Jan.1960 as date for 1960:Q1, as opposed to 31.3.1960. I prefer the latter as also in the US.data.out file
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
date.sequence <- seq(from=(as.Date(ti(shiftQuarter(data.start,-1),  'quarterly'))+1), to =
                          (as.Date(ti(shiftQuarter(data.end,-1),tif='quarterly'))+1), by = 'quarter')

# ADD DATES TO DATA.OUT TIS FRAME
EA.data.withdates       <- data.frame(matrix(NA, dim(EA.data.out)[1], dim(EA.data.out)[2]+1))
EA.data.withdates[,1]   <- date.sequence
EA.data.withdates[,2:(dim(EA.data.out)[2]+1)] <- EA.data.out
# check the series
# head(EA.data.withdates$X1)

# # SET COLUMN NAMES FOR CSV OUTPUT
output.col.names <- c("Date",colnames(EA.data.out))
colnames(EA.data.withdates) <- output.col.names
print(head2tail(EA.data.withdates))
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# NOW SAVE TO CSV FILE IF save.data.2.csv == 1
if (save.data.2.csv) {
  # check if store.data.dir exists, if not make it
  if (!dir.exists(store.data.dir)){dir.create(store.data.dir)}
	write.table(EA.data.withdates, save.file.name, col.names = output.col.names,
              quote=FALSE, row.names=FALSE, sep = ',', na = '')
  print( paste0("DONE! all data saved to ", store.data.dir) )
}

cat("ALL DONE")




#
# #------------------------------------------------------------------------------#
# # Output Data
# #------------------------------------------------------------------------------#
# # add dates
# date.sequence <- seq(from=(as.Date(ti(shiftQuarter(data.start,-1),'quarterly'))+1), to =
#                        (as.Date(ti(shiftQuarter(data.end,-1),tif='quarterly'))+1), by = 'quarter')
#
# data.out <- window(cbind(gdp.log, inflation, inflation.expectations, interest),start = data.start, end = data.end)
#
# data.withdates       <- data.frame(matrix(NA,dim(data.out)[1],dim(data.out)[2]+1))
# data.withdates[,1]   <- date.sequence
# data.withdates[,2:(dim(data.out)[2]+1)] <- data.out
#
# # # SET COLUMN NAMES FOR CSV OUTPUT
# output.col.names <- c("Date",colnames(data.out))
# colnames(data.withdates) <- output.col.names;
#
# # now write to csv file
# write.table(data.withdates, file = 'rstar.data.ea.csv', sep = ',',
#             col.names = TRUE, quote = FALSE, na = '.', row.names = FALSE)
#
# print(data.withdates)







# download.file('https://www.ons.gov.uk/generator?format=csv&uri=/economy/grossdomesticproductgdp/timeseries/abmi/pn2&series=&fromQuarter=Q1&fromYear=1955&toQuarter=Q4&toYear=2019&frequency=quarters', method = 'curl', destfile = 'here.csv')









