#------------------------------------------------------------------------------#
# File:        get.UK.data.R
###############################################################################################################
# Description: This file gets and prepares the data used in HLW for the UK.
# PACKAGES OECD, QUANDL are used. # URL LINK FOR BOE DATA ARE:
# https://www.bankofengland.co.uk/monetary-policy/the-interest-rate-bank-rate for bank rate changes data
# https://www.bankofengland.co.uk/-/media/boe/files/monetary-policy/baserate.xls?la=en&hash=EEB8729ABFFF4B947B85C328340AE5155A99AD0F
# GO TO THIS SITE AND MANNUALY DOWNLOAD HISTORY OF QUARTERLY AVERAGE OF HISTORICAL BANK RATE IUQABEDR:
# https://www.bankofengland.co.uk/boeapps/database/fromshowcolumns.asp?Travel=NIxAZxSUx&FromSeries=1&ToSeries=50&DAT=RNG&FD=1&FM=Jan&FY=1963&TD=31&TM=Dec&TY=2025&FNY=Y&CSVF=TT&html.x=66&html.y=26&SeriesCodes=IUQABEDR&UsingCodes=Y&Filter=N&title=IUQABEDR&VPD=Y#
###############################################################################################################
rm(list = ls())
# CLEAR THE CONSOLE
cat("\014")
# load helper functions
path.2.functions = "./R.local.Functions/"
source(paste0(path.2.functions,"./modified.utilities.R"))

# LOAD REQUIRED LIBRARIES ----------------------------------------------------------------------------------------------
# List of required packages
pkgs <- c('zoo','tis','seasonal','ecb','RCurl','OECD','Quandl')
# get packages that are not installed
get.pkgs <- pkgs[!pkgs %in% installed.packages()]
if(!length(get.pkgs)==0){cat('Getting packages\n'); for(lib in get.pkgs) {install.packages(lib, dependencies = TRUE) } }
sapply( get.pkgs, require, character = TRUE )
# LOAD PACKAGES  
cat('Loading packages\n')
for(ii in pkgs) library(ii, character.only = TRUE)

# SET TO ONE TO SAVE DATA TO CSV FILE
save.data.2.csv     <- 0
update.FRED.data    <- 0
update.ONS.gdp.data <- 0
download.data.vintage = '2020.Oct.5'
#------------------------------------------------------------------------------#
# Get Raw Data
#------------------------------------------------------------------------------#
# SET THE START AND END DATES OF THE DATA USED IN THE ESTIMATION
data.start      <- c(1955,1)
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
save.file.name  <- paste0(store.data.dir, 'UK.data.csv')

#------------------------------------------------------------------------------#
# they seem to have this data also as FRED, so not sure why not do a getFREDnoWget?
if (update.FRED.data) {
  # gdp.fred = gdp ONS x 1000;
  gdp.fred        <- getFREDnoWget('https://fred.stlouisfed.org/data/GBRRGDPQDSNAQ.txt'  , paste0(raw.source.dir, '/UK_GBRRGDPQDSNAQ_FRED.txt'))
  # cpi data from FRED
  cpi.fred        <- getFREDnoWget('https://fred.stlouisfed.org/data/GBRCPIALLQINMEI.txt', paste0(raw.source.dir, '/UK_GBRCPIALLQINMEI_FRED.txt'))
  core.cpi.fred   <- getFREDnoWget('https://fred.stlouisfed.org/data/GBRCPICORQINMEI.txt', paste0(raw.source.dir, '/UK_GBRCPICORQINMEI_FRED.txt'))
  # get quarterly BOE bank rate from fred, BOERUKQ
  boe.rate.fred.q     <- getFREDnoWget('https://fred.stlouisfed.org/data/BOERUKQ.txt', paste0(raw.source.dir, '/UK_BOERUKQ_FRED.txt'), convert2tif = FALSE)
  boe.rate.fred.q[,1] <- as.Date(boe.rate.fred.q[,1], format="%Y-%m-%d") # class(bbr)
} else {
  # LOAD FROM FILES DOWNLOADED TO raw.source.dir
  raw.data.dir <- raw.source.dir
  gdp.fred        <- readFREDrawData(paste0(raw.data.dir, 'UK_GBRRGDPQDSNAQ_FRED.txt'))
  cpi.fred        <- readFREDrawData(paste0(raw.data.dir, 'UK_GBRCPIALLQINMEI_FRED.txt'))
  core.cpi.fred   <- readFREDrawData(paste0(raw.data.dir, 'UK_GBRCPICORQINMEI_FRED.txt'))
  # get quarterly BOE bank rate from fred, BOERUKQ
  boe.rate.fred.q     <- readFREDrawData(paste0(raw.source.dir, '/UK_BOERUKQ_FRED.txt'), convert2tif = FALSE)
  boe.rate.fred.q[,1] <- as.Date(boe.rate.fred.q[,1], format="%Y-%m-%d") # class(bbr)
}

# DOWNLOAD DATA FROM ONS.GOV.UK AND BOE
# NOTE: HLW wrote: We seasonally adjust data rather than using HEXSA, HICPSA because of availability
gdp.store.file      <- paste0(raw.source.dir, 'uk_gdp_ons_abmi.csv')

if (update.ONS.gdp.data) { # change the date at the end of the http address
  download.file(url <- 'https://www.ons.gov.uk/generator?format=csv&uri=/economy/grossdomesticproductgdp/timeseries/abmi/pn2&series=&fromQuarter=Q1&fromYear=1955&toQuarter=Q4&toYear=2059&frequency=quarters', 
                        method='curl', destfile = gdp.store.file)
  print(" Updating GDP data from ONS.gov.uk  --> DONE")
}

#------------------------------------------------------------------------------#
# GDP DATA READ/GET FROM CSV FILE OR ONS.GOV.UK
#------------------------------------------------------------------------------#
if (!file.exists(gdp.store.file)){
  download.file(url <- 'https://www.ons.gov.uk/generator?format=csv&uri=/economy/grossdomesticproductgdp/timeseries/abmi/pn2&series=&fromquarter=q1&fromyear=1955&toquarter=q4&toyear=2059&frequency=quarters',
                        method='curl', destfile = gdp.store.file)
  print(" Downloading GDP data from ONS.gov.uk  --> DONE")
}

gdp.data  <- read.table(file = gdp.store.file, skip = 6, header = FALSE,
                        sep = ',', stringsAsFactors = FALSE)
gdp       <- tis(gdp.data$V2, start = data.start, tif='quarterly')
# Take log of real GDP
gdp.log   <- log(gdp)

#------------------------------------------------------------------------------#
# IMPORT CORE CPI AND CPI DATA
#------------------------------------------------------------------------------#
# PRIOR STEPS:
#     1. Download data from the OECD website:
#        Source: Consumer Prices (MEI); Series: Consumer prices - all items;
#                Consumer prices - all items non-food, non-energy
#     2. Save both data series in one CSV and specify the file name in cpi.file
core.cpi.start <- c(1970,1)
cpi.start      <- c(1955,1)

# THIS USES THE OECD R API.
# https://stats.oecd.org/restsdmx/sdmx.ashx/GetData/PRICES_CPI/GBR.CPGRLE01.IXOB.Q/all?startTime=1949-Q3&endTime=2020-Q1
df.core.cpi <- get_dataset("PRICES_CPI",filter = "GBR.CPGRLE01.IXOB.Q",start_time = 1947, end_time = data.end[1])
df.cpi      <- get_dataset("PRICES_CPI",filter = "GBR.CPALTT01.IXOB.Q",start_time = 1947, end_time = data.end[1])
# make TIS
core.cpi.nsa  <- tis(df.core.cpi$obsValue, start = core.cpi.start, tif = 'quarterly')
cpi.nsa       <- tis(df.cpi$obsValue, start = cpi.start, tif = 'quarterly')
# same as GBRCPIALLQINMEI from fred, but starts in 1960:Q1.

# Seasonally adjust CPI and core CPI data and re-format as tis series
cpi       <- final(seas(as.ts(naWindow(cpi.nsa),freq=4)))
cpi       <- as.tis(cpi,start=cpi.start,tif='quarterly')
core.cpi  <- final(seas(as.ts(naWindow(core.cpi.nsa),freq=4)))
core.cpi  <- as.tis(core.cpi,start=core.cpi.start,tif='quarterly')

# Create annualized core inflation and inflation series using the price indices
core.inflation.q  <- 400*log(core.cpi/Lag(core.cpi, k=1))
all.inflation.q   <- 400*log(cpi/Lag(cpi, k=1))
# MERGE Final inflation series: CPI series is used prior to 1970q2; # thereafter, use core CPI series
inflation  <- mergeSeries(window(all.inflation.q, end = core.cpi.start), window(core.inflation.q, start = shiftQuarter(core.cpi.start, 1)))

# Inflation expectations measure: 4-quarter moving average of past inflation
inflation.expectations <- (inflation + Lag(inflation, k=1) + Lag(inflation, k=2) + Lag(inflation, k=3))/4

# compare.inflation = cbind(all.inflation.q, core.inflation.q, union = TRUE)
# # # check with a plot if needed
# tisPlot( compare.inflation,
#          color = c("red", "blue","green"),
#          lineType = c(1,2,2),labelLeftTicks = TRUE,
#          yearLabels = TRUE,
#          lineWidth = c(6,5,4),tck = .01,
#          xMinorTickFreq = "annual", nberShade = TRUE)

#------------------------------------------------------------------------------#
# GET/READ INTEREST RATE DATA
#------------------------------------------------------------------------------#
boe.interest.rate.history.file <- paste0(raw.source.dir, 'IUQABEDR Bank of England Database.csv')
if (!file.exists(boe.interest.rate.history.file)) {
  print(" You need to manually donwload the CSV file from the website that will pop-up")
  print(" Click on the CSV button below Export source data, and put the csv file in")
  print(raw.source.dir)
  browseURL("https://www.bankofengland.co.uk/boeapps/database/fromshowcolumns.asp?Travel=NIxAZxSUx&FromSeries=1&ToSeries=50&DAT=RNG&FD=1&FM=Jan&FY=1963&TD=31&TM=Dec&TY=2025&FNY=Y&CSVF=TT&html.x=66&html.y=26&SeriesCodes=IUQABEDR&UsingCodes=Y&Filter=N&title=IUQABEDR&VPD=Y#")
}
ir.data.tmp <- read.table(file = boe.interest.rate.history.file,
                          skip = 1, header = FALSE,sep = ',', stringsAsFactors = FALSE)

interest.data.Q.recent <- tis(rev(ir.data.tmp$V2), start = c(1975,1), tif = 'quarterly')
# head(interest.data.Q.recent);tail(interest.data.Q.recent)

# trim to shorter size to be able to convert to tis
interest.data.Q.historical <- boe.rate.fred.q[boe.rate.fred.q[["DATE"]] >= "1900-01-01", ]
# now merge to one long tis object
ir.historical   <- tis(interest.data.Q.historical["VALUE"], start = c(1900,1), tif='quarterly')

interest.splice.date <- c(1975,1)
interest.spliced <- mergeSeries( window(ir.historical, end = shiftQuarter(interest.splice.date,-1) ),
                                 window(interest.data.Q.recent, start  = interest.splice.date))
# cbind(interest.spliced,interest.data.Q.recent,ir.historical,union = TRUE)

# Express interest rate data on a 365-day basis
interest.rate <- 100*((1+interest.spliced/36000)^365 -1)
real.rate     <- interest.rate - inflation.expectations

# ### interest rate ------------------------------------------------------------------
# ir.mean.q = Quandl("BOE/IUQABEDR")
# ir.eofp.q = Quandl("BOE/IUQLBEDR")
# daily is Quandl("BOE/IUDBEDR")


#------------------------------------------------------------------------------#
# Output Data
#------------------------------------------------------------------------------#
# COMBINE THE DATA. THE CBIND FUNCTION ALREADY OPERATES ON A TIS OBJECT, SO JUST ADD UNION TRUE AND IT WILL COMBINE WITH NANS
UK.data <- cbind(gdp.log,
                 real.rate,
                 interest.rate,
                 inflation,
                 inflation.expectations,
                 core.inflation.q,
            union = TRUE)

# NOW TRIM TO DESIRED LENGTH: START.DATE TO END.DATE
UK.data.out  <- window( UK.data,  start = data.start, end = data.end)
# head(UK.data.out); tail(UK.data.out)

###############################################################################################
# WRITE DATA WITH DATES.
# this format below takes 1.Jan.1960 as date for 1960:Q1, as opposed to 31.3.1960. I prefer the latter as also in the US.data.out file
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
date.sequence <- seq(from=(as.Date(ti(shiftQuarter(data.start,-1),  'quarterly'))+1), to =
                          (as.Date(ti(shiftQuarter(data.end,-1),tif='quarterly'))+1), by = 'quarter')

# ADD DATES TO DATA.OUT TIS FRAME
UK.data.withdates       <- data.frame(matrix(NA, dim(UK.data.out)[1], dim(UK.data.out)[2]+1))
UK.data.withdates[,1]   <- date.sequence
UK.data.withdates[,2:(dim(UK.data.out)[2]+1)] <- UK.data.out
# check the series
head(UK.data.withdates$X1)

# # SET COLUMN NAMES FOR CSV OUTPUT
output.col.names <- c("Date",colnames(UK.data.out))
colnames(UK.data.withdates) <- output.col.names
print(head2tail(UK.data.withdates))

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# NOW SAVE TO CSV FILE IF save.data.2.csv == 1
if (save.data.2.csv) {
  if (!dir.exists(store.data.dir)){dir.create(store.data.dir)}
	write.table(UK.data.withdates, save.file.name, col.names = output.col.names,
              quote=FALSE, row.names=FALSE, sep = ',', na = '')
  print( paste0("DONE! all data saved to ", store.data.dir) )
}
print("ALL DONE")


































































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


















# https://stats.oecd.org/restsdmx/sdmx.ashx/GetDataStructure/PRICES_CPI
# United Kingdom	CPALTT01	CPI: 01-12 - All items	IXOB	Index	Q	Quarterly	1955-Q1
# United Kingdom	CPALTT01	CPI: 01-12 - All items	IXOB	Index	Q	Quarterly	1955-Q2
# United Kingdom	CPALTT01	CPI: 01-12 - All items	IXOB	Index	Q	Quarterly	1955-Q3
# FRED HAS
# Organization for Economic Co-operation and Development, Consumer Price Index of All Items in the United Kingdom
# [GBRCPIALLQINMEI], retrieved from FRED, Federal Reserve Bank of St. Louis; https://fred.stlouisfed.org/series/GBRCPIALLQINMEI,
# March 7, 2020.
# Source: Organization for Economic Co-operation and Development
#
# Release: Main Economic Indicators
#
# Units:  Growth Rate Same Period Previous Year, Not Seasonally Adjusted
#
# Frequency:  Quarterly
#
# OECD descriptor ID: CPGRLE01
# OECD unit ID: GY
# OECD country ID: GBR
#  https://stats.oecd.org/restsdmx/sdmx.ashx/GetDataStructure/PRICES_CPI
# cli_filters <- list(
#     "CPGRLE01.IXOB",
#     c("GBR"),
#     "Q"
# )
# cli_raw <- get_dataset("PRICES_CPI", filter= "GBR")

