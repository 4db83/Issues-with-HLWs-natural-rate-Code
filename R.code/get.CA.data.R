#------------------------------------------------------------------------------#
# File:        get.CA.data.R
###############################################################################################################
# Description: This file compiles and prepares the data used in HLW for Canada.
# NOTE: They seem to download a GDP index from IMFs-IFS, rather than use a chain based measure. As of March 2020,
# the data for this series is not available for Q4-2019. Data is available also from OECD and FRED.
# SEE: https://data.oecd.org/api/sdmx-json-documentation/ for API documentation
###############################################################################################################


rm(list = ls())
# CLEAR THE CONSOLE
cat("\014")
# load helper functions
path.2.functions = "./R.local.Functions/"
source(paste0(path.2.functions,"./modified.utilities.R"))

# LOAD REQUIRED LIBRARIES ----------------------------------------------------------------------------------------------
# List of required packages
pkgs <- c('tis','seasonal','IMFData','cansim','xts','OECD','dplyr')
# get packages that are not installed
get.pkgs <- pkgs[!pkgs %in% installed.packages()]
if(!length(get.pkgs)==0){cat('Getting packages\n'); for(lib in get.pkgs) {install.packages(lib, dependencies = TRUE) } }
sapply( get.pkgs, require, character = TRUE )
# LOAD PACKAGES  
cat('Loading packages\n')
for(ii in pkgs) library(ii, character.only = TRUE)

# SET TO ONE TO SHOW PLOTS
plot_series_to_check = 0
# SET TO ONE TO SAVE DATA TO CSV FILE
save.data.2.csv     <- 0
get.new.data        <- 0
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
save.file.name  <- paste0(store.data.dir, 'CA.data.csv')

#------------------------------------------------------------------------------#
# ifs.gdp.file.name <-  'IMF_IFS_GDP_data_CA.RData'
store.QNA_ARCHIVE.OECD  <- paste0(raw.source.dir, 'QNA_ARCHIVE.OECD.df.RData')
store.QNA.OECD          <- paste0(raw.source.dir, 'QNA.OECD.df.RData')
store.IMF.IFS           <- paste0(raw.source.dir, 'IMF_IFS_GDP_data_CA.RData')

# MAKE SIMPLE LOG-DIFF FUNCTION
lndf = function(x) diff(log(x))

# ----------------------------------------------------------------------------------------------------------------------
# GET VARIOUS GDP DATA SERIES, SOME FROM FRED, SOME FROM IMF'S IFS SOME FROM OECD
# ----------------------------------------------------------------------------------------------------------------------
# NAEXKP01CAQ189S:  Gross Domestic Product by Expenditure in Constant Prices: Total Gross Domestic Product for Canada (NAEXKP01CAQ189S) National Currency, Seasonally Adjusted
# NAEXKP02CAQ661S:  Gross Domestic Product by Expenditure in Constant Prices: Private Final Consumption Expenditure for Canada (NAEXKP02CAQ661S) Index 2015=100, Seasonally Adjusted
# CANRGDPQDSNAQ:    Gross Domestic Product in Canada (DISCONTINUED) (CANRGDPQDSNAQ)	Millions of Chained 2002 Canadian Dollars, Seasonally Adjusted
# ----------------------------------------------------------------------------------------------------------------------
# GET GDP DATA FROM IMF IFS R-API CANADA NGDP_R_K_SA_IX: THIS IS WHAT THEY SAY THEY USE, BUT IT CANNOT BE THE CASE, ONLY UDATED TO 2019:Q2
# ----------------------------------------------------------------------------------------------------------------------
# 1. Download data from the International Financial Statistics (IFS) Database on the IMF website:
#    Series: "Gross Domestic Product, Real, Seasonally adjusted, Index"
# ----------------------------------------------------------------------------------------------------------------------

if ((!get.new.data) && (file.exists(store.IMF.IFS)) && (file.exists(store.QNA.OECD)) && (file.exists(store.QNA_ARCHIVE.OECD))) {
  # load from files downloaded to raw.source.dir
  raw.data.dir  <- raw.source.dir
  # IMF.IFS
  load(file = store.IMF.IFS )

  # OECD
  load(file = store.QNA.OECD)
  load(file = store.QNA_ARCHIVE.OECD)

  # FRED
  # Gross Domestic Product by Expenditure in Constant Prices: Total Gross Domestic Product for Canada (NAEXKP01CAQ189S) National Currency, Seasonally Adjusted
  gdp.fred.exp.tgdp <- readFREDrawData( paste0(raw.source.dir,  'CA_NAEXKP01CAQ189S_FRED.txt'))
  # Gross Domestic Product by Expenditure in Constant Prices: Private Final Consumption Expenditure for Canada (NAEXKP02CAQ661S) Index 2015=100, Seasonally Adjusted
  gdp.fred.exp.pfc  <- readFREDrawData( paste0(raw.source.dir,  'CA_NAEXKP02CAQ661S_FRED.txt'))
  # Gross Domestic Product in Canada (DISCONTINUED) (CANRGDPQDSNAQ)	Millions of Chained 2002 Canadian Dollars, Seasonally Adjusted
  gdp.fred.chained  <- readFREDrawData( paste0(raw.source.dir,  'CA_CANRGDPQDSNAQ_FRED.txt'))

} else {
  print(" ... Downloading GDP data from IMF's IFS, OECD and FRED")

  # IMF.IFS
  queryfilter <- list(CL_FREA = "Q", CL_AREA_IFS = "CA", CL_INDICATOR_IFS = "NGDP_R_K_SA_IX")
  Q.CA.NGDP.query <- CompactDataMethod("IFS", queryfilter, "1947-01-01", "2025-12-31", FALSE)  # NY.GDP.MKTP.KD.ZG
  print( paste0(" Done with IMF! saving data to:  ", store.IMF.IFS , ". -> Getting OECD data"))
  save(Q.CA.NGDP.query, file = store.IMF.IFS )

  # OECD
  measure <- c("CUR" ,"CQRSA" ,"DNBSA" ,"GYSA" ,"LNBQRSA" ,"PERSA" ,"CARSA" ,"DOBSA" ,"GPSA" ,"LNBARSA" ,"VOL" ,"PER" ,"CQR" ,"IND" ,"HRSSA" ,"VIXNB" ,"VNBQRSA" ,"CAR" ,"GRW" ,"HRS" ,"VIXNBSA" ,"VNBARSA" ,"CTQRGPSA" ,"CPCARSA" ,"POP" ,"JOBSA" ,"VIXOBSA" ,"VOBARSA" ,"JOB" ,"LNBQR" ,"VNBQR" ,"HCPCARSA" ,"VNBAR" ,"CD" ,"VPVOBARSA" ,"HVPVOBARSA")
  subject <- c("GDP", "B1_GE", "B1_GS1", "B1_GA", "B1_GAI3","B1_GI")
  filter_list <-  list( c("CAN"), subject, measure, c("Q"))
  QNA.OECD.df           <- get_dataset("QNA",filter = filter_list, start_time = 1947, end_time = data.end[1])
  QNA_ARCHIVE.OECD.df   <- get_dataset("QNA_ARCHIVE", filter = filter_list, start_time = 1947, end_time = data.end[1])
  save(QNA.OECD.df,         file = store.QNA.OECD)
  save(QNA_ARCHIVE.OECD.df, file = store.QNA_ARCHIVE.OECD)
  print( paste0(" Done with OECD! saving data to:  ", store.QNA.OECD, store.QNA_ARCHIVE.OECD , ". -> Now Getting FRED data"))

  # FRED
  # Gross Domestic Product by Expenditure in Constant Prices: Total Gross Domestic Product for Canada (NAEXKP01CAQ189S) National Currency, Seasonally Adjusted
  gdp.fred.exp.tgdp <- getFREDnoWget('https://fred.stlouisfed.org/data/NAEXKP01CAQ189S.txt', paste0(raw.source.dir,  'CA_NAEXKP01CAQ189S_FRED.txt'))
  # Gross Domestic Product by Expenditure in Constant Prices: Private Final Consumption Expenditure for Canada (NAEXKP02CAQ661S) Index 2015=100, Seasonally Adjusted
  gdp.fred.exp.pfc  <- getFREDnoWget('https://fred.stlouisfed.org/data/NAEXKP02CAQ661S.txt', paste0(raw.source.dir,  'CA_NAEXKP02CAQ661S_FRED.txt'))
  # Gross Domestic Product in Canada (DISCONTINUED) (CANRGDPQDSNAQ)	Millions of Chained 2002 Canadian Dollars, Seasonally Adjusted
  gdp.fred.chained  <- getFREDnoWget('https://fred.stlouisfed.org/data/CANRGDPQDSNAQ.txt',   paste0(raw.source.dir,    'CA_CANRGDPQDSNAQ_FRED.txt'))
  print( paste0(" Done with FRED! Now reading in the series from IMF and OECD"))
}

# IMF.IFS CONVERT SERIES TO TIS OBJECT
gdp.ifs <- tis(data.matrix(Q.CA.NGDP.query$Obs[[1]][2]), start = data.start, tif = 'quarterly')
colnames(gdp.ifs) ="IMF.volindex"

# gdp from OECD data expenditiure approach, Volume index, SA.
gdp.oecd.ge.volindex <- tis((QNA.OECD.df %>% filter(SUBJECT == "B1_GE" & MEASURE == "VIXOBSA"))$obsValue, start = c(1961, 1), tif = 'quarterly')
# gdp from OECD data expenditiure approach, Volume index, SA. VNBQR
# LNBQRSA	Nationalcurrency,chainedvolumeestimates,nationalreferenceyear,quarterlylevels,seasonallyadjusted
gdp.oecd.ge.chain    <- tis((QNA.OECD.df %>% filter(SUBJECT == "B1_GE" & MEASURE == "LNBQRSA"))$obsValue, start = c(1961, 1), tif = 'quarterly')

# the 4 series below give nearly identical growth rates
gdp.comp = cbind( gdp.ifs,
                  gdp.oecd.ge.volindex,
                  gdp.oecd.ge.chain,
                  gdp.fred.exp.tgdp,
                  union = TRUE)
# write.table(lndf(gdp.comp), 'can.gdp.growth.csv', col.names = TRUE, quote=FALSE, row.names=FALSE, sep = ',', na = '')
# print((gdp.comp))
# head2tail(lndf(gdp.comp))

# series.2.plot = window(lndf(gdp.comp), start = c(1948,1) , end = c(2019,4))
series.2.plot = window((gdp.comp), start = c(1948,1) , end = data.end)

# check with a plot if needed
if (plot_series_to_check) {
tisPlot( series.2.plot,
         color = c("red", "blue","green","black","slateblue4","burlywood"),
         lineType = c(1,2,3,2,2,2),labelLeftTicks = TRUE,
         yearLabels = TRUE,
         lineWidth = c(4,6,4,4,6,4),tck = .01,
         xMinorTickFreq = "annual", nberShade = TRUE)

tisPlot( window(lndf(gdp.comp), start = c(1948,1) , end = c(1978,4)),
         color = c("red", "blue","green","black","slateblue4","burlywood"),
         lineType = c(1,2,3,2,2,2),labelLeftTicks = TRUE,
         yearLabels = TRUE,
         lineWidth = c(4,6,4,4,6,4),tck = .01,
         xMinorTickFreq = "annual", nberShade = TRUE)
}

# BEAUSE OF THE MASSIVE KINK IN THE IMFIFS SERIES, MAKE A NEW SERIES FROM SPLICED GROWTH RATES OF IMF SERIES AND OECD VOLINDEX
lndf.gdp.merged <- mergeSeries(window(lndf(gdp.ifs),          end   = c(1961,1)),
                          window(lndf(gdp.oecd.ge.volindex),  start = c(1961,2)))

# now make log.gdp series by taking log of initial condition and then cumsum lndf.gdg.merged
gdp.log <- 1*log(gdp.ifs[1]) + cumsum(lndf.gdp.merged)

series.2.plot <- cbind(exp(gdp.log), gdp.ifs, gdp.oecd.ge.volindex, union = TRUE)

if (plot_series_to_check) {
# check with a plot if needed
tisPlot( lndf(series.2.plot),
         color = c("red", "blue","green","black","slateblue4","burlywood"),
         lineType = c(1,2,3,2,2,2),labelLeftTicks = TRUE,
         yearLabels = TRUE,
         lineWidth = c(4,6,4,4,6,4),tck = .01,
         xMinorTickFreq = "annual", nberShade = TRUE)
}

### ----------------------------------------------------------------------------------------------------------------------
### there are no other useful series
### ----------------------------------------------------------------------------------------------------------------------
### test/search the series of interest here. use info.out = get_data_structure("QNA") to get Dataframe info on QNA.
### gdp.OECD.tmp = QNA.OECD.df %>% filter(SUBJECT == "B1_GE" & MEASURE == "LNBQRSA")
##gdp.OECD.tmp <- QNA.OECD.df %>% filter(SUBJECT == "B1_GI")
### head2tail(gdp.OECD.tmp)
##unique(gdp.OECD.tmp$MEASURE)
### ----------------------------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------------#
# Import monthly core CPI and CPI data (two series - the second CPI series will be used to extend data back)
#---------------------------------------------------------------------------------------------------------#
core.cpi.start <- c(1984,1)
cpi.start      <- c(1992,1)
cpi.back.start <- c(1959,1)
# GET DATA FROM CANSIM (STATISTICS CANADA)
core.cpi.data <- get_cansim_vector("v112593705","1984-01-01")
cpi.data      <- get_cansim_vector("v41690914", "1992-01-01")
cpi.back.data <- get_cansim_vector("v41690973", paste0(cpi.back.start[1], "-01-01"))
cpi.long.data <- get_cansim_vector("v41690973", "1947-01-01")
# FORMAT DATA AS TIME SERIES WITH MONTHLY FREQUENCY
#---------------------------------------------------------------------------------------------------------#
core.cpi.m     <- tis(core.cpi.data$VALUE,  start = core.cpi.start, tif = 'monthly')
cpi.m          <- tis(cpi.data$VALUE,       start = cpi.start,      tif = 'monthly')
cpi.back.m.nsa <- tis(cpi.back.data$VALUE,  start = cpi.back.start, tif = 'monthly')
# seasonally adjust CPI back series
cpi.back.m <- final(seas(as.ts(naWindow(cpi.back.m.nsa),freq=12)))
cpi.back.m <- as.tis(cpi.back.m,start=cpi.back.start,tif='monthly')
# get the very long series
cpi.long.start = c(1947,1)
cpi.long.m.nsa <- tis(cpi.long.data$VALUE,  start = cpi.long.start , tif = 'monthly')
# seasonally adjsut the long cpi series
cpi.long.m <- final(seas(as.ts(naWindow(cpi.long.m.nsa),freq=12), x11 = ""))
cpi.long.m <- as.tis(cpi.long.m, start = cpi.long.start , tif='monthly')
# head2tail(cpi.long.m)

# CPI series: cpi.m is used from Jan 1992 to present and extended back to
# Jan 1959 using growth rates with splice() function in utilities.R
cpi.series.m <- splice(cpi.back.m, cpi.m, cpi.start, freq = 'monthly')

# Aggregate core CPI and CPI series to quarterly frequency by taking the average
core.cpi.q <- aggregate(core.cpi.m,   nf = 4, FUN = mean)
cpi.q      <- aggregate(cpi.series.m, nf = 4, FUN = mean)
cpi.long.q <- aggregate(cpi.long.m,   nf = 4, FUN = mean)

# Create annualized core inflation and inflation series using the price indices
core.inflation.q  <- 400*log(core.cpi.q/Lag(core.cpi.q, k=1))
inflation.q       <- 400*log(cpi.q/Lag(cpi.q, k=1))
inflation.long    <- 400*log(cpi.long.q/Lag(cpi.long.q, k=1))

# Final inflation series: CPI series is used pUS.data.out4q2;
# thereafter, use core CPI series
inflation <- mergeSeries(window(inflation.q, end = core.cpi.start),window(core.inflation.q, start = shiftQuarter(core.cpi.start,1)))

# Inflation expectations measure: 4-quarter moving average of past inflation
inflation.expectations <- (inflation + Lag(inflation, k=1) + Lag(inflation, k=2) + Lag(inflation, k=3))/4
# PLOTTING -----
# # data.2.plot = cbind(cpi.series.m, core.cpi.m, cpi.m, cpi.long, union = TRUE)
# data.2.plot = cbind(inflation.long, inflation.q, union = TRUE)
# # check with a plot if needed
# tisPlot( data.2.plot,
#          color = c("red", "blue","green","black"),
#          lineType = c(1,2,2,2),labelLeftTicks = TRUE,
#          yearLabels = TRUE,
#          lineWidth = c(4,4,4,4),tck = .01,
#          xMinorTickFreq = "annual", nberShade = TRUE)

#---------------------------------------------------------------------------------------------------------#
# Import daily bank rate data
#---------------------------------------------------------------------------------------------------------#
# 1. Download data from the Statistics Canada (CANSIM) website: Series: v122530
# 2. Download data from the Statistics Canada (CANSIM) website: Series: v39079
#---------------------------------------------------------------------------------------------------------#
# GET DATA FROM CANSIM (STATISTICS CANADA)
bank.rate.start   <- "1927-01-01"
bank.rate.data    <- get_cansim_vector("v122530", bank.rate.start)
# format data as time series with monthly frequency
bank.rate.m       <- tis(bank.rate.data$VALUE, start = bank.rate.data$REF_DATE[1], tif='monthly')

# IMPORT DAILY TARGET RATE DATA AND TAKE EOP MONTHLY VALUES
target.rate.start <- c(2001,1)
cat("NOTE: Canadian data only available AFTER Canadian market has OPENED, ie. AFTER 15:30CET\n")
cat("THUS, best to run script in the afternoon European time!\n")
target.rate.data  <- get_cansim_vector("v39079", "1947-01-01")
# format data as time series with monthly frequency
target.rate.d <- as.xts(target.rate.data$VALUE,
                 order.by = as.Date(target.rate.data$REF_DATE) , frequency = 365)
# remove data with zero values (these indicate non-business days)
target.rate.d <- target.rate.d[!(target.rate.d==0.00)]
# Aggregate data to monthly frequency by taking end-of-period values
target.rate.m <- apply.monthly(target.rate.d, tail, n=1)
# Format data as a tis time series with monthly frequency
target.rate.m <- as.tis(data.frame(target.rate.m)$target.rate.m, start = c(1992,12), tif = 'monthly')
# # Bank rate is used prior to May 2001; thereafter, use the target rate
interest.m    <- mergeSeries(window(bank.rate.m, end = c(2001,4)),window(target.rate.m, start = c(2001,5)))
# Aggregate interest rate data to quarterly frequency by taking the average
interest.q    <- aggregate(interest.m, nf = 4, FUN = mean)
# Express interest rate data on a 365-day basis
interest.rate <- 100*((1+interest.q/36000)^365 -1)
# make the real rate
real.rate     <- interest.rate - inflation.expectations

# PLOTTING -----
# # data.2.plot = cbind(cpi.series.m, core.cpi.m, cpi.m, cpi.long, union = TRUE)
# data.2.plot = cbind(target.rate.m, bank.rate.m, interest.m, union = TRUE)
# # check with a plot if needed
# tisPlot( data.2.plot,
#          color = c("red", "blue","green","black"),
#          lineType = c(1,2,3,2),labelLeftTicks = TRUE,
#          yearLabels = TRUE,
#          lineWidth = c(4,4,4,4),tck = .01,
#          xMinorTickFreq = "annual", nberShade = TRUE)

#------------------------------------------------------------------------------#
# Output Data
#------------------------------------------------------------------------------#
# COMBINE THE DATA. THE CBIND FUNCTION ALREADY OPERATES ON A TIS OBJECT, SO JUST ADD UNION TRUE AND IT WILL COMBINE WITH NANS
CA.data <- cbind(gdp.log,
                 real.rate,
                 interest.rate,
                 inflation,
                 inflation.expectations,
                 inflation.long,
            union = TRUE)
# colnames(CA.data)[1] = "gdp.log"
# head2tail(CA.data)

# NOW TRIM TO DESIRED LENGTH: START.DATE TO END.DATE
CA.data.out  <- window( CA.data,  start = data.start, end = data.end)
# head2tail(CA.data.out)

###############################################################################################
# WRITE DATA WITH DATES.
# this format below takes 1.Jan.1960 as date for 1960:Q1, as opposed to 31.3.1960. I prefer the latter as also in the US.data.out file
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
date.sequence <- seq(from=(as.Date(ti(shiftQuarter(data.start,-1),  'quarterly'))+1), to =
                          (as.Date(ti(shiftQuarter(data.end,-1),tif='quarterly'))+1), by = 'quarter')

# ADD DATES TO DATA.OUT TIS FRAME
CA.data.withdates       <- data.frame(matrix(NA, dim(CA.data.out)[1], dim(CA.data.out)[2]+1))
CA.data.withdates[,1]   <- date.sequence
CA.data.withdates[,2:(dim(CA.data.out)[2]+1)] <- CA.data.out
# check the series
# head(CA.data.withdates$X1)

# # SET COLUMN NAMES FOR CSV OUTPUT
output.col.names <- c("Date",colnames(CA.data.out))
colnames(CA.data.withdates) <- output.col.names
print(head2tail(CA.data.withdates))
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# NOW SAVE TO CSV FILE IF save.data.2.csv == 1
if (save.data.2.csv) {
  # check if store.data.dir exists, if not make it
  if (!dir.exists(store.data.dir)){dir.create(store.data.dir)}
	write.table(CA.data.withdates, save.file.name, col.names = output.col.names,
              quote=FALSE, row.names=FALSE, sep = ',', na = '')
  print( paste0("DONE! all data saved to ", store.data.dir) )
}
print("ALL DONE")















# 2            GDP                                     Gross domestic product
# 3          B1_GE              Gross domestic product - expenditure approach
# 4         B1_GS1                                     Gross domestic product
# 52         B1_GI                   Gross domestic product - income approach
# 88         B1_GA  Gross domestic product at market prices - output approach
# 104      B1_GAI3  Gross domestic product at market prices - output approach
# $MEASURE
#        CUR 	1  CurrentpricesCurrentprices
#       CQRSA	2  Nationalcurrency,currentprices,quarterlylevels,seasonallyadjusted
#       DNBSA	3  Deflator,nationalbase/referenceyear,seasonallyadjusted
#        GYSA	4  Growthratecomparedtothesamequarterofpreviousyear,seasonallyadjusted
#     LNBQRSA	5  Nationalcurrency,chainedvolumeestimates,nationalreferenceyear,quarterlylevels,seasonallyadjusted
#       PERSA	6  Persons,seasonallyadjusted
#       CARSA	7  Nationalcurrency,currentprices,annuallevels,seasonallyadjusted
#       DOBSA	8  Deflator,OECDreferenceyear,seasonallyadjusted
#        GPSA	9  Growthratecomparedtopreviousquarter,seasonallyadjusted
#     LNBARSA	10 Nationalcurrency,chainedvolumeestimates,nationalreferenceyear,annuallevels,seasonallyadjusted
#         VOL	11 Volumes
#         PER	12 Persons
#         CQR	13 Nationalcurrency,currentprices,quarterlylevels
#         IND	14 Volumeandpriceindices
#       HRSSA	15 Hoursworked,seasonallyadjusted
#       VIXNB	16 Volumeindex,nationalbase/referenceyear
#     VNBQRSA	17 Nationalcurrency,constantprices,nationalbaseyear,quarterlylevels,seasonallyadjusted
#         CAR	18 Nationalcurrency,currentprices,annuallevels
#         GRW	19 Growthrates
#         HRS	20 Hoursworked
#     VIXNBSA	21 Volumeindex,nationalbase/referenceyear,seasonallyadjusted
#     VNBARSA	22 Nationalcurrency,constantprices,nationalbaseyear,annuallevels,seasonallyadjusted
#    CTQRGPSA	23 ContributionstoQ-o-QGDPgrowth,seasonallyadjusted
#     CPCARSA	24 USdollars,currentprices,currentPPPs,annuallevels,seasonallyadjusted
#         POP	25 Populationandemploymentmeasures
#       JOBSA	26 Jobs,seasonallyadjusted
#     VIXOBSA	27 Volumeindex,OECDreferenceyear,seasonallyadjusted
#     VOBARSA	28 Nationalcurrency,volumeestimates,OECDreferenceyear,annuallevels,seasonallyadjusted
#       LNBQR	30 Nationalcurrency,chainedvolumeestimates,nationalreferenceyear,quarterlylevels
#       VNBQR	31 Nationalcurrency,constantprices,nationalbaseyear,quarterlylevels
#    HCPCARSA	32 PerHead,US$,currentprices,currentPPPs,seasonallyadjusted
#       VNBAR	33 Nationalcurrency,constantprices,nationalbaseyear,annuallevels
#          CD	34 NationalcurrencyperUSDollar
#   VPVOBARSA	35 USdollars,volumeestimates,fixedPPPs,OECDreferenceyear,annuallevels,seasonallyadjusted
#  HVPVOBARSA	36 PerHead,US$,constantprices,fixedPPPs,OECDreferenceyear,seasonallyadjusted
























# get_data_structure("QNA") %>%
#   pluck("SUBJECT") %>%
#   filter(grepl("Gross domestic product", label)) %>%
#   {if (is_html_output()) print_table(.) else .}
# aa = get_dataset("QNA", filter = list(c("CAN"), c("B1_GE","GDP","B1_GA","B1_GI")))
# library("tibble")
# library("tidyverse")
# get_data_structure("QNA") %>%
#   pluck("SUBJECT") %>%
#   filter(grepl("Gross domestic product", label)) %>%



#core.cpi.file  <- "rawData/cansim_core_cpi_v41690926.csv"
# core.cpi.file  <- "rawData/cansim_core_cpi_v112593705 .csv"
# cpi.file       <- "rawData/cansim_cpi_v41690914.csv"
# cpi.back.file  <- "rawData/cansim_cpi_back_v41690973.csv"

# core.cpi.data  <- read.table(core.cpi.file, skip = 2, header = TRUE, sep = ',', stringsAsFactors = FALSE)
# cpi.data       <- read.table(cpi.file, skip = 2, header = TRUE, sep = ',', stringsAsFactors = FALSE)
# cpi.back.data  <- read.table(cpi.back.file, skip = 2, header = TRUE, sep = ',', stringsAsFactors = FALSE)

# TT = length(core.cpi.data$VALUE) # dim data
# T0 = min(which(!is.na(core.cpi.data$VALUE)))  # find first non-nan index
# core.cpi.data$VALUE[T0:TT]

# #----
# aout = get_data_structure("SNA_TABLE1")
# bout = get_data_structure("QNA")
# cout = get_data_structure("SNA_TABLE1_ARCHIVE")
#
# gdp.oecd <- get_dataset("QNA",filter = "CAN+ZAF.B1_GE+P7.DOBSA+VIXOBSA.Q",start_time = 1947, end_time = 2019)
# head(gdp.oecd,n=100)
#
# str(dstruc, max.level = 1)
#
# filter_list <- list("CAN", "B1_GA", c("V","VIXOB"))
# # GET Data from OECD from Quartely National accounts seasonally adjusted Index measure
# # ------------------------------------------------------------------------------------------------------------------------
# # filter_list = "AUS+CAN.GDP+B1_GE+B1_GA.CUR+IND+VIXOBSA.Q"
# # filter_list <-  list( c("CAN"), c("GDP", "B1_GE", "B1_GA"), c("VOL", "VIXOBSA"), c("Q"))
# # filter_list <-  list( c("CAN"), c("GDP", "B1_GE", "B1_GA"), measure, c("Q"))
#
# measure <- c("CUR" ,"CQRSA" ,"DNBSA" ,"GYSA" ,"LNBQRSA" ,"PERSA" ,"CARSA" ,"DOBSA" ,"GPSA" ,"LNBARSA" ,"VOL" ,"PER" ,"CQR" ,"IND" ,"HRSSA" ,"VIXNB" ,"VNBQRSA" ,"CAR" ,"GRW" ,"HRS" ,"VIXNBSA" ,"VNBARSA" ,"CTQRGPSA" ,"CPCARSA" ,"POP" ,"JOBSA" ,"VIXOBSA" ,"VOBARSA" ,"JOB" ,"LNBQR" ,"VNBQR" ,"HCPCARSA" ,"VNBAR" ,"CD" ,"VPVOBARSA" ,"HVPVOBARSA")
# subject <- c("GDP", "B1_GE", "B1_GS1", "B1_GA", "B1_GAI3","B1_GI")
# filter_list <-  list( c("CAN"), subject, measure, c("Q"))
#
# gdp.from.OECD.df  <- get_dataset("QNA",filter = filter_list, start_time = 1947, end_time = 2019)
# head2tail((gdp.from.OECD.df),100)
# From SNA_TABLE1_ARCHIVE
# subject <- c("B1_GA", "B1_GAI3","B1_GI")
# measure <- c("VIXOBSA")
# filter_list <-  list( c("CAN"), subject, measure)
# QNA_ARCHIVE.df <- get_dataset("QNA_ARCHIVE", filter = filter_list, start_time = 1947, end_time = 2019)
# print(aout)
#
# flist = list(c("CAN"), "VIXOB")
# df <- get_dataset("SNA_TABLE1", filter = flist)
# head(df)

# browse_metadata("SNA_TABLE1")

#----
# CAN.QGDP.VOLIDX.IDX.Q
# "GBR.CPGRLE01.IXOB.Q"
# gdp.oecd <- get_dataset("SNA_TABLE1",filter = "",start_time = 1947, end_time = 2019)

# get data from Statsitics Canada Final domestic demand = v62305753
# Chained (2012) dollars 	Seasonally adjusted at annual rates
# gdp.final.dd <- get_cansim_vector("v112593705","1984-01-01")
# if (!file.exists(store.IMF.IFS))  {
#   # print(" ... Downloading GDP data from IMF's IFS")
#   Q.CA.NGDP.query <- CompactDataMethod("IFS", queryfilter, "1947-01-01", "2019-12-31", FALSE)  # NY.GDP.MKTP.KD.ZG
#
#   print( paste0(" Done! saving data to:  ", paste0(raw.source.dir, ifs.gdp.file.name)))
#   # save locally
#   save(Q.CA.NGDP.query, file = paste0(raw.source.dir, ifs.gdp.file.name))
# } else {
#   # if exists, load from file
#   load(paste0(raw.source.dir, ifs.gdp.file.name))
#   print(" Loading/Reading GDP data")
# }
# gdps = cbind(gdp.data.ifs, gdp.fred.expenditure,  union = TRUE)
# print(gdps)  # FROM FRED
# Source: Organization for Economic Co-operation and Development  Release: Main Economic Indicators
# Units:  Index 2015=100, Seasonally Adjusted
# OECD descriptor ID: NAEXKP02# OECD unit ID: IXOBSA# OECD country ID: CAN

