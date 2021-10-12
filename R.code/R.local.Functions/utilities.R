#------------------------------------------------------------------------------#
# File:        utilities.R
#
# Description: This file contains basic functions that will be used throughout
#              HLW.
#------------------------------------------------------------------------------#

shiftQuarter <- function(original.start,shift){
#################################################################
# This function takes in a (year,quarter) date in time series format
# and a shift number, and returns the (year,quarter) date corresponding
# to the shift. Positive values of shift produce leads and negative values
# of shift produce lags.
# For example, entering 2014q1 with a shift of -1 would return 2013q4.
# Entering 2014q1 with a shift of 1 would return 2014q2.
# In each case, the first argument of the function must be entered as
# a two-element vector, where the first element corresponds to the year
# and the second element corresponds to the quarter.
# For example, Q12014 must be entered as "c(2014,1)".
################################################################    

# Leads (positive values of shift)
    if (shift > 0) {
        new.start = c(0,0)
        sum = original.start[2] + shift
    
        # Get the year value
        if (sum <= 4) {
            new.start[1] = original.start[1]
        }
        else {
            new.start[1] = original.start[1] + ceiling(sum/4) - 1
        }

        # Get the quarter value
        if (sum %% 4 > 0) {
            new.start[2] = sum %% 4
        }
        else {
            new.start[2] = sum %% 4 + 4
        }
    }

# Lags (negative values of shift)
    else {
        new.start = c(0,0)
        diff = original.start[2] - abs(shift)
    
        # Get the year value
        if (diff > 0) {
            new.start[1] = original.start[1]
        }
        else {
            new.start[1] = original.start[1] - (1 + floor(abs(diff)/4))
        }

        # Get the quarter value
        if (diff %% 4 > 0) {
            new.start[2] = diff %% 4
        }
        else {
            new.start[2] = diff %% 4 + 4
        }
    }
        
return(new.start)}


shiftMonth <- function(original.start,shift){
#################################################################
# This function takes in a (year,month) date in time series format
# and a shift number, and returns the (year,month) date corresponding
# to the shift. Positive values of shift produce leads and negative values
# of shift produce lags.
# For example, entering 2014m1 with a shift of -1 would return 2013m12.
# Entering 2014m1 with a shift of 1 would return 2014m2.
# In each case, the first argument of the function must be entered as
# a two-element vector, where the first element corresponds to the year
# and the second element corresponds to the month.
# This function is analogous to shiftQuarter().
################################################################    

# Leads (positive values of shift)
    if (shift > 0) {
        new.start = c(0,0)
        sum = original.start[2] + shift
    
        # Get the year value
        if (sum <= 12) {
            new.start[1] = original.start[1]
        }
        else {
            new.start[1] = original.start[1] + ceiling(sum/12) - 1
        }

        # Get the month value
        if (sum %% 12 > 0) {
            new.start[2] = sum %% 12
        }
        else {
            new.start[2] = sum %% 12 + 12
        }
    }

# Lags (negative values of shift)
    else {
        new.start = c(0,0)
        diff = original.start[2] - abs(shift)
    
        # Get the year value
        if (diff > 0) {
            new.start[1] = original.start[1]
        }
        else {
            new.start[1] = original.start[1] - (1 + floor(abs(diff)/12))
        }

        # Get the month value
        if (diff %% 12 > 0) {
            new.start[2] = diff %% 12
        }
        else {
            new.start[2] = diff %% 12 + 12
        }
    }
        
return(new.start)}


getFRED <- function(url, freq = "Quarterly") {
##########################################################################################
# This function downloads data from FRED. It returns quarterly data.
# User must provide the FRED url.
########################################################################################### 
    # Download the data from FRED
    
    #download.file(url, destfile = 'FREDtemp.txt', method = "wget")
    #FREDraw <- readLines('FREDtemp.txt')
    
    txt.file.name <- paste0("rawData/",substr(url, regexpr('[a-zA-z0-9]*.txt',url),1000))
    if (!file.exists(txt.file.name)){
        # Download the data from FRED
        #download.file(url, destfile = 'FREDtemp.txt', method = "wget")
        system(paste0('wget --no-check-certificate "', url, '"'))
        system(paste('mv',substr(url, regexpr('[a-zA-z0-9]*.txt',url),1000),txt.file.name))
    }
    FREDraw <- readLines(txt.file.name) 

    # Frequency
    freq.FRED <- gsub(' ', '',substr(FREDraw[which(regexpr('Frequency', FREDraw)==1)],
                                     (nchar('Frequency')+2),100))    

    # Where does the data start
    datastart = which(gsub(' ', '',FREDraw)=='DATEVALUE') - 2

    #data <- read.table('FREDtemp.txt', skip = datastart, header = TRUE)
    data <- read.table(txt.file.name, skip = datastart, header = TRUE)

    first.year  <- as.numeric(format(as.Date(data$DATE[1]),'%Y'))
    first.month <- as.numeric(format(as.Date(data$DATE[1]),'%m'))
    
    # Adjust frequency
    if (freq.FRED == 'Quarterly'){
        first.q  <- (first.month-1)/3 + 1
        data.tis <- tis(data$VALUE, start = c(first.year, first.q), tif = 'quarterly')
    } else if (freq.FRED == 'Monthly') {
        data.tis <- tis(data$VALUE, start = c(first.year, first.month), tif = 'monthly')
    }

    # Convert frequency
    if (freq.FRED == 'Monthly' & freq == 'Quarterly') {
        data.tis <- convert(data.tis, tif = 'quarterly', method = 'constant', observed. = 'averaged')
    }

    return(data.tis)
} 


splice <- function(s1, s2, splice.date, freq) {
##########################################################################################
# This function splices two series, with the series s2 beginning at splice.date
# and extended back using the growth rate at the splice.date times series s1
# The freq argument accepts two values - 'quarterly' and 'monthly' -
# but it could be modified to take more.
##########################################################################################    
    t <- splice.date #renaming for convenience
    if (freq == "quarterly" | freq == "Quarterly") {
        t.minus.1 <- shiftQuarter(t,-1)
    }
    else if (freq == "monthly" | freq == "Monthly") {
        t.minus.1 <- shiftMonth(t,-1)
    }
    else { stop("You must enter 'quarterly' or 'monthly' for freq.") }
    ratio <- as.numeric(window(s2,start = t, end = t)/
                        window(s1,start = t, end = t))

    return(mergeSeries(ratio*window(s1,end = t.minus.1),window(s2, start = t)))
}


gradient <- function(f, x, delta = x * 0 + 1.0e-5) {
##########################################################################################
# This function computes the gradient of a function f given a vector input x.
##########################################################################################   
    g <- x * 0
    for (i in 1:length(x)) {
        x1 <- x
        x1[i] <- x1[i] + delta[i]
        f1 <- f(x1)
        x2 <- x
        x2[i] <- x2[i] - delta[i]
        f2 <- f(x2)
        g[i] <- (f1 - f2) / delta[i] / 2
    }
    return(g)
}
