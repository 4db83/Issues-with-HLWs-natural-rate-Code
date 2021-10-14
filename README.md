### Overview
<!-- Replication files for Buncic, D. (2021) "On a standard Method for Measuring the Natural Rate of Interest" -->

This repo contains estimation files in **Matlab** and **R** to replicate the results in Buncic (2021) "*On a standard Method for Measuring the Natural Rate of Interest*". All estimates reported in the paper use the *'earlier'* vintage data contained in the directory **'R.data.for.estimation.2020.May.28'**. Make sure to use that vintage of data to exactly replicate the numerical results that are reported in the paper. 

The provided **R** code can also be used to automatically update the empirical data for all four estimates, but be aware that it may be difficult to estimate the model satisfactorily because of the sharp drops in GDP due to the pandemic. The input data can be updated by using any of the **'get.XX.data.R'** R scripts, which will automatically download the "raw" source data and then perform the required transformations and splicing of data to produce the final data set use in the analysis.

<!-- The repo also provides data files, with the latest vintage of data stored in the directory **'R.data.for.estimation.2020.Oct.5'**. The R code can be used to automatically update the empirical data, but be aware that it may not be possible to successfully estimate the model due  -->
<!-- **NOTE:** All estimates reported in the paper use the *'earlier'* vintage data contained in the directory **'R.data.for.estimation.2020.May.28'**. Make sure to use that vintage of data to exactly replicate the results that are reported. If the input data is updated by using any of the **'get.XX.data.R'** R scripts, the results will be quantitatively different.  -->

### Matlab Users
If you are a **Matlab** user, go to the **matlab.code** directory and familiarize yourself with the file **'estimate_HLW.m'**. This is the main estimation script for all four countries and produces all Stage 2 and Stage 3 output, as well as the factors, plots etc. There are a number of options for printing and plotting, which need to be enabled/disabled by setting an indicator to 1 or 0 (see PRINT_FACTORS_EXCEL = 0, for instance).

#### Mex
I have added a mex file for the main KF loop using Armadillo C++ for the computation of the log-likelihood in the numerical optimisation. This is compiled for a Windows machine. I have not compiled it for Mac as yet. 

If you want your code to run somewhat faster, you can got to the **'kalman_filter_loglike.m'** file in the **'./utility.Functions/'** folder and uncomment line 189 to use **kalman_filter_loop_mex** instead of **kalman_filter_loop_fast** on line 188.

### R Users
If you are an **R** user, go to the **R.code** directory and familiarize yourself with the file **'fit.HLW.R'**. This is the main file for estimation in R. The estimation results are stored in various **'data/R.HLW.results/'** files for each country of interest and most results are also printed to screen. 


### The directories structure is as follows:
1. R.code
	- R.local.Functions: directory with modified R functions from HLWs R code needed for the estimation of the models
	- fit.HLW.R: main R script that replicates the results of the paper
	- get.CA.data.R: gets/makes the Canadian CA.data.csv data and saves the source data in source data dir
	- get.EA.data.R: gets/makes the Euro Area EA.data.csv data and saves the source data in source data dir
	- get.UK.data.R: gets/makes the UK UK.data.csv data and saves the source data in source data dir
	- get.US.data.R: gets/makes the US US.data.csv data and saves the source data in source data dir
	- full_estimation_output.txt: diary/log file with all the estimation output from R

2. data
	- R.data.for.estimation.2020.May.28: data directory created by R in csv format, used in/by Matlab for estimation 
	- R.HLW.results: various output files from running the R-code
	- R.source.data.2020.May.28: source data directory
	- US_trend_growth_1947Q1-2019Q4.mat: Mat file with various US trend growth estimates, not really needed, only for plotting purposes if an alternative benchmark is needed

3. matlab.code
	- _csminwel_output: directory csmin optimisation routine output files
	- local.Functions: directory with main local functions used in the estimate_HLW.m script
	- utility.Functions: directory with various Matlab utility functions that are called
	- estimate_HLW.m: main script that does the estimation and produces the output
	- get_HLW_factors.m: script that downloads/gets HLWs factors for various vintages from FRBNY website, mainly needed for plotting. 

4. matlab.output
	- CA_factors.xls: Canadian filtered and smoothed factors in xls file  
	- EA_factors.xls: Euro Area filtered and smoothed factors in xls file
	- UK_factors.xls: UK filtered and smoothed factors in xls file
	- US_factors.xls: US filtered and smoothed factors in xls file
	- CA_results.txt: all Matlab estimation results in stored in text file
	- EA_results.txt: all Matlab estimation results in stored in text file
	- UK_results.txt: all Matlab estimation results in stored in text file
	- US_results.txt: all Matlab estimation results in stored in text file


Stockholm, 14.10.2021
