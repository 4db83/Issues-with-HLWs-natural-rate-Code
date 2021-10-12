This directory contains estimation files in Matlab as well as R to replicate the results in the paper "On Measuring the Natural Rate of Interest" by Daniel Buncic. Data files are also provided, with the latest vintage data in the directory 'R.data.for.estimation.2020.Oct.5'. 

NOTE: All estimates my paper use the earlier vintage data 'R.data.for.estimation.2020.May.28' to exactly replicate the R-File results from HLWs website at the FRBNY. 

If you are a Matlab users, go to the matlab.code directory and familiarize yourself with the file 'estimate_HLW.m'. This is the main estimation script for all four countries and produces all Stage 2 and Stage 3 output as well as factors, plots etc. There are a number of options for printing and plotting. 
 
If you are a R users, go to the R.code directory and familiarize yourself with the file 'fit.HLW.R'. This is the main file for estimation in R. The estimation results are stored in various 'data/R.HLW.results/' files for each country of interest and most results are also printed to screen. 

Let me know if there are any issues. 
---------------------------------------
Daniel Buncic, Stockholm.
Last updated: 06.10.2020.
---------------------------------------


The directories are as follows:
1) data
		- R.data.for.estimation.2020.May.28 -> data created by R in csv format, used in/by Matlab for estimation 
		- R.data.for.estimation.2020.Oct.5 	-> data created by R in csv format, used in/by Matlab for estimation
		- R.HLW.results -> various output files from running the R-code
		- R.source.data.2020.May.28 -> source data files
		- R.source.data.2020.Oct.5	-> source data files
		- US_trend_growth_1947Q1-2019Q4.mat -> Mat file with various US trend growth estimates, not really needed, only for plotting purposes if an alternative benchmark is needed
2) matlab.code
		- _csminwel_output -> csmin optimisation routine output files
		- local.Functions	-> directory with main local functions used in the estimate_HLW.m script
		- utility.Functions -> various Matlab utility functions that are called
		- estimate_HLW.m -> main script that does the estimation and produces the output
		- get_HLW_factors.m -> scritp that gets HLWs factors for various vintages from FRBNY website, mainly needed for plotting. 
3) R.code
		- R.local.Functions -> directory with modified R functions from HLWs R code needed for the estimation of the models
		- fit.HLW.R -> main R script that replicates the results of the paper
		- get.CA.data.R -> gets/makes the Canadian CA.data.csv data and saves the source data in source data dir
		- get.EA.data.R -> gets/makes the Euro Area EA.data.csv data and saves the source data in source data dir
		- get.UK.data.R -> gets/makes the UK UK.data.csv data and saves the source data in source data dir
		- get.US.data.R -> gets/makes the US US.data.csv data and saves the source data in source data dir
		- full_estimation_output.txt -> diary/log file with all the estimation output from R
4) matlab.output
		- CA_factors.xls -> Canadian filtered and smoothed factors in xls file  
		- EA_factors.xls -> Euro Area filtered and smoothed factors in xls file
		- UK_factors.xls -> UK filtered and smoothed factors in xls file
		- US_factors.xls -> US filtered and smoothed factors in xls file
		- CA_results.txt -> all Matlab estimation results in stored in text file
		- EA_results.txt -> all Matlab estimation results in stored in text file
		- UK_results.txt -> all Matlab estimation results in stored in text file
		- US_results.txt -> all Matlab estimation results in stored in text file


