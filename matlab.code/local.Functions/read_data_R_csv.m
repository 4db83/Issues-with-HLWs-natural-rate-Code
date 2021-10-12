function data = read_data_R_csv(CSV_FILE_NAME, CSV_DIR_INPUT, GET_NEW_DATA)
% READ HLW MACRO DATA: INFLATION, REAL RATE AND INFLATION EXPECATIONS AND GDP GROWTH
% ---------------------------------------------------------------------------------------------
% CSV.DATA [] = country
% 	- '[].data.csv'
% created with 
% 	- 'get.[]data.R'.
% in the DIRECTORY:
% 				'./inputDatafromR'
% =============================================================================================

% clc; clear;
SetDefaultValue(3,'GET_NEW_DATA',0)
mat_file_name = [CSV_DIR_INPUT CSV_FILE_NAME '.mat'];

% if file does not exist or get new data
if GET_NEW_DATA || ~(exist(mat_file_name,'file')==2)
	% READ IN THE DATA FROM CSV FILE
	file2read			= [CSV_DIR_INPUT CSV_FILE_NAME '.csv'];
	[DAT0, TXT0]	= xlsread(file2read);

	% RE-ARRANGE DATA
	varNames	= TXT0(1,2:end);
	dates	= datetime( TXT0(2:end,1),'InputFormat','dd.MM.yyyy' );
	% %% hlw_data	= timetable(gdp, infl, infl_exp, int_rate, real_rate, 'RowTimes', dates);
	data	= array2timetable(DAT0,'RowTimes', dates, ...
							'VariableNames',strrep(varNames,'.','_')); 
	% add some user information to the dateime object
	data.Properties.UserData = {'Data is from R-File: get.rstar.data.2019.Q4.R'};
%  	head2tail(data)
	
	% SAVE THE DATA
	disp(['saving data to ' CSV_DIR_INPUT])
	save([CSV_DIR_INPUT CSV_FILE_NAME '.mat'], 'data');
else 
	% READ FROM FILE
	load([CSV_DIR_INPUT CSV_FILE_NAME '.mat'], 'data');
end

	
% EOF