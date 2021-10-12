% Get HLW Real time and current estimates from NYFRB website and read into matlab for plotting.
clear; clc; tic;
set(groot,'defaultLineLineWidth',2); % sets the default linewidth to 1.5;
% ADD LOCAL FUNCTION PATH WITHOUT SUBFOLDERS (IE. '_OLD') (only add what is needed)
addpath('local.Functions','utility.Functions');
GET_FROM_FRBNY		= 1;
SAVE_HLW_FACTORS	= 0;

% file locations
xls_store_dir = '../data/R.HLW.results/';
FRBNY_website	= 'https://www.newyorkfed.org/medialibrary/media/research/economists/williams/data/';
current_HLW		= 'Holston_Laubach_Williams_current_estimates.xlsx';
realtime_HLW	= 'Holston_Laubach_Williams_real_time_estimates.xlsx';

% to get the data from web
if GET_FROM_FRBNY == 1
	% get the xls files from the web and store localy
	websave([xls_store_dir current_HLW],	[FRBNY_website current_HLW]	);
	websave([xls_store_dir realtime_HLW],	[FRBNY_website realtime_HLW]);
	fprintf(' Data stored in directory: %s \n', xls_store_dir)
end

%% now read in the data. 
[dat0,tmp0] = xlsread([xls_store_dir current_HLW],'HLW Estimates');
dates0	= tmp0(7:end,1);
dates		= datetime( dates0,'InputFormat','dd.MM.yyyy' );
TT = numel(dates);
% iS = [1 6 11 16];
varNames = {'rstar';'g';'z';'cycle'};
iS = [16 6 11 1]; % varnams ordering
% date files
US = dat0(:,iS);
CA = dat0(:,iS+1);
EA = dat0(:,iS+2);
UK = dat0(:,iS+3);
cntry	= {'US';'CA';'EA';'UK'};

% loop through countries
for ii = 1:length(cntry)
	hlw_current.(cntry{ii}) = array2timetable(dat0(:,iS+ii-1),...
		'RowTimes',dates,'VariableNames',varNames);
end

% save current factors
if SAVE_HLW_FACTORS
	save([xls_store_dir 'hlw_current_FRBNY.mat'],'hlw_current')
	fprintf(' HLW current estimates saved to: %s \n', xls_store_dir);
end

%% Getting the real time estimats
[~,sheet_name] = xlsfinfo([xls_store_dir realtime_HLW]);
% for k=1:numel(sheet_name)
ns = numel(sheet_name);
us_mat = nan(TT,ns,4);
ca_mat = nan(TT,ns,4);
ea_mat = nan(TT,ns,4);
uk_mat = nan(TT,ns,4);

for k = ns:-1:1
	[dat0,tmp0]	= xlsread([xls_store_dir realtime_HLW],sheet_name{k});
	us_mat(:,k,:) = [dat0(:,iS);	nan(ns-k,4)];
	ca_mat(:,k,:) = [dat0(:,iS+1);nan(ns-k,4)];
	ea_mat(:,k,:) = [dat0(:,iS+2);nan(ns-k,4)];
	uk_mat(:,k,:) = [dat0(:,iS+3);nan(ns-k,4)];
end
%% now make a structure with 
hlw_realtime.('US') = us_mat;
hlw_realtime.('CA') = ca_mat;
hlw_realtime.('EA') = ea_mat;
hlw_realtime.('UK') = uk_mat;
hlw_realtime.info		= {	'Dimensions are: [T  Vintage Factor]'; ...
												'Factors are: [rstar, g, z, cycle]'; ...
												'Vintage dimension is from 2015:Q4 - now'};
hlw_realtime.dates		= dates;
hlw_realtime.vintages = sheet_name;
% save real time factors
if SAVE_HLW_FACTORS
	save([xls_store_dir 'hlw_realtime_FRBNY.mat'], 'hlw_realtime');
	fprintf(' HLW real-time estimates saved to: %s \n', xls_store_dir);
end 

%% some plotting
% series to plot
tp  = [1 9 17]; % time periods to plot
legNames = hlw_realtime.vintages(tp);
C = 1;	% countries:	{'US';'CA';'EA';'UK'};
F = 3;	% factor:			{'rstar';'g';'z';'cycle'};

s2p	= hlw_realtime.(cntry{C})(:,:,F);
clf;
f1 = figure(1); clf; set(f1,'WindowState','maximized','Position',[1441 641 1200 1843]);
	plot(s2p(:,tp))
	grid on; box on;	
	setplot(.22, 14,1)
	setdateticks(hlw_realtime.dates)
	hline(0);	set(gca,'GridLineStyle',':','GridAlpha',1/3);
	setoutsideTicks; add2yaxislabel;
	addlegend(legNames); 
	subtitle([cntry{C} ': factor = ' varNames{F}],-1.18)
	
	
	






















	
	
	
	
	
	
	
	
	





% EOF

