% ESTIMATES THE HLW model with both correct and incorrct Stage 2 procedures, as well as
% \sigma(g) and \sigma(z) directly.
% ---------------------------------------------------------------------------------------------------
% Reminder: their SSF is
% ---------------------------------------------------------------------------------------------------
% 	Y(t) = A*X(t) + H*Xi(t) + v(t), Var[v(t)] = R.
%  Xi(t) = F*Xi(t-1) + S*e(t),			Var[e(t)] = Q.
% ---------------------------------------------------------------------------------------------------
% MY STATE SPACE FORM IS:
% 		Observed:	Y(t)			= D(t) + M*alpha(t)			+ e(t);		Var[e(t)] = H.
% 		State:		alpha(t)	= C(t) + Phi*alpha(t-1)	+ S*n(t);	Var[n(t)] = Q.
% ===================================================================================================
clear; clc; close all; tic;
set(groot,'defaultLineLineWidth',2); % sets the default linewidth to 2 in plots;
%% ADD LOCAL FUNCTION PATH WITHOUT SUBFOLDERS (IE. '_OLD') (only add what is needed)
addpath('local.Functions','utility.Functions');
cntr_ = {'US', 'EA', 'UK', 'CA'};

% MAKE OUTPUT/GRAPHICS DIRECTORIES IF THEY DO NOT EXIST
matlab_output_dir		= '../matlab.output/';
latex_graphics_dir	= '../../graphics/';
latex_table_dir			= '../../table.input/';
% making dirs
if ~exist(matlab_output_dir,	'dir');	mkdir(matlab_output_dir);		end
if ~exist(latex_graphics_dir,	'dir');	mkdir(latex_graphics_dir);	end
if ~exist(latex_table_dir,		'dir');	mkdir(latex_table_dir);			end

% DEFINDE WHICH VINTAGE OF DATA TO BE USE IN ESTIMATION. THESE ARE STORED IN DIFFERENT DIRECTORIES, 
DATA_DIR_INPUT	= '../data/R.data.for.estimation.2020.May.28/';		% DATA ENDS IN Q4-2019

for CI = 1
% for CI = 1:4
	% DEFINE COUNTRY, IE 
	COUNTRY		= cntr_{CI};
	% DEFINE SAMPLE-END. NOTE IF DATA IN CSV/MAT FILE WHICH IS READ IS LESS THAN SMPL_END, DATA ENDS THERE.
	SMPL_END	= 'Q4-2019'; 
	% SET VARIOUS PLOTTING AND DATA READING/WRITING 
	CSV_READ			= 0;				% SET TO 1 TO READ NEW DATA FROM CSV FILE, OTHERWISE LOAD THE .MAT CONVERTED FILE
	PLOT_ON				= 0;
	PLOT_F_STATS	= 1;
	PLOT_GDP_RR		= 0;
	PRINT_PLOT_TO_PDF		= 0;
	PRINT_FACTORS_EXCEL = 0;
	PRINT_RESULTS_LATEX = 0;
	% KEEP A LOG FILE USING THE DIARY FUNCTION
	DIARY_ON = 0;

	% WRITE RESULTS TO LOG FILE
	if DIARY_ON
		% save as textfile
    %   diary_file_name = ['./_diary_dir/' char(cntr_(CI)) '_results.txt'];
		% save to _output folder as well
		diary_file_name = [matlab_output_dir char(cntr_(CI)) '_results.txt'];
		if exist(diary_file_name,'file'); delete(diary_file_name); end
		diary(diary_file_name)
	end
	
	%% OPTIMISATION SETTINGS FOR FMINUNC AND FMINCON
	% ***************************************************************************************************
	% NOTE: UK ESTIMATES SENSITIVE TO NUMERICAL ALGORITHM. USE SQL IN
	% options_con FOR fmincon, RATHER THAN THE SUGGESTED 'interior-point'.
	% ---------------------------------------------------------------------------------------------
	% FMINUNC HAS FOLLOWING ALGORITHMS: 
	% ---------------------------------------------------------------------------------------------
	% (1) 'quasi-newton' (default); (2) 'trust-region'.
	% ---------------------------------------------------------------------------------------------
	% Use optimoptions to set the Algorithm option at the command line.Recommendations: If your 
	% objective function includes a gradient, use 'Algorithm' = 'trust-region', and set 
	% the SpecifyObjectiveGradient 
	% option to true. Otherwise, use 'Algorithm' = 'quasi-newton'.
	% ---------------------------------------------------------------------------------------------
	options_unc = optimoptions(@fminunc, ...
								'Algorithm'								,	'quasi-newton' , ... 
								'Display'									, 'none'	, ...
								'MaxIterations'						, 5e5			, ... 
								'MaxFunctionEvaluations'	, 5e5			, ... 
								'OptimalityTolerance'			,	1e-8		, ...
								'FunctionTolerance'				, 1e-8		, ...
								'StepTolerance'						, 1e-8		, ...
								'FiniteDifferenceType'		, 'central');
	% 						'HessianApproximation'		, 'bfgs'	, ...						
	% 						'Algorithm'		, 'trust-region', ...
	% 						'Algorithm'		, 'quasi-newton', ...
	% ---------------------------------------------------------------------------------------------
	% FMINCON HAS FOLLOWING ALGORITHMS: 
	% ---------------------------------------------------------------------------------------------
	% (1) 'interior-point' (default); (2) 'trust-region-reflective'; (3) 'sqp'; (4) 'sqp-legacy';
	% (5) 'active-set'
	% ---------------------------------------------------------------------------------------------
	% Use the 'interior-point' algorithm first.
	% For help if minimization fails, see When the Solver Fails or When the Solver Might Succeed. 
	% Try 'sqp' next, and 'active-set' last. Use 'trust-region-reflective' when applicable.	
	% ---------------------------------------------------------------------------------------------
	% Set Algorithm to 'active-set' to get their estimates, 
	% changing this to 'interior-point', to get the higher loglike value
	% ---------------------------------------------------------------------------------------------
	options_con = optimoptions(@fmincon, ... 
								'Algorithm'								,	'sqp'   , ... 
								'Display'									, 'none'	, ...
								'MaxIterations'						, 5e5			, ... 
								'MaxFunctionEvaluations'	, 5e5			, ... 
								'HessianApproximation'		, 'bfgs'	, ...
								'OptimalityTolerance'			,	1e-8		, ...
								'FunctionTolerance'				, 1e-8		, ...
								'StepTolerance'						, 1e-8		, ...
								'FiniteDifferenceType'		, 'central');
	% 						'Algorithm'		, 'trust-region-dogleg', ...
	% 						'Algorithm'		, 'active-set', ...												
	% 						'Algorithm'		, 'quasi-newton', ...
	% ***************************************************************************************************
	% PARAMETER NAMES
	par_names3 = {'a_y1'    ; ...
								'a_y2'    ; ...
								'a_r'     ; ...
								'b_pi'    ; ...
								'b_y'     ; ...
								'sigma_y~'; ...
								'sigma_pi'; ...
								'sigma_y*'; ...
								'sigma_g' ; ...
								'sigma_z' ; ...
								'Log-Like'; ...
								'Lambda_g'; ...
								'Lambda_z'};
	% ***************************************************************************************************

	%% LOAD THE HWL INPUT DATA TABLE
	% ---------------------------------------------------------------------------------------------------
	if strcmp(COUNTRY,'EA') 
		hlw_read_sample = '1972Q1.to.2019Q4';	
	else 
		hlw_read_sample = '1961Q1.to.2019Q4'; 
	end
	
	DATA_DIR_NAME	= [COUNTRY '.data'];
	%	READ IN DATA FROM MATLAB OR CSV FILE
	data = read_data_R_csv(DATA_DIR_NAME, DATA_DIR_INPUT, CSV_READ);
	% head2tail(data)
	% MAKE ANONYMOUS FUNCTION TO GET SIGMA/LAMBDA
	S3g	= @(theta,L_g) abs(theta(8)*L_g);
	S2g	= @(theta,L_g) abs(theta(10)*L_g);
	L2g	= @(theta) theta(11)/theta(10);
	L3g	= @(theta) theta(9)/theta(8);
	L3z	= @(theta) abs( theta(10)*theta(3)/theta(6) );
	S3z	= @(theta,L_z) abs(theta(6)*L_z/theta(3));
	% now for readign some data
	Sa00 = @(x) (['Stage' num2str(x) '.xi00.']	);		% 'x00.';
	SP00 = @(x) (['Stage' num2str(x) '.P00.']	);			% 'P00.';
	Sstr = @(x) (['Stage' num2str(x) '.theta.']);			% 'theta.';
	S2Lz = @(x) (['Stage' num2str(x) '.Lambda.z.']);	% 'Lambda.z.';
	% HLW R-code data output of parameters, some data etc, which is used as an input.
	HLW_INPUT_DATA_DIR = '../data/R.HLW.results/';

	% CHANGE SAMPLE SIZE IF NEEDED (use timerange function to shorten sample).
	data = data(timerange('Q1-1960',SMPL_END,'closed'),:);
% 	data = data(timerange('Q1-1960','Q4-2007','closed'),:); % does not matter if beg-sample is less than in data for EA.

	% ---------------------------------------------------------------------------------------------------
	% HLW: LOAD PARAMETER ESTIMATES FROM THEIR R FILES.
	% ---------------------------------------------------------------------------------------------------
	% function handle to load the HLW R output paramters etc
	getS = @(x) ([HLW_INPUT_DATA_DIR COUNTRY '/' x '.csv']);
	HLW_Stage3_R_File		= xlsread( getS([Sstr(3) hlw_read_sample]) );
	% THESE ARE THE INTIAL VALUES USED BY HLW (FIRST COLUMN OF HLW2017_THETA)
	HLW_initvals_R_File	= HLW_Stage3_R_File(1:8,1);
	% these are their estimates
	HLW_theta_R_File		= HLW_Stage3_R_File(1:8,2);

	% SET OTHER OPTIONAL PARAMETERS NEEDED IN PROCEDURE (MW,EW,QLR, Lambda_z Col(1) MUE.stat Col(2))
	% HLW_S2_Lz_	= dlmread([HLW_input_data_dir 'Stage2.Lambda_z.' SMPL_END '.csv'],',',1,1);
	Lambda_g = HLW_Stage3_R_File(end-1,2);
	Lambda_z = HLW_Stage3_R_File(end  ,2);

	% Stage 3 other paramters to be parsed to the function
	other_pars_S3.a00 = dlmread( getS([Sa00(3) hlw_read_sample]) ,',',0,0);
	other_pars_S3.P00 = dlmread( getS([SP00(3) hlw_read_sample]) ,',',0,0);
	% SIMPLE 0.2*I AS IN THEIR INITIALISATIONS
	% other_paras.P00 = 0.2*eye(size(other_paras.a00,1));
	other_pars_S3.Lambda_g	= Lambda_g;			%%% some reference values 0.0538690377714457
	other_pars_S3.Lambda_z	= Lambda_z;			%%% some reference values 0.0302172212554928
	% stage 2 mue Lambdazs and test statistics, order MW, EW, QLR
	S2MUEs = dlmread( getS([S2Lz(2) hlw_read_sample]) ,',',1,1);

	%% ARANANGE DATA INTO APPROPRIATE Y AND X
	% ---------------------------------------------------------------------------------------------------
	% prepare data for the matrices
	lnGDP	= 100*data.gdp_log;			        % ln(GDP) times 100: NOTE GDP already logged in HLW's R-File
	INFL	=			data.inflation;				    % INFLATION SERIES UNTRANSFORMED
	RR		=			data.real_rate;						% REAL RATE: 

	% defined lagged variables to be used in stage 3 SSF
	GPD_2_lags	= mlag(lnGDP,2);	        % [y_{t-1} y_{t-2}]
	RR_2_lags		= mlag(RR,2);			        % [r_{t-1} r_{t-2}]
	INFL_4_lags = mlag(INFL,4);		        % [\pi_{t-1} \pi_{t-2} \pi_{t-3} \pi_{t-4}]

	% MEASUREMENT MATRIX (2xT) TO BE PASSED TO KF ROUTINE. NOTE THE DATA TRIMMED TO (5:END)
	% (SO FROM 1961:Q1 ONWARDS)
	YY	= [lnGDP INFL]; 
	Y		=	YY(5:end,:);
	% EXOGENEOUS VECTOR X NEEDED FOR DT = A*X (7xT) 
	% [GDP(t-1) GDP(t-2) r(t-1) r(t-2) pi(t-1) mean(pi(t-2:t-4))]
	XX	= [ GPD_2_lags RR_2_lags INFL_4_lags(:,1) mean(INFL_4_lags(:,2:4),2) ];
	X		= XX(5:end,:);   
	% sample size
	T		= length(Y);

	% PRINT DATES OF TIME PERIOD CONSIDERED (5:end) is 1961:Q1.
	Dates	= datenum(data.Time(5:end));
	sep('=')
	fprintf( '		Country is %s: ', cntr_{CI}) 
	fprintf( ' Sample period is: %s', datestr(Dates(1),'yyyy:qq'));
	fprintf( ' to %s ', datestr(Dates(end),'yyyy:qq'));
	fprintf( '| Vintage of Data is %s: \n', DATA_DIR_INPUT(31:end-1) );
	sep('=')

	%% SET INITIAL VALUES OF PARAMETERS FOR STAGE 2 & 3 MODELS (based on full sample from 1960:Q1)
	% ---------------------------------------------------------------------------------------------------
	% run HP Filter on Y to get initial crude trend and cycle decomposition
	[HP_cycle, HP_trend] = hp_filter(YY(:,1), 36000); % Lambda = 36000 is in HLW paper
	inflt = YY(:,2) ; 
	% RUN OLS ON HP YCYCL AND RT_1 TO GET INTIAL ESTIMATES A_ AND B_ AND THEIR STANDARD ERRORS
	% y~ equation
	a_y		= fullols(HP_cycle, [ mlag(HP_cycle,2) mean(RR_2_lags,2) ], 1);
	% pi equation
	b_pi	= fullols( inflt - mean(INFL_4_lags(:,2:4),2) , ...
				 [ ( lag(inflt) - mean(INFL_4_lags(:,2:4),2) ) lag(HP_cycle)], 1);
	% Stage 2 output gap equation is a(L)y_t = a0 + ar/2(r(t_1)+r(t-2)) + ag*g(t-1) + e(t)		 
	a_y2	= fullols(HP_cycle, [ mlag(HP_cycle,2) mean(RR_2_lags,2) lag(delta(HP_trend))], 0);

	% SET INITIAL VALUES TO BE USED (same order to match initVals of HLW)
	% Theta2 vector is: [ a_y1,a_y2,a_r,a_0,a_g, b_pi,b_y, sig_y~,sig_pi,sig_y* ]
	%										[   1 ,  2 , 3 , 4 , 5 ,  6  , 7 ,   8   ,  9   ,  10   ] 
	initVals2 = [	a_y2.bhat([2:4 1 end]); 
								b_pi.bhat; 
								sqrt(a_y.sig2); 
								sqrt(b_pi.sig2);
								std( diff(HP_trend) ) ];

	% Theta3 vector is: [a_y1,a_y2,a_r,b_pi,b_y, s_y~,s_pi,s_y*]
	% 									[   1,   2,  3,   4,  5,    6,   7,  8 ] 
	initVals3 = [	a_y.bhat; 
								b_pi.bhat; 
								sqrt(a_y.sig2); 
								sqrt(b_pi.sig2);
								std( diff(HP_trend) ) ];
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%																			STAGE 2        
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%		Theta2 vec.:	[ a_y1,a_y2,a_r,a_0,a_g, b_pi,b_y, sig_y~,sig_pi,sig_y* ]
	%									[   1 ,  2 , 3 , 4 , 5 ,  6  , 7 ,   8   ,  9   ,  10   ] 
	% ---------------------------------------------------------------------------------------------------
	% make stage 2 paramete names
	par_names2	= [	par_names3([1 2 3]);{'a_0';'a_g'}	; 
									par_names3([4 5 6 7 8 9 11 12])  ];
	% make data with intercept in last entry
	data_in_C		= [Y X ones(T,1)]';		
	% ---------------------------------------------------------------------------------------------------
	HLW_Stage2_R_File	= xlsread( getS([Sstr(2) hlw_read_sample]));

	% read from Stage 2 R-Code output
	other_pars_S2.a00 = dlmread( getS([Sa00(2) hlw_read_sample]) ,',',0,0);
	other_pars_S2.P00 = dlmread( getS([SP00(2) hlw_read_sample]) ,',',0,0);
	other_pars_S2.Lambda_g	= Lambda_g;		
	% S2 for the other_paras input for the correct Stage 2 model M0.g
	other_pars_S2_M0g.a00	= other_pars_S3.a00(1:end-2);
	other_pars_S2_M0g.P00	= other_pars_S3.P00(1:end-2,1:end-2);
	other_pars_S2_M0g.Lambda_g	= Lambda_g;

	% STAGE 2: SET BOUNDS ON THE PARAMETER SPACE IF NEEDED
	% ---------------------------------------------------------------------------------------------------
	% In HLW, (b_y) b2.constraint = 0.025, (a_r) a3.constraint = -0.0025
	% ---------------------------------------------------------------------------------------------------
	A_	= [	1 1; -1 1]; bb	= [.99 .99];
	S2k	= size(initVals2,1);
	LB2	= struct; UB2	= struct; AA2	= struct;
	LB2.bl	= -Inf*ones(S2k,1);			LB2.bl(7)		=	 0.0250; 	% (b_y) b2.constraint >  0.0250	
	UB2.bl	=  Inf*ones(S2k,1);			UB2.bl(3)		= -0.0025; 	% (a_r) a3.constraint < -0.0025
	AA2.bl	= [A_ zeros(2,length(LB2.bl)-2)];
	% TO FORCE POSITIVE SIGMA (not needed, leads to marginally worse numerical results)
	% LB2.bl(end-2:end)		= 0;															% [sigma_y~ sigma_pi sigma_y*] > 0 

	%	MLE SIGMA_G 
	LB2.g		= -Inf*ones(S2k+1,1);		LB2.g(7)		=	 0.0250;	% (b_y) b2.constraint >  0.0250		
	UB2.g		=  Inf*ones(S2k+1,1);		UB2.g(3)		= -0.0025;	% (a_r) a3.constraint < -0.0025                                               
	AA2.g		= [A_ zeros(2,length(LB2.g)-2)];
	% TO FORCE POSITIVE SIGMA (not needed, leads to marginally worse numerical results)
	% LB2.g(end-3:end)		= 0;															% [sigma_y~ sigma_pi sigma_y* sigma_g] > 0

	% MLE M(0).(g) Correct Stage 2 model
	LB2.M0g = -Inf*ones(S2k-1,1); 	LB2.M0g(5)	=  0.0250;	% (b_y) b2.constraint >  0.0250		
	UB2.M0g =  Inf*ones(S2k-1,1); 	UB2.M0g(3)	= -0.0025;	% (a_r) a3.constraint < -0.0025                                               
	AA2.M0g	=	[A_ zeros(2,length(LB2.M0g)-2)];
	% TO FORCE POSITIVE SIGMA (not needed, leads to marginally worse numerical results)
	% LB2.M0g(end-3:end)	= 0;															% [sigma_y~ sigma_pi sigma_y* sigma_g] > 0

	% allocate structure space 
	LL2			= struct();		gradOpt2 = struct(); 
	theta2	= struct();		hessOpt2 = struct();

	% STAGE 2 MLE ESTIMATION
	% ---------------------------------------------------------------------------------------------------
	% LL at initival values
	LL2.initvals	= LogLike_Stage2_HLW_SSF(initVals2, data_in_C, other_pars_S2);
	% using csminwell
	[theta2.CS, LL2.CS, ~,~, gradOpt2.CS, hessOpt2.CS] = csminwel('LogLike_Stage2_HLW_SSF', initVals2, eye(length(initVals2)), [], 1e-8, 1e5, data_in_C,	other_pars_S2);
	% using fminunc (without constraints on parameters)
	[theta2.O,	LL2.O,  ~,~, gradOpt2.O, hessOpt2.O]	 = fminunc(@LogLike_Stage2_HLW_SSF,	initVals2, options_unc, data_in_C,	other_pars_S2);
	% S2 baseline model
	[theta2.bl,		LL2.bl]		= fmincon(@LogLike_Stage2_HLW_SSF,	theta2.O,							AA2.bl,	bb	,[],[],LB2.bl,	UB2.bl	,[],options_con,	data_in_C,						other_pars_S2);
	% estimate sigma_g directly by MLE in S2 baseline model
	[theta2.g,		LL2.g]		= fmincon(@LogLike_Stage2_HLW_SSF_g,[theta2.bl; .03],			AA2.g,	bb	,[],[],LB2.g,	UB2.g			,[],options_con,	data_in_C,						other_pars_S2);
	% MLE M(0).g Correct Stage 2 model
	[theta2.M0g,	LL2.M0g]	= fmincon(@LogLike_Stage2_HLW_M0g,	theta2.g([1:3 6:end]),AA2.M0g,bb	,[],[],LB2.M0g,	UB2.M0g ,[],options_con,	data_in_C(1:end-1,:),	other_pars_S2_M0g);
	% unconstrained optimisation M(0).g model
	% [theta2.M0ga,	LL2.M0ga]	= fminunc(@LogLike_Stage2_HLW_M0g,	theta2.g([1:3 6:end]),	options_unc, data_in,		other_paras_S2_M0g);
  
  % -----------------------------------------------------------------------------------------------------------
  % ADD THE MLE M(0).gC Correct Stage 2 model Plus an intercept term (a0) to show that it is not important. 
  % -----------------------------------------------------------------------------------------------------------
  %% 
%   clc
  [theta2.M0g_C,	LL2.M0g_C]	= fminunc(@LogLike_Stage2_HLW_M0g_C, [-.2; theta2.M0g], options_unc, data_in_C, other_pars_S2_M0g);
 
	% Combine for printing the STAGE 2 output
	sep(120)
	S2_initVals		= [initVals2; NaN; -LL2.initvals; NaN];
	S2_Rfile			=  HLW_Stage2_R_File(:,2);
	S2_replicated	= [theta2.bl; S2g( theta2.bl , Lambda_g ); -LL2.bl; Lambda_g];
	S2_MLE_g			= [theta2.g;		-LL2.g;			L2g(theta2.g)];
	S2_MLE_M0g		= [plugin(theta2.M0g,NaN(2,1),4);  -LL2.M0g;  L3g(theta2.M0g)];
  S2_MLE_M0g_C	= [plugin(theta2.M0g_C(2:end),[theta2.M0g_C(1);NaN(1,1)],4);  -LL2.M0g_C;  L3g(theta2.M0g_C(2:end))];
	
	% make Table
	tabS2 = table(par_names2, S2_initVals, S2_Rfile, S2_replicated, S2_MLE_g, S2_MLE_M0g, S2_MLE_M0g_C);% , S2_replicated, S2_MLE_g);
	fprintf('														Stage 2 Results           \n')
	sep(120)
	print2screen(tabS2, '%16.10f')

	%% **********************************************************************************************************
	%																	MUE FOR STAGE 2 begins here 
	%************************************************************************************************************
	% MAKE THE STAGE 2 KFS OBJECTs for MUE.S2(z) NOW
	[~,S2_KFS.bl]	 = LogLike_Stage2_HLW_SSF  ( theta2.bl,		data_in_C,						other_pars_S2,     1);
	[~,S2_KFS.g]	 = LogLike_Stage2_HLW_SSF_g( theta2.g,		data_in_C,						other_pars_S2,     1);
	[~,S2_KFS.M0g] = LogLike_Stage2_HLW_M0g	 ( theta2.M0g,	data_in_C(1:end-1,:),	other_pars_S2_M0g, 1);   % [11]

	% CALL THE MUE Stage 2 functions
	[L2z_bl,	Chow2_bl,		eout_S2_bl]		= call_MUE_HLW_stage2    ( theta2.bl,		S2_KFS.bl,	YY, XX); 
	[L2z_g,		Chow2_g	,		eout_S2_g ]		= call_MUE_HLW_stage2		 ( theta2.g,		S2_KFS.g ,	YY, XX); 
	[L2z_M0g,	Chow2_M0g,	eout_S2_M0g]	= call_MUE_HLW_stage2_M0g( theta2.M0g,	S2_KFS.M0g, YY, XX);
	% ---------------------------------------------------------------------------------------------------
	% NOTE ABOVE: FOR THE CORRECT STAGE 2 MODEL MLE(SIGMA_G).M0, WE NEED TO CALL DIFFERENT MUE FUNCTION 
	% BECAUSE THE WAY GY(t) IS CONSTRUCTED IS DIFFERENT: GY(t) AY(L)Y~ -AR(L)[R(t)-4G(t)] in sM0g Model
	% ---------------------------------------------------------------------------------------------------
	
	%% PRINT NICE MUE2 OUTPUT IF NEEDED
	%************************************************************************************************************
	sep(140); fprintf('								TIME VARYING (TV) PHI VERSION AS IN HLW \n'); sep(140);
	%************************************************************************************************************
	% NOW MAKE LOWER CI INTERVALS (Need to devide by sample size T)
	low_TV_bl		=	structfun(@(x) (x/eout_S2_bl.T), eout_S2_bl.Lambda_low_CI	, 'UniformOutput', false);
	low_TV_g		=	structfun(@(x) (x/eout_S2_bl.T), eout_S2_g.Lambda_low_CI	, 'UniformOutput', false);
	low_TV_M0g	= structfun(@(x) (x/eout_S2_bl.T), eout_S2_M0g.Lambda_low_CI, 'UniformOutput', false);

	% UPPER CI INTERVALS									
	up_TV_bl	=	structfun(@(x) (x/eout_S2_bl.T), eout_S2_bl.Lambda_up_CI	, 'UniformOutput', false);
	up_TV_g		=	structfun(@(x) (x/eout_S2_bl.T), eout_S2_g.Lambda_up_CI		, 'UniformOutput', false);
	up_TV_M0g	= structfun(@(x) (x/eout_S2_bl.T), eout_S2_M0g.Lambda_up_CI , 'UniformOutput', false);

	% make struture form R files output
	L2z_RFile	= struct('L',NaN,'MW',S2MUEs(1,1),'EW',S2MUEs(2,1),'QLR',S2MUEs(3,1));
	printstructs(	L2z_RFile , ...
								L2z_bl    ,	low_TV_bl , up_TV_bl , ...
								L2z_g     , low_TV_g  , up_TV_g  , ...
								L2z_M0g   , low_TV_M0g, up_TV_M0g, '12.6f');						

	sep(140);  fprintf('		CORRESPONDING (TV) STRUCTURAL BREAK F-STATISTICS\n'); sep(140); 
	bstat2_TV_bl	= eout_S2_bl.stats_ols;
	bstat2_TV_g		= eout_S2_g.stats_ols;
	bstat2_TV_M0g	= eout_S2_M0g.stats_ols;
	
	% get the p-values on the F-stats from eout_S2_bl.
	pvals_TV_bl		= eout_S2_bl.Fstat_pvals;
	pvals_TV_g		= eout_S2_g.Fstat_pvals;
	pvals_TV_M0g	= eout_S2_M0g.Fstat_pvals;
	pvals_RFile		= struct('L',NaN,'MW',NaN,'EW',NaN,'QLR',NaN);
	
	% make struture form R files output
	bstat2_RFile = struct('L',NaN,'MW',S2MUEs(1,2),'EW',S2MUEs(2,2),'QLR',S2MUEs(3,2));
	printstructs(	bstat2_RFile,		pvals_RFile	, ...
								bstat2_TV_bl,		pvals_TV_bl	, ...
								bstat2_TV_g,		pvals_TV_g 	,	...
								bstat2_TV_M0g,	pvals_TV_M0g , 6);							
	sep(140);fprintf('\n')
	
	% get the p-values on the F-stats from eout_S2_bl.
	pvals_TV_bl		= eout_S2_bl.Fstat_pvals;
	pvals_TV_g		= eout_S2_g.Fstat_pvals;
	pvals_TV_M0g	= eout_S2_M0g.Fstat_pvals;

	%************************************************************************************************************
	sep(140); fprintf('			CONSTANT (C) PHI (GY COMPUTED ONCE OUTSIDE THE LOOP)\n'); sep(140);
	%************************************************************************************************************
	% GET CIS / P-VALUES FOR CONSTANT PHI VERSION (NOTE: MAKE_MUE_LAMBDA_FROM_BREAK_STATS DIVIDES BY T)
	% --------------------------------------------------------------------------------------------------------------
	[~,pvals_C_bl , low_C_bl , up_C_bl ] = make_MUE_Lambda_from_break_stats(eout_S2_bl.stats_chow , eout_S2_bl.T );
	[~,pvals_C_g  , low_C_g  , up_C_g  ] = make_MUE_Lambda_from_break_stats(eout_S2_g.stats_chow	, eout_S2_g.T  );
	[~,pvals_C_M0g, low_C_M0g, up_C_M0g] = make_MUE_Lambda_from_break_stats(eout_S2_M0g.stats_chow, eout_S2_M0g.T);
	% --------------------------------------------------------------------------------------------------------------

	% PRINT NICE MUE2 OUTPUT IF NEEDED
	printstructs(	Chow2_bl,	low_C_bl , up_C_bl	, ...
								Chow2_g,	low_C_g  , up_C_g		,	...
								Chow2_M0g,low_C_M0g, up_C_M0g , '12.6f');

	sep(140); fprintf('		CORRESPONDING (C) STRUCTURAL BREAK F-STATISTICS\n'); 
	bstat2_C_bl		= eout_S2_bl.stats_chow;
	bstat2_C_g		= eout_S2_g.stats_chow;
	bstat2_C_M0g	= eout_S2_M0g.stats_chow;

	sep(140)					
	printstructs(	bstat2_C_bl , pvals_C_bl	, ...
								bstat2_C_g  ,	pvals_C_g		, ...
								bstat2_C_M0g, pvals_C_M0g , 8);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%																			STAGE 3        
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%		Theta3 vec.:	[a_y1,a_y2,a_r,b_pi,b_y, s_y~,s_pi,s_y*]
	%									[   1,   2,  3,   4,  5,    6,   7,  8 ] 
	% ---------------------------------------------------------------------------------------------------
	% make data input file
	data_in		= [Y X]';

	% ADD THE CORRECT LAMBDA_Z FROM CORRECT STAGE 2 MUE: -> choose which Lambda to use
	% ---------------------------------------------------------------------------------------------------
	other_pars_S3_M0 = other_pars_S3;
	other_pars_S3_M0.Lambda_z = Chow2_M0g.EW;		% constant Phi version (EW as used by HLW)
	% other_pars_S3_M0.Lambda_z = L2z_M0g.EW;		% time-varying Phi version (EW as used by HLW)
	% other_pars_S3_M0.Lambda_z = 0;						% set to zero if needed for testing.

	% STAGE 3: SET BOUNDS ON THE PARAMETER SPACE IF NEEDED
	% ---------------------------------------------------------------------------------------------------
	% In HLW, (b_y) b2.constraint = 0.025, (a_r) a3.constraint = -0.0025
	% ---------------------------------------------------------------------------------------------------
	S3k	= size(initVals3,1);
	LB	= struct; UB	= struct; AA	= struct; 
	
	LB.bl	= -Inf*ones(S3k,1);		LB.bl(5)	=	 0.0250;	% (b_y) b2.constraint >  0.0250
	UB.bl	=  Inf*ones(S3k,1);		UB.bl(3)	= -0.0025;	% (a_r) a3.constraint < -0.0025
	AA.bl	= [A_ zeros(2,length(LB.bl)-2)];
	% TO FORCE POSITIVE SIGMA (not needed, leads to marginally worse numerical results)
	% LB.bl(end-2:end)	= 0;													% [sigma_y~ sigma_pi sigma_y*] > 0 
	
	%	MLE SIGMA_G 
	LB.g 	= -Inf*ones(S3k+1,1);  LB.g(5)	=	 0.0250;  % (b_y) b2.constraint >  0.0250
	UB.g 	=  Inf*ones(S3k+1,1);  UB.g(3)	= -0.0025;  % (a_r) a3.constraint < -0.0025
	AA.g 	=	[A_ zeros(2,length(LB.g)-2)]; 
	% TO FORCE POSITIVE SIGMA (not needed, leads to marginally worse numerical results)
	% LB.g(end-3:end)		= 0;													% [sigma_y~ sigma_pi sigma_y* sigma_g] > 0
	
	% MLE SIGMA_G,SIGMA_Z
	LB.gz = -Inf*ones(S3k+2,1); LB.gz(5)	=	 0.0250;	% (b_y) b2.constraint >  0.0250
	UB.gz =  Inf*ones(S3k+2,1); UB.gz(3)	= -0.0025;	% (a_r) a3.constraint < -0.0025
	AA.gz	=	[A_ zeros(2,length(LB.gz)-2)];
	% TO FORCE POSITIVE SIGMA (not needed, leads to marginally worse numerical results)
	% LB.gz(end-4:end)	= 0;													% [sigma_y~ sigma_pi sigma_y* sigma_g sigma_z] > 0

	% allocate structure space 
	LL3			= struct();		gradOpt3 = struct();	
	theta3	= struct();		hessOpt3 = struct();

	% STAGE 3 MLE ESTIMATION
	% ---------------------------------------------------------------------------------------------------
	% NOTE: I use the initVals for the baseline model, and then use the baseline models estimates
	% as initial values in the estimations that follow for the paramters that are the same.
	% ---------------------------------------------------------------------------------------------------
	% LL at initival values
	LL3.initvals	= LogLike_Stage3_HLW_SSF(initVals3, data_in, other_pars_S3);
	% unconstraint estimates
	% using CSMINWELL 
% 	[theta3.CS, LL3.CS, retCode,~, gradOpt3.CS, hessOpt3.CS] = csminwel('LogLike_Stage3_HLW_SSF', initVals3, eye(length(initVals3)), [], 1e-8, 1e5, data_in,	other_pars_S3);
	% using FMINUNC 
	[theta3.O,	LL3.O,~,~, gradOpt3.O, hessOpt3.O]	= fminunc(@LogLike_Stage3_HLW_SSF,	initVals3, options_unc, data_in,	other_pars_S3);

	% BASELINE MODEL
	[theta3.bl, LL3.bl]	= fmincon(@LogLike_Stage3_HLW_SSF,					theta3.O,					AA.bl,	bb	,[],[],LB.bl,	UB.bl ,[],	options_con, data_in, other_pars_S3);
	% \SIGMA_G FREELY ESTIMATED BY MLE WITH \Lambda_z(HLW)
	[theta3.g,	LL3.g]	= fmincon(@LogLike_Stage3_HLW_SSF_sigma_g,	[theta3.bl; .03],	AA.g,		bb	,[],[],LB.g,	UB.g	,[],	options_con, data_in, other_pars_S3);
	% \SIGMA_G FREELY ESTIMATED BY MLE WITH CORRECT STAGE2 \Lambda_z(M0)
	[theta3.M0g,LL3.M0g]= fmincon(@LogLike_Stage3_HLW_SSF_sigma_g,	[theta3.bl; .03],	AA.g,		bb	,[],[],LB.g,	UB.g	,[],	options_con, data_in, other_pars_S3_M0);
	% \SIGMA_G AND \SIGMA_Z FREELY ESTIMATED BY MLE 
	[theta3.gz, LL3.gz]	= fmincon(@LogLike_Stage3_HLW_SSF_sigma_gz,	[theta3.g;	.01],	AA.gz,	bb	,[],[],LB.gz,	UB.gz ,[],	options_con, data_in, other_pars_S3);
	% CSMINWELL on .gz model (without constraints on parameters)
% [theta3.gzCS, LL3.gzCS, retCodegz,~, gradOpt3.gzCS, hessOpt3.gzCS] = csminwel('LogLike_Stage3_HLW_SSF_sigma_gz', theta3.gz, 0.001*eye(length(theta3.gz)), [], 1e-8, 1e5, data_in,	other_pars_S3);
	% UNCONSTRAINED FMINUNC on .gz model (without constraints on parameters) for extra numerical precision for very small changes in LogLikelihood function
	[theta3.gzO,	LL3.gzO,~,~, gradOpt3.gzO, hessOpt3.gzO] = fminunc(@LogLike_Stage3_HLW_SSF_sigma_gz, theta3.gz, options_unc, data_in,	other_pars_S3);

	%% GROUP FOR PRINTING TO SCREEN
	% take absolute vaues of standard deviations to have postive values, -> does not matter for estimation
	theta3.bl(end-2:end)	= abs(theta3.bl (end-2:end));
	theta3.g(end-3:end)		= abs(theta3.g  (end-3:end));
	theta3.M0g(end-3:end) = abs(theta3.M0g(end-3:end));
	theta3.gz(end-4:end)	= abs(theta3.gz (end-4:end));
	% ---------------------------------------------------------------------------------------------------
	% HLW_initvals	= [HLW_initvals_R_File; nan(2,1); -LL0_HLW_R_File; nan(2,1)];
	S3_initVals		= [initVals3; NaN(2,1); -LL3.initvals; NaN(2,1)];
	S3_HLW_Rfile	= [HLW_Stage3_R_File(1:11,2); other_pars_S3.Lambda_g; other_pars_S3.Lambda_z ];	
	S3_replicated	= [theta3.bl; S3g(theta3.bl,other_pars_S3.Lambda_g);S3z(theta3.bl,other_pars_S3.Lambda_z); -LL3.bl ; other_pars_S3.Lambda_g; other_pars_S3.Lambda_z];
	S3_MLE_g			= [theta3.g; S3z(theta3.g,other_pars_S3.Lambda_z); -LL3.g; L3g(theta3.g); other_pars_S3.Lambda_z];
	S3_MLE_M0g		= [theta3.M0g; S3z(theta3.M0g,other_pars_S3_M0.Lambda_z); -LL3.M0g; L3g(theta3.M0g); other_pars_S3_M0.Lambda_z];
	S3_MLE_gz			= [theta3.gz; -LL3.gz; L3g(theta3.gz); L3z(theta3.gz)];
	S3_MLE_gzO		= [theta3.gzO; -LL3.gzO; L3g(theta3.gzO); L3z(theta3.gzO)];
	
	% PRINT RESULTS TO SCREEN / THIS MAKES A MATLABE TABLE OBJECT
% 	Cols2Print = 1:5;
	tabS3 = table(par_names3, S3_initVals, S3_HLW_Rfile, S3_replicated, S3_MLE_g, S3_MLE_M0g, S3_MLE_gz);
%   tabS3 = table(par_names3, S3_HLW_Rfile, S3_replicated, S3_MLE_M0g, S3_MLE_gz);
	sep(120)
	fprintf('																Stage 3 Results           \n')
	sep(120)
% 	print2screen(tabS3([1:10 12:13 11],[1 3:end]), '%16.10f')
	print2screen(tabS3, '%16.10f'); sep(120)

	%	GET KFS OUTPUT, MAKE HLW FACTORS [R*, trend growth g,  z, and cycle], STORE IN STRUCTURE HLW
	% ---------------------------------------------------------------------------------------------------
	% allocate structure space
	HLW = struct();
	% BASELINE REPLICATION
	[~,S3_KFS.db]		= LogLike_Stage3_HLW_SSF(				theta3.bl,		data_in, other_pars_S3, 1);
	HLW.bl					= make_HLW_factors(S3_KFS.db); 

	% SIGMA_G ESTIMATED BY MLE \LAMBDA_Z FROM HLW S2 MODEL
	[~,S3_KFS.g]		= LogLike_Stage3_HLW_SSF_sigma_g(theta3.g,		data_in, other_pars_S3, 1);
	HLW.g						= make_HLW_factors(S3_KFS.g);

	% \SIGMA_G FREELY ESTIMATED BY MLE WITH CORRECT STAGE2 \Lambda_z(M0)
	[~,S3_KFS.M0g]	= LogLike_Stage3_HLW_SSF_sigma_g(theta3.M0g,	data_in, other_pars_S3_M0, 1);
	HLW.M0g					= make_HLW_factors(S3_KFS.M0g);

	% SIGMA_G and \SIGMA_Z ESTIMATED BY MLE
	[~,S3_KFS.gz]		= LogLike_Stage3_HLW_SSF_sigma_gz(theta3.gz,	data_in, other_pars_S3, 1);
	HLW.gz					= make_HLW_factors(S3_KFS.gz);

	% % SIGMA_G ESTIMATED BY MLE \LAMBDA_Z FROM CORRECT S2 MODEL (M0)
	% [~,S3_sigma_g_zM0]	= LogLike_Stage3_HLW_SSF_sigma_g(theta_sig_g_zM0, data_in, other_paras_zM0, 1);
	% HLW_sigma_g_zM0			= make_HLW_factors(S3_sigma_g_zM0);

	% GET HLW FACTORS FROM R-CODE FOR PLOTTING
	HLW.Rsmoothed = dlmread(getS(['Stage3.Smoothed.' hlw_read_sample]),',',1,1);
	HLW.Rfiltered = dlmread(getS(['Stage3.Filtered.' hlw_read_sample]),',',1,1);
	% NOW

	%% NOW DO THE PLOTTING % ----------------------------------------------------------------------------
	% FILTERED OR SMOOTHED FOR PLOTTING
	FLTRD = 1;		% set to 1 for filtered, smoothed otherwise
	% DEFINE THE TIME RANGE (TR) FOR THE TIMETABLE VARIABLES IN THE PLOTS STARTING IN 1961:Q1
	% TR = timerange(datestr(data.Time(5),'qq-yyyy'), SMPL_END,'quarters');
	TR = timerange(datestr(data.Time(5),'qq-yyyy'), SMPL_END,'closed');
	% ---------------------------------------------------------------------------------------------------------
	% LOAD IN THE EXCEL FILE BASED HLW ESTIMATES REPORTED ON WEB AT: 'HTTPS://WWW.NEWYORKFED.ORG/MEDIALIBRARY/
	% ---------------------------------------------------------------------------------------------------------
	% load real gdp data and other trend growth estimates. 
	load(['../data/' 'US_trend_growth_1947Q1-2019Q4.mat'])
	USdat = US_trend_growth;
	% media/research/economists/williams/data/Holston_Laubach_Williams_current_estimates.xlsx';
	load([HLW_INPUT_DATA_DIR 'hlw_current_FRBNY.mat'])
	% ---------------------------------------------------------------------------------------------------------	
	% recession indicator
	RecI_clr = .9*ones(3,1); 
	
	% ---------------------------------------------------------------------------------------------------------	
	% SOME PLOTTING COMMANDS
	% ---------------------------------------------------------------------------------------------------------	
	set(groot,'defaultLineLineWidth',2.5); 
	% PLOTTING DIMENSIONS FOR SUB-FIGURES ---------------------------------------------------------------
	p.fs = 15;			% FONT SIZE                                       
	p.df = 21;			% DATE FREQUENCY                                  
	p.aj = -1.21;	  % SUBTITLE POSITION ADJUSTMENT                  
	% p.dm = @(x)([.8-(x-1)*.215 .86 .165]);% FIG DIMENSION      
	p.dm = @(x)([.785-(x-1)*.235 .86 .185]);% FIG DIMENSION      
	% call as: setplot(p.dm(2),p.fs,1)
	% ---------------------------------------------------------------------------------------------------
	subtitles_ = {'(a) Natural rate ($r_t^\ast$)';				% % plot_names = {'(a) Natural real rate $r^{\ast}_t=4g_t+z_t$';
								'(b) Trend growth ($g_t$)';							% % 							'(b) Trend growth $4g_t$ (annualized)';       
								'(c) Other factor ($z_t$)';							% % 							'(c) Other factor $z_t$'};                    
								'(d) Output gap ($\tilde{y}_t$)'};			% % 							'(d) Cycle $\tilde{y}_t$'};                    

	legNames_= {'HLW R-Files';
							'Replicated';
							'MLE($\sigma_g|\lambda_z^{\mathcal{HLW}})$';
							'MLE($\sigma_g|\hat\lambda_z^{\mathrm{Correct}})$';
							'MLE($\sigma_g,\sigma_z)$'};
	
	legNames = legNames_([1 4 5]);

	% SOME CONTROLS
	legI	= [];		% Legend index
	numS	= 4;		% number of subplots to show in figure
	hlw_xls = hlw_current.(cntr_{CI});
	
	%% PLOTTING 
	if PLOT_ON
		f1 = figure(1); clf; set(f1,'WindowState','maximized','Position',[1441 641 1200 1843]);
		for jj = 1:numS
		% ---------------------------------------------------------------------------------------------------
			subplot(numS,1,jj); 
		%	THIS ADDS US RECESSION INDICATORS/COMMENT OUT IF NOT NEEDED
		if CI == 1
			% US Recessions: exclude for the other countries if not needed/wanted 
				if (jj == 1 | jj == 2)
					bar( USdat{TR,{'RecI'}}*6.15,1,'FaceColor',RecI_clr,'EdgeColor',RecI_clr,'ShowBaseLine','off'); hold on;
				else
					bar( USdat{TR,{'RecI'}}*6,1,'FaceColor',RecI_clr,'EdgeColor',RecI_clr,'ShowBaseLine','off'); hold on;
				end
				
				if (jj == 4 )
					bar(-USdat{TR,{'RecI'}}*9.5,1,'FaceColor',RecI_clr,'EdgeColor',RecI_clr,'ShowBaseLine','on' );
				else
					bar(-USdat{TR,{'RecI'}}*2.0,1,'FaceColor',RecI_clr,'EdgeColor',RecI_clr,'ShowBaseLine','on' );
				end
		end

			if FLTRD == 1
				legI(1) =	plot(HLW.Rfiltered(:,jj),	'Color', clr(3)); hold on
				legI(2) =	plot(HLW.M0g.att(:,jj)	,	'Color', clr(2), 'LineStyle','-'); 				
				legI(3) =	plot(HLW.gz.att(:,jj)		,	'Color', clr(1), 'LineStyle','--'); 				
			else
				legI(1) =	plot(HLW.Rsmoothed(:,jj),	'Color', clr(3)); hold on;
				legI(2) =	plot(HLW.M0g.atT(:,jj)	,	'Color', clr(2), 'LineStyle','-'); 				
				legI(3) =	plot(HLW.gz.atT(:,jj)		,	'Color', clr(1), 'LineStyle','--'); 				
			end	
			
			hold off; 
			grid on; box on;
			setplot(p.dm(jj),p.fs,1)
			setdateticks(Dates, p.df)
 			if jj == 1; setyticklabels( -0	: 1	: 6	, 0); hline(0); ylim([-0.25 6.2]);end
 			if jj == 2; setyticklabels( -0	: 1	: 6	, 0);						ylim([-0		6.2]);end
			if jj == 3; setyticklabels(-1.5	:.5 : 1	, 1); hline(0); ylim([-1.8  1.2]);end
			if jj == 4; setyticklabels( -8	: 2	: 6	, 0); hline(0); ylim([-9.5	6	 ]);end
			set(gca,'GridLineStyle',':','GridAlpha',1/3);
			setoutsideTicks; add2yaxislabel

% 			if (jj == 3 | jj == 3); addlegend(legI, legNames); end % addlegend(aG, legNames); end;
			if (jj == 2); addlegend(legI, legNames); end % addlegend(aG, legNames); end;
			subtitle(subtitles_{jj},p.aj,p.fs-0,1);
		end

		diff_name = 'HLW_prior';
		if FLTRD == 1; KFS_name	= 'filtered'; else; KFS_name = 'smoothed';	end

		% % print graphics
		figout_name1 = [ COUNTRY '_estimates_' KFS_name '_' diff_name '_' SMPL_END ];
 		if PRINT_PLOT_TO_PDF == 1; print2pdf(figout_name1, latex_graphics_dir); end
	end

	diary off 
	toc

  %% PRINT ALL THE 'FACTORS' TO EXCEL FILE WITH DATES
  if PRINT_FACTORS_EXCEL == 1
  	warning('off','MATLAB:xlswrite:AddSheet');
  	factor_names_short = {'r(t)*','g(t)','z(t)','y~(t)'};
  	factors_xls_name = ['..\matlab.output\' char(cntr_(CI)) '_factors.xls'];
  	PWD_ = pwd;
  
  	%%%  FILTERED
  	% dates
  	xlswrite(factors_xls_name, cellstr('Dates')									, 'Filtered', 'A2')
  	xlswrite(factors_xls_name, cellstr(datestr(Dates,'YYYY-QQ')), 'Filtered', 'A3')
  
  	% filtered factors HLW-R Files
  	xlswrite(factors_xls_name, cellstr('HLW-R Files')			,				'Filtered', 'B1')
  	xlswrite(factors_xls_name, cellstr(factor_names_short), 			'Filtered', 'B2')
  	xlswrite(factors_xls_name, HLW.Rfiltered							, 			'Filtered', 'B3') 

  	% filtered factors MLE.M0(g) MLE(?(g)|?(z)^Correct))
  	xlswrite(factors_xls_name, cellstr('MLE(sigma(g)|lambda(z)^Correct))'),'Filtered', 'G1')
  	xlswrite(factors_xls_name, cellstr(factor_names_short), 			'Filtered', 'G2')
  	xlswrite(factors_xls_name, HLW.M0g.att								, 			'Filtered', 'G3') 
  
  	% filtered factors MLE(g)
  	xlswrite(factors_xls_name, cellstr('MLE.(g,z)')				, 			'Filtered', 'L1')
  	xlswrite(factors_xls_name, cellstr(factor_names_short), 			'Filtered', 'L2')
  	xlswrite(factors_xls_name, HLW.gz.att									, 			'Filtered', 'L3') 

% THESE ARE NOT REPORTED IN THE PAPER  
%   	% filtered factors Baseline replication
%   	xlswrite(factors_xls_name, cellstr('Baseline Replicated'),		'Filtered', 'G1')
%   	xlswrite(factors_xls_name, cellstr(factor_names_short), 			'Filtered', 'G2')
%   	xlswrite(factors_xls_name, HLW.bl.att									, 			'Filtered', 'G3') 
  
%   	% filtered factors MLE(g)
%   	xlswrite(factors_xls_name, cellstr('MLE(g|Lambda_g^HLW)'), 		'Filtered', 'L1')
%   	xlswrite(factors_xls_name, cellstr(factor_names_short), 			'Filtered', 'L2')
%   	xlswrite(factors_xls_name, HLW.g.att									, 			'Filtered', 'L3') 
  
  	%%%  SMOOTHED
  	% dates
  	xlswrite(factors_xls_name, cellstr('Dates')									, 'Smoothed', 'A2')
  	xlswrite(factors_xls_name, cellstr(datestr(Dates,'YYYY-QQ')), 'Smoothed', 'A3')
  
  	% smoothed factors HLW-R Files
  	xlswrite(factors_xls_name, cellstr('HLW-R Files')			, 			'Smoothed', 'B1')
  	xlswrite(factors_xls_name, cellstr(factor_names_short), 			'Smoothed', 'B2')
  	xlswrite(factors_xls_name, HLW.Rsmoothed							, 			'Smoothed', 'B3') 
  
  	% smoothed factors MLE.M0(g)
  	xlswrite(factors_xls_name, cellstr('MLE(sigma(g)|lambda(z)^Correct))'),'Smoothed', 'G1')
  	xlswrite(factors_xls_name, cellstr(factor_names_short), 			'Smoothed', 'G2')
  	xlswrite(factors_xls_name, HLW.M0g.atT								, 			'Smoothed', 'G3') 
  
  	% smoothed factors MLE(g)
  	xlswrite(factors_xls_name, cellstr('MLE.(g,z)')				, 			'Smoothed', 'L1')
  	xlswrite(factors_xls_name, cellstr(factor_names_short), 			'Smoothed', 'L2')
  	xlswrite(factors_xls_name, HLW.gz.atT									, 			'Smoothed', 'L3') 

% THESE ARE NOT REPORTED IN THE PAPER
%   	% smoothed factors Baseline replication
%   	xlswrite(factors_xls_name, cellstr('Baseline Replicated'),		'Smoothed', 'Q1')
%   	xlswrite(factors_xls_name, cellstr(factor_names_short), 			'Smoothed', 'Q2')
%   	xlswrite(factors_xls_name, HLW.bl.atT									, 			'Smoothed', 'Q3') 
%   
%   	% smoothed factors MLE(g)
%   	xlswrite(factors_xls_name, cellstr('MLE(g|Lambda_g^HLW)'), 		'Smoothed', 'V1')
%   	xlswrite(factors_xls_name, cellstr(factor_names_short), 			'Smoothed', 'V2')
%   	xlswrite(factors_xls_name, HLW.g.atT									, 			'Smoothed', 'V3') 

  
  	disp([' Factors printed to excel file ' factors_xls_name])
  	
  	% delete the Sheet1 excel sheet, needs absolute path input
  	delete_excel_sheets( [PWD_(1:end-12)  factors_xls_name(3:end)],'Sheet1');
  end

  % PRINT RESULTS TO LATEX FILES TO BE READ-IN							
	if PRINT_RESULTS_LATEX == 1
		% ---------------------------------------------------------------------------------------------------
		% S2 MLE results
		% ---------------------------------------------------------------------------------------------------
		rowNames = {'$\hsp[3]a_{y,1}  $              ','$\hsp[3]a_{y,2}  $              ','$\hsp[3]a_{r}    $              ','$\hsp[3]a_{0}    $              ','$\hsp[3]a_{g}    $              ','$\hsp[3]b_{\pi } $              ','$\hsp[3]b_{y}    $              ','$\hsp[3]\sigma _{\tilde{y}}$    ','$\hsp[3]\sigma _{\pi }     $    ','$\hsp[3]\sigma _{y^{\ast }}$    ','$\hsp[3]\sigma _{g}$ {(implied)}','$\hsp[3]\lambda_g  $ {(implied)}'};				
% 		TS2		= tabS2{:, [3:end]}
		TS2		= tabS2{:, [3 5 6 7]};
		TS2a	= TS2;
		TS2a(end-1,:) = [];
		TS2b	= TS2(end-1,:);

 		LT_S21a = latexmat(rowNames, [TS2a], '{\hsp[7]---}', 'nomath', '% 4.8f');
 		fid	= fopen([latex_table_dir 'Table_S2a_' COUNTRY '.tex'],'wt'); fprintf(fid, '%s\n', LT_S21a); fclose(fid);
% 
 		LT_S21b = latexmat({'{Log-likelihood}'}, [TS2b], [], 'nomath', '% 4.8f');
 		fid	= fopen([latex_table_dir 'Table_S2b_' COUNTRY '.tex'],'wt'); fprintf(fid, '%s\n', LT_S21b); fclose(fid);

		% ---------------------------------------------------------------------------------------------------
		% now S2 MUE results
		% ---------------------------------------------------------------------------------------------------
		rowNames = {'$L$' ,'MW ' ,'EW ' ,'QLR'};				
		% Lambda z
		Lz_TV = [	struct2array(L2z_RFile    )' , ...
							struct2array(L2z_g        )' , struct2array(low_TV_g)'	, struct2array(up_TV_g)'		, ...
							struct2array(L2z_M0g      )' , struct2array(low_TV_M0g)', struct2array(up_TV_M0g)'	];

		Lz_C  = [	struct2array(Chow2_bl     )' , struct2array(low_C_bl)'	, struct2array(up_C_bl)'		, ...
							struct2array(Chow2_g      )' , struct2array(low_C_g)'		, struct2array(up_C_g)'			, ...
							struct2array(Chow2_M0g    )' , struct2array(low_C_M0g)'	, struct2array(up_C_M0g)'		];

 		LT_1a = latexmat(rowNames, [Lz_TV inf(4,1) Lz_C],'{\hsp[3]---}', 'nomath', '% 4.6f');
 		fid		= fopen([latex_table_dir 'Table_MUEa_' COUNTRY '.tex'],'wt'); fprintf(fid, '%s\n', LT_1a); fclose(fid);

		% Fstatistics
		F_TV  = [	struct2array(bstat2_RFile  )' , ...
							struct2array(bstat2_TV_g   )' , struct2array(pvals_TV_g		)' , ...
							struct2array(bstat2_TV_M0g )' , struct2array(pvals_TV_M0g )' ];
						

		F_C  = [	struct2array(bstat2_C_bl  )'  , struct2array(pvals_C_bl	)' , ...
							struct2array(bstat2_C_g   )'  , struct2array(pvals_C_g		)' , ...
							struct2array(bstat2_C_M0g )'  , struct2array(pvals_C_M0g )' ];

		LT_1b = latexmat(rowNames, [F_TV inf(4,1) F_C],'{\hsp[3]---}', 'nomath', '% 4.6f');
 		fid		= fopen([latex_table_dir 'Table_MUEb_' COUNTRY '.tex'],'wt'); fprintf(fid, '%s\n', LT_1b); fclose(fid);

		% ---------------------------------------------------------------------------------------------------
		% S3 MLE results
		% ---------------------------------------------------------------------------------------------------
		rowNames = {'$\hsp[3]a_{y,1}  $               ','$\hsp[3]a_{y,2}  $               ','$\hsp[3]a_{r}    $               ','$\hsp[3]b_{\pi } $               ','$\hsp[3]b_{y}    $               ','$\hsp[3]\sigma _{\tilde{y}}$     ','$\hsp[3]\sigma _{\pi }     $     ','$\hsp[3]\sigma _{y^{\ast }}$     ','$\hsp[3]\sigma _{g}$ {(implied)} ','$\hsp[3]\sigma _{z}$ {(implied)} ','$\hsp[3]\lambda_g  $ {(implied)} ','$\hsp[3]\lambda_z  $ {(implied)} '};

% 		TS3		= tabS3{:, 3:end};
		TS3		= tabS3{:, [3 end-1:end]};
		TS3a	= TS3;
		TS3a(end-2,:) = [];
		TS3b	= TS3(end-2,:);

		LT_S31a = latexmat(rowNames, [TS3a], '{\hsp[8]---}', 'nomath', '% 4.8f');
		fid	= fopen([latex_table_dir 'Table_S3a_' COUNTRY '.tex'],'wt'); fprintf(fid, '%s\n', LT_S31a); fclose(fid);

		LT_S31b = latexmat({'{Log-likelihood}'}, [TS3b], [], 'nomath', '% 4.8f');
		fid	= fopen([latex_table_dir 'Table_S3b_' COUNTRY '.tex'],'wt'); fprintf(fid, '%s\n', LT_S31b); fclose(fid);

		disp('  --- >   Results written to latex tables to table.input')
	end

	%% PLOT THE SEQUENCE OF FSTATISTICS
	if PLOT_F_STATS == 1
		f2 = figure(2); clf; set(f2,'WindowState','maximized','Position',[1441 641 1200 1843]);
		set(groot,'defaultLineLineWidth',1.75); 
 		% PLOTTING DIMENSIONS FOR SUB-FIGURES ---------------------------------------------------------------
 		p.fs = 15;			% FONT SIZE                                       
 		p.df = 21;			% DATE FREQUENCY                                  
 		p.aj = -1.225;	% SUBTITLE POSITION ADJUSTMENT                  
 		p.dm = @(x)([.785-(x-1)*.255 .86 .195]);% FIG DIMENSION      
 		% call as: setplot(p.dm(2),p.fs,1)
 		% ---------------------------------------------------------------------------------------------------

		subplot(2,1,1)
		sw = 3;
		hold on;
			plot(addnans( eout_S2_bl.Fstat_ols		,3,4 ),'-'	,'Color',clr(2)); 
			plot(addnans( eout_S2_bl.Fstat_chow		,3,4 ),'-'	,'Color',clr(1)); 
			hold off;
		% hline(0)
		box on; grid on;
		setplot(p.dm(2),p.fs,1)
		setdateticks(Dates, p.df,[],p.fs)
 		setyticklabels(-0		:2	:12	, 0); ylim([0 12.5]);
		set(gca,'GridLineStyle',':','GridAlpha',1/3)
		setoutsideTicks
		add2yaxislabel
% 		subtitle('(a) Misspecified: $GY_t = a_{y}(L)\tilde{y}_{t}-a_{r}(L)r_{t}-a_{g}g_{t-1}-a_{0}$', p.aj, p.fs+1, 1)
		subtitle('(a) HLW''s (misspecified) Stage 2 model (right column block of equation 6)', p.aj, p.fs+1, 0)    
							% subtitle('Demeand series', -1.22, 14, 1)					
% 		legendflex({'Time varying \bf\phi';
% 								'Constant \bf\phi'}, ...
% 								'Interpreter','Tex','anchor', 3*ones(1,2))
 		legendflex({'HLW';
								'SW'}, ...
								'anchor', 3*ones(1,2))

		subplot(2,1,2)
		hold on;
			plot(addnans( eout_S2_M0g.Fstat_ols		,3,4 ),'-','Color',	clr(2)); 
			plot(addnans( eout_S2_M0g.Fstat_chow	,3,4 ),'-','Color',	clr(1)); 
		hold off;
		box on; grid on;
		setplot(p.dm(3),p.fs,1)
		setdateticks(Dates, p.df,[],p.fs)
 		setyticklabels(-0		:2	:12	, 0); ylim([0 12.5]);
		set(gca,'GridLineStyle',':','GridAlpha',1/3)
		setoutsideTicks
		add2yaxislabel
% 		subtitle('(b) Correctly specified: $GY_{t}=a_{y}(L)\tilde{y}_{t}-a_{r}(L)[r_{t}-4g_{t}]$', p.aj, p.fs+1, 1)
		subtitle('(b) Correct Stage 2 model (left column block in equation 6)', p.aj, p.fs+1)    

% 		legendflex({'Time varying \bf\phi';
% 								'Constant \bf\phi'}, ...
% 								'Interpreter','Tex','anchor', 3*ones(1,2))

 		legendflex({'HLW';
								'SW'}, ...
								'anchor', 3*ones(1,2))

		fstat_fig_outname = [ COUNTRY '_Fstat_MUE2_' SMPL_END ];
		if PRINT_PLOT_TO_PDF == 1; print2pdf(fstat_fig_outname, latex_graphics_dir); end
%     print2pdf(fstat_fig_outname, latex_graphics_dir)
	end
	
	%% PLOT GDP GROWTH AND REAL RATE
	if PLOT_GDP_RR == 1
		f3 = figure(3); clf; set(f3,'WindowState','maximized','Position',[1441 641 1200 1843]);
	% 	set(get(gcf),'WindowState','maximized','Position',[1441 641 1200*ff 1843*ff]);
		%	make GDP growth (annualized) from 1961:Q1 onwards --> drop first 4 periods
		gdp_growth	= drop(delta(data.gdp_log,1)*400, 4);
		real_rate		= drop(data.real_rate, 4);
		gdp_rr = gdp_growth - real_rate;
		zt = HLW.Rfiltered(:,3);

		subplot(2,1,1)
		if CI == 1
		% US Recessions: exclude for the other countries if not needed/wanted 
			bar( USdat{TR,{'RecI'}}*16 ,1,'FaceColor',RecI_clr,'EdgeColor',RecI_clr,'ShowBaseLine','off'); hold on;
			bar(-USdat{TR,{'RecI'}}*10 ,1,'FaceColor',RecI_clr,'EdgeColor',RecI_clr,'ShowBaseLine','off'); 
		end
		hold on;
			Lg(1) = plot(gdp_growth, 'Color',.0*ones(1,3),'Marker','.','MarkerSize',15,'LineStyle',':','LineWidth',3/2);
			Lg(2) = plot(real_rate , 'Color', clr(2), 'LineStyle','-','LineWidth',2);

			setplot(p.dm(1),p.fs,0)
			setdateticks(Dates,p.df);
			box on; grid on;
			set(gca,'GridLineStyle',':','GridAlpha',1/3)
			hline(0)
			% ADJUST THE YTICK LABELS 
			setyticklabels(-8:4:16,0);	ylim([-10 16])
			if CI == 2; setyticklabels(-12:4:8,0); ylim([-14 10]); end
			if CI == 3; setyticklabels(-12:4:20,0); ylim([-12 21]); end
		hold off;	
		% MAKE OUTSIDE TICKS FOR TWO AXIS Y LABELS 
		set(gca,'GridLineStyle',':','GridAlpha',1/3)
		setoutsideTicks
		add2yaxislabel
		addlegend(Lg,{'GDP growth','Real rate'}, [], p.fs-3)
		subtitle('(a) Time series evolution of real GDP growth and real interest rate', p.aj+.02, p.fs+1, 1)

		subplot(2,1,2)
		if CI == 1
		% US Recessions: exclude for the other countries if not needed/wanted 
			bar( USdat{TR,{'RecI'}}*16 ,1,'FaceColor'	,RecI_clr,'EdgeColor',RecI_clr,'ShowBaseLine','off'); hold on;
			bar(-USdat{TR,{'RecI'}}*16 ,1,'FaceColor'	,RecI_clr,'EdgeColor',RecI_clr,'ShowBaseLine','off'); 
		end
		hold on;
			plot(gdp_rr)
	%  		plot(-zt)
			setplot(p.dm(2),p.fs,0)
			setdateticks(Dates,p.df);
			box on; grid on;
			set(gca,'GridLineStyle',':','GridAlpha',1/3)
			hline(0)
			% ADJUST THE YTICK LABELS 
			setyticklabels(-16:4:16,0);% 	
			if CI == 2; setyticklabels(-12:4:8,0);	ylim([-14 9]); end
			if CI == 3; setyticklabels(-16:4:20,0); ylim([-16 20]); end
			% 	ylim([-15 16])
		hold off;	
		% MAKE OUTSIDE TICKS FOR TWO AXIS Y LABELS 
		set(gca,'GridLineStyle',':','GridAlpha',1/3)
		setoutsideTicks
		add2yaxislabel
		subtitle('(b) Time series evolution (real GDP growth -- real interest rate)', p.aj+.02, p.fs+1, 1)

		gdp_rr_fig_outname = [ COUNTRY '_GDP_RR_TS_PLOT_' SMPL_END ];
		if PRINT_PLOT_TO_PDF == 1; print2pdf(gdp_rr_fig_outname, latex_graphics_dir); end

	% PLOT ACF OF gdp-rr
		f4 = figure(4); clf; set(f4,'WindowState','maximized','Position',[1441 641 1200 1843]);
		plotacf(gdp_rr);
		acfplot_fig_outname = [ COUNTRY '_GDP_RR_ACF_PLOT_' SMPL_END ];
		if PRINT_PLOT_TO_PDF == 1; print2pdf(acfplot_fig_outname, latex_graphics_dir); end
	end
	
	
	
	
end % end of country loop

%% ESTIMATE AN AR model
% AR_names = @(AR_order) (cellstr([repmat('AR(' ,AR_order,1) num2str([1:AR_order]') repmat(')',AR_order,1)]));
% AR_order	= 2;
% AR_out		= ols(gdp_rr,mlag(gdp_rr,AR_order),0,AR_names(AR_order)');
% AR_mu			= AR_out.bhat(1)/(1-sum(AR_out.bhat(2:end)));
% disp(AR_mu)


































% EOF end of file 
	% ---------------------------------------------------------------------------------------------
	% 	figure('WindowState','maximized','Position',[1441 641 1200 1843]);
	% 	figure('WindowState','maximized','Position',[1 641 1200 1843]);
	%		figure('WindowState','maximized','Position',[-1441 641 1200 1843]);
	% ---------------------------------------------------------------------------------------------


