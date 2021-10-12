function [Lambdas_out, p_values, Lambda_low_CI, Lambda_up_CI, Lambda_alt] = make_MUE_Lambda_from_break_stats(Stats,T)
	% NOW GET/LOAD THE CRITICAL LAMBDA VALUES FROM THE COARSE GRID AS THE TABLE IN SW1998 AND AS USE IN HLW2017.
	Lambda_table = (0:30)';
	% L				MW				EW				QLR            % corresponding Lambda values from the SW GRID 
	CritVals = [                                                                                
	0.118	 	 0.689		 0.426	 	 3.198           %     0                                        
	0.127	 	 0.757		 0.476		 3.416           %     1                                        
	0.137	 	 0.806		 0.516		 3.594           %     2                                        
	0.169	 	 1.015		 0.661		 4.106           %     3                                        
	0.205	 	 1.234		 0.826		 4.848           %     4                                        
	0.266	 	 1.632		 1.111		 5.689           %     5                                        
	0.327	 	 2.018		 1.419		 6.682           %     6                                        
	0.387	 	 2.390		 1.762		 7.626           %     7                                        
	0.490	 	 3.081		 2.355		 9.160           %     8                                        
	0.593	 	 3.699		 2.910		10.660           %     9                                        
	0.670	 	 4.222		 3.413		11.841           %    10                                        
	0.768	 	 4.776		 3.868		13.098           %    11                                        
	0.908	 	 5.767		 4.925		15.451           %    12                                        
	1.036	 	 6.586		 5.684		17.094           %    13                                        
	1.214	 	 7.703		 6.670		19.423           %    14                                        
	1.360	 	 8.683		 7.690		21.682           %    15                                        
	1.471	 	 9.467		 8.477		23.342           %    16                                        
	1.576	 	10.101		 9.191		24.920           %    17                                        
	1.799	 	11.639		10.693		28.174           %    18                                        
	2.016	 	13.039		12.024		30.736           %    19                                        
	2.127	 	13.900		13.089		33.313           %    20                                        
	2.327	 	15.214		14.440		36.109           %    21                                        
	2.569	 	16.806		16.191		39.673           %    22                                        
	2.785	 	18.330		17.332		41.955           %    23                                        
	2.899	 	19.020		18.699		45.056           %    24                                        
	3.108	 	20.562		20.464		48.647           %    25                                        
	3.278	 	21.837		21.667		50.983           %    26                                        
	3.652	 	24.350		23.851		55.514           %    27                                        
	3.910	 	26.248		25.538		59.278           %    28                                        
	4.015	 	27.089		26.762		61.311           %    29                                        
	4.120	 	27.758		27.874		64.016           %    30                                        
	];
	
	%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	%					ALTERNATIVELY: LOAD THE FINER GRID FROM ORIGINAL SW FILES AS USED IN SW98
	% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	% 	SWdatadir = 'D:\matlab.tools\db.toolbox\SW98.inputs';
	% 	SWdatadir = '';
	SWdatadir = './local.Functions/';
	load([SWdatadir 'SW1998_MUE_lookup_table.mat']);
	load([SWdatadir 'SW1998_MUE_lookup_TSTCDF.mat']);
	% -----------------------------------------------------------------------------------------
	%			UNCOMMENT TO USE THE FINER GRIDED SW1998_ORIGINAL_MUE_lookup_table.MAT VALUES
	%			USING THE TABLE ABOVE EXACTLY REPLICATES THE INTERPOLATION OVER THE COARSER GRID.
	% -----------------------------------------------------------------------------------------
	% 	MUE_vals_tmp	= lookup_table.Variables;
	% 	Lambda_table	= MUE_vals_tmp(:,1);
	% 	CritVals			= MUE_vals_tmp(:,2:5); % order is ['L';'MW';'EW';'QLR']; % SW = QLR
	% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	%% Final calculations for clean output
	testStats = struct2array(Stats);
	% testNames = {'L';'MW';'EW';'QLR'}; % SW = QLR
	testNames = fieldnames(Stats);
	% fprintf(' Test is: %s \n', char(testNames(ii)))

	IL_MAX = length(Lambda_table);

	% loop through the 4 main tests
	for ii = 1:4
		test_stat = testStats(ii);
		crit_val	= CritVals(:,ii);
		% get indicator
		isGreater = test_stat > crit_val;
		IL1 = find(isGreater,	1, 'Last');
		IL0 = find(isGreater);
		IL	= max(IL1, max(IL0));

		if(isempty(IL))
			Lambda_hat(ii,1) = 0;
		elseif(IL==(IL_MAX)) % if we are at the last table entry, just return that without interpolation
			Lambda_hat(ii,1) = Lambda_table(IL_MAX) / T; % if IL is greater than that largest lambda value of 30, return 30/T
		else
			% INTERPOLATE TWO ADJACENT OBSERVATIONS AND DIVIDE BY T. THIS IS A DOUBLE
			Lambda_hat(ii,1) = interp1(crit_val(IL:IL+1), Lambda_table(IL:IL+1), test_stat) / T;
		end
%		% MAKE A STRUCTURE WITH THE FIELDS EQUAL TO THE TESTNAMES 
%		Lambdas_out.(testNames{ii}) = Lambda_hat(ii);
	end
	
	% THIS IS NEEEDED FOR STAGE 2 MUE ONLY FOR QLR STAT REQUIRING Y-XBETA TYPE OF COMPUTATION. 
	% if there are any nans, find them
	Inan = find(anynan(testStats'));
	% now set those to nan
	if ~isempty(Inan)
		Lambda_hat(Inan) = nan;
	end
	% add to the output structure
	for ii = 1:4	
			% now make a structure with the fields equal to the testNames 
		Lambdas_out.(testNames{ii}) = Lambda_hat(ii);
	end
	
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% the code below computes some additional output based on the lookup_table_mue.m function and the        %% 
  %% 'SW1998_ORIGINAL_MUE_lookup_table.mat']); and 'SW1998_ORIGINAL_MUE_lookup_TSTCDF.mat' critical values  %% 
  %% matfiles that were created from SW GAUSS code and used in the EXCACT replication of their results in		%%  
  %% the file SW1998_MUE_replication.m in the D:\_research\_current\ecb\HLW_matlab directory.               %%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
	% compute other lambdas from lookup_table_mue function
	Lambda_alt.L	 = lookup_table_mue( testStats(1) , lookup_table.Lambda, lookup_table.L  );
	Lambda_alt.MW	 = lookup_table_mue( testStats(2) , lookup_table.Lambda, lookup_table.MW );
	Lambda_alt.EW	 = lookup_table_mue( testStats(3) , lookup_table.Lambda, lookup_table.EW );
	Lambda_alt.QLR = lookup_table_mue( testStats(4) , lookup_table.Lambda, lookup_table.SW );
	
	% COMPUTE THE P-VALUES for the structural break stats
	% inline function to find p-value
	get_p_value = @(x) ( 1-TSTCDF.pdvec(x == min(x)) );
	p_values.L		=	get_p_value( abs( ( TSTCDF.L (1,:) - testStats(1) ) ) );  % L
	p_values.MW		= get_p_value( abs( ( TSTCDF.MW(1,:) - testStats(2) ) ) );  % MW
	p_values.EW		= get_p_value( abs( ( TSTCDF.EW(1,:) - testStats(3) ) ) );  % EW
	p_values.QLR	= get_p_value( abs( ( TSTCDF.SW(1,:) - testStats(4) ) ) );  % QLR

	% % COMPUTE THE LOWER CI BOUND ON LABMDA 90% CI
	CI 	= .90;
	TAR = (1-CI)/2;
	% find columns corresponding to this target 
	[~,Iup ] = min(abs( TSTCDF.pdvec - TAR ));
	[~,Ilow] = min(abs( TSTCDF.pdvec - 1-TAR ));

	% lower CI
	Lambda_low_CI.L		= lookup_table_mue( testStats(1) , lookup_table.Lambda, TSTCDF.L (:,Ilow) );	% L
	Lambda_low_CI.MW	= lookup_table_mue( testStats(2) , lookup_table.Lambda, TSTCDF.MW(:,Ilow) );	% MW
	Lambda_low_CI.EW	= lookup_table_mue( testStats(3) , lookup_table.Lambda, TSTCDF.EW(:,Ilow) );	% EW
	Lambda_low_CI.QLR	= lookup_table_mue( testStats(4) , lookup_table.Lambda, TSTCDF.SW(:,Ilow) );	% QLR
	% upper CI
	Lambda_up_CI.L		= lookup_table_mue( testStats(1) , lookup_table.Lambda, TSTCDF.L (:,Iup) );	% L
	Lambda_up_CI.MW		= lookup_table_mue( testStats(2) , lookup_table.Lambda, TSTCDF.MW(:,Iup) );	% MW
	Lambda_up_CI.EW		= lookup_table_mue( testStats(3) , lookup_table.Lambda, TSTCDF.EW(:,Iup) );	% EW
	Lambda_up_CI.QLR	= lookup_table_mue( testStats(4) , lookup_table.Lambda, TSTCDF.SW(:,Iup) );	% QLR
	
	% divide them by T to get the HLW \lambda_z = \lambda/T
	Lambda_low_CI = structfun(@(x) (x/T), Lambda_low_CI , 'UniformOutput', false);
	Lambda_up_CI	= structfun(@(x) (x/T), Lambda_up_CI  , 'UniformOutput', false);
	
end

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%  THESE are the smallest structural break statistics values that lead to a non-zero lower 90%CI
% -------------------------------------------------------------------------------------------------
%        Break-stat        Pvals			 Lam_low_CI     Lam_up_CI 
% -------------------------------------------------------------------------------------------------
% L      0.84570000     0.00500000     0.00005466    44.89719682 
% MW     5.26970000     0.00500000     0.00008120    45.60976097 
% EW     3.83020000     0.00500000     0.00006895    38.13258002 
% QLR   13.18150000     0.00500000     0.00007674    36.43378650 
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

























% EOF 

















































