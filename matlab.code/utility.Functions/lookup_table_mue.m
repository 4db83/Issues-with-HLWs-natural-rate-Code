function Lambda_hat = lookup_table_mue(test_stat, lookup_table_Lamda, lookup_table_column)
% FUNCTION: look-up table to get the lambda value for a given value of the test statistic. 
% call as: Lambda_hat = lookup_table_mue(statistic, lookup_table, whichTest)
% 
% where the inputs are:
%
%		test_stat			= the test statistic of interst
%		lookup_table	= the matlab Tabel with all the values for lambda, L, EW etc as computed by SW(1998)
%	  whichTest			= string, either L, MW, EW, SW for the 4 main tests. have a look at the SW(1998) 
%										what these tests arepaper. NOTE SW == QLR in table.
% ------------------------------------------------------------------------------------------------------

% lookup_table_Lamda = lookup_table;

% lookup_table_column 

% % switch whichTest
% % 	case 'L'
% % 		lookup_table_column = lookup_table.(2);
% % 	case 'MW'
% % 		lookup_table_column = lookup_table.(3);
% % 	case 'EW'
% % 		lookup_table_column = lookup_table.(4);
% % 	case 'SW'
% % 		lookup_table_column = lookup_table.(5);
% % 	otherwise
% % 			disp('No other choices')
% % end;
			
% now use the interpolation procedure suggested by SW to find Lambda MU
isGreater = test_stat > lookup_table_column;

% Note: the lookup table is not monotonic. 
% IL_1 finds the last value where test_stat > tab_statistic, while IL_2 gives the first value
% so we have to choose one to take and can take either, depending on what is wanted.

IL_1 = find(isGreater,	1, 'Last');
IL_2 = find(~isGreater,	1, 'First') -1;

ILL = [IL_1 IL_2];
% choose the one that gives the lower or larger lambda, ie the first time time or last time test_stat > tab_statistic
IL = min(ILL);

% disp([ILL IL])

if (isempty(IL) || IL==0)
	Lambda_hat = 0;
else
	if IL >= 201
		Lambda_hat = lookup_table_Lamda(end);
	else
		% by hand interpolation is Y(2) = Y(1) + (X(2)-X(1))*DY/DX = yhat = b0 + X*b1
% 		DX = diff(lookup_table_column(IL:IL+1));
% 		DY = diff(lookup_table_Lamda(IL:IL+1));
% 		Lambda_hat_byhand = lookup_table_Lamda(IL) + (test_stat-lookup_table_column(IL))*DY/DX;
		
		% NOW DO THE INTERPOLATION YNEW = INTERP1(X, Y, XNEW):
		% NOTE: if the values in the X variabe are not monotonically increasing, it will return an NA.
		Lambda_hat = interp1(lookup_table_column(IL:IL+1)	, lookup_table_Lamda(IL:IL+1) , test_stat);

% 		disp([Lambda_hat Lambda_hat_byhand])
	end
end








% ------------------------------------------------------------------------------------------------------
% now use the interpolation procedure suggested by SW to find Lambda MU
% IL	= find( Lstat				> lookup_table.L,		1, 'Last');
% Imw	= find( mean_Fstat	> lookup_table.MW,	1, 'Last');
% Iew	= find( exp_Wald		> lookup_table.EW,	1, 'Last');
% Isw	= find( sup_Fstat		> lookup_table.QLR, 1, 'Last');
% ------------------------------------------------------------------------------------------------------
% THIS IS WHAT IS WRITTEN IN TABLE 3 OF SW, BUT THIS IS NOT WHAT THEY DO IN THE GAUSS CODE
% ------------------------------------------------------------------------------------------------------
% Lambda_ew = lookup_table.Lamda(Iew) + (exp_Wald - lookup_table.EW(Iew))/...
% 						(lookup_table.EW(Iew+1)-lookup_table.EW(Iew))
% ------------------------------------------------------------------------------------------------------
% NOW DO THE INTERPOLATION YNEW = INTERP1(X, Y, XNEW):
% NOTE: in the paper what is written in the footnote is wrong, this is exctly matched to GAUSS code results
% ------------------------------------------------------------------------------------------------------
% if Lstat < lookup_table.L(1)
% 	disp(' ERROR L_stat < MUE LOOK-UP TABLE VALUE');
% end
% clc
% Lambdas.L	 	= interp1(lookup_table.L(IL:IL+1)		 , lookup_table.Lamda(IL:IL+1)  , Lstat);
% Lambdas.MW  = interp1(lookup_table.MW(Imw:Imw+1) , lookup_table.Lamda(Imw:Imw+1), mean_Fstat);
% Lambdas.EW  = interp1(lookup_table.EW(Iew:Iew+1) , lookup_table.Lamda(Iew:Iew+1), exp_Wald );
% Lambdas.SW  = interp1(lookup_table.QLR(Isw:Isw+1), lookup_table.Lamda(Isw:Isw+1), sup_Fstat );									
% ------------------------------------------------------------------------------------------------------





















% Imw	= find( mean_Fstat	> lookup_table.MW,	1, 'Last');
% Iew	= find( exp_Wald		> lookup_table.EW,	1, 'Last');
% Isw	= find( sup_Fstat		> lookup_table.QLR, 1, 'Last');
% ------------------------------------------------------------------------------------------------------
% THIS IS WHAT IS WRITTEN IN TABLE 3 OF SW, BUT THIS IS NOT WHAT THEY DO IN THE GAUSS CODE
% ------------------------------------------------------------------------------------------------------
% Lambda_ew = lookup_table.Lamda(Iew) + (exp_Wald - lookup_table.EW(Iew))/...
% 						(lookup_table.EW(Iew+1)-lookup_table.EW(Iew))
% ------------------------------------------------------------------------------------------------------
% NOW DO THE INTERPOLATION YNEW = INTERP1(X, Y, XNEW):
% NOTE: in the paper what is written in the footnote is wrong, this is exctly matched to GAUSS code results
% ------------------------------------------------------------------------------------------------------


% if statistic < lookup_table.L(1)
% 	disp(' ERROR L_stat < MUE LOOK-UP TABLE VALUE');
% end
% Lambdas.L	 	= interp1(lookup_table.L(IL:IL+1)		 , lookup_table.Lamda(IL:IL+1)  , statistic);
% Lambdas.MW  = interp1(lookup_table.MW(Imw:Imw+1) , lookup_table.Lamda(Imw:Imw+1), mean_Fstat);
% Lambdas.EW  = interp1(lookup_table.EW(Iew:Iew+1) , lookup_table.Lamda(Iew:Iew+1), exp_Wald );
% Lambdas.SW  = interp1(lookup_table.QLR(Isw:Isw+1), lookup_table.Lamda(Isw:Isw+1), sup_Fstat );									
