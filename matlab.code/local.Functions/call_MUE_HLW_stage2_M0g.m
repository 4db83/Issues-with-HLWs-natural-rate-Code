function [LAMBDA_Z_db, Lambda_chow, extra_output] = call_MUE_HLW_stage2_M0g(theta, S2object, YY, XX)
% this is a modified version of the call_MUE_HLW_stage2.m file because we now need to construct
% GY(t) = y~(t)-a_1y~(t)-a_1y~(t) - ar[ (r(t-1)+r(t-2))/2 - 4*(g(t-1)+g(t-2))/2 ];
% get the theta parameter input vector with right dimension
theta4_MUE2 = theta(1:3); % get the order right [ay1,ay2,ar]

% MAKE THE Y AND X INPUT VARIABLES AS IN HLW(2017)
% YY	  = [lnGDP INFL]; 
% [GDP(t-1) GDP(t-2) r(t-1) r(t-2) pi(t-1) mean(pi(t-2:t-4))]
% XX	= [ GPD_2_lags RR_2_lags INFL_4_lags(:,1) mean(INFL_4_lags(:,2:4),2) ];

YX_db = make_yx_4_stage2_M0g(YY, XX, S2object);

yy = YX_db.tT(:,1);			% y~(t) (T x 1) (y~(t) = cycle variable)
xx = YX_db.tT(:,2:end); % [y~(t-1) y~(t-2) rr_2-gg_2] (T x 3)

% GY = yy - xx*theta4_MUE2;
% MUE_HLW_Stage1(GY)

% CALL TO MUE PROCEDURE TO COMPUTE LAMBDA_Z
[LAMBDA_Z_db, Lambda_chow, Stats_ols, Stats_chow, extra_output]	= MUE_HLW_Stage2_M0g(yy, xx, theta4_MUE2);

% adding y and x as well as parameter input vector to extra_output
extra_output.yy = yy;
extra_output.xx = xx;
extra_output.YX = YX_db;
extra_output.stats_ols		= Stats_ols;
extra_output.stats_chow		= Stats_chow;
extra_output.theta4_MUE2	= theta4_MUE2;

