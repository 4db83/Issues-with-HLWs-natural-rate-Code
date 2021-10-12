function [LAMBDA_Z_db, Lambda_chow, extra_output] = call_MUE_HLW_stage2(theta, S2object, YY, XX)
% GET THE THETA PARAMETER INPUT VECTOR WITH RIGHT ORDER AND DIMENSION.
% re-order to [ay1,ay2,ar,ag,a0], since theta_hat is orderd as: [ay1,ay2,ar,a0,ag]
theta4_MUE2 = theta([1:3 5 4]); 

% MAKE THE Y AND X INPUT VARIABLES AS IN HLW(2017)
YX_db = make_yx_4_stage2(YY, XX, S2object);

yy = YX_db.tT(:,1);			% y~(t) (T x 1) (y~(t) = cycle variable)
xx = YX_db.tT(:,2:end); % [y~(t-1) y~(t-2) rr_2-gg_2] (T x 3)

% CALL TO MUE PROCEDURE TO COMPUTE LAMBDA_Z
[LAMBDA_Z_db, Lambda_chow, Stats_ols, Stats_chow, extra_output]	= MUE_HLW_Stage2(yy, xx, theta4_MUE2);

% adding y and x as well as parameter input vector to extra_output
extra_output.yy = yy;
extra_output.xx = xx;
extra_output.YX = YX_db;
extra_output.stats_ols		= Stats_ols;
extra_output.stats_chow		= Stats_chow;
extra_output.theta4_MUE2	= theta4_MUE2;

