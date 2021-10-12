function YX_Stage2 = make_yx_4_stage2_M0g(YY, XX, Stage2_KFS_object)
% Function: makes the Y and X variables needed for HLWs Stage 2 MUE inputs. 
% ---------------------------------------------------------------------------------------------
% INPUT: - YY, XX data matrices from main file, and 
% 			 - Stage2_KFS_object obtained from running:
% [~,S2.C0_g]	= LogLike_Stage2_HLW_C0_g(theta_C0_g,	data_in_no_intercept, other_paras_g2lags, 1);   
% NOTE:	This is the SSF and KFS object of the correct version of the Stage 2 model that I list in 
%       equation (49, GYcorr) on page 30. that is: GY = ay(L)y~ - ar(L)[r(t) - 4g(t)] of the 
%       MLE(sigma.g)M0 model. SSF of the model that estimates sigma.g or not is the same.
% OUTPUT: 
% YX_Stage2 structure with tt and tT elements referring to fitlered and smoothed
% estimates.
% db
% --------------------------------------------------------------------------------------------- 

RR_2_lags = XX(:,3:4);
% T = length(Stage2_KFS_object.D);

% MAKE THE STAGE 2 Y AND X VARIABLES FROM SMOOTHED KFS OUTPUT OF THE CORRECT MODEL.
% NOTE: xi(t) = [y*(t) y*(t-1) y*(t-2) g(t-1) g(t-2)], 
% so take the first two row entries in y*(t-2) and stack on top of y*(t) 
ystar_tT	= [Stage2_KFS_object.KFS.atT(1:2,3); Stage2_KFS_object.KFS.atT(:,1)];

% smoothed output gap (or cycle) no make the cycle and then lag two times. This exactly reproduces HLW
cycle_tT  = YY(3:end,1)	-	ystar_tT;

% lag the cycle by two periods
cycle_2_lags = mlag(cycle_tT,2);

% Y = [ y~(t) ] output gap from 1960:Q3 (1961:Q1 is the start of the in-sample
YStage2 = cycle_tT(3:end);

% pre-compute [r(t-1)+r(t-2)]/2 and [4*g(t-1)+4*g(t-2)]/2
rr_2 = mean(RR_2_lags(5:end,:),2);									% [r(t-1)+r(t-2)]/2 		NO NEED TO DIVIDE AR BY 2 LATER.
gg_2 = mean(4*Stage2_KFS_object.KFS.atT(:,4:5),2);  % [4*g(t-1)+4*g(t-2)]/2	THIS IS ANNUALIZED,IE., 4*G(T).

% X = [ y~(t-1) y~(t-2) (rr_2-gg_2) ]
XStage2 = [	cycle_2_lags(3:end,:), ...      % y~(t-1) y~(t-2) 
            (rr_2-gg_2)	                ];	% rr_2-gg_2 ( [r(t-1)+r(t-2)]/2 - [4*g(t-1)+4*g(t-2)]/2 );

YX_Stage2.tT = [YStage2 XStage2];					


% MAKE THE STAGE 2 Y AND X VARIABLES FROM FILTERED KFS OUTPUT
% filtered potential or y*(t)
ystar_tt	= [Stage2_KFS_object.KFS.att(1:2,3); Stage2_KFS_object.KFS.att(:,1)];
% filtered output gap (or cycle)
cycle_tt  = YY(3:end,1)	-	ystar_tt;

% lag the cycle by two periods
cycle_2_lags = mlag(cycle_tt,2);

% Y = [ y~(t) ] output gap from 1960:Q3 (1961:Q1 is the start of the in-sample
YStage2_tt = cycle_tt(3:end);

% CHANGE gg_2 TO NOW USE FILTERED SERIES, rr_2 IS FROM DATA, SO NOT INFLUENCED BY THAT.
gg_2 = mean(4*Stage2_KFS_object.KFS.att(:,4:5),2);		% 0.5*( 4*g_(t-1)+4*g_(t-2) )	this is not annualized in any way.

% X = [ y~(t-1) y~(t-2) (rr_2-gg_2) ]
XStage2_tt = [ cycle_2_lags(3:end,:), ...     % y~(t-1) y~(t-2)                                         
               (rr_2-gg_2)                ];	% rr_2-gg_2 ( [r(t-1)+r(t-2)]/2 - [4*g(t-1)+4*g(t-2)]/2 );

YX_Stage2.tt = [YStage2_tt XStage2_tt];										