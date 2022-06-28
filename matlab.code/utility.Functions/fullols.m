function olsout = fullols(y,x,no_const,L,INCLUDE_PW,Recession_indicator)
% F: fit beta_hat and se(beta_hat) only
% Call as:		olsout = fullols(y,x)
% or 
% FULL CALL:	olsout = fullols(y,x,no_const,L,INCLUDE_PW,Recession_indicator)
% Input: 
%		y 
%		x 
%		no_constant:	set to 1 to exclude constant from regression, default is 0.
%		L:						Pre-whitening lag
% ------------------------------------------------------------------------------------------------------
% Print results in nice format with print_fullols(olsout);
% ------------------------------------------------------------------------------------------------------
% NOTE: I have changed this arround. now constant indicator is third entry.
% previous calls were:
% function olsout = fullols(y,x,L,INCLUDE_PW,Recession_indicator,no_const)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE: ORDER OF VARIABLS TO BE CALL HAS CHANGED: DB JUNE 2018
% ------------------------------------------------------------------------------------------------------

SetDefaultValue(3,'no_const',0)			% default is to include a constant
SetDefaultValue(4,'L',[]);
SetDefaultValue(5,'INCLUDE_PW',0);	% set to 1 to include pre-whitening in LRV computation
SetDefaultValue(6,'Recession_indicator',[]);

% remove any missing values
Inan	= anynans([y x]);
y = y(~Inan,:); 
x = x(~Inan,:);

[T,K0] = size(x);				

% check if first column is unit vector for constant
x1 = x(:,1);
Ic = (sum(x1==1)==T); % Ic == 1 if constant is included in x

if no_const==0
% check if constant already in, if not, add a constant
	if ~Ic
		x = [ones(T,1) x];
	end
else
	if Ic
		disp(' You are trying to exclude a constant in the regression model, but x includes a vector of ones')
		disp(' The constant (vector of ones) was excluded from the regression model')
		x = x(:,2:end);
	end
end

% NEWEY AND WEST 1994 Rule of Thumb truncation lag L. for BARTLET
L_NW_rot	= ceil(4*(T/100)^(2/9));		% use ceil instead of floor to be more conservative

if isempty(L)
	L = L_NW_rot;
end

[T,K] = size(x);				

xpx 	= x'*x;
invxpx= xpx\eye(K);
beta	= invxpx*(x'*y);
% make the hat matrix X*inv(X'X)*X'(y) = yhat = X*bhat
hatmatrix = x*invxpx*x';
% beta  = x\y; is much much slower!!!
u			= y-x*beta;
sig2	= u'*u/(T-K);							% this is the OLS and not MLE variance
VCV		= invxpx*sig2;						% standard OLS VCV NOT HAC!
se		= sqrt(diag(VCV));				% standard OLS SE NOT HAC!
R2		= 1 - sum(u.^2)/sum((y-mean(y)).^2);
SSE		= u'*u;
SSY		= T*var(y,1);
u_		= u(2:end)-u(1:end-1);
DW		= u_'*u_/SSE;
% Rbar_squared or adjusted R2
R2adj = 1 - (T-1)/(T-K)*(1-R2);
% (concentrated) log-likelihood function
LL		= -T/2*( log(2*pi) + log(SSE/T) + 1 );

%% compute recession/expasion based R2
if ~isempty(Recession_indicator)
	[R2_exp,R2_rec] = Rsquared_rec(u,y,Recession_indicator);
	
% 	y_	= demean(y);
% 	% compute the 1-sum(I.*uhat.^2)/sum(I.*(r-rbar).^2) ratio
% 	R2_exp	= 1 - mean(u(US_rec==0).^2)/mean(y_(US_rec==0).^2);
% 	R2_rec	= 1 - mean(u(US_rec==1).^2)/mean(y_(US_rec==1).^2);
end

olsout.bhat 	= beta;     	% fitted regression parameters.
olsout.uhat 	= u;        	% fitted residuals
olsout.sse  	= SSE;				% Sum of Squared errors
olsout.se   	= se;       	% NON-HAC standard error of beta_hat
olsout.tstat	= beta./se; 	% tstat of beta_hat
% two sided p-value thus times 2! 
olsout.pval		= (1-normcdf(abs(beta./se)))*2; 
% Both Eviews and Stata also report two sided p-values. To get one sided p-value, devide by 2.
olsout.N			= T;					% sample size
olsout.k			= K;					% number of variables plus constant.	
olsout.xpx		= xpx;				% x'x
olsout.y			= y;					% return y
olsout.ssy		= SSY;				% return Sum of Squared Y
olsout.sig2		= sig2;				% variance of Regression (error)
olsout.yhat		= x*beta;			% fitted values
olsout.DW			= DW;					% Durbin Watson stat.
olsout.R2			= R2;					% R-squared.
olsout.R2_bar	= R2adj;			% Rbar-squared or adjusted R-squared.
olsout.LL			= LL;					% log-likelihood function
olsout.DF			= T-K;				% degrees of freeedom, ie Sample size minus constant minus number of regressors. N-k-1
														% here K = k+1, so only need to compute T-K.

if ~isempty(Recession_indicator)
	olsout.R2_exp	= R2_exp;		% R-squared expansion.
	olsout.R2_rec	= R2_rec;		% R-squared recession.
end
olsout.VCV			= VCV;			% Simple Variance/Covariance Matrix

% INFORMATION CRITERIA (log(L) = -(T/2)*log(SSE/T) for GAUSSIAN LOG LIKELIHOOD)
% ALL INFOCRITERIA ARE OF THE FORM -2log(L) + P, WHERE P IS A PENALTY
MSE_	= SSE/T;
IC_AIC 	= T*log(MSE_) + 2*K;		
IC_BIC	= T*log(MSE_) + K*log(T);
IC_AICc = IC_AIC			+ 2*K*(K+1)/(T-K-1);
IC_HQ		= T*log(MSE_) + 2*K*log(log(T));

% DEFLATE BY SAMPLE SIZE T (or N) TO BE CONSISTENT WITH STANDARD SOFTWARE OUTPUT.
olsout.aic	= IC_AIC/T;
olsout.aicc	= IC_AICc/T;
olsout.bic	= IC_BIC/T;
olsout.hq		= IC_HQ/T;

% HAC VCV and se baseed on simple Newy and West 1994 4*T/100^2/9 Rule of Thumb (ROT)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% moment condition x.*uhat = x.*(y-xbhat) at Optimum
Xu	= bsxfun(@times,x,u);

% % fit to each Xu an AR(1) and an ARMA(1,1) model TO BE DONE LATER


% add to output variable
olsout.Xu	= Xu;					% cross product of regressors times uhat.														

% Compute HAC VCV
S				= NW_var(Xu,L);				% L = is the truncation lag for HAC, if L == 0, we get White HC.
D_inv		= -T*invxpx;					
HAC_VCV = (D_inv'*S*D_inv)/T; % HAC VAR is 1/T*(D_inv*S*D_inv), D_inv = -X*X/T; S=Var(X*u).

% STORE THE RESULTS 
olsout.HAC_VCV		= HAC_VCV;
olsout.hac_se			= sqrt(diag(HAC_VCV));
olsout.hac_tstat	= beta./olsout.hac_se;	% standard error of beta_hat
olsout.hac_pval		= (1-normcdf(abs(beta./olsout.hac_se)))*2; % HAC two sided p-value

if INCLUDE_PW == 1
% matlabs HAC function with AR1OLS and ARMA11 prewhitening
	% NOTE: ALWAYS EXCLUDE THE INTERCEPT TERM => SET 'INTERCEPT',FALSE, AS IT IS ALREADY ADDED BY
	% DEFAULT TO X ABOVE, SEE LINES 17 TO 35.
	HAC_VCV_pw_AR1			=	HAC(x,y,'bandwidth','AR1OLS','weights','QS','whiten',1,'smallT',true,'intercept',false,'display','off');
	HAC_VCV_pw_ARMA11		= HAC(x,y,'bandwidth','ARMA11','weights','QS','whiten',1,'smallT',true,'intercept',false,'display','off');	
	% store these results
	olsout.HAC_VCV_pw_AR1		= HAC_VCV_pw_AR1;
	olsout.HAC_VCV_pw_ARMA11= HAC_VCV_pw_ARMA11;
end

% parse no constant info to outupt
olsout.no_const = no_const;	% equal to 1 if constant is surpressed, ie. not included.

% add some joint cofficient restrition option of the from R*bhat=r which will have distribution
q		= K0;								% number restrictions, starting from last variable.
RR	= [zeros(q,K-q) eye(q)];
rr	= zeros(q,1);
d		= (RR*beta-rr);
had_var_d	= (RR*HAC_VCV*RR');
var_d	= (RR*VCV*RR');
	
% stats are
hac_F_stat		= d'*inv(had_var_d)*d/q;
hac_Chi2_stat	= d'*inv(had_var_d)*d;
	
% p-values are
hac_Chi2_stat_pval	= chi2cdf(hac_Chi2_stat,q,'upper');
hac_F_stat_pval			= fcdf(hac_F_stat,q,T,'upper');

% NON-HAC
F_stat		= d'*inv(var_d)*d/q;
Chi2_stat	= d'*inv(var_d)*d;
	
% p-values are
Chi2_stat_pval	= chi2cdf(Chi2_stat,q,'upper');
F_stat_pval			= fcdf(F_stat,q,T,'upper');

olsout.hac_Fstat			= hac_F_stat;
olsout.hac_pval_Fstat	= hac_F_stat_pval;
olsout.Fstat					= F_stat;
olsout.pval_Fstat			= F_stat_pval;

olsout.hac_Chi2_stat			= hac_Chi2_stat;
olsout.hac_pval_Chi2_stat	= hac_Chi2_stat_pval;
olsout.Chi2_stat					= Chi2_stat;
olsout.pval_Chi2_stat			= Chi2_stat_pval;
olsout.L									= L;						% HAC truncation lag
olsout.INCLUDE_PW					= INCLUDE_PW;		% indicator if prewhitening is included
olsout.hatmatrix					= hatmatrix;
olsout.invxpxxp						= invxpx*(x');
olsout.x									= x;


% Fstat = (R*bhat-r)'inv([R*HAC_VCV*R'])(R*bhat-r)/q, q = No. of restrictions
% % RR = [zeros(4,K-4) eye(4)];
% % rr = zeros(4,1);
% % 
% % BB = (RR*beta-rr)*inv(RR*HAC_VCV*RR')*(RR*beta-rr)/4


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Newey West VAR/COVAR nwvar.m in db.toolbox is old way relying on JPLesage code. 
% This is much neater and faster
function S = NW_var(g,m)
% NWFn    Calculates covariance matrix of sqrt(T)*sample average.

T = size(g,1);						% g is Txq
m = min( [m;(T-1)] );				% number of lags, truncation lag L.
g	= g - repmat(mean(g),T,1);      % Normalizing to Eg=0 (NOT NEEDED REALLY), beacuse FOC should be 0

% Omega0 = g'*g/T;                  %(qxT)*(Txq)
Omega0 = cov(g);

S = Omega0;
for s = 1:m
  Omega_s = g(s+1:T,:)'*g(1:T-s,:)/T;   % same as Sum[g(t)*g(t-s)',t=s+1,T]
  S				= S  +  ( 1 - s/(m+1) ) * (Omega_s + Omega_s');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%