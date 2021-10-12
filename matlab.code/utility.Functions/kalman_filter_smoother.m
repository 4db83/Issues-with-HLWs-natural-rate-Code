function KFS = kalman_filter_smoother(y, varargin) 
% ******************************************************************************************************
% FUNCTION: General Kalman Filter and Smoother routine (FAST) 
% ******************************************************************************************************
% NOTE: y input has to be (Ky x T) --> Time loop is over columns for speed (not over rows)
% ------------------------------------------------------------------------------------------------------
% CALL AS: 
% 		KFS = kalman_filter_smoother(y, D, M, H, C, Phi, Q, S, a00, P00, DK)  
% OR SIMPLEY AS: 
% 		KFS = kalman_filter_smoother(y, PARS_IN) 
% WHERE:
% 		PARS_IN = structure containing the relevant D, M, H, C, Phi, Q, S, a00, P00, DK parameters.  
% ------------------------------------------------------------------------------------------------------
% KALMAN FILTER for multivariate states (alphas) general notation as in HARVEY 1992 and/or DK Book
% STATE SPACE MODEL:
% 
% 		Observed:	y(t)		 = D(t) + M*alpha(t)			+ e(t);		Var[e(t)] = H.
% 		State:		alpha(t) = C(t) + Phi*alpha(t-1)	+ S*n(t);	Var[n(t)] = Q.
% 
% where [e(t); n(t)]  ~ MNorm([0; 0],[H 0; 0 Q] and 
% S a selection matrix (rather than R as in Harvey).
% ------------------------------------------------------------------------------------------------------
% INPUT: 
% ------------------------------------------------------------------------------------------------------
% Observed/measurement equation:
% 	y		  = (Ky x T) vector of observed data, where Ky = dim(y(t)) and T is the time dimension.
% 	D		  = (Ky x T) or (Ky x 1) matrix of exogenous variables and/or intercepts,
% 	M		  = (Ky x Ka) matrix of state to observed mapping parameters, Ka = dim(alpha(t)).
% 	H		  = (Ky x Ky) Variance/Covariance matrix.
% State equation:
% 	C			= (Ka x T) or (Ka x 1) matrix of exogenous variables and/or intercepts in state vector alpha.
% 	Phi		= (Ka x Ka) matrix of state to observed mapping parameters.
% 	S			= Selection matrix to make up the right dimension for Variance/covariance matrix. I use
% 					S (for selection) instead of R as in the textbook treatment of Harvey and others.
% 	Q			= (Ka x Ka) Variance/Covariance matrix of state vector alpha(t).
% ------------------------------------------------------------------------------------------------------
% OUTPUT:
% ------------------------------------------------------------------------------------------------------
% KFS        = a structure with Filtered and Smoothed states, as well as correspdoning variances and 
%              log-likelihood functions.
% KFS.att 	 = filtered states
% KFS.atT 	 = smoothed states
% KFS.Ptt 	 = filtered states Var/Cov
% KFS.LL_sum = sum Log-Likelihood
% KFS.LL_t 	 = full Log-Likelihood
% KFS.Ft 		 = (Ky x Ky x T) Variance of the prediction error
% KFS.uhat	 = (T x 1) Prediction error
% ------------------------------------------------------------------------------------------------------
% INITIALIZATION:
% ------------------------------------------------------------------------------------------------------
% a00 and P00 are the initial values of the states and its Variance/Covariance matrix. 
% If states are  % stationary, can initize at uncondtional mean of alpha_t, that is, AR(1) unconditonal 
% mean of inv(I-Phi)*C and its variance inv(I-PhixPhi)vec(R*QR') or vec(Q).
% ------------------------------------------------------------------------------------------------------
% SOME NOTES:
% ------------------------------------------------------------------------------------------------------
% Durbin-Koopman use a(t) = a(t|t-1) and P(t) = a(t|t-1), but for the initialization
% they define a(1) = E[alpha(1)] and P(1) = Var[alpha(1)], ie., unconditionally. For the first time period, 
% ie., at t=1 iteration the following are used:
% 													a(t|t_1) = a1; 
% 													P(t|t_1) = P1;
% 
% In Koop and Korobilis (and most other treatments, see the DMA/DMS code), there is an extra updating 
% step, so that the following are used at time t=1:
% 													a(t|t_1) = C + Phi*a(1);
% 													P(t|t_1) = Phi*P1*Phi' + RQR; 
% 
% THIS IS JUST AN INITIALIZATION, thus not important per se for the Filter as it dissipates quickly, but 
% good to keep in mind when comparing output to be exactly the same. Ultimately, it is a question of 
% where the initialization happens, ie., either t=1 (DK) or at t=0, so that for t=1 we use already the 
% predicted att_1 values (C + Phi*a1). 
% ------------------------------------------------------------------------------------------------------
% REFERENCES: 
% ------------------------------------------------------------------------------------------------------
% I use a mix of notation from, mainly from Harvey's books, but I replace Z in the state vector with M.
% In the state vector, I use Phi instead of T for the state dynamics and S instead of R for the selection
% matrix in case of identities.
% The Smoothing recursions are most cleanly stated in Kim and Nelson (1999, page 37), so I use those.
% ------------------------------------------------------------------------------------------------------
% db (06.03.2018)
% UPDATED from previous Kalmanfilter routine in the ECB folder.
% db (09.03.2018)
% db (05.04.2018) cleaned up
% db (12.06.2018) made faster
% ******************************************************************************************************

% SetDefaultValue(11,'DK',0);

% passing a structure argument for easier KFS output generation from given model
if nargin == 2
	% get the fieldnames first
	fnames = fieldnames(varargin{1});
	I_DK = sum(strcmpi(fnames,'DK'));
	
	D   = varargin{1}.D  ;
	M   = varargin{1}.M  ;
	H   = varargin{1}.H  ;
	C   = varargin{1}.C  ;
	Phi = varargin{1}.Phi;
	Q   = varargin{1}.Q  ;
	S   = varargin{1}.S  ;
	a00 = varargin{1}.a00;
	P00 = varargin{1}.P00;
	% if DK does not exist, set it to 0 PER DEFAULT
	if ~sum(I_DK)
		DK = 0;
	else
		DK = varargin{1}.DK;
	end
else
	D   = varargin{1};
	M   = varargin{2};
	H   = varargin{3};
	C   = varargin{4};	
	Phi = varargin{5};
	Q   = varargin{6};
	S   = varargin{7};
	a00 = varargin{8};
	P00 = varargin{9};
	if length(varargin) < 10
	% if DK does not exist, set it to 0 PER DEFAULT
		DK = 0;
	else
		DK	= varargin{10};
	end
end

% PRECOMPUTE log(2*pi)
log_2_pi = log(2*pi);

% DIMENSION OF y(t) DIM(y(t)) AND SAMPLE SIZE T
[Ky, T] = size(y);
% CHECK INPUT SIZE OF OBSERVABLE VECTOR Y WHICH SHOULD BE (Ky x T)
if Ky > T
	error('Y has to be a (Ky x T) matrix, to loop over the columns');
end
% DIMENSION OF STATE VECTOR ALPHA dim(alpha(t))
Ka = size(Phi,1);		

% Get the dimension of the 'constants' C(t)
if size(C,2) == 1; C = repmat(C,1,T); end
% and D(t) and make them of the papropriate size
if size(D,2) == 1; D = repmat(D,1,T); end

% % SET S TO IDENTITY IF NOT SUPPLIED.
% if nargin < 8; S = eye(Ka); end;

% PRE-COMPUTE S*Q*S'
SQS = S*Q*S';

% ------------------------------------------------------------------------------------------------------	
% KALMAN FILTER RECURSIONS (HARVEY 1989, STARTS WITH PREDICTION STEP, SEE CHAPTER 4)
% ------------------------------------------------------------------------------------------------------	
% INITIALISE THE KALMAN FILTER RECURSIONS 
% ------------------------------------------------------------------------------------------------------	
% Koop & Korobilis start at a00 and P00 (-> initialisation at t=0) and then predict a(t|t-1), and 
% continue from there. this is also what the ssm toolbox in matlab uses.
% ------------------------------------------------------------------------------------------------------	
att_1 = C(:,1) + Phi*a00;				
Ptt_1 = Phi*P00*Phi' + SQS;
% ------------------------------------------------------------------------------------------------------	
% Durbin & Koopman initialise at a(1)=a(1|0) and P(1)=P(1|0) (-> init at t=1) (2012, see p.32 and p.124) 
% ------------------------------------------------------------------------------------------------------	
if DK
	att_1 = a00;	
	Ptt_1 = P00;
end

% INITIALIZE LOG-LIKELIHOOD SUM AT 0.
LL_sum = 0;
LLt	= zeros(1,T);
I		= eye(Ka); % Identiy matrix

% SPACE ALLOCATION FOR THINGS THAT ARE STORED FOR LATER ON [Ka = length(alpha), Ky = length(y)]
at	= zeros(Ka,T);		  uhat_store	= zeros(Ky,T);		
Pt	= zeros(Ka,Ka,T);		Ft_store		= zeros(Ky,Ky,T);		
Gt	= zeros(Ka,Ky,T);
PHI	= zeros(Ka,Ka,T);			% storage 1 for part of filter weights
PSI	= zeros(Ka,Ka,T);			% storage 2 for part of filter weights
Gybar	= zeros(Ka,T);			% G*ybar = G*(y - Ax)
ybar	= (y - D);					% (y - Ax) without the state vector.
Psi0	= I;								% initialisation

% MAIN FORWARD LOOP (LOOPING OVER COLUMNS IN MATLAB IS FASTER THAN OVER ROWS)
for t = 1:T
	% FORECAST ERROR AND ITS MSE
	uhat	= y(:,t) - M*att_1 - D(:,t);
	Ft = M*Ptt_1*M' + H;
	% PRECOMPUTE PM_invF = P(t|t-1)*M'inv(F(t))
%  	P_MT_invF = Ptt_1*M'/Ft;
	P_MT_invF = Ptt_1*M'*pinv(Ft);
	% UPDATING/FILTERING a(t|t) and P(t|t)	
	att		= att_1 + P_MT_invF*uhat;
	Ptt		= Ptt_1 - P_MT_invF*M*Ptt_1;	
	% PREDICTION/FORECASTING a(t|t-1) and P(t|t-1)
	att_1	= Phi*att + C(:,t);
	Ptt_1	= Phi*Ptt*Phi' + SQS;
	% ----------------------------------------------------------------------------------------------------		
	% STORE THE FILTERED STATES/MSEs as well AS Ft and uhat (prediction error and its variance for LogLike)
	% ----------------------------------------------------------------------------------------------------	
 	uhat_store(:,t)	= uhat;		
 	Ft_store(:,:,t)	= Ft;
	at(:,t)	  = att;
	Pt(:,:,t) = Ptt;
	
	% now compute and Store the G(t), BigPhi(t), and PSI(t) quantaties for the filter weights
	Gt(:,:,t) = P_MT_invF;
	PHI(:,:,t) = (I-Gt(:,:,t)*M)*Phi;
	% I need reverse ordering of matrices, unlike cumprod, so do not use it, make in filter loop
	Psi0 = PHI(:,:,t)*Psi0;
	PSI(:,:,t) = Psi0;
	Gybar(:,t) = Gt(:,:,t)*ybar(:,t);

	% ----------------------------------------------------------------------------------------------------	
	% LOG-LIKELIHOOD 
	% ----------------------------------------------------------------------------------------------------	
	LL_t		= -0.5*( Ky*log_2_pi + log(det(Ft)) + uhat'/Ft*uhat );
	LL_sum	= LL_sum + LL_t;
	LLt(1,t)= LL_t;
end

% ALLOCATE SPACE for backward recursions AND initialize at a(T|T) and P(T|T)
atT = at;  PtT = Pt;

% BACKWARDS RECURSION ON a(t|T) and P(t|T) see page 37 in Kim and Nelson (1999) [at Pt]
for n = T-1:-1:1
  % FIRST PRE-COMPUTE SOME QUANTITIES
  % P(t+1|t) = Phi*P(t|t)Phi' + SQS
  Pt1_t = Phi*Pt(:,:,n)*Phi' + SQS;
  % KG(t) = P(t|t)*Phi'inv(P(t+1|t))
% 	KGt = Pt(:,:,n)*Phi'/Pt1_t; 
	% HARVEY SUGGESTS TO USE GENERALIZED INVERSE HERE PINV
	KGt = Pt(:,:,n)*Phi'*pinv(Pt1_t);
  atT(:,n)    = at(:,n)   + KGt*( atT(:,n+1) - C(:,t) - Phi*at(:,n) );
  PtT(:,:,n)  = Pt(:,:,n) + KGt*( PtT(:,:,n+1) - Pt1_t )*KGt'; 
end

% NOW RETURN BACK TO FUNCTION AS STRUCTURE
KFS.att		= at'; 		      % Filtered states
KFS.atT		= atT';					% Smoothed states
KFS.Ptt		= Pt;			      % Filtered states Var/Cov
KFS.PtT		= PtT;			    % Smoothed states Var/Cov
KFS.uhat	= uhat_store';	% (Ky x T) Prediction error
KFS.Ft		= Ft_store;			% (Ky x Ky x T) Variance of the prediction error
KFS.LL_t	= LLt';					% Full Log-Likelihood
KFS.LLsum = LL_sum;				% Sum Log-Likelihood
KFS.Y			= y';						% retunr the input observed data as well
KFS.a00		= a00;					% as well as the priors
KFS.P00		= P00;					% 
KFS.Gt		= Gt;						% matrix infront of uhat
KFS.ybar	= ybar;					% (y-Ax)
KFS.PHI		= PHI;					% big PHI Weights.
KFS.PSI		= PSI;					% PSI weights


% % RETURN LOG-LIKELIHOOD ONLY
% logLikeout = LL_sum;






































% EOF