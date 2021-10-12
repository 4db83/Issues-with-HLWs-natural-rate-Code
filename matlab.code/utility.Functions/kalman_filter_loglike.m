function logLikeout = kalman_filter_loglike(y, D, M, H, C, Phi, Q, S, a00, P00, DK) 
% ******************************************************************************************************
% FUNCTION: General Kalman Filter and Smoother routine. 
% ******************************************************************************************************
% NOTE: y input has to be (Ky x T) --> Time loop is over columns for speed (not over rows)
% ------------------------------------------------------------------------------------------------------
% CALL AS: 
% 		logLikeout = kalman_filter_loglike(y, D, M, H, C, Phi, Q, S, a00, P00, DK) 
% 
% 
% WHERE:
% 		the relevant parameters are D, M, H, C, Phi, Q, S, a00, P00, DK parameters.  
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
% logLikeout = sum of the Log-Likelihood function L(t). Full (Tx1) can be 
%              can be computed by un-commenting line (200) and parsing it out.
% 
% 
% 
% 
% 
% 
% 
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

SetDefaultValue(11,'DK',0);

%% PASSING A STRUCTURE ARGUMENT FOR EASIER KFS OUTPUT GENERATION FROM GIVEN MODEL
% % if nargin == 2
% % 	% get the fieldnames first
% % 	fnames = fieldnames(varargin{1});
% % 	I_DK = sum(strcmpi(fnames,'DK'));
% % 	
% % 	D   = varargin{1}.D  ;
% % 	M   = varargin{1}.M  ;
% % 	H   = varargin{1}.H  ;
% % 	C   = varargin{1}.C  ;
% % 	Phi = varargin{1}.Phi;
% % 	Q   = varargin{1}.Q  ;
% % 	S   = varargin{1}.S  ;
% % 	a00 = varargin{1}.a00;
% % 	P00 = varargin{1}.P00;
% % 	% if DK does not exist, set it to 0 PER DEFAULT
% % 	if ~sum(I_DK)
% % 		DK = 0;
% % 	else
% % 		DK = varargin{1}.DK;
% % 	end
% % else
% % 	D   = varargin{1};
% % 	M   = varargin{2};
% % 	H   = varargin{3};
% % 	C   = varargin{4};	
% % 	Phi = varargin{5};
% % 	Q   = varargin{6};
% % 	S   = varargin{7};
% % 	a00 = varargin{8};
% % 	P00 = varargin{9};
% % 	if length(varargin) < 10
% % 	% if DK does not exist, set it to 0 PER DEFAULT
% % 		DK = 0;
% % 	else
% % 		DK	= varargin{10};
% % 	end;
% % end

%% PRECOMPUTE log(2*pi)
% log_2_pi = log(2*pi); % 
% Computed with R package library("Rmpfr"): ( Pi <- Const("pi", 1000 *log2(10)) )
% pi_ = 3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798214808651328230664709384460955058223172535940812848111745028410270193852110555964462294895493038196442881097566593344612847564823378678316527120190914564856692346034861045432664821339360726024914127372458700660631558817488152092096282925409171536436789259036001133053054882046652138414695194151160943305727036575959195309218611738193261179310511854807446237996274956735188575272489122793818301194912983367336244065664308602139494639522473719070217986094370277053921717629317675238467481846766940513200056812714526356082778577134275778960917363717872146844090122495343014654958537105079227968925892354201995611212902196086403441815981362977477130996051870721134999999837297804995105973173281609631859502445945534690830264252230825334468503526193118817101000313783875288658753320838142061717766914730359825349042875546873115956286388235378759375195778185778053217122680661300192787661119590921642019;
% Computed with R package library("Rmpfr"): 
% ( log2_ = Const("log2", 1000*log2(10) ) ); 
% ( Pi <- Const("pi", 1000 *log2(10)) );
% ( log_pi = mpfr( log(Pi), 1000*log2(10) ) ); 
% ( log2_pi = mpfr( log2_ + log_pi, 1000*log2(10) ) );
% log_2_pi = 1.83787706640934548356065947281123527972279494727556682563430308096553139185452079538948659727190839524401129324926867489273372576368158714431175183044536278720712148509471733809279181198276161126032646974618925474925103650338990895482019171870278396322319626114801069539077212991798446242791138554869994220056703919663898506278854129259137294882312495242609747363056899875868876466079702589530931456386347597570617137884627256430794616720529505853098298007871119999920741269437051440471524307006872475920543169750097227190768496265835824853999227536792803027895754591002020664176839367123881595143325254117505076497245186050590421609903624039361045196009176107714976706588822781361565555347544450762667651879014828040523867874263374089441371189156869826552081590826015367960940350517749618771749114464650668778489385596557499370542251617516233174875058017696896618350778815259190881989693579607832426181446570287357290751247594207086908526347557529234407222834527535937679132380540148826095822827999;
% OR FROM PROJECT GUTENBERG AT:
% log_2_pi = 1.83787706640934548356065947281123527972279494727556682563430308096553139185452079538948659727190839524401129324926867489273372576368158714431175183044536278720712148509471733809279181198276161126032646974618925474925103650338990895482019171870278396322319626114801069539077212991798446242791138554869994220056703919663898506278854129259137294882312495242609747363056899875868876466079702589530931456386347597570617137884627256430794616720529505853098298007871119999920741269437051440471524307006872475920543169750097227190768496265835824853999227536792803027895754591002020664176839367123881595143325254117505076497245186050590421609903624039361045196009176107714976706588822781361565555347544450762667651879014828040523867874263374089441371189156869826552081590826015367960940350517749618771749114464650668778489385596557499370542251617516233174875058017696896618350778815259190881989693579607832426181446570287357290751247594207086908526347557529234407222834527535937679132380540148826095822827999;
% http://www.gutenberg.org/files/634/634.txt

% DIMENSION OF y(t) DIM(y(t)) AND SAMPLE SIZE T
[Ky, T] = size(y);
% CHECK INPUT SIZE OF OBSERVABLE VECTOR Y WHICH SHOULD BE (Ky x T)
if Ky > T
	error('Y has to be a (Ky x T) matrix, to loop over the columns');
end
% DIMENSION OF STATE VECTOR ALPHA dim(alpha(t))
% Ka = size(Phi,1);		

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
% continue from there. this is also what the ssm toolbox in matlab does.
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
sumLL = 0;
% LLt	= zeros(1,T);
% SPACE ALLOCATION FOR THINGS THAT ARE STORED FOR LATER ON [Ka = length(alpha), Ky = length(y)]
% at = zeros(Ka,T);		  uhat_store	= zeros(Ky,T);		
% Pt = zeros(Ka,Ka,T);  Ft_store		= zeros(Ky,Ky,T);		

% MAIN FORWARD LOOP (LOOPING OVER COLUMNS IN MATLAB IS FASTER THAN OVER ROWS)
sumLL = kalman_filter_loop_fast(y, D, M, H, C, Phi, SQS, att_1, Ptt_1, T, Ky);
% sumLL = kalman_filter_loop_mex(y, D, M, H, C, Phi, SQS, att_1, Ptt_1, T, Ky);

% for t = 1:T
% 	% FORECAST ERROR AND ITS MSE
% 	uhat	= y(:,t) - M*att_1 - D(:,t);
% 	Ft		= M*Ptt_1*M' + H;
% 	% PRECOMPUTE PM_invF = P(t|t-1)*M'inv(F(t))
%  	P_MT_invF = Ptt_1*M'/Ft;
% 	% UPDATING/FILTERING a(t|t) and P(t|t)	
% 	att		= att_1 + P_MT_invF*uhat;
% 	Ptt		= Ptt_1 - P_MT_invF*M*Ptt_1;
% 	% PREDICTION/FORECASTING a(t|t-1) and P(t|t-1)
% 	att_1	= Phi*att + C(:,t);
% 	Ptt_1	= Phi*Ptt*Phi' + SQS;
% 	% ----------------------------------------------------------------------------------------------------		
% 	% STORE THE FILTERED STATES/MSEs as well AS Ft and uhat (prediction error and its variance for LogLike)
% 	% ----------------------------------------------------------------------------------------------------	
% %  	uhat_store(:,t)	= uhat;		
% %  	Ft_store(:,:,t)	= Ft;
% % 	at(:,t)	  = att;
% % 	Pt(:,:,t) = Ptt;
% 	% ----------------------------------------------------------------------------------------------------	
% 	% LOG-LIKELIHOOD 
% 	% ----------------------------------------------------------------------------------------------------	
% 	LL_t	= -0.5*( Ky*log_2_pi + log(det(Ft)) + uhat'/Ft*uhat );
% 	sumLL	= sumLL + LL_t;
% % 	LLt(1,t)= LL_t;
% end

% % ALLOCATE SPACE for backward recursions AND initialize at a(T|T) and P(T|T)
% atT = at;  PtT = Pt;
% 
% % BACKWARDS RECURSION ON a(t|T) and P(t|T) see page 37 in Kim and Nelson (1999)
% for n = T-1:-1:1
%   % FIRST PRE-COMPUTE SOME QUANTITIES
%   % P(t+1|t) = Phi*P(t|t)Phi' + SQS
%   Pt1_t = Phi*Pt(:,:,n)*Phi' + SQS;
%   % KG(t) = P(t|t)*Phi'inv(P(t+1|t))
%   KGt = Pt(:,:,n)*Phi'/Pt1_t;
%   atT(:,n) = at(:,n) + KGt*( atT(:,n+1) - C(:,t) - Phi*at(:,n) );
%   PtT(:,:,n) = Pt(:,:,n) + KGt*( PtT(:,n+1) - Pt1_t )*KGt'; 
% end
% 
% % NOW RETURN BACK TO FUNCTION AS STRUCTURE
% KFS.att  = at'; 		      % Filtered states
% KFS.atT  = atT';					% Smoothed states
% KFS.Ptt  = Pt;			      % Filtered states Var/Cov
% KFS.PtT  = PtT;			      % Smoothed states Var/Cov
% KFS.uhat = uhat_store';		% (Ky x T) Prediction error
% KFS.Ft	 = Ft_store;			% (Ky x Ky x T) Variance of the prediction error
% KFS.LL_t = LLt';					% Full Log-Likelihood
% KFS.LL_sum = LL_sum;			% Sum Log-Likelihood

% RETURN LOG-LIKELIHOOD ONLY
logLikeout = sumLL;
end 












































% EOF