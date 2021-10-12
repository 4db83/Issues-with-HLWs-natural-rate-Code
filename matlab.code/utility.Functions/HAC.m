function [EstCoeffCov,se,coeff] = HAC(Data,varargin)
%HAC Heteroscedasticity and autocorrelation consistent covariance estimators
%
% Syntax:
%
%   [EstCoeffCov,se,coeff] = hac(X,y)
%   [EstCoeffCov,se,coeff] = hac(DS)
%   [EstCoeffCov,se,coeff] = hac(Mdl)
%   [EstCoeffCov,se,coeff] = hac(...,param,val,...)
%
% Description:
%
%   HAC computes robust covariance estimates for ordinary least squares
%   (OLS) coefficient estimates of multiple linear regression models
%
%   y = X*coeff + e
%
%   under general forms of heteroscedasticity and autocorrelation in the
%   innovations process e. Estimates are of the form:
%               
%   PhiHat = X'*OmegaHat*X
%   EstCoeffCov = c*inv(X'*X)*PhiHat*inv(X'*X)
%
%   where OmegaHat is an estimate of the innovations covariance and c is an
%   optional small-sample correction.
%
% Input Arguments:  
%
%   X - numObs-by-numPreds matrix of predictor data for a multiple linear
%       regression model.
%
%   y - numObs-by-1 vector of response data for a multiple linear
%       regression model.
%
%   DS - numObs-by-numPreds+1 dataset array of data for a multiple linear
%       regression model, with predictor data X in the first numPreds
%       columns and response data y in the last column.
%
%   Mdl - Linear model of the form Mdl = LinearModel.fit(X,y).
%
%   Observations with missing (NaN) values in the predictors or the
%   response are removed from the data.
%
% Optional Input Parameter Name/Value Pairs:
%
%   NAME        VALUE
%
%   'varNames'  Cell vector of variable name strings of length numPreds to
%               be used in displays of the results. Names should include
%               the intercept term (for example, 'Const') and higher-order
%               terms (for example, 'x1^2' or 'x1:x2'), if present. The
%               default for matrix X is {'x1','x2',...}. The default for
%               dataset array DS is DS.Properties.VarNames. The default for
%               linear models Mdl is Mdl.CoefficientNames.
%
%   'intercept' Logical value indicating whether or not to add an intercept
%               when fitting the model. The default is true. When the input
%               is already a fitted model, Mdl, 'intercept' is ignored.
%
%   'type' 		Type of estimator. Values are 'HC' and 'HAC'. 'HC' computes
%               heteroscedasticity-consistent estimates, as described in
%               [9], [7], and [3]. 'HAC' computes heteroscedasticity-and-
%               autocorrelation-consistent estimates, as described in [8],
%               [5], [2], and [1]. The default is 'HAC'.
%
%   'weights' 	Weights used to compute PhiHat = X'*OmegaHat*X. Values are
%               a numerical vector w of length numObs or a string
%               indicating a method for computing w from the data.
%
%               When 'type' is 'HC', OmegaHat = diag(w). Elements w(i)
%               estimate the innovations variance at each observation i.
%               Data-driven w are computed from model residuals u(i), their
%               leverages h(i), and their number of degrees of freedom dfe.
%               Values are:
%
%               o 'CLM' w(i) = (1/dfe)*sum(u(i)^2). This is the Classical
%                       Linear Model estimator, used in the absence of
%                       heteroscedasticity.
%
%               o 'HC0' w(i) = u(i)^2. This is the White estimator in [9].
%                       This is the default when 'type' is 'HC'.
%
%               o 'HC1' w(i) = (numObs/dfe)*u(i)^2. This is the HC1
%                       estimator in [7].
%
%               o 'HC2' w(i) = (u(i)^2)/(1-h(i)). This is the HC2
%                       estimator in [7].
%
%               o 'HC3' w(i) = (u(i)^2)/((1-h(i))^2). This is the HC3
%                       estimator in [7].
%
%               o 'HC4' w(i) = (u(i)^2)/((1-h(i))^d(i)), where d(i) is
%                       min(4,h(i)/mean(h)). This is the estimator in [3].
%
%               o w     Vector of user-specified weights, of length numObs.
%                       If data contain a missing value at observation i,
%                       w(i) is removed during estimation.
%
%               When 'type' is 'HAC', component products forming PhiHat,
%               X(i,:)'*(u(i)*u(j))*X(j,:), are weighted by w(l), where 
%               l = abs(i-j) is the lag length between observations.
%               Elements w(l) represent the strength of autocorrelations at
%               each lag length. Data-driven w(l) = k(l/b) are computed by
%               a kernel density estimator k, with b specified by the
%               'bandwidth' parameter. Values are:
%
%               o 'TR' Truncated kernel:
%
%                             | 1  for abs(x) <= 1
%                      k(x) = |
%                             | 0  otherwise
%
%                      This produces the White estimator in [10, p. 152].
%
%               o 'BT' Bartlett kernel:
%
%                             | 1-abs(x)  for abs(x) <= 1
%                      k(x) = |
%                             | 0         otherwise
%
%                      This produces the Newey-West estimator in [8]. This
%                      is the default when 'type' is 'HAC'.
%
%               o 'PZ' Parzen kernel:
%
%                             | 1-6*x^2+6*abs(x)^3  for 0 <= abs(x) <= 1/2
%                      k(x) = | 2*(1-abs(x))^3      for 1/2 <= abs(x) <= 1
%                             | 0                   otherwise
%
%                      This produces the Gallant estimator in [5, p. 533].
%
%               o 'TH' Tukey-Hanning kernel:
%
%                             | (1+cos(pi*x))/2  for abs(x) <= 1
%                      k(x) = |
%                             | 0                otherwise
%
%                      This estimator was introduced by Andrews in [1].
%
%               o 'QS' Quadratic spectral kernel:
%
%                      k(x) = [25/(12*pi^2*x^2)]*...
%                             [sin(6*pi*x/5)/(6*pi*x/5)-cos(6*pi*x/5)]
%
%                      This estimator was introduced by Andrews in [1].
%
%               o w    Vector of user-specified weights, of length numObs.
%                      If data contain a missing value at observation i,
%                      w(i) is removed during estimation.
%
%   'bandwidth' Scalar bandwidth parameter b used by the kernel estimator
%               in 'weights' when 'type' is 'HAC'. If 'type' is 'HC',
%               'bandwidth' is ignored. The default is data-driven, as in
%               [1]. Values of 'AR1' or 'ARMA11' specify the model used for
%               data-driven bandwidth selection. The default model is
%               'AR1'. The method used to estimate the 'AR1' model can be
%               specified by values of 'AR1OLS', to use OLS, or 'AR1MLE',
%               to use maximum likelihood. The default method is 'AR1MLE'.
%
%   'smallT'    Logical flag indicating whether or not to apply the small-
%               sample correction c = T/dfe, where T is the effective
%               sample size and dfe is the number of degrees of freedom of
%               the model residuals. Values are true and false. The default
%               is false when 'type' is 'HC' and true when 'type' is 'HAC'.
%
%   'whiten'    Nonnegative integer specifying the lag order of a vector
%               autoregressive model used as a prewhitening filter when
%               'type' is 'HAC', as in [2]. If 'type' is 'HC', 'whiten' is
%               ignored. The default is 0, which bypasses the filter.
%
%   'display'   String to control the display of results to the command
%               window. The display shows outputs in tabular form. Values
%               are 'cov' to display only EstCoeffCov, 'full' to display
%               all of coeff, se, and EstCoeffCov, or 'off' to turn off the
%               display. The default is 'cov'.
% 
% Output Arguments:
%
%   EstCoeffCov - numPreds-by-numPreds matrix of coefficient covariance
%                 estimates.
%
%   se - numPreds-by-1 vector of coefficient standard error estimates. This
%        is sqrt(diag(EstCoeffCov)).
%
%   coeff - numPreds-by-1 vector of OLS coefficient estimates.
%
% Notes:
%
%   o Estimates formed by HAC are often called "sandwich estimates," with
%     X'*OmegaHat*X the "meat" and inv(X'*X) the "bread."
%
%   o For kernels with unit-interval support, the bandwidth parameter b is
%     often called the "lag-truncation parameter," since w(l) = k(l/b) = 0
%     for lags l > b.
%
%   o HAC with default settings computes a Newey-West estimator using a
%     Bartlett kernel with a data-driven bandwidth. To obtain the standard
%     Newey-West estimate described in [8], set 'bandwidth' to maxLag+1,
%     with maxLag = floor(4*(T/100)^(2/9)).
%
%   o The original White HC estimator, HC0, is justified asymptotically.
%     HC1, HC2, HC3, and HC4 are meant to improve small-sample performance.
%     [6] and [3] recommend HC3 and HC4, respectively, in the presence of
%     influential observations.
%
%   o HAC estimators formed with the truncated kernel are not guaranteed to
%     be positive semidefinite in finite samples. [8] proposes the Bartlett
%     kernel as a remedy, but the resulting estimator is suboptimal in
%     terms of its rate of consistency. The quadratic spectral kernel
%     achieves an optimal rate of consistency.
%
%   o The default estimation method for HAC bandwidth selection, 'AR1MLE',
%     is generally more accurate, but slower, than the 'AR1' alternative,
%     'AR1OLS'. The 'ARMA11' model must be estimated by maximum likelihood.
%     Bandwidth-selection models may exhibit sensitivity to the relative
%     scale of the predictors in X.
%
%   o [2] recommends prewhitening for HAC estimators to reduce bias. The
%     procedure tends to increase estimator variance and mean-squared
%     error, but can improve confidence interval coverage probabilities and
%     reduce the over-rejection of t statistics.
%
% Example:
%
%   % Compute OLS coefficients and Newey-West standard errors for a
%   % regression of nominal GNP (GNPN) on the consumer price index (CPI),
%   % real wages (WR), and the money stock (MS):
%
%   load('Data_NelsonPlosser.mat')
%   DS = Dataset(:,[8,10,11,2]);
%   T = sum(~any(isnan(double(DS)),2));
%   maxLag = floor(4*(T/100)^(2/9));
%   EstCoeffCov = hac(DS,'bandwidth',maxLag+1,'display','full');
%
% References:
% 
%   [1] Andrews, D. W. K. "Heteroskedasticity and Autocorrelation
%       Consistent Covariance Matrix Estimation." Econometrica. v. 59,
%       1991, pp. 817-858.
%
%   [2] Andrews, D. W. K., and J. C. Monohan. "An Improved
%       Heteroskedasticity and Autocorrelation Consistent Covariance Matrix
%       Estimator." Econometrica. v. 60, 1992, pp. 953-966.
%
%   [3] Cribari-Neto, F. "Asymptotic Inference Under Heteroskedasticity of
%       Unknown Form." Computational Statistics & Data Analysis. v. 45,
%       2004, pp. 215-233.
%
%   [4] den Haan, W. J., and A. Levin. "A Practitioner's Guide to Robust
%       Covariance Matrix Estimation." In Handbook of Statistics. Edited by
%       G. S. Maddala and C. R. Rao. Amsterdam: Elsevier, 1997.
%
%   [5] Gallant, A. R. Nonlinear Statistical Models. Hoboken, NJ: John
%       Wiley & Sons, Inc., 1987.
%
%   [6] Long, J. S., and L. H. Ervin. "Using Heteroscedasticity-Consistent
%       Standard Errors in the Linear Regression Model." The American
%       Statistician. v. 54, 2000, pp. 217-224.
%
%   [7] MacKinnon, J. G., and H. White. "Some Heteroskedasticity-Consistent
%       Covariance Matrix Estimators with Improved Finite Sample
%       Properties." Journal of Econometrics. v. 29, 1985, pp. 305-325.
%
%   [8] Newey, W. K., and K. D. West. "A Simple, Positive-Definite,
%       Heteroskedasticity and Autocorrelation Consistent Covariance
%       Matrix." Econometrica. v. 55, 1987, pp. 703-708.
%
%   [9] White, H. "A Heteroskedasticity-Consistent Covariance Matrix and a
%       Direct Test for Heteroskedasticity." Econometrica. v. 48, 1980, pp.
%       817-838.
%
%  [10] White, H. Asymptotic Theory for Econometricians. New York: Academic
%       Press, 1984.
%
% See also LinearModel.fit, lscov.

% Copyright 2013 The MathWorks, Inc.

% Parse inputs and set defaults:

parseObj = inputParser;
parseObj.addRequired('Data',@DataCheck);
parseObj.addOptional('y',[],@yCheck);
parseObj.addParamValue('varNames',[],@varNamesCheck);
parseObj.addParamValue('intercept',true,@interceptCheck);
parseObj.addParamValue('type','HAC',@typeCheck);
parseObj.addParamValue('weights',[],@weightsCheck);
parseObj.addParamValue('bandwidth',[],@bandwidthCheck);
parseObj.addParamValue('smallT',[],@smallTCheck);
parseObj.addParamValue('whiten',0,@whitenCheck);
parseObj.addParamValue('display','cov',@displayCheck);

parseObj.parse(Data,varargin{:});

Data = parseObj.Results.Data;
y = parseObj.Results.y;
varNames = parseObj.Results.varNames;
iFlag = parseObj.Results.intercept;
estType = upper(parseObj.Results.type);
weights = upper(parseObj.Results.weights);
b = parseObj.Results.bandwidth;
ssFlag = parseObj.Results.smallT;
p = parseObj.Results.whiten;
display = lower(parseObj.Results.display);

% Convert inputs to LinearModel object:

if isnumeric(Data)
    
    if isempty(y)
        
        error(message('econ:hac:EmptyResponseVector'))
          
    else
    
        try
            M = LinearModel.fit(Data,y,'Intercept',iFlag);
        catch exception
            throw(exception)
        end
        
    end
    
elseif isa(Data,'dataset')
    
    try   
        M = LinearModel.fit(Data,'Intercept',iFlag);    
    catch exception  
        throw(exception)
    end
    
else
    
    M = Data;
    iFlag = M.Formula.HasIntercept;
    
end

% Extract relevant model information:

missingData = M.ObservationInfo.Missing;
T = M.NumObservations;
coeff = M.Coefficients.Estimate;
numPreds = M.NumCoefficients;
u = M.Residuals.Raw;
u(missingData) = [];
dfe = M.DFE;
s = warning('off','MATLAB:structOnObject');
MStruct = struct(M);
X = MStruct.Design;
X(missingData,:) = [];
clear MStruct
warning(s)

% Specify default weights:

if isempty(weights)
    if strcmp(estType,'HC')
        weights = 'HC0';
    else
        weights = 'BT';
    end
end

% Check user-specified weights:

userWeights = isnumeric(weights);

if userWeights
    
    try
        w = userWeightsCheck(weights,T,missingData);
    catch exception
        throw(exception)
    end
    
end

% Compute PhiHat:
        
switch estType

    case 'HC'
        
        if ~userWeights

            u2 = u.^2;
            h = M.Diagnostics.Leverage;
            h(missingData) = [];

            switch weights

                case 'CLM'

                    sse = sum(u2);
                    w = repmat(sse/dfe,size(u));

                case 'HC0'

                    w = u2;

                case 'HC1'

                    w = (T/dfe)*u2;

                case 'HC2'

                    w = u2./(1-h);

                case 'HC3'

                    w = u2./((1-h).^2);

                case 'HC4'

                    d = min(4,h/mean(h));
                    w = u2./((1-h).^d);

                otherwise

                    error(message('econ:hac:HCWeightsInvalid'))

            end
        
        end
        
        PhiHat = X'*spdiags(w,0,T,T)*X;

    case 'HAC'

        % Scores matrix:

        V = bsxfun(@times,X,u);

        % Prewhiten scores:

        if p > 0

            if p > T-1

                error(message('econ:hac:WhiteningOrderTooBig'))

            end

            VAR = vgxset('n',numPreds,'nAR',p);
            V0 = V(1:p,:);
            V1 = V(p+1:end,:);
            [VARfit,~,~,V] = vgxvarx(VAR,V1,[],V0);
            VARcoeffs = vgxget(VARfit,'AR');

        end
        
        if ~userWeights

            % Compute bandwidth, if unspecified:

            if isempty(b) || ischar(b)

                % Specify method:
                if isempty(b)
                    model = 'AR1MLE';
                elseif ischar(b)                
                    model = b;  
                end

                try
                    b = getBW(V,weights,model,iFlag);
                catch exception
                    throw(exception)
                end

            end

            % Compute kernel weights:

            lags = 0:T-1;
            w = zeros(T,1);
            x = lags/b;           

            switch weights

                case 'TR'

                    TR = (abs(x) <= 1);
                    w(TR) = 1;

                case 'BT'

                    BT = (abs(x) <= 1);                
                    w(BT) = 1-abs(x(BT));

                case 'PZ'

                    PZ1 = (abs(x) >= 0) & (abs(x) <= 1/2);
                    PZ2 = (abs(x) >= 1/2) & (abs(x) <= 1);  
                    w(PZ1) = 1-6*x(PZ1).^2+6*abs(x(PZ1)).^3;
                    w(PZ2) = 2*(1-abs(x(PZ2))).^3;

                case 'TH'

                    TH = (abs(x) <= 1);
                    w(TH) = (1+cos(pi*x(TH)))/2;

                case 'QS'

                    argQS = 6*pi*x/5;
                    w1 = 3./(argQS.^2);
                    w2 = (sin(argQS)./argQS)-cos(argQS);
                    w = w1.*w2;
                    w(x == 0) = 1;

                otherwise

                    error(message('econ:hac:HACWeightsInvalid'))

            end
        
        end
        
        PhiHat = w(1)*(V'*V);
        for i = 1:T-p-1
            LagCov = (V(1+i:T-p,:)'*V(1:T-p-i,:));
            PhiHat = PhiHat+w(i+1)*(LagCov+LagCov');
        end

end

% Recolor prewhitened PhiHat:

if strcmp(estType,'HAC') && (p > 0)

    C = eye(numPreds)-sum(reshape([VARcoeffs{:}],numPreds,numPreds,p),3);
    PhiHat = (C\PhiHat)/C';

end

% Compute the sandwich estimate:

% EstCoeffCov = inv(X'*X)*PhiHat*inv(X'*X)  
%             = inv(R'*Q'*Q*R)*Z'*Z*inv(R'*Q'*Q*R)
%             = ((R'*R)\Z')*(Z/(R'*R))
%             = (S\Z')*(Z/S)

try    
    Z = chol(PhiHat);    
catch exception % 'TR' or w fails to produce a positive definite matrix   
    throw(exception)    
end

[~,R] = qr(X,0);
S = R'*R;
EstCoeffCov = (S\Z')*(Z/S);

% Specify the default small-sample correction:

if isempty(ssFlag)
    if strcmp(estType,'HC')
        ssFlag = false;
    else
        ssFlag = true;
    end
end

% Apply the small-sample correction:

if ssFlag
    EstCoeffCov = (T/dfe)*EstCoeffCov;
end

% Compute standard errors:

se = sqrt(diag(EstCoeffCov));

% Display results:

if ~strcmp(display,'off')
  
    if isempty(varNames)
        
        varNames = M.CoefficientNames;
        
        if iFlag
            
            varNames{1} = 'Const';
            
        end
        
    else
        
        if length(varNames) < numPreds

            error(message('econ:hac:VarNamesTooFew'))

        elseif length(varNames) > numPreds

            error(message('econ:hac:VarNamesTooMany'))

        end
        
    end

    % Print output:
    
    fprintf('\nEstimator type: %s',estType)
    if userWeights        
        fprintf('\nEstimation method: user-defined weights')        
    else        
        fprintf('\nEstimation method: %s',weights)
        if strcmp(estType,'HAC')
            fprintf('\nBandwidth: %.4f',b)
            fprintf('\nWhitening order: %u',p)
        end        
    end
    fprintf('\nEffective sample size: %u',T)
    if ssFlag
        fprintf('\nSmall sample correction: on\n')
    else
        fprintf('\nSmall sample correction: off\n')
    end

    if strcmp(display,'full')        
        fprintf('\nCoefficient Estimates:\n\n')
        internal.econ.tableprint([coeff,se],...
                                 'colNames',{'Coeff','SE'},...
                                 'rowNames',varNames)
    end
    
    fprintf('\nCoefficient Covariances:\n\n')
    internal.econ.tableprint(EstCoeffCov,...
                             'colNames',varNames,...
                             'rowNames',varNames)
    
end

%-------------------------------------------------------------------------
% Check input Data
function OK = DataCheck(Data)

if isempty(Data)

    error(message('econ:hac:DataUnspecified'))

elseif ~isnumeric(Data) && ~isa(Data,'dataset') && ~isa(Data,'LinearModel')

    error(message('econ:hac:DataFormatUnsupported'))

else
    
    OK = true;

end

%-------------------------------------------------------------------------
% Check response vector
function OK = yCheck(y)

if ~isnumeric(y)

    error(message('econ:hac:ResponseDataNonNumeric'))

elseif ~isvector(y)

    error(message('econ:hac:ResponseDataNonVector'))

else

    OK = true;

end

%-------------------------------------------------------------------------
% Check value of 'varNames' parameter
function OK = varNamesCheck(varNames)

if ~isvector(varNames)

    error(message('econ:hac:VarNamesNonVector'))

elseif isnumeric(varNames) || (iscell(varNames) && any(cellfun(@isnumeric,varNames)))

    error(message('econ:hac:VarNamesNumeric'))

else

    OK = true;

end

%-------------------------------------------------------------------------
% Check value of 'intercept' parameter
function OK = interceptCheck(intercept)

if ~islogical(intercept)

    error(message('econ:hac:InterceptNonBoolean'))

elseif ~isscalar(intercept)

    error(message('econ:hac:InterceptNonScalar'))

else

    OK = true;

end

%-------------------------------------------------------------------------
% Check value of 'type' parameter
function OK = typeCheck(estType)
    
if ~ischar(estType) || ~isvector(estType)

    error(message('econ:hac:EstTypeNonString'))

elseif ~ismember(upper(estType),{'HC','HAC'})

    error(message('econ:hac:EstTypeInvalid'))

else

    OK = true;

end

%-------------------------------------------------------------------------
% Check value of 'weights' parameter
function OK = weightsCheck(weights)
    
if ~isvector(weights)

    error(message('econ:hac:WeightsNonVector'))

elseif ischar(weights) && ...
       ~ismember(upper(weights),...
       {'CLM','HC0','HC1','HC2','HC3','HC4','TR','BT','PZ','TH','QS'})

   error(message('econ:hac:WeightsInvalid'))

else

    OK = true;

end

%-------------------------------------------------------------------------
% Check value of 'bandwidth' parameter
function OK = bandwidthCheck(b)

if ~isvector(b)

    error(message('econ:hac:BandwidthNonVector'))
      
elseif isnumeric(b) && ~isscalar(b)

    error(message('econ:hac:BandwidthValueNonScalar'))
      
elseif ischar(b) && ~ismember(upper(b),{'AR1','AR1OLS','AR1MLE','ARMA11'})
    
    error(message('econ:hac:BandwidthEstimationMethodInvalid'))

else

    OK = true;

end

%-------------------------------------------------------------------------
% Check value of 'smallT' parameter
function OK = smallTCheck(ssFlag)
    
if ~islogical(ssFlag)

    error(message('econ:hac:SmallSampleFlagInvalid'))

else

    OK = true;

end

%-------------------------------------------------------------------------
% Check value of 'whiten' parameter
function OK = whitenCheck(p)

if ~isnumeric(p) || ~isscalar(p) || mod(p,1) ~= 0 || p < 0

    error(message('econ:hac:WhiteningOrderInvalid'))

else

    OK = true;

end
    
%-------------------------------------------------------------------------
% Check value of 'display' parameter
function OK = displayCheck(display)
    
if ~isvector(display)

    error(message('econ:hac:DisplayFlagNonVector'))

elseif isnumeric(display) || (iscell(display) && any(cellfun(@isnumeric,display)))

    error(message('econ:hac:DisplayFlagNumeric'))

elseif ~all(ismember(lower(display),{'cov','full','off'}))

    error(message('econ:hac:DisplayFlagInvalid'))

else

    OK = true;

end
    
%-------------------------------------------------------------------------
% Check user-defined weights
function w = userWeightsCheck(weights,numObs,missingData)

if length(weights) ~= numObs

    error(message('econ:hac:UserWeightsDimensionMismatch'))

else

    weights(missingData) = [];
    w = weights;

end

%-------------------------------------------------------------------------
% Compute data-driven bandwidth
function b = getBW(V,weights,model,iFlag)

[T,n] = size(V);

% Use the scale-invariant alpha weights in [1]. See [4] for alternatives.

wAlpha = ones(n,1); % Unit weight for slope coefficients
if iFlag
    wAlpha(1) = 0; % Zero weight for intercept coefficient
end

% Set MLE options:

options = optimset('fmincon');
options = optimset(options,'display','off','diagnostics','off',...
                           'algorithm','sqp','tolCon',1e-7);
                               
switch model
    
    case {'AR1','AR1OLS','AR1MLE'}

        rho = zeros(n,1);     % AR coefficient estimates
        sigmaSq = zeros(n,1); % Residual variances
        
        switch model

            case {'AR1','AR1MLE'}

                AR1 = arima(1,0,0); % AR(1) model

                for j = 1:n
                    
                    ARfit = estimate(AR1,V(:,j),'print',false,'options',options);
                    rho(j) = ARfit.AR{1};
                    sigmaSq(j) = ARfit.Variance;
                
                end

            case 'AR1OLS'

                for j = 1:n
                    
                    v = V(:,j);
                    vt = v(2:T);
                    vtLag1 = v(1:T-1);
                    rhoj = vtLag1\vt;
                    res = vt-rhoj*vtLag1;
                    rho(j) = rhoj;
                    sigmaSq(j) = res'*res/(T-1);
                
                end

        end

        num0 = 4*(rho.^2).*(sigmaSq.^2);
        den = (sigmaSq.^2)./((1-rho).^4);
        
    case 'ARMA11'

        rho = zeros(n,1);      % AR coefficient estimates
        psi = zeros(n,1);      % MA coefficient estimates
        sigmaSq = zeros(n,1);  % Residual variances
        
        ARMA11 = arima(1,0,1); % ARMA(1,1) model
        
        for j = 1:n
            
            v = V(:,j);
            ARMAfit = estimate(ARMA11,v,'print',false,'options',options);
            rho(j) = ARMAfit.AR{1};
            psi(j) = ARMAfit.MA{1};
            sigmaSq(j) = ARMAfit.Variance;
            
        end
        
        num0 = 4*((1+rho.*psi).^2).*((rho+psi).^2).*(sigmaSq.^2);
        den = (((1+psi).^4).*(sigmaSq.^2))./((1-rho).^4);
        
end

% Compute alphas:

num1 = ((1-rho).^6).*((1+rho).^2);
num2 = (1-rho).^8;
denSum = sum(wAlpha.*den);
alpha1 = sum(wAlpha.*(num0./num1))/denSum;
alpha2 = sum(wAlpha.*(num0./num2))/denSum;

% Compute optimal bandwidth:

switch weights
    
    case 'TR'
        
       b = (0.6611)*(alpha2*T)^(1/5); 
        
    case 'BT'
        
        b = (1.1447)*(alpha1*T)^(1/3);
                
    case 'PZ'
        
        b = (2.6614)*(alpha2*T)^(1/5);
                
    case 'TH'
        
        b = (1.7462)*(alpha2*T)^(1/5);
                
    case 'QS'
      
        b = (1.3221)*(alpha2*T)^(1/5);
        
    otherwise
        
        error(message('econ:hac:HACWeightsInvalid'))
        
end