function [Inan, varargout] = anynans(X,Y)
% find any rows with NANs in regressor matrix X and also in Y if inputed and return indicator Inan
% useful for workign with time series data that has been lagged by mlag or differenced by delta.
% can be used to trim out the require entries for OLS etc. 
% then return rows with anynans.

if strcmp(class(X),'fints')
	anyX = any(isnan(fts2mat(X)),2);
	anyXall = isnan(fts2mat(X));
else
	anyX = any(isnan(X),2);
	anyXall = isnan(X);
end

if nargin < 2
  Inan = anyX;
%  disp('This')
else
  anyY = any(isnan(Y),2);
  Inan = anyX|anyY;
end;

if nargout > 1
	varargout{1} = anyXall;
end;
	
%%Inan = Inan==0;
