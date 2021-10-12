function x = demean(y)
% Demeans an N-dimensional vector
% 
% USAGE:
%   DEMEANE_DATA = demean(DATA)
% 
% INPUTS:
%   DATA     - An (TxK)-dimensional matrix
%
% OUTPUTS:
%   DEMEANED - A mean 0 matrix with the same size as DATA  (column-by-column)
%
% COMMENTS:
%
% See also STANDARDIZE


mu	= nanmean(y,1);
x		= bsxfun(@minus,y,mu);

