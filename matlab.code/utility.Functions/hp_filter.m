function [cycle,trend,m] = hp_filter(data,lambda)
%{F: fast and simple HP filter implementation.
%===============================================================================
% very fast HP filter using sparse matrix functions
%-------------------------------------------------------------------------------
% 	USGAGE		[cycle,trend,m] = hp_filter(data,lambda)
%-------------------------------------------------------------------------------
% 	INPUT  
%	  data:			(Tx1) data vector to be HP filtered.
%		lambad:		scalar, smoothing parameter. 
%
% 	OUTPUT       
%	  cycle:		(Tx1) vector of cyclical component.
%		trend:		(Tx1) vector of permanent/trend component.
%		m:				(TxT) matrix that shows the moving average weights of the HP filter.
%===============================================================================
% 	NOTES :   none.
%-------------------------------------------------------------------------------
% Created :		long time ago.
% Modified:		24.01.2014.
% Copyleft:		Daniel Buncic.
%------------------------------------------------------------------------------%}

T		= size(data,1);
I		= speye(T);							% sparse identiy matrix.
D		= diff(I,2);						% Difference maker matrix.
m		= (I + lambda*(D'*D));	% shows the moving average weight matrix if needed.
trend	= m\data;							% permanent component.
cycle	= data - trend;				% transitory component. 
