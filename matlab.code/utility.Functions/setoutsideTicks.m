function Hout = setoutsideTicks(tickSize, currentHandle)
%{F: Sets Y axis ticks on the outside when Yaxis labels are printed on boths sides
% ---------------------------------------------------------------------------------------------
% CALL AS: 
%		setoutsideTicks										or as
%		setoutsideTicks(tickSize)
% 	setoutsideTicks(tickSize, currentHandle)
%										
%	INPUT:	
%		ticksize = ticksize adjustment shrinks the ticksize to 3/4 of original length.
%		currentHandle = handle to current figure (Defaults to gca)
% OUTPUT: 
%		Hout = handle to axis objects. if not called not returned.
% ---------------------------------------------------------------------------------------------
% Created :		07.01.2019.
% Modified:		07.01.2019.
% Copyleft:		Daniel Buncic.
% ---------------------------------------------------------------------------------------------}

SetDefaultValue(1,'tickSize'			,3/4);
SetDefaultValue(2,'currentHandle'	,gca);

h1 = currentHandle;
h2 = copyobj(h1, gcf);
% set(h1, 'TickDir', 'out', 'YTick', []);
set(h1, 'TickDir', 'out', 'XTick', []);
set(h2, 'TickDir', 'in'	, 'YTickLabel', []);

tickshrink(tickSize, h2);

if length(tickSize)>1
	tickshrink(tickSize(2), h1);
	tickshrink(tickSize(1), h2);
else
	tickshrink(tickSize, h2);
end


if nargout
	Hout = [h1 h2];
end