function setxticklabels(xticks,fig_handle)
% set or adjust the yticklables and ticks at desired location vector ticks.
% =======================================================================================
% input
% ticks:			(kx1), location vector of ticks which will also be the lables.
% fig_handle:	scalar, figure handle, default is gca.
%
% call as:		setyticklabel([-.1:.1.1]) to set ticks and ticklables as the vector.
%							[-.1:.1.1] in the y-axis ticks.
% =======================================================================================
% db 14.06.2012
% modfified on 10.07.2013
% =======================================================================================

if nargin < 2
	fig_handle = gca;
end

xlim([xticks(1) xticks(end)])

xticksLabls = strtrim(cellstr(num2str(xticks')));
% regexprep(cellstr(num2str(xticks')),' ','')

set(fig_handle,'XTick'			,xticks);
set(fig_handle,'XTicklabel'	,xticksLabls)

end
