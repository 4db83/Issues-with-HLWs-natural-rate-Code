function varargout= setyticklabels(ticks,d,FNS,fig_handle)
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


FNS_0 = get(gca,'FontSize');
FNS0	= FNS_0 - 0;

SetDefaultValue(2, 'd'	, 2);
SetDefaultValue(3, 'FNS', FNS0);
SetDefaultValue(4, 'fig_handle' , gca);

% % if nargin < 3;
% % 	fig_handle = gca;
% % end;

if ticks(1)<ticks(end)
	ylim(fig_handle,[ticks(1) ticks(end)]);
end

% find the zero crossing
% f0		= find(ticks==0);           % find 0 value.
f0		= find(abs(ticks)<5*eps);
s22f	= ['%2.' num2str(d) 'f'];  %

% this is the new y-ticke lable with formatting s22f
ffyy = cellstr(num2str(ticks',s22f));

% add normal 0 lable for 0 axis
if ~isempty(f0)
	ffyy{f0} = '0';
end

set(fig_handle,'YTick'			,ticks);
set(fig_handle,'YTicklabel'	,ffyy)
set(fig_handle,'FontSize'		,FNS);

% % if nargin == 3;
% setytick(d,fig_handle);
% % end;

if nargout > 0
	varargout{1} = fig_handle;
end


end
