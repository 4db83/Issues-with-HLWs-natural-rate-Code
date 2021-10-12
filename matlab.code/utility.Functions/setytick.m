function [] = setytick(d,FNS,fig_handle)
% set the yticklabel to d digits. 
% =======================================================================================
% input
%			d:				scalar, specifies the number of digits, default is 2.
%	fig_handle:		scalar, figure handle, default is gca
%
% call as: 
%								setytick
%								or pass in d to specify number of digits and fig_handle to specify which Figure handle.
% =======================================================================================
% db 14.06.2012
% modfified on 10.07.2013
% modfified on 10.06.2016 (for R2015b)
% =======================================================================================

% SetDefaultValue(2, 'd'	, 2);
SetDefaultValue(2, 'FNS', 10);
SetDefaultValue(3, 'fig_handle' , gca);

% % fig_handle = gca; d = 2; plot(jmp_w/5- 0.01);

ytck		= get(fig_handle,'YTick');
ytcklbl = get(fig_handle,'YTickLabel');

f0      = find(ytck==0);           % find 0 value.
s22f    = ['%2.' num2str(d) 'f'];  %


% this is the new y-ticke lable with formatting s22f
ffyy = cellstr(num2str(ytck',s22f));

% add normal 0 lable for 0 axis
if ~isempty(f0)
	ffyy{f0} = '0';
end;

Matversion = version;
if str2double(Matversion(13:16)) < 2015
	%set(fig_handle,'YTickLabel',ffyy,'YTick',ytck,'FontName','Times New Roman');
	set(fig_handle,'YTickLabel',ffyy,'YTick',ytck,'FontName','Palatino');
else 
	set(fig_handle,'YTickLabel',ffyy,'YTick',ytck,'FontName','Times New Roman');
end

set(fig_handle,'FontSize'		,FNS);

