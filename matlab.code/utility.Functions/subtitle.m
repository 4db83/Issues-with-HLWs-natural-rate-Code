function varargout = subtitle(titlename,adjst,FontSize,useLaTex)
% call as: subtitle(Name,postition,FN,1)
%F: puts tile (especially in subfigure) below the subfigure like in latex.
% Usage  subtitle('Name of subplot',ymove_dwon_by_howmuch) -1.1 is default, try -1.2 to move down
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get title handle 
title_handle	= title(titlename);
FNS0 = get(title_handle,'FontSize');

SetDefaultValue(2 ,'adjst'		, -1.1);
% set default font size 
SetDefaultValue(3 ,'FontSize'	, FNS0);
% use Latex in subtilte
SetDefaultValue(4 ,'useLaTex'	, 0);

% use Latex if useLaTex == 1
if useLaTex==1
	set(title_handle,'Interpreter','Latex');
elseif useLaTex==2
	set(title_handle,'Interpreter','Tex');
end

% Use normalised units, works better as it does not depend on the y-axis scale.
set(title_handle  ,'Units','normalized');

current_pos		= get(title_handle,'Position');

if length(adjst)>1
	new_pos				= current_pos + [adjst(2) adjst(1) 0];
else
	new_pos				= current_pos + [0 adjst(1) 0];
end

set(title_handle,'Position',new_pos);
set(title_handle,'FontSize',FontSize);
set(title_handle,'FontWeight','Normal','Units','normalized','FontName','Times New Roman');

if nargout > 0
	varargout{1} = title_handle;
end


end

% % if nargin < 3
% % 	FontSize = 9;
% % end;
% % 
% % % title_handle					= title(titlename,'FontSize',FontSize);
% % % remove boldface
% % set(title_handle,'Fontweight','Normal');
% % 
% % 
% % position_handle				= get(title_handle,'Position');
% % new_position_handle		= [position_handle(1) position_handle(2)-ymove position_handle(3)];
% % set(title_handle,'Position',new_position_handle);

