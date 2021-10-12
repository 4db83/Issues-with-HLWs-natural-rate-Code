function setplot(figdim,fntsize,y_digits,axLineWidth,fighandle) 
%{F: specifies plots dimensions, fontsize and fontname
%===============================================================================
% Makes it easy to view strings in excel which are actually numbers without 
% converting them to numbers`
%-------------------------------------------------------------------------------
% 	USGAGE:	(1) setplot([height width], Fontsize)
%						(2) setplot([leftpos height width], Fontsize) (rightpos at 0)
%						(2) setplot([leftpos rightpos height width], Fontsize) 
%           [leftpos rightpos] is optional and will only be called if one 
%           needs to center the plot. Default font size is 10, fontname is Palation.
%-------------------------------------------------------------------------------
% 	INPUT : 
%	  figdim		=  (1x2) matrix with [height width] dimension of plot
%   fntsize		=  scalar, fontsize. default fontsize is 10.
%		fighandle =	 sclar, figure handle. 
% 	OUTPUT:       
%	  zero arguments: adjusted plot.
%===============================================================================
% 	NOTES :   Always call it at the end of the plotting commands, after hold off;
%-------------------------------------------------------------------------------
% Created :		27.07.2013.
% Modified:		27.07.2013.
% Copyleft:		Daniel Buncic.
%------------------------------------------------------------------------------%}
widh_def = .86;
SetDefaultValue(1, 'figdim'			, [widh_def .2]);

FNS_0 = get(gca,'FontSize');
FNS0	= FNS_0 - 0;
% since 2021a, matlab seesm to rotate the Xtick labels. -> remove that 
set(gca,'XTickLabelRotation',0);

xd = size(figdim);

if ~(max(xd) == 1)
	Big1 = figdim > 1;
	if any(Big1) > 0
		figdim = figdim./10;
	end
end

% SetDefaultValue(3, 'BoxLineWidth', 1.25);

SetDefaultValue(2, 'fntsize'		, FNS0);
SetDefaultValue(3, 'y_digits'		, 1);
SetDefaultValue(4, 'axLineWidth', 6/5);
SetDefaultValue(5, 'fighandle'	, gca);

% ff = get(fighandle,'Position');

% added the font name conversion now to print in TNR as opposed to Palation!
FName = 'Times New Roman';
top_	= .45;
left_ = .07;
% left_ = .05;

% if nargin < 2
  % if only one input argument is given
  if max(xd) == 2
    % sets width and height, takes the centered location [.1 .2]
    set(fighandle,'Position',[[left_ top_] figdim],'FontSize',fntsize,'FontName',FName);
  elseif max(xd) == 3
    % sets width and height as well as the location
    set(fighandle,'Position',[left_ figdim(1:3)],'FontSize',fntsize,'FontName',FName);
	elseif max(xd) == 4
    % sets width and height as well as the location
    set(fighandle,'Position',figdim,'FontSize',fntsize,'FontName',FName);
  elseif max(xd) == 1
    % only scalar input means that only the font size is to be set.
%     fntsize = figdim;
%     set(fighandle,'FontSize',fntsize,'FontName',FName);
		set(fighandle,'Position',[[left_ top_ widh_def] figdim],'FontSize',fntsize,'FontName',FName);
  else 
    disp('Error in plot diminsion');
	end

ytcks = get(gca,'YTick');
	
% set the ylim digits 
if ~isempty(y_digits)
	setyticklabels(ytcks)
	setytick(y_digits,fntsize,fighandle)
% 	setyticklabels(13:.5:16);
end
	
% get axes handle and increase the thickness of lnes and choose grid value for color
	ax = gca; 
% 	ax.GridAlpha = .15;  
 	ax.LineWidth = axLineWidth;
	
	% else
  % if two input arguments are given then same as above except for last bit
% %   pst0 = get(fighandle,'Position');
% % 	pst0(1:2) = [.05 .55];
% %   if max(xd) == 2
% %     set(fighandle,'Position',[pst0(1:2) figdim],'FontSize',fntsize,'FontName',FName);
% %   elseif max(xd) == 3
% %     set(fighandle,'Position',[figdim(1) ff(2) figdim(2:3)],'FontSize',fntsize,'FontName',FName);
% % 	elseif max(xd) == 4
% %     set(fighandle,'Position',figdim,'FontSize',fntsize,'FontName',FName);    
% %   else 
% %     disp('Error in plot diminsion');
% %   end
% % end

% set axis laver over the top of the plots.
% set(gca, 'Layer','top')	


%set(gca,'Position',[[.1 .2] figdim],'FontSize',fntsize,'FontName',FName);






% % % if nargin < 2
% % %   % if only one input argument is given
% % %   if max(xd) == 2
% % %     figdim = figdim;
% % %     fntsize = 10;     % fontsize is 10 as default value
% % %     % sets width and height, takes the centered location [.1 .2]
% % %     set(fighandle,'Position',[[.055 .2] figdim],'FontSize',fntsize,'FontName',FName);
% % %   elseif max(xd) == 3
% % %     figdim = figdim;
% % %     fntsize = 10;     % fontsize is 10 as default value
% % %     % sets width and height as well as the location
% % %     set(fighandle,'Position',[figdim(1) 0 figdim(2:3)],'FontSize',fntsize,'FontName',FName);
% % % 	elseif max(xd) == 4
% % %     figdim = figdim;
% % %     fntsize = 10;     % fontsize is 10 as default value
% % %     % sets width and height as well as the location
% % %     set(fighandle,'Position',figdim,'FontSize',fntsize,'FontName',FName);
% % %   elseif max(xd) == 1
% % %     % only scalar input means that only the font size is to be set.
% % %     fntsize = figdim;
% % %     set(fighandle,'FontSize',fntsize,'FontName',FName);
% % %   else 
% % %     disp('Error in plot diminsion');
% % %   end
% % % else 
% % %   % if two input arguments are given then same as above except for last bit
% % %   pst0 = get(fighandle,'Position');
% % % 	pst0(1:2) = [.05 .55];
% % %   if max(xd) == 2
% % %     set(fighandle,'Position',[pst0(1:2) figdim],'FontSize',fntsize,'FontName',FName);
% % %   elseif max(xd) == 3
% % %     set(fighandle,'Position',[figdim(1) ff(2) figdim(2:3)],'FontSize',fntsize,'FontName',FName);
% % % 	elseif max(xd) == 4
% % %     set(fighandle,'Position',figdim,'FontSize',fntsize,'FontName',FName);    
% % %   else 
% % %     disp('Error in plot diminsion');
% % %   end
% % % end
