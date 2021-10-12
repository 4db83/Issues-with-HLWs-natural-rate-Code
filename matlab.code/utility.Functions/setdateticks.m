% function setdateticks(dates,Width,dateformat,date_space,fig_handle)
function step_sizes_out = setdateticks(dates,Width,dateformat,FNTS,date_space,fig_handle)
% SET OR ADJUST the yticklables and ticks at desired location vector ticks to dates.
% ----------------------------------------------------------------------------------------------
% CALL AS:		setdateticks(dates,Width,dateformat,FNTS,date_space,fig_handle)
% 
% SET WIDTH equal to (length(dates)-1)/x = 35 to 45, with x as close to an
% integer as possible (with date_space = {0,3} to show first and last date % entry).
% 35 to 45 are the number of ticks that are shown
% ----------------------------------------------------------------------------------------------
% db 14.06.2012
% modfified on 10.07.2013
% ----------------------------------------------------------------------------------------------

% get defaulf font size
fonts0 = get(gca,'FontSize');

SetDefaultValue(2, 'Width', 15);
SetDefaultValue(3, 'dateformat', 0);
SetDefaultValue(4, 'FNTS', fonts0);
SetDefaultValue(5, 'date_space', 0);
SetDefaultValue(6, 'fig_handle', gca);

% Font name
FName = 'Times New Roman';

% if datetime convert to datenum
if isdatetime(dates)
	dates = datenum(dates);
end

% mod the format of the date string showing
if ~isstring(dateformat)
	if dateformat == 0
		dateformat = 'yyyy:QQ';
	elseif dateformat == 1
		dateformat = 'mmm-yyyy';
	elseif dateformat == 2
		dateformat = 'dd.mm.yyyy';
	end
end
	
% get the x-axis (DATES) increments
x_ticks0 = get(gca,'XTick');
% check if the minimum is less than 1000 (dates are usuall from 700000 onwards
Imin_x = min(x_ticks0) < datenum('01-Jan-1600');

if Imin_x
	% this is the case in a normal setdatesticks call
	TT = size(dates,1);
	T0 = 1;
else 
	% for calls of setdatesticks with Financial time series objects (fts).
	TT	= x_ticks0(end);
	T0	= x_ticks0(1);
end

% DETERMIN INTEGER DATE_STEP_WIDTH;
step_range				= (2:TT-2)';
integer_date_with = (length(dates)-1)./step_range;
% suggested integer step width
date_step_width = step_range(integer_date_with == round(integer_date_with));

tmp_ = (length(dates)-1)./Width;
is_integer = tmp_ == round(tmp_);

if date_space == 0
	if is_integer
		DT	= (T0:Width:TT);
	else
		DT	= floor(linspace(T0,TT,TT/Width));
	end
elseif date_space == 2
	DT		= (TT:-Width:T0)';
	DT		= flipud(DT);
elseif date_space == 1
	DT		= (T0:Width:TT)';
elseif date_space == 3
	DT		= (T0:Width:TT)';
%	DT		= flipud(DT);
%	DT(end) = T;				% take out the last element
% DT	= [DT; TT];			% use 1 as the first entry here.
%	DT	= [(T:-W:1)'];
end

% diff(DT)

% if we put in a char vector of dates
if ~isfloat(dates)
	set(fig_handle,'XTick'			, DT);
	set(fig_handle,'XTicklabel'	, dates(DT,:) );
	set(fig_handle,'FontName'		, FName);
else
	% else do the normal plotting
	set(fig_handle,'XTick'			, DT);
	if Imin_x
		set(fig_handle,'XTicklabel'	, datestr(dates(DT), dateformat))
	else
		set(fig_handle,'XTicklabel'	, datestr(DT, dateformat))
	end
	set(fig_handle,'FontName'		, FName)
end

set(gca,'FontName',FName,'FontSize',FNTS);     

tickshrink(.5);
datelim(2);

if nargout > 0
	step_sizes_out = date_step_width';
end


%%%%%%%%%%%%%O%%%%%%%%L%%%%%%%D%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THIS IS THE OLD SETDATETICKS FUNCTION
%%%%%%%%%%%%%O%%%%%%%%L%%%%%%%D%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % SetDefaultValue(2, 'dateformat', 0);
% % % SetDefaultValue(3, 'W', 40);
% % % SetDefaultValue(4, 'date_space', 0);
% % % SetDefaultValue(5, 'fig_handle', gca);
% % % 
% % % % mod the format of the date string showing
% % % if ~isstring(dateformat)
% % % 	if dateformat == 0
% % % 		dateformat = 'yyyy:QQ';
% % % 	elseif dateformat == 1
% % % 		dateformat = 'mmm-yyyy';
% % % 	elseif dateformat == 2
% % % 		dateformat = 'dd.mm.yyyy';
% % % 	end
% % % end
% % % 	
% % % T = size(dates,1);
% % % FName = 'Times New Roman';
% % % 
% % % if date_space == 0
% % % 	DT		= round(linspace(1,T,W));
% % % elseif date_space == 2
% % % 	DT		= [(T:-W:1)'];
% % % 	DT		= flipud(DT);
% % % elseif date_space == 1
% % % 	DT		= [(1:W:T)'];
% % % elseif date_space == 3
% % % %	DT		= [(T:-W:1)'];
% % % 	DT		= [(1:W:T)'];
% % % %	DT		= flipud(DT);
% % % %	DT(end) = T;						% take out the last element
% % % 	DT			= [DT; T];			% use 1 as the first entry here.
% % % end;
% % % 
% % % % if we put in a char vector of dates
% % % if ischar(dates)
% % % 	set(fig_handle,'XTick'			, DT);
% % % 	set(fig_handle,'XTicklabel'	, dates(DT,:) );
% % % 	set(fig_handle,'FontName'		, FName);
% % % 
% % % % else do the normal plotting
% % % else
% % % 	%DT = round(linspace(1,T,W))
% % % 	%diff(DT)
% % % 	% size(DT)
% % % 
% % % 	set(fig_handle,'XTick'			, DT);
% % % 	set(fig_handle,'XTicklabel'	, datestr(dates(DT), dateformat))
% % % 	set(fig_handle,'FontName'		, FName)
% % % end
%%%%%%%%%%%%%O%%%%%%%%L%%%%%%%D%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%