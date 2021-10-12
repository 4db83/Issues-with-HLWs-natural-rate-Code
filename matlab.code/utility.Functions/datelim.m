function [] = datelim(date_limits)
% date axis limits padding funciton
% date_end is padding towards at the end.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if length(date_limits) < 2
	date_beg = date_limits;
	date_end = date_limits;
else 
	date_beg = date_limits(1);
	date_end = date_limits(2);
end

get_date_limits = get(gca,'XTick');

xlim([get_date_limits(1)-date_beg get_date_limits(end) + date_end])