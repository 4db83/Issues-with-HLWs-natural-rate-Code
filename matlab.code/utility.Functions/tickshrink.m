function [] = tickshrink(fct,fig_hndl)
%% F: reduces the lenght of the ticks to small/larger values by specifying a shrink/scale factor fct.
%	 USAGE: 
%					tickshrink(fct,figh_in)
%					
%	 INPUT: fct				= shrink/scale factor to increase reduce the scale of the ticks.
%					fig_hndl	= figure handle.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%}

if nargin < 2
	fig_hndl = gca;
end

current_tick_length = get(fig_hndl,'TickLength');

set(fig_hndl,'TickLength',fct.*current_tick_length);


