function [] = add2yaxislabel
% function create a second y-axis label.

ylim_	= get(gca,'YLim');
ylim_ticks = get(gca,'YTick');
ylim_lables_0 = get(gca,'YTickLabel');

% make the zere entry alligned to the right
% ylim_lables(find(ylim_ticks==0)) = {' aa'}

ylim_lables = strrep(ylim_lables_0,' ','');

yyaxis(gca,'right');
ylim(ylim_);
set(gca,'YColor',[0 0 0],'YTickLabel',ylim_lables,'YTick',ylim_ticks);

% % make the zere entry alligned to the right
% ylim_lables(find(ylim_ticks==0)) = {' aa'}
% 
% set(gca,'YColor',[0 0 0],'YTickLabel',ylim_lables,'YTick',ylim_ticks);
% 
% a = 1

