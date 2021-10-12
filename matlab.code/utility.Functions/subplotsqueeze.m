function subplotsqueeze(hFig, nF)
%{F Stretch width and height (and column gap) of all subplots in a figure window
% subplotsqueeze(H, F) will stretch subplots in figure with handle H by the
% proportional factor F.
% Input: 
% 				hFig : figure handle
% 				nF   : (3x1) vector of [width height column_gap]; all in [0 1] interval
% 
% Examples (1 input only)
%  subplotsqueeze(gcf, 1.2) will expand all axes by 20%
%  subplotsqueeze(gcf, 0.8) will contract all axes by 20%
% 
% Examples (2 inputs)
% subplotsqueeze(gcf, [0.5 1.1]) narrower by 50%, higher by 10%
% subplotsqueeze(gcf, [1.2 0.5]) wider by 20%		, shorter by 50%
% 
% Examples (3 inputs)
% subplotsqueeze(gcf, [0.5 1.1 0.5]) narrower by 50%, higher by 10% and reduce column gap by 0.5.
% this works best for 6,2 subplots or so.
% ========================================================================================
% 	NOTES :   Expansion and contraction is equal in both directions and axes remain
%             centered on their current locations.                                 
% ----------------------------------------------------------------------------------------
% Created :		12.08.2014.
% Copyleft:		Daniel Buncic.
% ----------------------------------------------------------------------------------------%}


dim = length(nF);

if dim == 1
	nF = [nF; nF];
end;

hAx = findobj(hFig, 'type', 'axes');
% % for h = 1:length(hAx)
% %     vCurrPos = get(hAx(h), 'position'); % current position
% %     set(hAx(h), 'position', ...
% % 			 (vCurrPos.*[1 1 nF(1) nF(2)]) ... 
% % 			-[0 vCurrPos(4)*(nF(2)-1)/2 0 0]);
% % end

for h = 1:length(hAx)
  vCurrPos = get(hAx(h), 'position'); % current position
  set(hAx(h), 'position', (vCurrPos.*[1 1 nF(1) nF(2)]) - [vCurrPos(3)*(nF(1)-1)/2 vCurrPos(4)*(nF(2)-1)/2 0 0]);
	if dim == 3
		set(hAx(h), 'position', vCurrPos.*[nF(3) 1 nF(1) nF(2)]);		
	end;
	
end

return

% for h = 1:length(hAx)
%   vCurrPos = get(hAx(h), 'position'); % current position
%   set(hAx(h), 'position', vCurrPos.*[nF(1) 1 nF(2) nF(3)]);
% end
% return 