function SetDefaultValue(position, argName, defaultValue)
% Call as:			SetDefaultValue(1, 'x', 10);
% Initialise a missing or empty value in the caller function.
% 
% SETDEFAULTVALUE(POSITION, ARGNAME, DEFAULTVALUE) checks to see if the
% argument named ARGNAME in position POSITION of the caller function is
% missing or empty, and if so, assigns it the value DEFAULTVALUE.
% 
% Example:
% function x = TheCaller(x)
% SetDefaultValue(1, 'x', 10);
% end
% TheCaller()    % 10
% TheCaller([])  % 10
% TheCaller(99)  % 99
% 
% $Author: Richie Cotton $  $Date: 2010/03/23 $
% MODIFIED BY DANIEL TO ALLOW FOR USE WITH MATLAB R2015b ACCORDING TO THE HINT GIVEN HERE:
% http://www.mathworks.com/matlabcentral/fileexchange/27056-set-default-values. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In order to set default values for variables, I find the start of my functions littered with
% 
% if nargin < 1 || isempty(x) 
% x = 1; 
% end
% 
% if nargin < 2 || isempty(y) 
% y = 3; 
% end 
% etc.
% 
% This is pretty ugly, so I've created a wrapper to prettify it. Honestly, it's so simple that 
% I nearly didn't upload this, but it does make your functions cleaner. 
% Now the above is transformed to
% 
% SetDefaultValue(1, 'x', 1); 
% SetDefaultValue(2, 'y', 3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% if evalin('caller', 'nargin') < position || ...
if ~evalin('caller', ['exist(''' argName ''', ''var'')']) || ...
      isempty(evalin('caller', argName))
   assignin('caller', argName, defaultValue);
end

end
