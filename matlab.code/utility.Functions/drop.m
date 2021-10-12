function Xdrop= drop(x,k1,k2)
SetDefaultValue(3,'k2',0);
% drops the first k values
Xdrop		= trimr(x,k1,k2);
%Xdrop = x(k+1:end,:);
