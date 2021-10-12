function aout = plugin(A,b,i)
% plug-in vector b after positoin i in Matrix A;
% output Aout = [	A(1:i-1,:);
% 								b;
%									A(1:i-1,:)];
% ----------------------------------------------------------------------------------------

[ar,ac] = size(A);
[br,bc] = size(b);

if ac > bc
	b = repmat(b,br,ac);
end

aout = [ A(1:i-1,:);...
				 b;...
				 A(i:end,:)];
end