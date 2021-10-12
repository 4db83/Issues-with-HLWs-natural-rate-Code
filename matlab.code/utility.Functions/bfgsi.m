function H = bfgsi(H0,dg,dx,PRINT_ITERS_)
% H = bfgsi(H0,dg,dx)
% dg is previous change in gradient; dx is previous change in x;
% 6/8/93 version that updates inverse hessian instead of hessian
% itself.
% Copyright by Christopher Sims 1996.  This material may be freely
% reproduced and modified.

SetDefaultValue(4, 'PRINT_ITERS_', 1);

if size(dg,2)>1
   dg=dg';
end
if size(dx,2)>1
   dx=dx';
end
Hdg = H0*dg;
dgdx = dg'*dx;
if (abs(dgdx) >1e-12)
	H = H0 + (1+(dg'*Hdg)/dgdx)*(dx*dx')/dgdx - (dx*Hdg'+Hdg*dx')/dgdx;
else
	if PRINT_ITERS_
		disp('bfgs update failed.')
		disp(['|dg| = ' num2str(sqrt(dg'*dg)) '|dx| = ' num2str(sqrt(dx'*dx))]);
		disp(['dg''*dx = ' num2str(dgdx)])
		disp(['|H*dg| = ' num2str(Hdg'*Hdg)])
	end
  H = H0;
end

dir_store = './_csminwel_output/';

% save H.dat H  %% why dat and not mat???
save([dir_store 'H.dat'], 'H');

