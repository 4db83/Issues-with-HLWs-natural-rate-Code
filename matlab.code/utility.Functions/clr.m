function clr_out = clr(indx)

Nindx = length(indx);

if Nindx == 1

switch indx
	case 'b' ; indx = 1;	% blue
	case 'r' ; indx = 2;	% red
	case 'o' ; indx = 3;	% orange
	case 'p' ; indx = 4;	% purple
	case 'g' ; indx = 5;	% green
	case 'c' ; indx = 6;	% cyan
	case 'm' ; indx = 7;	% morone
end

end

	
c0=[	 	 0    0.4470    0.7500;
    0.8500    0.2250    0.0580;							% 0.8500    0.3250    0.0980;
    0.9290    0.6940    0.1250;
    0.4940    0.1840    0.5560;
    0.4660    0.6740    0.1880;
    0.3010    0.7450    0.9330;
    0.6350    0.0780    0.1840;
         0    0.7500    0.7500;
    0.7500         0    0.7500;
    0.7500    0.7500         0;
		1							 0         0;
		1							 1         1;
    0.2500    0.2500    0.2500;];
    
    
clr_out = c0(indx,:);

