function HLW_factors = make_HLW_factors(KFS_structure)
% F: make the HLW factors. Returns a strucuture with filtered (.att) and smoothed (.atT)
%														[1		2								3		4				]
% HLW factors, in the order [R*, trend growth g,  z, and cycle]
% ---------------------------------------------------------------------------------------------

% MAKE THE 4 SET OF FACTORS FOR PLOTTING [r*=4*g+z 4*g z] NO NEED TO PLOT THE CYCLE HERE
% FILTERED
HLW_factors.att = [ 4*KFS_structure.KFS.att(:,4) + KFS_structure.KFS.att(:,6) ,...	% r*
										4*KFS_structure.KFS.att(:,4)	, ...															% trend growth
											KFS_structure.KFS.att(:,6)	, ...															% z
											KFS_structure.KFS.Y(:,1) - KFS_structure.KFS.att(:,1)		,...	% cycle
									];

% SMOOTHED								
HLW_factors.atT = [ 4*KFS_structure.KFS.atT(:,4) + KFS_structure.KFS.atT(:,6) ,...	% r*
										4*KFS_structure.KFS.atT(:,4)	, ...															% trend growth
											KFS_structure.KFS.atT(:,6)	, ...															% z
											KFS_structure.KFS.Y(:,1) - KFS_structure.KFS.atT(:,1)		,...	% cycle
									];

%



























%EOF
%%	