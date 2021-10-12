%%Function: Converts structure to matrix.
%{==================================================================
% USAGE:	[a]=struct2mat(structure);
%-------------------------------------------------------------------
% INPUT:    
%			structure.
%        	        
% OUTPUT:   
%			a = matrix.
%	
% NOTES:
%			None.
%			
%-------------------------------------------------------------------
% created by Daniel Buncic
% Date: 26/07/2006.
%-------------------------------------------------------------------%}
function [a]=struct2mat(structure);
a=cell2mat((struct2cell(structure))');
