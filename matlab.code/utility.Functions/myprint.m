function myprint(y,RowNames,ColNames,xlsout,FMT,WIDTH,SEP)
% PURPOSE: print an (nobs x nvar) matrix in formatted form
% CALL AS: function myprint(y,RowNames,ColNames,xlsout,FMT,WIDTH,SEP)
% 
% ---------------------------------------------------------------------------------------------------
%  USAGE:    mprint(y,RowNames,ColNames) 
% 
% 	y					= (Nrow x Ncols) matrix (or vector) to be printed
% 	RowNames 	= Cell of dimension Nrow or (Nrow + 1), ie.
% 	RowNames 	= {'R1';'R2'; ... } or RowNames = {'RowHeader';'R1';'R2'; ... }.
% 	ColNames 	= Cell of dimension Ncol, ie.
% 	ColNames 	= {'C1';'C2'; ... }.
% 	FMT				= Optional format of number printing, default is '%14.6f'.
% 	WIDTH			= Optional width of lines to be printed. default is 200.
% 	SEP				= Optional seperation symbol, default is '-'.
% 
% ---------------------------------------------------------------------------------------------------
% DB: 11.08.2015.
% NOTES: - this is a modified version of the orginal print file from LP toolbox
% ---------------------------------------------------------------------------------------------------

[Nr,Nc] = size(y);

% set default values
ColName_default	= repmat({' '},1,Nc);

SetDefaultValue(3,'ColNames',ColName_default);
SetDefaultValue(4,'xlsout'	,0);
SetDefaultValue(5,'FMT'		,'%14.10f');
SetDefaultValue(6,'WIDTH'	,200);
SetDefaultValue(7,'SEP'		,'-');

if nargin < 2
% 	info0.fmt = '%14.6f';
	info0.fmt = FMT;
	myprint_struct(y,info0);	
else 
	% if second input is a structure, then 
	if isstruct(RowNames)
		RowNames.width	= WIDTH;
		myprint_struct(y,RowNames);	
	else
		if ~isempty(RowNames)
			if ~iscell(RowNames)
				RowNames = cellstr(RowNames);
			end
		end
% 		ColNames = [repmat({' '},1,Nc-length(ColNames)); ColNames]
		info_in.rnames	= RowNames;
		if isempty(RowNames) 
			info_in.rnames	= cellstr(repmat({''},Nr,1)); 
		end
		info_in.cnames	= [repmat({' '},1,Nc-length(ColNames)); ColNames];
		info_in.fmt			= FMT;
		info_in.width		= WIDTH;
		info_in.sep			= SEP;
		myprint_struct(y,info_in);
	end
end

if nargin > 1
	% make sure ColNames is 1xNc vector and not Ncx1
	[Cr,Cc] = size(ColNames);
	if Cr > Cc
		ColNames = ColNames';
	end
	% make sure RowNames is Nrx1 vector and not 1xNr 
	[Rr,Rc] = size(RowNames);
	if Rc > Rr
		RowNames = RowNames';
	end
end

if (xlsout~=0)
	if ischar(xlsout)
		xlsout_name = xlsout;
	else
		xlsout_name = 'myprint_out.xls';
	end
	% now write to xls.
	xlswrite(xlsout_name,ColNames,	'Sheet1','B1');
	if length(RowNames) > Nr
		xlswrite(xlsout_name,RowNames,'Sheet1','A1');
	else
		xlswrite(xlsout_name,RowNames,'Sheet1','A2');
	end
	xlswrite(xlsout_name,y			 ,	'Sheet1','B2');
end


































