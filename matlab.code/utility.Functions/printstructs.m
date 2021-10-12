function varargout = printstructs(varargin)
% Function: print output from mulitple input structures to screen
% CALL AS: printstructs(struct1,stuctur2,...,structEND)
% varargin input:
% struct1, stuctur2, etc all having the same fields.
% optional is the printing format styel, default is '%10.6f';
% this function supersedes the printstruct function for single struture input
% db
% 20.05.2015
% --------------------------------------------------------------------------------------------------
% SetDefaultValue(2,'fmt','%8.8f');

% get name of input as string
Nin = length(varargin);

% CHECK IF LAST INPUT ARGUMENT IS A STRUCTURE OR NOT
last_varargin = varargin{Nin};
isAstructure	= isa(last_varargin,'struct');

% old1 = isstr(last_varargin);
% old1 = ischar(last_varargin);

XLS_OUTPUT_NAME = 0;

% check if we inlude the format as an input
if isAstructure 
	frmt = '%14.8f';
else
	lastIn = varargin{Nin};
	if ischar(lastIn)
		if contains(lastIn,'.xls')
			XLS_OUTPUT_NAME = lastIn;
			frmt = '%18.14f';
		else
			frmt = lastIn;
			XLS_OUTPUT_NAME = 0;
		end
	else
		frmt = ['%14.' num2str(lastIn) 'f'];
	end
end

% allocate some space
XX = [];

% strct_count = 0;
% loop through input arguments that are structures
for ii = 1:(Nin - ~isAstructure)
	names_str{ii} = inputname(ii);
	is_struct = isa(varargin{ii},'struct');
	if is_struct
		XX = [XX struct2array(varargin{ii})'];
	else
		XX = [XX varargin{ii}];
	end
	I_struct(ii) = is_struct;
end

row_names = fieldnames(varargin{find(I_struct,1,'first')});


% xlsout		= 0;
sepWidth	= 97;
sepWidth	= 140;
% disp(lastIn)

% myprint(XX,row_names,names_str,xlsout,frmt,sepWidth)
% xlsout = Xls_output_name;
myprint(XX,row_names,names_str, XLS_OUTPUT_NAME, frmt, sepWidth)

% myprint(XX,in)

if nargout > 0
	varargout{1}  = XX;
end	



% in.fmt		= fmt;
% in.rnames = row_names;
% in.cnames = names_str;
% in.width	= 90;

% % % get structure inputs.
% % m1 = fieldnames(varargin);
% % a1 = struct2array(varargin)';
% % 
% % xout = [ repmat('  ',size(a1)) char(m1) repmat(' :   ',size(a1)) num2str(a1,fmt)];
% % 
% % [Xr,Xc] = size(xout);
% % 
% % C = ' ';
% % str_out = [repmat(C,1,Xc-size(names_str,2)-1) names_str];
% % disp(str_out)
% % disp(xout)




