function s = latexmat(rowNames, M, NaN_replace, varargin)
% function s = latexmat(rowNames, M, NaN_replace, remove_last_, varargin)
% ---------------------------------------------------------------------------------------------------
%latexmat   Print a matrix in LaTeX tabular format.
% 		% ---------------------------------------------------------------------------------------------------
% 		% now S2 MUE results
% 		% ---------------------------------------------------------------------------------------------------
% 		rowNames = {'$L$' ,'MW ' ,'EW ' ,'QLR'};				
% 		% Lambda z
% 		Lz_TV = [	struct2array(L2z_RFile    )' , ...
% 							struct2array(L2z_bl       )' , ...
% 							struct2array(L2z_g        )'	, ...
% 							struct2array(L2z_M0g      )' ];
% 
% 		Lz_C  = [	struct2array(Chow2_bl     )' , ...
% 							struct2array(Chow2_g      )'	, ...
% 							struct2array(Chow2_M0g    )' ];
% 
% 		LT_1a = latexmat(rowNames, [Lz_TV inf(4,1) Lz_C],[], 'nomath', '% 4.8f');
% 		fid		= fopen(['../../table.input/Table_MUEa_' COUNTRY '.tex'],'wt'); fprintf(fid, '%s\n', LT_1a); fclose(fid);
% 
% 		% Fstatistics
% 		F_TV = [	struct2array(bstat2_RFile )' , ...
% 							struct2array(bstat2_TV_bl )' , ...
% 							struct2array(bstat2_TV_g  )'	, ...
% 							struct2array(bstat2_TV_M0g)' ];
% 
% 		F_C  = [	struct2array(bstat2_C_bl  )' , ...
% 							struct2array(bstat2_C_g   )'	, ...
% 							struct2array(bstat2_C_M0g )' ];
% 
% 		LT_1b = latexmat(rowNames, [F_TV inf(4,1) F_C],[], 'nomath', '% 4.8f');
% 		fid		= fopen(['../../table.input/Table_MUEb_' COUNTRY '.tex'],'wt'); fprintf(fid, '%s\n', LT_1b); fclose(fid);
% 
% 		% ---------------------------------------------------------------------------------------------------
% 		% S3 MLE results
% 		% ---------------------------------------------------------------------------------------------------
% 		rowNames = {'$\hsp[3]a_{y,1}  $               ','$\hsp[3]a_{y,2}  $               ','$\hsp[3]a_{r}    $               ','$\hsp[3]b_{\pi } $               ','$\hsp[3]b_{y}    $               ','$\hsp[3]\sigma _{\tilde{y}}$     ','$\hsp[3]\sigma _{\pi }     $     ','$\hsp[3]\sigma _{y^{\ast }}$     ','$\hsp[3]\sigma _{g}$ {(implied)} ','$\hsp[3]\sigma _{z}$ {(implied)} ','$\hsp[3]\lambda_g  $ {(implied)} ','$\hsp[3]\lambda_z  $ {(implied)} '};
% 
% 		TS3		= tabS3{:, 3:end};
% 		TS3a	= TS3;
% 		TS3a(end-2,:) = [];
% 		TS3b	= TS3(end-2,:);
% 
% 		LT_S31a = latexmat(rowNames, [TS3a], '{\hsp[8]---}', 'nomath', '% 4.10f');
% 		fid	= fopen(['../../table.input/Table_S3a_' COUNTRY '.tex'],'wt'); fprintf(fid, '%s\n', LT_S31a); fclose(fid);
% 
% 		LT_S31b = latexmat({'{Log-likelihood}'}, [TS3b], [], 'nomath', '% 4.10f');
% 		fid	= fopen(['../../table.input/Table_S3b_' COUNTRY '.tex'],'wt'); fprintf(fid, '%s\n', LT_S31b); fclose(fid);
% ---------------------------------------------------------------------------------------------------
%   latexmat(M) prints out the numeric matrix M in a LaTeX tabular
%   format. The '&' character appears between entries in a row, '\\'
%   is appended to the ends of rows, and each entry is set in math
%   mode. Complex numbers are understood, and exponentials will be
%   converted to a suitable format.
%
%   LATEX(M,'nomath') does not include the $$ needed to put each 
%   entry in math mode (e.g., for use with the amsmath matrix modes).
%   
%   LATEX(M,FMT) uses a format specifier FMT of the SPRINTF type for
%   each entry.
%   
%   LATEX(M,FMT1,FMT2,...) works through the given format specifiers
%   on each row of M. If fewer are given than the column size of M,
%   the last is used repeatedly for the rest of the row.
%   
%   S = LATEX(M,...) does not display output but returns a character
%   array S.
%   
%   Examples:
%     latexmat( magic(4) )
%     latexmat( magic(4), '%i', 'nomath' )
%     latexmat( magic(4), '%i', '%.2f' )
%   
%   See also SPRINTF, SYM/LATEX.

%   Copyright 2002 by Toby Driscoll. Last updated 12/06/02.

if ~isa(M,'double')
  error('Works only for arrays of numbers.')
elseif ndims(M) > 2
  error('Works only for 2D arrays.')
end

if nargin < 2
  fmt = {'%#.4g'};
  mathstr = '$';
% 	veryLastentry = 'yep';
else
  fmt = varargin;
  idx = strmatch('nomath',fmt);
% 	veryLastentry = char(fmt(end));
	
  if isempty(idx)
    mathstr = '$';
  else  
    mathstr = '';
    fmt = fmt([1:idx-1 idx+1:end]);
    if isempty(fmt), fmt = {'%#.8g'}; end
	end 
% 	veryLastentry = 'AOUT';
end

% Extend the format specifiers.
[m,n] = size(M);
if n > length(fmt)
  [fmt{end:n}] = deal(fmt{end});
end

% create \begin tabular string
begt = ['\n\\begin{tabular}{' repmat('c',1,size(M,2)) '}\n\n'];
begt = ['\n'];

% Create one format for a row.
rowfmt = ' & ';
for p = 1:n
  if p == 1
		%fprintf('\n\\begin{tabular} \n\n');
		fprintf(begt);
	end
	
	% Remove blanks.
  thisfmt = deblank(fmt{p});
	
  % Add on imaginary part if needed.
  if ~isreal(M(:,p)) 
    % Use the same format as for the real part, but force a + sign for
    % positive numbers. 
    ifmt = thisfmt;
    j = findstr(ifmt,'%');
    if ~any(strcmp(ifmt(j+1),['-';'+';' ';'#']))
      ifmt = [ifmt(1:j) '+' ifmt(j+1:end)];
    else
      ifmt(j+1) = '+';
    end
    ifmt = [ifmt 'i'];
    thisfmt = [thisfmt ifmt];
  end

  % Add to row.
%   rowfmt = [rowfmt mathstr thisfmt mathstr ' & '];
	rowfmt = [ rowfmt mathstr thisfmt mathstr ' & '];
end

% After last column, remove column separator and put in newline.
rowfmt(end-1:end) = [];
rowfmt = [rowfmt '\\\\\n'];

% Use it.
% A = M.';
% A = M(1,:)';
% if isreal(M)
%   S = sprintf(rowfmt,A);
% else
%   S = sprintf(rowfmt,[real(A(:)) imag(A(:))].');
% end

S = [];
for jj = 1:m
	S = [S sprintf(['%s' rowfmt], rowNames{jj}, M(jj,:)')];
end

% Remove extraneous imaginary part for real entries.
if ~isreal(M)
  zi = sprintf(ifmt,0);
  S = strrep(S,zi,blanks(length(zi)));
end

% Remove NaNs.
%S = strrep(S,'$NaN$','--');        %% commented out by Daniel Buncic

if isempty(NaN_replace)
	S = strrep(S,'NaN', '{---}');	
else
	S = strrep(S,'NaN', NaN_replace);		
end

S = strrep(S,'Inf', '');							

% Convert 'e' exponents to LaTeX form. This is probably really slow, but
% what can you do without regular expressions?
% S = strrep(S,'e','E');
% ex = min(findstr(S,'E'));
% while ~isempty(ex)
%   % Find first non-digit character. Where is ISDIGIT?
%   j = ex+2;
%   while ~isempty(str2num(S(j))) & ~strcmp(S(j),'i')
%     j = j+1;
%   end
% 
%   % This strips off leading '+' and zeros.
%   num = sprintf('%i',str2num(S(ex+1:j-1)));
%   
%   ee = ['\times 10^{' num '}'];
%   S = [S(1:ex-1) ee S(j:end)];
%   
%   ex = ex + min(findstr(S(ex+1:end),'E'));
% end

% --------------------------------------------------------------------------------------------------
% For good form, remove that last '\\'. end of column
% --------------------------------------------------------------------------------------------------
% S(end-2:end-1) = '';              %% commented out by Daniel Buncic
% if remove_last_ == 1
% 	S(end-2:end-1) = '';              %% commented out by Daniel Buncic
% end																	
																	
% S(end-2:end-1) = veryLastentry;
% B = [S veryLastentry];
% disp(B)
% C = S;
% C(end-2:end-1) = 'asdfasfa'
																		
% Display or output?
if nargout==0
  disp(S)
else
  s = S;
end
% fprintf('\\end{tabular} \n');
fprintf('');

