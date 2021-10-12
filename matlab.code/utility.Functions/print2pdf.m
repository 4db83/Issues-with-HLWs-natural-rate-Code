function [] = print2pdf(filename, dirname, make_eps_figure, make_emf_figure)
% F: 
% USAGE: 
% 						print2pdf(filename,dirname,make_emf_figure)
%		or 
% print2pdf(filename,dirname,make_emf_figure); 
% 
% to create an emf figure director as well, so that they can be included in WORD. If 
% make_emf_figure == 2, then no extra EMF directory is created.
% NOTE: dir name is now the second argument. 
% db 24.02.2017.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


SetDefaultValue(2 ,'dirname'				, './');
SetDefaultValue(3 ,'make_eps_figure', 0);
SetDefaultValue(4 ,'make_emf_figure', 0);

% if nargin < 4; 	make_emf_figure = 0; end
% if nargin < 3; 	make_eps_figure = 0; end

% if nargin < 2
% 	dirname = './';
% else
% 	dirname = ['./' dirname];
% end;

% if nargin < 2
% 	dirname = './';
% else
% 	dirname = [dirname];
% end

% GSversion = 9.14;             % set to the GSversion that is install in C:\Program Files\
% GSversion = 9.20;             % set to the GSversion that is install in C:\Program Files\
% GSversion = 9.26;             % set to the GSversion that is install in C:\Program Files\
% GSversion = 9.52;             % set to the GSversion that is install in C:\Program Files\
% MAKE IT A STRING NOW
GSversion = '9.53.3';           % set to the GSversion that is install in C:\Program Files\

% C:\Program Files\gs\gs9.53.3\bin

file_ 	= char(filename); % make sure it is not a cell input
% GS			= ['C:\Program Files\gs\gs' num2str(GSversion,'%2.2f') '\bin\gswin64c.exe'];
GS			= ['C:\Program Files\gs\gs' GSversion '\bin\gswin64c.exe'];
% NOT NEEDED ANY MORE. 
% % print('-depsc', file_);
% % h = gcf;
% % set(h,'PaperPositionMode','auto');         

% SINCE MATLAB R2016a, NEED DIFFERENT DEFAULT RENDERER (NOT SURE WHY THEY DID THIS, BUT THIS WORKS NOW)
% SET default renderer FROM 'zbuffer' TO 'painters'
set(gcf,'renderer','painters')

% NOW PRINT TO EPS FILE
% first make an eps directory
if make_eps_figure
	if ~(exist([dirname '/eps.figs/'], 'dir')==7) 
			 mkdir([dirname '/eps.figs/']); 	
	end
		print(gcf,'-depsc2', [dirname '/eps.figs/' file_]);
else
		print(gcf,'-depsc2', [dirname '/' file_]);
end

% [dirname '/' file_]

print(gcf,'-depsc2', [dirname '/' file_]);

set(gcf,'PaperPositionMode','auto');         

% PRINT ALSO TO EMW FILE TO BE INCLUDE IN WORD/POWER POINT TYPE DOCUEMNTS AT HIGH RESOLUTION
if make_emf_figure==1
	% make a storage directory first
	if ~(exist([dirname '/emf.figs/'], 'dir')==7)
		   mkdir([dirname '/emf.figs/'])
	end
	print(gcf,'-dmeta', [dirname '/emf.figs/' file_]);
elseif make_emf_figure==2
	print(gcf,'-dmeta', [dirname '/' file_]);
end

if exist(GS)
	% [rslt] = eps2pdf_funct([dirname '/eps.figs/' file_ '.eps'],GS,0);
	if make_eps_figure
		[rslt] = eps2pdf_funct([dirname '/eps.figs/' file_ '.eps'],GS,0);
	else 
		[rslt] = eps2pdf_funct([dirname '/' file_ '.eps'],GS,0);
	end
		
	if rslt== 0
		if ~make_eps_figure 
			% [file_ '.eps']
% 			[dirname '/' file_ '.eps']
			delete([dirname '/' file_ '.eps']);
% 			[dirname '/eps.figs/']
% 			rmdir([dirname './eps.figs/'])
		end
	else
		fprintf('-------------------------------------------------------------------------------\n')
		fprintf(['| Ghostscript ver. ' GSversion ' is NOT installed on your systerm in C:\\Program Files\\ |\n'])
		fprintf('| You can get it from https://www.ghostscript.com/download/                   |\n')
		fprintf('| Here the path is to the 64 bit install in C:\Program Files\gs\              |\n')
		fprintf('| No pdf file was created, only a .eps file                                   |\n');
		fprintf('-------------------------------------------------------------------------------\n')
	end
else
	fprintf('-------------------------------------------------------------------------------------\n');
	fprintf(['| Program is looking for Ghostscript Version ' GSversion ' in C:\\Program Files\\gs\\gs          |\n']);
	fprintf(['| But could not find it. You need to install Ghostscript Version ' GSversion ' on your system |\n']);
	fprintf('| No pdf file was created, only a .eps file which can be manually converted to pdf  |\n');
	fprintf('-------------------------------------------------------------------------------------\n');
end

if make_eps_figure 
	movefile([dirname './eps.figs/' file_ '.pdf'],[dirname '/' file_ '.pdf'])
	delete([dirname '/' file_ '.eps']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % SOME POSITIONING OPTIONS
% h=gcf;
% set(h,'PaperPositionMode','auto');         
% set(h,'PaperOrientation','landscape');
% set(h,'Position',[50 50 1200 800]);
% print(gcf, '-dpdf', 'test1.pdf')

% % MORE POSITIONING OPTIONS
% set(h,'PaperOrientation','landscape');
% set(h,'PaperPosition', [1 1 28 19]);
% 
% set(h,'PaperOrientation','landscape');
% set(h,'PaperUnits','normalized');
% set(h,'PaperPosition', [0 0 1 1]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [result,msg] = eps2pdf_funct(epsFile,fullGsPath,orientation)
% EPS2PDF Converts an eps file to a pdf file using GhostScript (GS)
% 
%   [result,msg] = eps2pdf(epsFile,fullGsPath,orientation)
% 
%   - epsFile:      eps file name to be converted to pdf file
%   - fullGsPath:   (optional) FULL GS path, including the file name, to
%                   the GS executable (on win32 it could be c:\program
%                   files\gs\gs8.14\bin\gswin32c.exe). The existence for
%                   fullGsPath will be checked for if given. On the other
%                   hand, if fullGsPath is not given or empty it defaults
%                   to 'gswin32c' for pcs and 'gs' for unix and the
%                   existence will not be checked for. But in this, latter
%                   case, GS's path must be in the system's path variable.
%   - orientation:  (optional) a flag that tells how the orientation tag in eps file should be treated
%                   just before the conversion (orientation tag is changed or even removed):
%                       0 -> no change (default)
%                       1 -> flip orientation
%                       2 -> remove orientation
%   
%   - result:       0 if ok, anything else otherwise
%   - msg:          resulting status on file being processed
%
%   NOTES: GhostScript is needed for this function to work. Orientation can
%   also be changed - use this only if you have problems with the orientation - in
%   such a case try with orientation=1 first and then orientation=2 if the first option is
%   not the right one.
%
%   EPS2PDF converts an existing EPS file to a PDF file using Ghostscript. EPS2PDF
%   reads an eps file, modifies the bounding box and creates a pdf file whose size
%   is determined by the bounding box and not by the paper size. This can not be
%   accomplished by using Ghostscript only. So, all that one needs is of course
%   Matlab and Ghostscript drivers.
% 
%   This tool is especially suited for LaTeX (TeX) users who want to create pdf
%   documents on the fly (by including pdf graphics and using either pdftex or
%   pdflatex). An example would be, if you are using LaTeX (TeX) to typeset
%   documents then the usual (simple) way to include graphics is to include eps
%   graphics using for example (if there exists myfigure.eps)
%   \begin{figure}
%       \centering
%       \includegraphics[scale=0.8]{myfigure}\\
%       \caption{some caption.}
%   \end{figure}
%   To use pdflatex (pdftex) you do not need to change anything but provide another
%   file myfigure.pdf in the same directory along with myfigure.eps. And this file,
%   of course, can be generated by EPS2PDF.
%
%   This function was tested on win32 system running Matlab R13sp1. It should work
%   on all platforms, if not, contact the author.
%
%   SOURCE:     The original idea came from the "eps-to-pdf" converter written in perl by Sebastian Rahtz
%
%   Primoz Cermelj, 24.08.2004
%   (c) Primoz Cermelj, primoz.cermelj@email.si
%
%   Version: 0.9.3
%   Last revision: 04.10.2004
%--------------------------------------------------------------------------
global epsFileContent

if ispc
    DEF_GS_PATH = 'gswin32c.exe';
else
    DEF_GS_PATH = 'gs';
end
GS_PARAMETERS = '-q -dNOPAUSE -dBATCH -dDOINTERPOLATE -dUseFlateCompression=true -sDEVICE=pdfwrite -r1200';

% error(nargchk(1,3,nargin));

if nargin < 2 || isempty(fullGsPath)
    fullGsPath = DEF_GS_PATH;
else
    if ~exist(fullGsPath)
        status = ['Ghostscript executable could not be found: ' fullGsPath];
        if nargout,      result = 1;    end;
        if nargout > 1,  msg = status;  else, disp(status);  end;
        return
    end
end

if nargin < 3 || isempty(orientation)
    orientation = 0;
end
orientation = abs(round(orientation));
orientation = orientation(1);
if orientation < 0 | orientation > 2
    orientation = 0;
end

epsFileContent = [];

%---------
% Get file name, path
%---------
source = epsFile;
[pathstr,sourceName,ext] = fileparts(source);
if isempty(pathstr)
    pathstr = cd;
    source = fullfile(pathstr,source);
end

targetName = [sourceName '.pdf'];
target = fullfile(pathstr,targetName);    % target - pdf file

tmpFileName = sourceName;
tmpFile = fullfile(pathstr,[tmpFileName ext '.eps2pdf~']);


% Create tmp file,...
[ok,errStr] = create_tmpepsfile(source,tmpFile,orientation);
if ~ok
    status = ['Error while creating temporary eps file: ' epsFile ' - ' errStr];
    if nargout,      result = 1;    end;
    if nargout > 1,  msg = status;  else, disp(status); end;
else
    % Run Ghostscript
    comandLine = ['"' fullGsPath '"' ' ' GS_PARAMETERS ' -sOutputFile=' '"' target '"' ' -f ' '"' tmpFile '"'];
    [stat, result] = system(comandLine);
    if stat
        status = ['pdf file not created - error running Ghostscript - check GS path: ' result];
        if nargout,      result = 1;    end;
        if nargout > 1,  msg = status;  else, disp(status);  end;  
    else
        status = 'pdf successfully created'; 
        status = '';
        if nargout,      result = 0;    end;
        if nargout > 1,  msg = status;  else, disp(status);  end; 
    end
end

% Delete tmp file
if exist(tmpFile)
    delete(tmpFile);
end






%/////////////////////////////////////////////////////////////////////
%                       SUBFUNCTIONS SECTION
%/////////////////////////////////////////////////////////////////////

%--------------------------------------------------------------------
function [ok,errStr] = create_tmpepsfile(epsFile,tmpFile,orientation)
% Creates tmp eps file - file with refined content
global epsFileContent

ok = 0;
errStr = [];
[ok,errStr] = read_epsfilecontent( epsFile );
if ~ok
    return
end
[ok,errStr] = update_epsfilecontent( epsFile,orientation );
if ~ok
    return
end
fh = fopen(tmpFile,'w');
if fh == -1
    errStr = ['Temporary file cannot be created. Check write permissions.'];
    return
end
try
    fwrite(fh,epsFileContent,'char');  % fwrite is faster than fprintf
    ok = 1;
catch
    errStr = ['Error writing temporary file. Check write permissions.'];
end
fclose(fh);
%--------------------------------------------------------------------


%--------------------------------------------------------------------
function [ok,errStr] = read_epsfilecontent( epsFile )
% Reads the content of the eps file into epsFileContent
global epsFileContent

ok = 0;
errStr = [];
fh = fopen(epsFile,'r');
if fh == -1
    errStr = ['File: ' epsFile ' cannot be accessed or does not exist'];
    return
end
try
    epsFileContent = fread(fh,'char=>char')';       % fread is faster than fscanf
    ok = 1;
catch
    errStr = lasterror;
end
fclose(fh);
%--------------------------------------------------------------------


%--------------------------------------------------------------------
function [ok,errStr] = update_epsfilecontent(epsFile,orientation)
% Updates eps file by adding some additional information into the header
% section concerning the bounding box (BB)
global epsFileContent

ok = 0;
errStr = [];

% Read current BB coordinates
ind = strfind( lower(epsFileContent), lower('%%BoundingBox:'));
if isempty(ind)
    errStr = ['Cannot find Bounding Box in file: ' epsFile];
    return
end
ii = ind(1) + 14;
fromBB = ii;
while ~((epsFileContent(ii) == sprintf('\n')) | (epsFileContent(ii) == sprintf('\r')) | (epsFileContent(ii) == '%'))
    ii = ii + 1;
end
toBB = ii - 1;
coordsStr = epsFileContent(fromBB:toBB);
coords = str2num( coordsStr );
if isempty(coords)
    errStr = ['Error reading BB coordinates from file: ' epsFile];
    return
end
NL = getnl;
w = abs(coords(3)-coords(1));
h = abs(coords(4)-coords(2));

% Change the orientation if requested
changeOrientation = 0;
if orientation ~= 0
    ind = strfind( lower(epsFileContent), lower('%%Orientation:'));
    if ~isempty(ind)
        ii = ind(1) + 14;
        fromOR = ii;
        while ~((epsFileContent(ii) == sprintf('\n')) | (epsFileContent(ii) == sprintf('\r')) | (epsFileContent(ii) == '%'))
            ii = ii + 1;
        end
        toOR = ii - 1;
        orientStr = strim(epsFileContent(fromOR:toOR));
        if ~isempty(orientStr) & orientation == 1           % flip
            if strfind(lower(orientStr),'landscape')
                changeOrientation = 1;
                orientStr = 'Portrait';
            elseif strfind(lower(orientStr),'portrait')
                changeOrientation = 1;
                orientStr = 'Landscape';                
            end            
        elseif  ~isempty(orientStr) & orientation == 2      % remove
            if strfind(lower(orientStr),'landscape') | strfind(lower(orientStr),'portrait')
                changeOrientation = 1;
                orientStr = ' ';
            end
        end
    end
end

% Refine the content - add additional information and even change the
% orientation
addBBContent = [' 0 0 ' int2str(w) ' ' int2str(h) ' ' NL...
                    '<< /PageSize [' int2str(w) ' ' int2str(h) '] >> setpagedevice' NL...
                    'gsave ' int2str(-coords(1)) ' ' int2str(-coords(2)) ' translate'];
if changeOrientation
    if fromOR > fromBB
        epsFileContent = [epsFileContent(1:fromBB-1) addBBContent epsFileContent(toBB+1:fromOR-1) orientStr epsFileContent(toOR+1:end)];
    else
        epsFileContent = [epsFileContent(1:fromOR-1) orientStr epsFileContent(toOR+1:fromBB-1) addBBContent epsFileContent(toBB+1:end)];
    end
else
    epsFileContent = [epsFileContent(1:fromBB-1) addBBContent  epsFileContent(toBB+1:end)];
end

ok = 1;
%--------------------------------------------------------------------

%--------------------------------------------------------------------
function NL = getnl
% Returns new-line string as found from first occurance from epsFileContent
global epsFileContent

NL = '\r\n';        % default (for Windows systems)
ii = 1;
len = length(epsFileContent);
while ~(epsFileContent(ii)==sprintf('\n') | epsFileContent(ii)==sprintf('\r') | ii<len)
    ii = ii + 1;
end
if epsFileContent(ii)==sprintf('\n')
    NL = '\n';
elseif epsFileContent(ii)==sprintf('\r')
    NL = '\r';
    if epsFileContent(ii+1)==sprintf('\n')
        NL = [NL '\n'];
    end
end
NL = sprintf(NL);
%--------------------------------------------------------------------


%--------------------------------------------------------------------
function outstr = strim(str)
% Removes leading and trailing spaces (spaces, tabs, endlines,...)
% from the str string.
if isnumeric(str);
    outstr = str;
    return
end
ind = find( ~isspace(str) );        % indices of the non-space characters in the str    
if isempty(ind)
    outstr = [];        
else
    outstr = str( ind(1):ind(end) );
end
%--------------------------------------------------------------------