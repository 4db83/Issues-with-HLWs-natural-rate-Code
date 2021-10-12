function [Status ,Error_message ] = delete_excel_sheets(file_name,Sheet_names)
%%**********************************************************************************************************
%   Name          : Delete_sheets_Excel
%   Author        : Pruthvi Raj G  ::  (9677066394 :: www.prudhvy.com )
%   Version       : Version 1.0 - 2011b Compactible
%   Description   : Deleting Excel Sheets required after writing data to Excel. 
%   Input         : File_Name with path included , Sheet_name / Sheet_names.
%   Date          : 22-April-2019
%  
%   Examples      : delete_excel_sheets('D:\Pruthvi\Test_file.xls',{'Sheet1','Sheet2'})  %To delete 2 sheets
%                   delete_excel_sheets('D:\Pruthvi\Test_file.xls','Sheet1')
%                   delete_excel_sheets('D:\Pruthvi\Test_file.xls') % Takes 'Sheet1' as Default
%************************************************************************************************************
if nargin < 2
Sheet_names = 'Sheet1';    
end
 Status = 1 ;
 Error_message = 'Succesfully Deleted the Sheets' ; 
% Open the Excel file (*.xls or *.xlsx).
objExcel = actxserver('Excel.Application');
objExcel.Workbooks.Open(file_name);
% Deleting  sheets .
try
if ischar(Sheet_names)
    objExcel.ActiveWorkbook.Worksheets.Item(Sheet_names).Delete;
else
    if size(Sheet_names,1) > size(Sheet_names,2)
        Size_sheets =  size(Sheet_names,1);
    else
        Size_sheets =  size(Sheet_names,2);
    end
    
    for i = 1:Size_sheets
        objExcel.ActiveWorkbook.Worksheets.Item(char(Sheet_names(i))).Delete;
    end
end
catch ME 
    Status = 0 ;
    
    if regexp(ME.message,'Invoke Error, Dispatch Exception: Invalid index.','once')
        Error_message = 'Check the Sheet Name Properly ! It doesn''t Exist';
    else
        Error_message = 'Check the Input Sheet Name';
    end
 end
% Save Excel Sheet
objExcel.ActiveWorkbook.Save;
% Close Excel Sheet
objExcel.ActiveWorkbook.Close;
% quit Excel Object
objExcel.Quit;
objExcel.delete;