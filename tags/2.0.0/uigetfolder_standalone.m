function folder = uigetfolder_standalone(title, initial_path)
%UIGETFOLDER_STANDALONE   Standard Windows browse for folder dialog box that
%                         runs in a standalone executable created with the MATLAB
%                         Compiler. The C/C++ Graphics Library is required for the
%                         errordlg function.
%
%   folder = uigetfolder_standalone(title, initial_path)
%
%   Output: folder       = selected folder (empty string if dialog cancelled)
%   Inputs: title        = title string (OPTIONAL)
%           initial_path = initial path (OPTIONAL, defaults to PWD)
%
%   Examples:   folder = uigetfolder_standalone                          - default title and initial path
%               folder = uigetfolder_standalone('Select results folder') - default initial path
%               folder = uigetfolder_standalone([], 'C:\Program Files')  - default title
%
%   See also UIGETFOLDER, UIGETFILE, UIPUTFILE

%-----------------------------------------------------------------------------------------------

if ~strcmp(computer, 'PCWIN')
   error_dialog_handle = errordlg(['The function ', upper(mfilename), ' only works on a MS-Windows PC'], ...
                                  mfilename, 'modal');
   folder = '';
else
   if nargin < 2
      initial_path = pwd;
   end
   
   if nargin < 1 | isempty(title)
      title = 'Select a folder';
   end
   
   % Error checking
   if ~ischar(title)
      error('The title must be a string')
   end
   if ~ischar(initial_path)
      error('The initial path must be a string')
   end
   % Following commented out to work in standalone executable
   %if ~exist(initial_path, 'dir')
   %   error(['The initial path: ', initial_path, ' does not exist!'])
   %end
   
    folder = uigetfolder_win32(title, initial_path);
% 	[t,r] = strtok(folder);
% 	if (~isempty(r))
%         warndlg(['If you had RTFM you should know that names with blanks '...
%                 '(wite spaces) are totaly FORBIDEN here.'],'ERROR')
%         folder = [];
% 	end
end
%-----------------------------------------------------------------------------------------------
