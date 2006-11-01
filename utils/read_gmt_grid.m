function [X,Y,Z,head] = read_gmt_grid(h_calling_fig, opt)
% This function loads a GMT grid and is meant to be used in GUIs that
% may want to read process directly the GMT grids.

X = [];     Y = [];     Z = [];     head = [];
if (nargin == 0)
    h_calling_fig = [];     opt = [];
elseif (nargin == 1)
    opt = [];
elseif (nargin > 2)
    errordlg('Error in read_gmt_grid. Wrong number of arguments.','Error')
    return
end

if (isempty(opt))       % Otherwise 'opt' already transmited the file name.
    if (~isempty(handles.h_calling_fig))                    % If we know the handle to the calling fig
        cfig_handles = guidata(handles.h_calling_fig);      % get handles of the calling fig
        last_dir = cfig_handles.last_dir;
        home = cfig_handles.home_dir;
    else
        last_dir = [];
    end

    if (~isempty(last_dir)),    cd(last_dir);   end
    [FileName,PathName] = uigetfile({'*.grd;*.GRD', 'Grid files (*.grd,*.GRD)';'*.*', 'All Files (*.*)'},'Select GMT grid');
    pause(0.01);
    if (~isempty(last_dir)),    cd(home);   end
    if isequal(FileName,0);     return;     end
    fname = [PathName FileName];
end

% Because GMT and Surfer share the .grd extension, find out which kind grid we are dealing with
[fid, msg] = fopen(fname, 'r');
if (fid < 0)
    errordlg([PathName FileName ': ' msg],'ERROR'); return
end
ID = fread(fid,4,'*char');
ID = strread(ID,'%s');
if strcmp(ID,'DSBB')
    fname = [fname '=6'];
elseif strcmp(ID,'DSAA')
    warndlg('I don''t know and do not intend to learn how to read ASCII Surfer grids.','Warning')
    return
end
[X,Y,Z,head] = grdread_m(fname);
