function out = read_textfile(file)
% READ_TEXTFILE(FILE) returns a cell array with the contents of the file FILE
%

fid = fopen(file, 'r');
if fid < 0
    errordlg(['Can''t open file:  ' file],'Error');
    return
end

out = {};   nl = 0;

while ~feof(fid)
    nl = nl + 1;
    tline = deblank(fgetl(fid));
    if ~ischar(tline) || strcmp(tline,''),    continue;     end
    out{nl,1} = tline;
end
fclose(fid);