function out_name = decompress(full_name, opt)
% Decompress files ZIP or GZIP compressed.
% FULL_NAME means what it says. That is, it must contain file's absolute path, extensions included
% OPT == 'warn' means that a warning fig will exist during the decompression.
% OUT_NAME is the name of the uncompressed file.

out_name = [];      do_warn = 0;
if ( (nargin == 2 && strcmp(opt,'warn')) )
    do_warn = 1;
end

[PATH,fname,EXT] = fileparts(full_name);
% Check that the file is zip or gzip compressed
if (~(strcmpi(EXT,'.zip') || strcmpi(EXT,'.gz')) )
    errordlg(['Error in decompress. File ' full_name ' is not compressed'],'Error')
    return
end

if (isempty(PATH))
    errordlg(['Error in decompress. File ' full_name ' must have the absolute path'],'Error');
end

str = ['gunzip -q -N -f -c ' full_name ' > ' [PATH filesep fname]];
if (do_warn)
    h_wait = waitfig(['Uncompressing ' full_name]);
end

if (isunix),    s = unix(str);
elseif ispc,   s = dos(str);
else            errordlg('Unknown platform.','Error');
end
if ~(isequal(s,0))                  % An error as occured
    errordlg(['Error decompressing file ' full_name],'Error');
    if (do_warn),   delete(h_wait);     end
    return
end
if (do_warn),   delete(h_wait);     end

out_name = [PATH filesep fname];