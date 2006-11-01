function [s,w] = run_and_report_error(str)
% [s,w] = run_and_error_report(str) is a replacement to the "dos" or "unix" command when
% used in codes that will be compiled. The Matlab original don't accept two
% outputs when they are compiled.
%
% It differs from mat_lyies because the latter screwed up things with the instruction:
% dos([str ' > text_stdout 2> text_stderr']). When gmt code outputs to stdout, the result would
% go into text_stdout and not to the output file name

w = '';   s = 0;
if isunix           % UNIX
    unix([str ' 2> text_stderr']);
elseif ispc         % Windows
    dos([str ' 2> text_stderr']);
else
    errordlg('Unknown platform.','Error');
end

% Now find out if it worked and get the output message into w

fid=fopen('text_stderr');
while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    if findstr(tline,'is not recognized')    % case windows error
        s = 1;        w = [str '. Command not found'];
    elseif findstr(tline,'command not found')    % case unix (or cygwin) error
        s = 1;        w = [str '. Command not found'];
    elseif findstr(tline,'ERROR')
        s = 1;        w = [w tline];
    elseif ~isempty(tline)          % Gave up. I cannot guess all error messages, so if any report it
        s = 1;        w = [w tline];
    end
end
fclose(fid);

delete text_stderr;
