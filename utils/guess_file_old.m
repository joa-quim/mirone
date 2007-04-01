function [bin,n_column,multi_seg,n_headers] = guess_file(fiche)
% [bin,n_column,multi_seg,n_headers] = guess_file(fiche) trys to guess if file "fiche"
% is ascii or binary. 
% If it detects that "fiche" is ascii this function tries to find out wether the multisegment
% symbol (">") is present, the number of columns in the file and if it has header lines.
% NOTE that this last tests may not allways give relyable results.

% Error testing
if nargin ~= 1
    errordlg('function guess_file: must give an input file name','File Error')
    bin = '';    multi_seg = '';  n_headers = '';  n_column = '';
    return
end

% Now remove leading white space in file name "fiche"
fiche = ddewhite(fiche);

type = exist(fiche,'file');
if type == 0
    errordlg(['function guess_file: file "' fiche '" does not exist'],'File Error')
    bin = '';    multi_seg = '';  n_headers = '';  n_column = '';
    return
end

if nargout == 0
    errordlg('function guess_file: must specify at least one output argument','Error')
    bin = '';    multi_seg = '';  n_headers = '';  n_column = '';
    return
end

% ----------------------------------------------------------------------------------------

global home_dir
global file_guess_on
bin = 0;    multi_seg = 0;  n_headers = 0;  n_column = 0;  nl_max = 30;

if isempty(home_dir)        % Case when this function was called directly
    home_dir = pwd;
    f_path = ['tmp' filesep];
else
    f_path = [home_dir filesep 'tmp' filesep];
end

fid =fopen(fiche);  A = fscanf(fid,'%c');  fclose(fid);
dumb = strread(A,'%s','delimiter','\n');
nl_max = min(nl_max,length(dumb));  clear dumb;
A = double(A);
sup = find(A > 126 & A < 192);
if ~isempty(sup)    % Binary files have bytes with values greater than 126 (but so is the ç char)
    bin = 1;
end

% -------------------------------------------------------------------------------------
% This is the skeleton of the future attempt to guess what's inside a binary file (to be continued)
%NaN_s = [255 255 255 127];  NaN_d = [255 255 255 255 255 255 255 127];
%i = 1;  n_nan =0;
%while i < 78
%    Bs = A(i:i+3);    Bd = A(i:i+7);
%    is_nan_s = find(Bs == NaN_s);    is_nan_d = find(Bd == NaN_d);
%    if isempty(is_nan_s) | length(is_nan_s) ~= 4;
%        break
%    else
%        n_nan = n_nan + 1;
%        i = i + 4;
%    end
%end
% -------------------------------------------------------------------------------------

% For ascii files make a crude test to find number of columns and number of headers
if bin == 0 
    fid=fopen(fiche);
    i = 1;  n_col(1:nl_max) = 0;    n_multi = 0;    t_line = cell(nl_max+1,1);
    t_line{1} = 1;                  % condition for the next loop may start
    while (t_line{i} > 0 & i <= nl_max)
        t_line{i} = fgetl(fid);                     % This is infameously slow, but it only reads at most 30 lines
        if (t_line{i} == -1),   break;      end;    % We are donne (maybe due to '>' header lines)
        if isempty(t_line{i});  continue;   end;    % Jump blank lines
        if findstr(t_line{i},'>')   % multisegmet line, so the rest of it is of no interest (but count them)
            n_multi = n_multi + 1;
            continue
        end
        t_line{i} = deblank(t_line{i});          % Blanks make a mess to the guessing code
        [tok,rem]=strtok(t_line{i});
        if ~isempty(rem);   n_col(i) = n_col(i) + 1;    end     % count first column
        while ~isempty(rem)
            [tok,rem]=strtok(rem);
            n_col(i) = n_col(i) + 1;
        end
        i = i + 1;
        t_line{i} = t_line{i-1};
    end
    i = i - 1;      % i was incremented 1 to much at the end of the previous loop
    fclose(fid);
    multi_seg = n_multi;

    % Now decide how many columns have the data lines. The easeast is to assume that the info is in the last line
    % However this may fail if last line contains, for example, the multisegment symbol (">").
    % So, do another test.
    m = i - n_multi;
    if (m > 1)
        n_c1 = n_col(m);    n_c2 = n_col(m-1);
        if (n_c1 ~= n_c2 && isempty(find(t_line{n_col(m-1)} > 57 & t_line{n_col(m-1)} < 127)))
            n_column = max(n_c1,n_c2);
        else
            n_column = n_c1;
        end
    else
        n_column = n_col;
    end
    
    % Well, this is a stupid patch for the case where the file has only one column.
    % In that case the above test failed
    if (n_column == 0),     n_column = 1;   end     % (unless the file is empty!!!)
    
    % Now test if header lines are present (ascii 65:122 contain upper and lower case letters)
    for j=1:m
        head = find((t_line{j} > 32 & t_line{j} < 43) | t_line{j} > 58);
        tmp = find(t_line{j} == 78);    % I'm searching for a NaN string (ascii 78,97,78)
        if ~isempty(tmp) && length(tmp) == 2
            tmp = find(t_line{j}(tmp(1)+1) == 97);
        else
            tmp = [];
        end
        if ~isempty(tmp)
            % NaN string found. But now we have a problem. If we reach here it's because a text string
            % was found. But that may equaly happen on a header line, or a data line with NaNs.
            % So how to decide which was the case? The way that occurs to me is to decide on basis
            % of the number of columns, but this is far from being 100% full proof.
            if n_col(j) == n_column
                head = [];              % NaNs on a data line, or in a header line with same n columns as data line
            else
                head = 1;               % NaNs on a header line
            end
            % However, exclude eventualy NaNs in a header line with the same n columns as a data line
            if isempty(head);
                w = strrep(t_line{j},'NaN','');
                head = find((w > 31 & w < 43) | w > 58);
            end
        end
        if ~isempty(head);   n_headers = n_headers + 1;  end
    end
end


% -------------------------------------------------------------------------------------
% Provide the user with the sumary of findings (if he wants so)
%if file_guess_on == 0;     return;     end     % Don't show the report of file findings

message = cell(1,6);
message{1} = ['Report on the findings about file: "' fiche '"'];
%message_intro = ['The following is what I guessed from your file:'];
if bin == 0 || (multi_seg ~= 0) || (n_headers ~= 0)
    if (multi_seg ~= 0)
        message{3} = '     - ASCII multisegment file.';
        message{4} = ' ';
    end
    if (n_headers ~= 0)
        message{5} = ['     - ASCII file with ' num2str(n_headers) ' header(s) lines.'];
        message{6} = ' ';
    end
    if ~isempty(message{3}) || ~isempty(message{5})
        %message{7} = ['NOTE: If these window warnings annoy you (like they do to me), you may disable ' ...
                %'them on the checkbox below. For getting them back (advised mode for beguinners), check ' ...
                %'"Guess what''s in data Files" checkbox in File -> Preferences of the main M_GMT window'];
        %wonder_file_guesser(message)
        %msgbox(message,'Wonder File Guesser')
    end
elseif bin == 1
    message{3} = [' - It''s an BINARY file but I cannot know if it''s single or double precision. ' ...
        'I''ll assume a double precision file but it''s up to you to decide if that''s correct.'];
        %message{4} = ' ';
        %message{5} = ['NOTE: If these window warnings annoy you (like they do to me), you may disable ' ...
                %'them on the checkbox below. For getting them back (advised mode for beguinners), check ' ...
                %'"Guess what''s in data Files" checkbox in File -> Preferences of the main M_GMT window'];
        %wonder_file_guesser(message)
        %msgbox(message,'Wonder File Guesser')
end
