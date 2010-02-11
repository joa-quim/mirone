function val = test_dms(s)
% Test if the input text string S is in the form dd:mm or dd:mm:ss
% val = test_dms(s) returns a cell array with the individual tokens
% if an error ocurrs or an empty imput is given, it returns an empty array
if isempty(s),    val = [];    return;   end
error = 0;
if (strcmp(s(1),':') | strcmp(s(end),':') | findstr(s,'::'))
    errordlg('One or more of the dd, mm or ss fields is empty','Error')
    error = error + 1;
end
[t,r]=strtok(s,':');
unit{1} = t;
if isnan(str2double(unit{1}))
    str = [unit{1} '   is not a valid number'];
    error = error + 1;
    errordlg(str,'Error')
end
i = 2;
while ~isempty(r)
    [t,r]=strtok(r,':');
    unit{i} = t;
    if isnan(str2double(unit{i}))
        str = [unit{i} '   is not a valid number'];
        error = error + 1;
        errordlg(str,'Error')
    end
    i = i + 1;
end
if (error), val = {[]};
else        val = unit;     end
