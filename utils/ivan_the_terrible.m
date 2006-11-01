function inc = ivan_the_terrible(dx,n,mod)
% This function was writen to deal with the TERRORIST dictatorship of GMT_grd_RI_verify.
% INC = IVAN_THE_TERRIBLE(DX,N,MODE)
% If MODE == 1, where DX is a real and N an integer (number of lines), returns a text string with INC
% that complies with GMT_grd_RI_verify.
% If MODE == 2, N now contains the candidate INC. The new INC is computed by a very simple algorithm
% (see below) that, on the contrary of a previous one (supposadly clever but very slow), works prety well.
%
% J. Luis   8/11/02

if nargin ~= 3
    errordlg('function ivan_the_terrible called with a wrong number of arguments','Error')
    return
end

if mod == 1
    inc = dx / (n-1);
elseif mod == 2
    t_inc = n;
    if (dx / t_inc) < 2
        msgbox(['You must be joking. With this increment your grid would have only two lines ' ...
            'along this dimension. I''ll give you the least idiot increment, which is the one that ' ...
            'alows three lines. But you better reconsider your choice.'],'Chico Clever');
        inc = dx / 2;
        return
    end
    nl = round(dx / t_inc);
    inc = dx / (nl - 1);
else
    errordlg('function ivan_the_terrible called in a wrong mode','Error')
end
