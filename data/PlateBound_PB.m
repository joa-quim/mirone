function PlateBound_PB()
fid = fopen('PB2002_steps.dat','r');
todos = fread(fid,'*char');
[sis pb_id lon1 lat1 lon2 lat2 len azim vel azim_vel div slip bat age class] = ...
    strread(todos,'%d %s %f %f %f %f %f %f %f %f %f %f %d %d %s');
fclose(fid);
clear todos

n_otf = 1;  np_otf = 0;
n_osr = 1;  np_osr = 0;
n_crb = 1;  np_crb = 0;
n_ctf = 1;  np_ctf = 0;
n_ccb = 1;  np_ccb = 0;
n_ocb = 1;  np_ocb = 0;
n_sub = 1;  np_sub = 0;

np = length(lat1);

OTF(n_otf).x_otf = [];      OTF(n_otf).y_otf = [];
OSR(n_osr).x_osr = [];      OSR(n_osr).y_osr = [];
CRB(n_crb).x_crb = [];      CRB(n_crb).y_crb = [];
CTF(n_ctf).x_ctf = [];      CTF(n_ctf).y_ctf = [];
CCB(n_ccb).x_ccb = [];      CCB(n_ccb).y_ccb = [];
OCB(n_ocb).x_ocb = [];      OCB(n_ocb).y_ocb = [];
SUB(n_sub).x_sub = [];      SUB(n_sub).y_sub = [];

vel_sum = 0;   azim_vel_sum = 0;

for i=1:np
    if (length(pb_id{i}) == 6), pb_id{i} = pb_id{i}(2:end);    end
end

for i=1:np
    if (sis(i) ~= -1)
        k = strfind(class{i},':');
        if isempty(k),  k = 1;
        else            k = 2;  end
        switch class{i}(k:k+2)
            case 'OTF'
                OTF(n_otf).x_otf = [OTF(n_otf).x_otf lon1(i)];  OTF(n_otf).y_otf = [OTF(n_otf).y_otf lat1(i)];
                np_otf = np_otf + 1;
                OTF(n_otf).x_otf = [OTF(n_otf).x_otf lon2(i)];  OTF(n_otf).y_otf = [OTF(n_otf).y_otf lat2(i)];
                vel_sum = vel_sum + vel(i);     azim_vel_sum = azim_vel_sum + azim_vel(i);
                np_otf = np_otf + 1;
            case 'OSR'
                OSR(n_osr).x_osr = [OSR(n_osr).x_osr lon1(i)];  OSR(n_osr).y_osr = [OSR(n_osr).y_osr lat1(i)];
                np_osr = np_osr + 1;
                OSR(n_osr).x_osr = [OSR(n_osr).x_osr lon2(i)];  OSR(n_osr).y_osr = [OSR(n_osr).y_osr lat2(i)];
                vel_sum = vel_sum + vel(i);     azim_vel_sum = azim_vel_sum + azim_vel(i);
                np_osr = np_osr + 1;
            case 'CRB'
                CRB(n_crb).x_crb = [CRB(n_crb).x_crb lon1(i)];  CRB(n_crb).y_crb = [CRB(n_crb).y_crb lat1(i)];
                np_crb = np_crb + 1;
                CRB(n_crb).x_crb = [CRB(n_crb).x_crb lon2(i)];  CRB(n_crb).y_crb = [CRB(n_crb).y_crb lat2(i)];
                vel_sum = vel_sum + vel(i);     azim_vel_sum = azim_vel_sum + azim_vel(i);
                np_crb = np_crb + 1;
            case 'CTF'
                CTF(n_ctf).x_ctf = [CTF(n_ctf).x_ctf lon1(i)];  CTF(n_ctf).y_ctf = [CTF(n_ctf).y_ctf lat1(i)];
                np_ctf = np_ctf + 1;
                CTF(n_ctf).x_ctf = [CTF(n_ctf).x_ctf lon2(i)];  CTF(n_ctf).y_ctf = [CTF(n_ctf).y_ctf lat2(i)];
                vel_sum = vel_sum + vel(i);     azim_vel_sum = azim_vel_sum + azim_vel(i);
                np_ctf = np_ctf + 1;
            case 'CCB'
                CCB(n_ccb).x_ccb = [CCB(n_ccb).x_ccb lon1(i)];  CCB(n_ccb).y_ccb = [CCB(n_ccb).y_ccb lat1(i)];
                np_ccb = np_ccb + 1;
                CCB(n_ccb).x_ccb = [CCB(n_ccb).x_ccb lon2(i)];  CCB(n_ccb).y_ccb = [CCB(n_ccb).y_ccb lat2(i)];
                vel_sum = vel_sum + vel(i);     azim_vel_sum = azim_vel_sum + azim_vel(i);
                np_ccb = np_ccb + 1;
            case 'OCB'
                OCB(n_ocb).x_ocb = [OCB(n_ocb).x_ocb lon1(i)];  OCB(n_ocb).y_ocb = [OCB(n_ocb).y_ocb lat1(i)];
                np_ocb = np_ocb + 1;
                OCB(n_ocb).x_ocb = [OCB(n_ocb).x_ocb lon2(i)];  OCB(n_ocb).y_ocb = [OCB(n_ocb).y_ocb lat2(i)];
                vel_sum = vel_sum + vel(i);     azim_vel_sum = azim_vel_sum + azim_vel(i);
                np_ocb = np_ocb + 1;
            case 'SUB'
                SUB(n_sub).x_sub = [SUB(n_sub).x_sub lon1(i)];  SUB(n_sub).y_sub = [SUB(n_sub).y_sub lat1(i)];
                np_sub = np_sub + 1;
                SUB(n_sub).x_sub = [SUB(n_sub).x_sub lon2(i)];  SUB(n_sub).y_sub = [SUB(n_sub).y_sub lat2(i)];
                vel_sum = vel_sum + vel(i);     azim_vel_sum = azim_vel_sum + azim_vel(i);
                np_sub = np_sub + 1;
        end
    else
        k = strfind(class{i-1},':');
        if isempty(k),  k = 1;
        else            k = 2;  end
        switch class{i-1}(k:k+2)
            case 'OTF'
                n = (np_otf / 2);
                OTF(n_otf).vel = vel_sum / n;       % vel media do segmento
                OTF(n_otf).azim_vel = azim_vel_sum / n;    % azim medio do segmento
                OTF(n_otf).pb_id = pb_id(i-1);      % Identificador das fronteiras
                OTF(n_otf).class = class{i-1}(k:k+2);      % Identificador de tipo de fronteira
                vel_sum = 0;   azim_vel_sum = 0;    % reinicializa estas porque vao ser reutilizadas
                n_otf = n_otf + 1;  np_otf = 0;     % incrementa o numero de segs e reinicializa o num de pts no seg
                OTF(n_otf).x_otf = [];      OTF(n_otf).y_otf = [];  % reinicializa
            case 'OSR'
                n = (np_osr / 2);
                OSR(n_osr).vel = vel_sum / n;       % vel media do segmento
                OSR(n_osr).azim_vel = azim_vel_sum / n;    % azim medio do segmento
                OSR(n_osr).pb_id = pb_id(i-1);      % Identificador das fronteiras
                OSR(n_osr).class = class{i-1}(k:k+2);      % Identificador de tipo de fronteira
                vel_sum = 0;   azim_vel_sum = 0;    % reinicializa estas porque vao ser reutilizadas
                n_osr = n_osr + 1;  np_osr = 0;
                OSR(n_osr).x_osr = [];      OSR(n_osr).y_osr = [];
            case 'CRB'
                n = (np_crb / 2);
                CRB(n_crb).vel = vel_sum / n;       % vel media do segmento
                CRB(n_crb).azim_vel = azim_vel_sum / n;    % azim medio do segmento
                CRB(n_crb).pb_id = pb_id(i-1);      % Identificador das fronteiras
                CRB(n_crb).class = class{i-1}(k:k+2);      % Identificador de tipo de fronteira
                vel_sum = 0;   azim_vel_sum = 0;    % reinicializa estas porque vao ser reutilizadas
                n_crb = n_crb + 1;  np_crb = 0;
                CRB(n_crb).x_crb = [];      CRB(n_crb).y_crb = [];
            case 'CTF'
                n = (np_ctf / 2);
                CTF(n_ctf).vel = vel_sum / n;       % vel media do segmento
                CTF(n_ctf).azim_vel = azim_vel_sum / n;    % azim medio do segmento
                CTF(n_ctf).pb_id = pb_id(i-1);      % Identificador das fronteiras
                CTF(n_ctf).class = class{i-1}(k:k+2);      % Identificador de tipo de fronteira
                vel_sum = 0;   azim_vel_sum = 0;    % reinicializa estas porque vao ser reutilizadas
                n_ctf = n_ctf + 1;  np_ctf = 0;
                CTF(n_ctf).x_ctf = [];      CTF(n_ctf).y_ctf = [];
            case 'CCB'
                n = (np_ccb / 2);
                CCB(n_ccb).vel = vel_sum / n;       % vel media do segmento
                CCB(n_ccb).azim_vel = azim_vel_sum / n;    % azim medio do segmento
                CCB(n_ccb).pb_id = pb_id(i-1);      % Identificador das fronteiras
                CCB(n_ccb).class = class{i-1}(k:k+2);      % Identificador de tipo de fronteira
                vel_sum = 0;   azim_vel_sum = 0;    % reinicializa estas porque vao ser reutilizadas
                n_ccb = n_ccb + 1;  np_ccb = 0;
                CCB(n_ccb).x_ccb = [];      CCB(n_ccb).y_ccb = [];
            case 'OCB'
                n = (np_ocb / 2);
                OCB(n_ocb).vel = vel_sum / n;       % vel media do segmento
                OCB(n_ocb).azim_vel = azim_vel_sum / n;    % azim medio do segmento
                OCB(n_ocb).pb_id = pb_id(i-1);      % Identificador das fronteiras
                OCB(n_ocb).class = class{i-1}(k:k+2);      % Identificador de tipo de fronteira
                vel_sum = 0;   azim_vel_sum = 0;    % reinicializa estas porque vao ser reutilizadas
                n_ocb = n_ocb + 1;  np_ocb = 0;
                OCB(n_ocb).x_ocb = [];      OCB(n_ocb).y_ocb = [];
            case 'SUB'
                n = (np_sub / 2);
                SUB(n_sub).vel = vel_sum / n;       % vel media do segmento
                SUB(n_sub).azim_vel = azim_vel_sum / n;    % azim medio do segmento
                SUB(n_sub).pb_id = pb_id(i-1);      % Identificador das fronteiras
                SUB(n_sub).class = class{i-1}(k:k+2);      % Identificador de tipo de fronteira
                vel_sum = 0;   azim_vel_sum = 0;    % reinicializa estas porque vao ser reutilizadas
                n_sub = n_sub + 1;  np_sub = 0;
                SUB(n_sub).x_sub = [];      SUB(n_sub).y_sub = [];
        end
    end
end

% Limpa os pts repetidos da OTF
for i=1:n_otf-1
    n = length(OTF(i).x_otf);   k = [];
    for j = 2:n
        if (OTF(i).x_otf(j) == OTF(i).x_otf(j-1))
            k = [k j];
        end
    end
    OTF(i).x_otf(k) = [];   OTF(i).y_otf(k) = [];
end;    OTF(n_otf) = [];

% Limpa os pts repetidos da OSR
for i=1:n_osr-1
    n = length(OSR(i).x_osr);   k = [];
    for j = 2:n
        if (OSR(i).x_osr(j) == OSR(i).x_osr(j-1))
            k = [k j];
        end
    end
    OSR(i).x_osr(k) = [];   OSR(i).y_osr(k) = [];
end;    OSR(n_osr) = [];

% Limpa os pts repetidos da CRB
for i=1:n_crb-1
    n = length(CRB(i).x_crb);   k = [];
    for j = 2:n
        if (CRB(i).x_crb(j) == CRB(i).x_crb(j-1))
            k = [k j];
        end
    end
    CRB(i).x_crb(k) = [];   CRB(i).y_crb(k) = [];
end;    CRB(n_crb) = [];

% Limpa os pts repetidos da CTF
for i=1:n_ctf-1
    n = length(CTF(i).x_ctf);   k = [];
    for j = 2:n
        if (CTF(i).x_ctf(j) == CTF(i).x_ctf(j-1))
            k = [k j];
        end
    end
    CTF(i).x_ctf(k) = [];   CTF(i).y_ctf(k) = [];
end;    CTF(n_ctf) = [];

% Limpa os pts repetidos da CCB
for i=1:n_ccb-1
    n = length(CCB(i).x_ccb);   k = [];
    for j = 2:n
        if (CCB(i).x_ccb(j) == CCB(i).x_ccb(j-1))
            k = [k j];
        end
    end
    CCB(i).x_ccb(k) = [];   CCB(i).y_ccb(k) = [];
end;    CCB(n_ccb) = [];

% Limpa os pts repetidos da OCB
for i=1:n_ocb-1
    n = length(OCB(i).x_ocb);   k = [];
    for j = 2:n
        if (OCB(i).x_ocb(j) == OCB(i).x_ocb(j-1))
            k = [k j];
        end
    end
    OCB(i).x_ocb(k) = [];   OCB(i).y_ocb(k) = [];
end;    OCB(n_ocb) = [];

% Limpa os pts repetidos da SUB
for i=1:n_sub-1
    n = length(SUB(i).x_sub);   k = [];
    for j = 2:n
        if (SUB(i).x_sub(j) == SUB(i).x_sub(j-1))
            k = [k j];
        end
    end
    SUB(i).x_sub(k) = [];   SUB(i).y_sub(k) = [];
end;    SUB(n_sub) = [];

clear sis pb_id lon1 lat1 lon2 lat2 len azim vel azim_vel div slip bat age class;
clear i j k n_otf  np_otf n_osr np_osr n_crb np_crb n_ctf np_ctf n_ccb  np_ccb n_ocb  np_ocb n_sub np_sub
clear fid n np vel_sum azim_vel_sum

save('PB_boundaries','OTF', 'OSR', 'CRB', 'CTF', 'CCB', 'OCB', 'SUB')

%xx=1;