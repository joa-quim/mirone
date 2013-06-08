function [comp,dilat] = patch_meca(str1,dip1,rake1,str2,dip2,rake2)
% Compute focal mechanisms on the equal area Schmidt projection
% STR1,DIP1 & RAKE1 are the strike, dip and rake of the first nodal plane
% Optionaly, give also STR2,DIP2 & RAKE2 for the second nodal plane. If not
% given, they will be computed from the parameters of the first nodal plane.
% COMP  -> Nx2 (x,y) var in the range [-1 1] with the compressive quadrants
% DILAT -> Nx2 (x,y) var in the range [-1 1] with the dilatation quadrants

%	Copyright (c) 2004-2012 by J. Luis
%
% 	This program is part of Mirone and is free software; you can redistribute
% 	it and/or modify it under the terms of the GNU Lesser General Public
% 	License as published by the Free Software Foundation; either
% 	version 2.1 of the License, or any later version.
% 
% 	This program is distributed in the hope that it will be useful,
% 	but WITHOUT ANY WARRANTY; without even the implied warranty of
% 	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% 	Lesser General Public License for more details.
%
%	Contact info: w3.ualg.pt/~jluis/mirone
% --------------------------------------------------------------------

np = 45;
D2R = pi/180;
TWO_PI = 2 * pi;

if (nargin ~= 3 && nargin ~= 6)
    errordlg('PATCH_MECA: wrong number of arguments','Error')
    return
end

if (abs(dip1-90) < 1e-5),    dip1 = 90 - 1e-5;   end
if (rake1 > 180),    rake1 = rake1 - 360;    end
if (abs(rake1)-90 < 1e-5 || abs(rake1)-180 < 1e-5)
    rake1 = rake1 - sign(rake1) * 1e-5;
end
if (nargin == 6)
	if (abs(dip2-90) < 1e-5),   dip2 = 90 - 1e-5;   end
    if (rake2 > 180),    rake2 = rake2 - 360;    end
	if (abs(rake2)-90 < 1e-5 || abs(rake1)-180 < 1e-5)
        rake2 = rake2 - sign(rake2) * 1e-5;
	end
end
NP1.str = str1;
NP1.dip = dip1;
NP1.rake = rake1;

if (nargin == 3)
    [str2,dip2,rake2] = define_second_plane(NP1);
end
NP2.str = str2;
NP2.dip = dip2;
NP2.rake = rake2;

[P,T,N] = dc_to_axe(NP1,NP2);

% ---------------- First nodal plane
rot_s1 = [cos(NP1.str*D2R) -sin(NP1.str*D2R);  sin(NP1.str*D2R) cos(NP1.str*D2R)];   % Strike1 rotation matrix

% Angles from 0 till null axis and from null axil till 180 
ang1 = (N.str - NP1.str) * D2R;
if (ang1 < 0),      ang1 = ang1 + TWO_PI;  end
if (ang1 > pi),     ang1 = ang1 - pi;          end

teta1 = [0:pi/np:ang1 ang1];                teta2 = [ang1:pi/np:pi pi];
s_teta1 = sin(teta1);                       s_teta2 = sin(teta2);

radip = atan(tan(NP1.dip*D2R) * s_teta1);   rproj = sqrt(2)*sin((pi/2 - radip)/2);
x1_P1 = rproj .* s_teta1;                   y1_P1 = rproj .* cos(teta1);

radip = atan(tan(NP1.dip*D2R) * s_teta2);   rproj = sqrt(2)*sin((pi/2 - radip)/2);
x2_P1 = rproj .* s_teta2;                   y2_P1 = rproj .* cos(teta2);

% Rotate to P1_strike
plan1_a = [x1_P1' y1_P1'] * rot_s1;         plan1_b = [x2_P1' y2_P1'] * rot_s1;

% --------------- Second nodal plane
rot_s2 = [cos(str2*D2R) -sin(str2*D2R);  sin(str2*D2R) cos(str2*D2R)];   % Strike2 rotation matrix
ang2 = (NP2.str - N.str) * D2R;
if (ang2 < 0),      ang2 = ang2 + TWO_PI;  end
if (ang2 > pi),     ang2 = ang2 - pi;          end
teta1 = [0:pi/np:pi-ang2 pi-ang2];          teta2 = [pi-ang2:pi/np:pi pi];
s_teta1 = sin(teta1);                       s_teta2 = sin(teta2);

radip = atan(tan(dip2*D2R) * s_teta1);      rproj = sqrt(2)*sin((pi/2 - radip)/2);
x1_P2 = rproj .* s_teta1;                   y1_P2 = rproj .* cos(teta1);

radip = atan(tan(dip2*D2R) * s_teta2);      rproj = sqrt(2)*sin((pi/2 - radip)/2);
x2_P2 = rproj .* s_teta2;                   y2_P2 = rproj .* cos(teta2);

% Rotate to P1_strike
plan2_a = [x1_P2' y1_P2'] * rot_s2;         plan2_b = [x2_P2' y2_P2'] * rot_s2;

%disp(['str2=' num2str(str2) '  strNul=' num2str(N.str)])
clear x1_P1 x2_P1 x1_P2 x2_P2 ang2 ang2 teta1 teta2 rot_s1 rot_s2 s_teta1 s_teta2 rproj radip;
clear y1_P1 y1_P2 y2_P1 y2_P2;

% ---------- Compute the arc-circles that close compressive/dilatational parts
if (NP1.rake >= 0)
	if (NP2.str < NP1.str)
        a1 = NP2.str * D2R + pi;                a2 = NP1.str * D2R;
        aa = sort([a1 a2]);    ang_g = aa(2);   ang_l = aa(1);
	else
        a1 = NP2.str * D2R - pi;                a2 = NP1.str * D2R;
        aa = sort([a1 a2]);    ang_g = aa(2);   ang_l = aa(1);
	end
else
    ang_g = NP1.str * D2R;             ang_l = NP2.str * D2R;
    if (ang_l > ang_g),   ang_g = ang_g + TWO_PI;   end
end
teta = pi/2 - [ang_l:pi/np:ang_g ang_g]';
arc1_x = cos(teta(end:-1:1));               arc1_y = sin(teta(end:-1:1));
arc3_x = -arc1_x(end:-1:1);                 arc3_y = -arc1_y(end:-1:1);

if (NP1.rake >= 0)
	if (NP2.str > NP1.str)
        ang_g = NP1.str * D2R + pi;         ang_l = NP2.str * D2R - pi;
	else
        ang_g = NP1.str * D2R + pi;         ang_l = NP2.str * D2R + pi;
	end
else
    ang_l = NP1.str * D2R;                  ang_g = NP2.str * D2R + pi;
    if (ang_g - ang_l > TWO_PI),          ang_g = ang_g - TWO_PI;   end
    if (ang_l > ang_g)
        ang_l = NP2.str * D2R;          ang_g = NP1.str * D2R - pi;
        tmp = plan1_a;  plan1_a = plan2_a;  plan2_a = tmp;
        tmp = plan1_b;  plan1_b = plan2_b;  plan2_b = tmp;
    end
% 	if (NP2.str > NP1.str)
%         ang_g = NP1.str * D2R + pi;             ang_l = NP2.str * D2R;
%         if (ang_g < ang_l)  ang_l = ang_l - TWO_PI; end
%     else
%         %ang_g = NP1.str * D2R - 0*pi;           ang_l = NP2.str * D2R + 0*pi;
%         ang_g = NP2.str * D2R + 1*pi;           ang_l = NP1.str * D2R + 0*pi;
%     end
end
teta = pi/2 - [ang_l:pi/np:ang_g ang_g]';
arc2_x = cos(teta);                     arc2_y = sin(teta);
arc4_x = -arc2_x(end:-1:1);             arc4_y = -arc2_y(end:-1:1);

if (NP1.rake >= 0)
	pat1_x = [plan1_a(:,1); plan1_b(:,1); arc3_x; plan2_a(:,1); plan2_b(:,1); arc1_x];
	pat1_y = [plan1_a(:,2); plan1_b(:,2); arc3_y; plan2_a(:,2); plan2_b(:,2); arc1_y];
	pat2_x = [plan2_b(:,1); arc2_x; plan1_b(end:-1:1,1); plan1_a(end:-1:1,1); arc4_x; plan2_a(:,1)];
	pat2_y = [plan2_b(:,2); arc2_y; plan1_b(end:-1:1,2); plan1_a(end:-1:1,2); arc4_y; plan2_a(:,2)];
    comp = [pat1_x pat1_y];
    dilat = [pat2_x pat2_y];
else
    pat1_x = [plan1_a(:,1); plan1_b(:,1); arc3_x(end:-1:1,1); plan2_b(end:-1:1,1); plan2_a(end:-1:1,1); arc1_x(end:-1:1,1)];
    pat1_y = [plan1_a(:,2); plan1_b(:,2); arc3_y(end:-1:1,1); plan2_b(end:-1:1,2); plan2_a(end:-1:1,2); arc1_y(end:-1:1,1)];
    pat2_x = [plan1_a(:,1); plan1_b(:,1); arc4_x(end:-1:1,1); plan2_a(:,1); plan2_b(:,1); arc2_x(end:-1:1,1)];
    pat2_y = [plan1_a(:,2); plan1_b(:,2); arc4_y(end:-1:1,1); plan2_a(:,2); plan2_b(:,2); arc2_y(end:-1:1,1)];
    comp = [pat1_x pat1_y];
    dilat = [pat2_x pat2_y];
end

% axis([-1 1 -1 1]);
% axis('square');
% axis off;
% hold on;

% patch('XData',comp(:,1), 'YData',comp(:,2),'FaceColor','None')
% patch('XData',dilat(:,1), 'YData',dilat(:,2),'FaceColor','None')
% patch(comp(:,1), comp(:,2), [0 0 0])
% patch(dilat(:,1), dilat(:,2), [1 1 1])
% -----------------

% X = [P.dip T.dip N.dip; P.str T.str N.str]';
% theta = (90-X(:,2)) * D2R;              % az converted to trignometric angle
% rho = sqrt(2)*sin(pi*(90-X(:,1))/360);  % projected distance from origin
% xp = rho .* cos(theta);
% yp = rho .* sin(theta);
% plot(xp(1),yp(1),'ok')
% plot(xp(2),yp(2),'ok','MarkerFaceColor','k')
% plot(xp(3),yp(3),'ok','MarkerFaceColor','r')

% axis('square');

% -----------------------------------------------------------------------
function [str2,dip2,rake2] = define_second_plane(NP1)
%     Compute strike, dip, slip for the second nodal plane
%     when are given strike, dip and rake for the first one.
%     Genevieve Patau

str2 = computed_strike1(NP1);
dip2 = computed_dip1(NP1);
rake2 = computed_rake1(NP1);

% -----------------------------------------------------------------------
function str2 = computed_strike1(NP1)
%    Compute the strike of the decond nodal plane when are given
%    strike, dip and rake for the first nodal plane with AKI & RICHARD's
%    convention. Angles are in degrees.

D2R = pi/180;
sr = sin(NP1.rake * D2R);     cr = cos(NP1.rake * D2R);
ss = sin(NP1.str * D2R);      cs = cos(NP1.str * D2R);
cd1 = cos(NP1.dip * D2R);
if (NP1.rake == 0),  am = 1;
else                am = NP1.rake / abs(NP1.rake);  end
if (cd1 < eps && abs(cr) < eps)
    str2 = NP1.str + 180;
else
    temp = cr * cs;
    temp = temp + sr * ss * cd1;
    sp2 = -am * temp;
    temp = ss * cr;
    temp = temp - sr *  cs * cd1;
    cp2 = am * temp;
    str2 = datan2(sp2, cp2);
    str2 = zero_360(str2);
end

% -----------------------------------------------------------------------
function dip2 = computed_dip1(NP1)
%    Compute second nodal plane dip when are given strike,
%    dip and rake for the first nodal plane with AKI & RICHARD's
%    convention. Angles are in degrees.
%    Genevieve Patau
if (NP1.rake == 0),  am = 1;
else                am = NP1.rake / abs(NP1.rake);  end
dip2 = acos(am * sin(NP1.rake*pi/180) * sin(NP1.dip*pi/180)) / (pi/180);
    
% -----------------------------------------------------------------------
function rake2 = computed_rake1(NP1)
%    Compute rake in the second nodal plane when strike ,dip
%    and rake are given for the first nodal plane with AKI & 
%    RICHARD's convention. Angles are in degrees.
%    Genevieve Patau

D2R = pi/180;
str2 = computed_strike1(NP1);
dip2 = computed_dip1(NP1);
if (NP1.rake == 0),  am = 1;
else                am = NP1.rake / abs(NP1.rake);  end
sd = sin(NP1.dip*D2R);          cd = cos(NP1.dip*D2R);
ss = sin((NP1.str-str2)*D2R);   cs = cos((NP1.str-str2)*D2R);

if(abs(dip2 - 90.) < eps)
    sinrake2 = am * cd;
else
    sinrake2 = -am * sd * cs / cd;
end
rake2 = datan2(sinrake2, -am * sd * ss);

% -----------------------------------------------------------------------
function str = zero_360(ang)
% put an angle between 0 and 360 degrees
if(ang >= 360.),    str = ang - 360.;
elseif(ang < 0.),    str = ang + 360.;
else               str = ang;
end

% -----------------------------------------------------------------------
function arctg = datan2(y,x)
%  compute arctg in degrees, between -180 et +180.

R2D = 180/pi;
if (abs(x) < eps)
    if(abs(y) < eps)
        arctg = 0.;
    elseif (y < 0)
        arctg = -90;
    else
        arctg = 90;
    end
elseif (x < 0)
    if (y < 0),  arctg = atan(y / x) * R2D - 180;
    else          arctg = atan(y / x) * R2D + 180;    end
else 
    arctg = atan(y / x) * R2D;
end

% -----------------------------------------------------------------------
function [P,T,N] = dc_to_axe(NP1,NP2)
% From FORTRAN routines of Anne Deschamps :
% compute azimuth and plungement of P-T axis
% from nodal plane strikes, dips and rakes.

D2R = pi / 180;
TWO_PI = 2 * pi;
M_PI_2 = pi / 2;
M_SQRT2 = sqrt(2);
pure_strike_slip = 0;

if (abs(sin(NP1.rake*D2R)) > eps)
    im = NP1.rake / abs(NP1.rake);
elseif (abs(sin(NP2.rake*D2R)) > eps)
    im = NP2.rake / abs(NP2.rake);
else
    pure_strike_slip = 1;
end

if (pure_strike_slip)
    if(cos(NP1.rake * D2R) < 0.)
        P.str = zero_360(NP1.str + 45);
        T.str = zero_360(NP1.str - 45);
    else
        P.str = zero_360(NP1.str - 45);
        T.str = zero_360(NP1.str + 45);
    end
    P.dip = 0;
    T.dip = 0;
else
    cd1 = cos(NP1.dip * D2R) *  M_SQRT2;
    sd1 = sin(NP1.dip * D2R) *  M_SQRT2;
    cd2 = cos(NP2.dip * D2R) *  M_SQRT2;
    sd2 = sin(NP2.dip * D2R) *  M_SQRT2;
    cp1 = - cos(NP1.str * D2R) * sd1;
    sp1 = sin(NP1.str * D2R) * sd1;
    cp2 = - cos(NP2.str * D2R) * sd2;
    sp2 = sin(NP2.str * D2R) * sd2;

    amz = - (cd1 + cd2);
    amx = - (sp1 + sp2);
    amy = cp1 + cp2;
    dx = atan2(sqrt(amx * amx + amy * amy), amz) - M_PI_2;
    px = atan2(amy, - amx);
    if (px < 0),    px = px + TWO_PI;  end

    amz = cd1 - cd2;
    amx = sp1 - sp2;
    amy = - cp1 + cp2;
    dy = atan2(sqrt(amx * amx + amy * amy), - abs(amz)) - M_PI_2;
    py = atan2(amy, - amx);
    if (amz > 0),   py = py - pi;   end
    if (py < 0),     py = py + TWO_PI;  end

    if (im == 1)
            P.dip = dy;     P.str = py;
            T.dip = dx;     T.str = px;
    else
            P.dip = dx;     P.str = px;
            T.dip = dy;     T.str = py;
    end
end

T.str = T.str / D2R; 
T.dip = T.dip / D2R;
P.str = P.str / D2R;
P.dip = P.dip / D2R;

N.dip = null_axis_dip(NP1.str, NP1.dip, NP2.str, NP2.dip);
N.str = null_axis_strike(NP1.str, NP1.dip, NP2.str, NP2.dip);

% -----------------------------------------------------------------------
function phn = null_axis_strike(str1,dip1,str2,dip2)
%    Compute null axis strike when strike and dip are given
%    for each nodal plane. Angles are in degrees.
%    Genevieve Patau

D2R = pi / 180;
sd1 = sin(dip1*D2R);    cd1 = cos(dip1*D2R);
sd2 = sin(dip2*D2R);    cd2 = cos(dip2*D2R);
ss1 = sin(str1*D2R);    cs1 = cos(str1*D2R);
ss2 = sin(str2*D2R);    cs2 = cos(str2*D2R);

cosphn = sd1 * cs1 * cd2 - sd2 * cs2 * cd1;
sinphn = sd1 * ss1 * cd2 - sd2 * ss2 * cd1;
if (sin((str1 - str2)*D2R) < 0)
    cosphn = -cosphn;
    sinphn = -sinphn;
end 
phn = datan2(sinphn, cosphn);
if (phn < 0),   phn = phn + 360;    end

% -----------------------------------------------------------------------
function den = null_axis_dip(str1,dip1,str2,dip2)
%    compute null axis dip when strike and dip are given
%    for each nodal plane. Angles are in degrees.
%    Genevieve Patau
  
D2R = pi / 180;
den = asin(sin(dip1*D2R) * sin(dip2*D2R) * sin( (str1 - str2)*D2R) ) / D2R; 
if (den < 0),    den = -den;    end
