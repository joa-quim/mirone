function [theta,ampfac]=nskew(yr,rlat,rlon,zobs,slin,sdec,sdip)

% NSKEW  Compute skewness parameter and amplitude factor
%  following Schouten (1971)
%  NEW VERSION THAT COMPUTES GEOCENTRIC DIPOLE UNLESS GIVEN
%  DEC AND DIP OF MAGNETIZATION      
% Usage: [theta,ampfac]=skew(yr,rlat,rlon,zobs,slin,sdec,sdip)
%  Input variables:
%   YR: decimal year of survey
%   RLAT: regional decimal latitude degrees
%   RLON: regional decimal longitude degrees
%   ZOBS: level of observation in km above sealevel
%   SLIN: strike of lineations normal to profile
%   SDEC: magnetization declination cw degrees from north
%   SDIP: magnetization inclination degrees
% Output variables
%   THETA: phase angle
%   AMPFAC: amplitude factor
% Calls <magfd>
% Maurice A. Tivey February 3, 1993
%    checked April 1996
%  Joaquim Luis         May   2004
%       Cleaned a bit of the code and replaced the igrf routine
%---------------------------------------------------------
D2R=pi/180;
% get unit vectors
%y = igrf(rlat, rlon, zobs, yr, 'igrf_1945_2005.dat');
y = igrf_m(rlon, rlat, zobs, yr);
bx = y(3);      by = y(4);      bz = y(5);      bh = y(2);
decl1 = y(6);   incl1 = y(7);   clear y;

% compute skewness parameter
% If a window for verbose exists
figdmsg = wfindobj('figure','tag','Wdmsgfig');
if (isempty(figdmsg))
    try,    message_win('create',' EARTH''S MAGNETIC FIELD DIRECTION:');    pause(0.001)
    catch,  fprintf(' EARTH''S MAGNETIC FIELD DIRECTION:\n');  end
else
    try,    message_win('add',' EARTH''S MAGNETIC FIELD DIRECTION:');
    catch,  fprintf(' EARTH''S MAGNETIC FIELD DIRECTION:\n');  end
end
try,    message_win('add',sprintf(' %10.3f = MAGNETIC DECLINATION ( STRIKE, CW FROM N )',decl1));
catch,  fprintf(' %10.3f = MAGNETIC DECLINATION ( STRIKE, CW FROM N )\n',decl1);    end
try,    message_win('add',sprintf(' %10.4f = MAGNETIC INCLINATION ( DIP, POS DOWN )',incl1));
catch,  fprintf(' %10.4f = MAGNETIC INCLINATION ( DIP, POS DOWN )\n',incl1);    end
if (nargin > 5)
    %  NOTE FOR GEOCENTRIC DIPOLE TAN(INC)=2*TAN(LAT)
    try,    message_win('add',' NON-GEOCENTRIC MAGNETIZATION VECTOR SPECIFIED:');
    catch,  fprintf(' NON-GEOCENTRIC MAGNETIZATION VECTOR SPECIFIED:\n');   end
    try,    message_win('add',sprintf(' %10.4f = DESIRED MAGNETIZATION DECLINATION (+CW FROM N)',sdec));
    catch,  fprintf(' %10.4f = DESIRED MAGNETIZATION DECLINATION (+CW FROM N)\n',sdec); end
    try,    message_win('add',sprintf(' %10.4f = DESIRED MAGNETIZATION INCLINATION (+DN)', sdip));
    catch,  fprintf(' %10.4f = DESIRED MAGNETIZATION INCLINATION (+DN)\n', sdip);   end
else %
    sdip = atan2( 2.*sin(rlat*D2R),cos(rlat*D2R) )/D2R;
    sdec = 0;
    try,    message_win('add',' GEOCENTRIC MAGNETIZATION VECTOR SPECIFIED:');
    catch,  fprintf(' GEOCENTRIC MAGNETIZATION VECTOR SPECIFIED:\n');   end
    try,    message_win('add',sprintf(' %10.4f = GEOCENTRIC DIPOLE INCLINATION',sdip));
    catch,  fprintf(' %10.4f = GEOCENTRIC DIPOLE INCLINATION \n',sdip); end
    try,    message_win('add',sprintf(' %10.3f = GEOCENTRIC DECLINATION ASSUMED\n',sdec));
    catch,  fprintf(' %10.3f = GEOCENTRIC DECLINATION ASSUMED\n',sdec); end
end
% compute phase and amplitude factors
ra1 = incl1*D2R;    rb1 = (decl1-slin)*D2R;
ra2 = sdip*D2R;     rb2 = (sdec-slin)*D2R;
% compute phase and amplitude factors
inclm = atan2(tan(ra2),sin(rb2));
inclf = atan2(tan(ra1),sin(rb1));
ampfac= ((sin(ra2))*(sin(ra1)))/((sin(inclm))*(sin(inclf)));
theta = (inclm/D2R)+(inclf/D2R)-180;
% compute unit vectors for a check
hatm(1) = cos(sdip*D2R)*sin((sdec-slin)*D2R);
hatm(2) = cos(sdip*D2R)*cos((sdec-slin)*D2R);
hatm(3) = -sin(sdip*D2R);
hatb(1) = cos(incl1*D2R)*sin((decl1-slin)*D2R);
hatb(2) = cos(incl1*D2R)*cos((decl1-slin)*D2R);
hatb(3) = -sin(incl1*D2R);
%
try,    message_win('add',sprintf('  %10.6f %10.6f %10.6f = MAGNETIZATION UNIT VECTOR',hatm(1),hatm(2),hatm(3)));
catch,  fprintf('  %10.6f %10.6f %10.6f = MAGNETIZATION UNIT VECTOR\n',hatm(1),hatm(2),hatm(3));    end
try,    message_win('add',sprintf('  %10.6f %10.6f %10.6f = AMBIENT FIELD UNIT VECTOR',hatb(1),hatb(2),hatb(3)));
catch,  fprintf('  %10.6f %10.6f %10.6f = AMBIENT FIELD UNIT VECTOR\n',hatb(1),hatb(2),hatb(3));    end
try,    message_win('add','  COMPONENTS ARE (X,Y,Z=ALONG, ACROSS PROFILE, AND UP');
catch,  fprintf('  COMPONENTS ARE (X,Y,Z=ALONG, ACROSS PROFILE, AND UP\n\n');   end
