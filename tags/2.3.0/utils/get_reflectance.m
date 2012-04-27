function R = get_reflectance(x,y,z,s,k,limits)
%GET_REFLECTANCE(X,Y,Z,S) 3-D shaded surface with lighting.
%
%   The shading is based on a combination of diffuse, specular and ambient lighting models.
%
%   The relative contributions due to ambient light, diffuse reflection, specular reflection,
%   and the specular spread coefficient can be set by using five arguments
%   GET_REFLECTANCE(X,Y,Z,S,K) where K=[ka,kd,ks,spread]. S is the three vector S = [Sx,Sy,Sz]
%   that specifies the direction of the light source. S can also be
%   specified in view coordinates, S = [AZ,EL].

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

if (nargin < 5)
    error('Incorrect number of input arguments.');
elseif (nargin == 5)    % Limits will be searched in gcf appdata
    limits = [];
elseif (nargin == 6)
    a = limits;
end

if (length(k)~=4)
  error('Weighting vector k must have four components.');
end

[msg,x,y,z] = xyzchk(x,y,z); if ~isempty(msg), error(msg); end
if any(size(z)<[3 3]), error('X, Y, and Z must be at least 3-by-3.'); end

vaz = 0;  vel = 90;     % Direction to the viewer

if (length(s)~=2) && (length(s)~=3)
    error('S must be specified using [AZ,EL] or [Sx,Sy,Sz].');
end

ms = length(s(:));
if ms==2, % Compute source direction from [AZ,EL]
    az = s(1)*pi/180; el = s(2)*pi/180; % Convert to radians
    s = zeros(1,3);
    s(1) =  sin(az)*cos(el);
    s(2) = -cos(az)*cos(el);
    s(3) =  sin(el);
end

% Determine plot scaling factors for a cube-like plot domain.
if (isempty(limits))    % gcf appdata must contain the following info
    xx = getappdata(gcf,'Xmin_max');  yy = getappdata(gcf,'Ymin_max');   zz = getappdata(gcf,'Zmin_max');
    a = [xx yy zz];
end
Sx = a(2)-a(1);     Sy = a(4)-a(3);     Sz = a(6)-a(5);
scale = max([Sx,Sy,Sz]);
Sx = Sx/scale; Sy = Sy/scale; Sz = Sz/scale;

% Compute surface normals.  Rely on ordering to define inside or outside.
x = x/Sx; y = y/Sy; z = double(z)/Sz;       % 'double' because z may well be of different type

% Horizontally tile the grids when they are large. This allows saving lots of memory.
m = size(z,1);      % I probably new this already, but it's beter play safe
[ind_s,ind] = tile(m,200,2);
if size(ind_s,1) > 1
    R = [];
    for i = 1:size(ind_s,1)
        tmp1 = (ind_s(i,1):ind_s(i,2));     % Indexes with overlapping zone
        tmp2 = ind(i,1):ind(i,2);           % Indexes of chunks without the overlaping zone
        [nx ny nz] = getnormals(x(tmp1,:),y(tmp1,:),z(tmp1,:));
        nx = nx(tmp2,:);    ny = ny(tmp2,:);    nz = nz(tmp2,:);
        tmp = (k(1)+k(2)*diffuse(nx,ny,nz,s)+k(3)*specular(nx,ny,nz,s,[vaz,vel],k(4)))/ sum(k(1:3));
        R = [R; tmp];
    end
    clear x y z nx ny nz;
else
    % Compute Lambertian shading + specular + ambient light
    [nx,ny,nz] = getnormals(x,y,z);     clear x y z;
    R = (k(1)+k(2)*diffuse(nx,ny,nz,s)+k(3)*specular(nx,ny,nz,s,[vaz,vel],k(4)))/ sum(k(1:3));
end

%--------------------------------------------------------------------------
function r = specular(nx,ny,nz,s,v,k)
%SPECULAR Specular reflectance.
%   R = SPECULAR(Nx,Ny,Nz,S,V) returns the reflectance of a surface with
%   normal vector components [Nx,Ny,Nz].  S and V specify the direction
%   to the light source and to the viewer, respectively.

D2R = pi / 180;     [m,n] = size(nx);
s = s(:)/norm(s);   % Normalize
v = [sin(v(1)*D2R)*cos(v(2)*D2R); -cos(v(1)*D2R)*cos(v(2)*D2R); sin(v(2)*D2R)];

% mag = sqrt(nx.*nx+ny.*ny+nz.*nz);   % Normalize normal vectors.
% d = find(mag==0); mag(d) = eps*ones(size(d));

%r = max(0,2*(s(1)*nx+s(2)*ny+s(3)*nz).*(v(1)*nx+v(2)*ny+v(3)*nz)./mag - (v'*s)*ones(m,n));
r = max(0,2*(s(1)*nx+s(2)*ny+s(3)*nz).*(v(1)*nx+v(2)*ny+v(3)*nz) - (v'*s)*ones(m,n));
r = r.^k;

%--------------------------------------------------------------------------
function r = diffuse(nx,ny,nz,s)
%DIFFUSE Diffuse reflectance.
%   R = DIFFUSE(Nx,Ny,Nz,S) returns the reflectance for a surface with
%   normal vector components [Nx,Ny,Nz].  S is a three vector that
%   defines the direction to the light source. S can also be a two vector
%   S = [Theta Phi] specifying the direction in spherical coordinates.
%
%   Lambert's Law: R = cos(PSI) where PSI is the angle between the
%   surface normal and light source.
%
s = s/norm(s); % Normalize
% mag = sqrt(nx.*nx+ny.*ny+nz.*nz);   % Normalize normal vectors.
% d = find(mag==0); mag(d) = eps*ones(size(d));
% r = max(0,(s(1)*nx + s(2)*ny + s(3)*nz)./mag);
r = max(0,(s(1)*nx + s(2)*ny + s(3)*nz));
