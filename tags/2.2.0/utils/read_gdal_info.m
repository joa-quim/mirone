function out = read_gdal_info(file)
% ...

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

out.projection = [];    out.ellipsoid = [];     out.datum = [];     % when info is not provided/found
out.CTable = [];        out.NoData = [];        out.Zmin = [];      out.Zmax = [];
gotBand = 0;

% I found a subtle bug that manifested when the corner coords where glued. In that case the delimiter
% which the ',' was not recognized by the strtok function. That is why I use now the following definition
delims = ['\t' ' ' ',' '\n'];
[tok,rem] = strtok(file{1});       % This one we know it's here
out.Driver = rem;
for i=2:length(file)
    if ~isempty(findstr(file{i},'Size is'));    % Get image's nx, ny from a line of type "Size is 515, 515"
        [tok,rem] = strtok(file{i});    [tok,rem] = strtok(rem);    [tok,rem] = strtok(rem);
        out.Dim_nx = str2double(tok(1:end-1));   % Because of the trailing ","
        out.Dim_ny = str2double(rem);
    elseif ~isempty(findstr(file{i},'Pixel Size'));    % Get image's pixel size
        par_left = findstr(file{i},'(');     par_right = findstr(file{i},')');
        [tok,rem] = strtok(file{i}(par_left(1)+1:par_right(1)-1),[',' ' ']);
        out.pixel_x = str2double(tok);        out.pixel_y = str2double(rem(2:end));   % Because of the leading "," (Bug ?)
    elseif ~isempty(findstr(file{i},'UNIT'))      % There are 2 UNIT lines and I have to find how to desinguish them
    elseif ~isempty(findstr(file{i},'Upper Left'));
        par_left = findstr(file{i},'(');     par_right = findstr(file{i},')');
        [tok,rem] = strtok(file{i}(par_left(1)+1:par_right(1)-1),delims);
        out.UL_prj_xmin = str2double(tok);    out.UL_prj_ymax = str2double(rem);
        if length(par_left) == 2
            [tok,rem] = strtok(file{i}(par_left(2)+1:par_right(2)-1),delims);
            out.UL_geo_xmin = strGeo2numGeo(tok);        out.UL_geo_ymax = strGeo2numGeo(rem);
        else    out.UL_geo_xmin = [];   out.UL_geo_ymax = [];
        end
    elseif ~isempty(findstr(file{i},'Upper Right'))
        par_left = findstr(file{i},'(');     par_right = findstr(file{i},')');
        [tok,rem] = strtok(file{i}(par_left(1)+1:par_right(1)-1),delims);
        out.UR_prj_xmax = str2double(tok);    out.UR_prj_ymax = str2double(rem);
        if length(par_left) == 2
            [tok,rem] = strtok(file{i}(par_left(2)+1:par_right(2)-1),delims);
            out.UR_geo_xmax = strGeo2numGeo(tok);        out.UR_geo_ymax = strGeo2numGeo(rem);
        else    out.UL_geo_xmin = [];   out.UL_geo_ymax = [];
        end
    elseif ~isempty(findstr(file{i},'Lower Left'))
        par_left = findstr(file{i},'(');     par_right = findstr(file{i},')');
        [tok,rem] = strtok(file{i}(par_left(1)+1:par_right(1)-1),delims);
        out.LL_prj_xmin = str2double(tok);    out.LL_prj_ymin = str2double(rem);
        if length(par_left) == 2
            [tok,rem] = strtok(file{i}(par_left(2)+1:par_right(2)-1),delims);
            out.LL_geo_xmin = strGeo2numGeo(tok);        out.LL_geo_ymin = strGeo2numGeo(rem);
        else    out.UL_geo_xmin = [];   out.UL_geo_ymax = [];
        end
    elseif ~isempty(findstr(file{i},'Lower Right'))
        par_left = findstr(file{i},'(');     par_right = findstr(file{i},')');
        [tok,rem] = strtok(file{i}(par_left(1)+1:par_right(1)-1),delims);
        out.LR_prj_xmax = str2double(tok);    out.LR_prj_ymin = str2double(rem);
        if length(par_left) == 2
            [tok,rem] = strtok(file{i}(par_left(2)+1:par_right(2)-1),delims);
            out.LR_geo_xmax = strGeo2numGeo(tok);        out.LR_geo_ymin = strGeo2numGeo(rem);
        else    out.UL_geo_xmin = [];   out.UL_geo_ymax = [];
        end
    elseif ~isempty(findstr(file{i},'GEOGCS'))      % Get datum
        xx = findstr(file{i},'"');
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'NAD27');    out.datum = 'NAD27';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'NAD83');    out.datum = 'NAD83';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'WGS 84');    out.datum = 'WGS 84';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'WGS 72');    out.datum = 'WGS 72';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'Lisbon');    out.datum = 'Lisbon';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'ED50');    out.datum = 'ED50';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'ED87');    out.datum = 'ED87';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'ETRF89');    out.datum = 'ETRF89';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'Datum 73');    out.datum = 'Datum 73';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'Madrid 1870');    out.datum = 'Madrid 1870';    end
    elseif ~isempty(findstr(file{i},'SPHEROID'))      % Get ellipsoid
        xx = findstr(file{i},'"');                % Extractd from ellipsoid.cvs in gdal (probably from proj4)
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'Airy 1830');    out.datum = 'Airy 1830';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'Airy Modified 1849');    out.datum = 'Airy Modified 1849';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'Australian National Spheroid');    out.datum = 'Australian National Spheroid';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'Bessel 1841');    out.datum = 'Bessel 1841';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'Bessel Modified');    out.datum = 'Bessel Modified';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'Bessel Namibia');    out.datum = 'Bessel Namibia';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'Clarke 1858');    out.datum = 'Clarke 1858';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'Clarke 1866');    out.ellipsoid = 'Clarke 1866';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'Clarke 1866 Michigan');    out.ellipsoid = 'Clarke 1866 Michigan';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'Clarke 1880');    out.ellipsoid = 'Clarke 1880';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'Clarke 1880 (IGN)');    out.ellipsoid = 'Clarke 1880 (IGN)';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'Clarke 1880 (RGS)');    out.ellipsoid = 'Clarke 1866 (RGS)';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'Clarke 1880 (Arc)');    out.ellipsoid = 'Clarke 1880 (Arc)';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'Clarke 1880 (SGA 1922)');    out.ellipsoid = 'Clarke 1880 (SGA 1922)';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'Everest 1830 (1937 Adjustment)');    out.ellipsoid = 'Everest 1830 (1937 Adjustment)';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'Everest 1830 (1967 Definition)');    out.ellipsoid = 'Everest 1830 (1967 Definition)';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'Everest 1830 Modified');    out.ellipsoid = 'Everest 1830 Modified';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'Everest 1830 (1975 Definition)');    out.ellipsoid = 'Everest 1830 (1975 Definition)';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'GRS 1980');    out.ellipsoid = 'GRS 1980';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'Helmert 1906');    out.ellipsoid = 'Helmert 1906';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'Indonesian National Spheroid');    out.ellipsoid = 'Indonesian National Spheroid';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'International 1924');    out.ellipsoid = 'International 1924';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'Krassowsky 1940');    out.ellipsoid = 'Krassowsky 1940';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'NWL 9D');    out.ellipsoid = 'NWL 9D';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'Plessis 1817');    out.ellipsoid = 'Plessis 1817';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'Struve 1860');    out.ellipsoid = 'Struve 1860';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'War Office');    out.ellipsoid = 'War Office';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'WGS 84');    out.ellipsoid = 'WGS 84';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'GEM 10C');    out.ellipsoid = 'GEM 10C';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'OSU86F');    out.ellipsoid = 'OSU86F';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'OSU91A');    out.ellipsoid = 'OSU91A';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'Sphere');    out.ellipsoid = 'Sphere';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'GRS 1967');    out.ellipsoid = 'GRS 1967';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'Average Terrestrial System 1977');    out.ellipsoid = 'Average Terrestrial System 1977';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'WGS 72');    out.ellipsoid = 'WGS 72';    end
    elseif ~isempty(findstr(file{i},'PROJECTION'))      % Get projection
        xx = findstr(file{i},'"');
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'Albers_Conic_Equal_Area');    out.projection = 'Albers_Conic_Equal_Area';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'Azimuthal_Equidistant');    out.projection = 'Azimuthal_Equidistant';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'Cassini_Soldner');    out.projection = 'Cassini_Soldner';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'Cylindrical_Equal_Area');    out.projection = 'Cylindrical_Equal_Area';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'Eckert_IV');    out.projection = 'Eckert_IV';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'Eckert_VI');    out.projection = 'Eckert_VI';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'Equidistant_Conic');    out.projection = 'Equidistant_Conic';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'Equirectangular');    out.projection = 'Equirectangular';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'Transverse_Mercator');    out.projection = 'Transverse_Mercator';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'Gall_Stereographic');    out.projection = 'Gall_Stereographic';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'Gnomonic');    out.projection = 'Gnomonic';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'hotine_oblique_mercator');    out.projection = 'hotine_oblique_mercator';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'Krovak');    out.projection = 'Krovak';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'Laborde_Oblique_Mercator');    out.projection = 'Laborde_Oblique_Mercator';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'Lambert_Azimuthal_Equal_Area');    out.projection = 'Lambert_Azimuthal_Equal_Area';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'Lambert_Conformal_Conic_1SP');    out.projection = 'Lambert_Conformal_Conic_1SP';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'Lambert_Conformal_Conic_2SP');    out.projection = 'Lambert_Conformal_Conic_2SP';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'Lambert_Conformal_Conic_2SP_Belgium');    out.projection = 'Lambert_Conformal_Conic_2SP_Belgium';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'Cylindrical_Equal_Area');    out.projection = 'Cylindrical_Equal_Area';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'Mercator_1SP');    out.projection = 'Mercator_1SP';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'Mercator_2SP');    out.projection = 'Mercator_2SP';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'Miller_Cylindrical');    out.projection = 'Miller_Cylindrical';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'Mollweide');    out.projection = 'Mollweide';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'New_Zealand_Map_Grid');    out.projection = 'New_Zealand_Map_Grid';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'Oblique_Mercator');    out.projection = 'Oblique_Mercator';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'Oblique_Stereographic');    out.projection = 'Oblique_Stereographic';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'Orthographic');    out.projection = 'Orthographic';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'Polar_Stereographic');    out.projection = 'Polar_Stereographic';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'Polyconic');    out.projection = 'Polyconic';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'Robinson');    out.projection = 'Robinson';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'CT_ObliqueMercator_Rosenmund');    out.projection = 'CT_ObliqueMercator_Rosenmund';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'Sinusoidal');    out.projection = 'Sinusoidal';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'Swiss_Oblique_Cylindrical');    out.projection = 'Swiss_Oblique_Cylindrical';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'Swiss_Oblique_Mercator');    out.projection = 'Swiss_Oblique_Mercator';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'Stereographic');    out.projection = 'Stereographic';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'Transverse_Mercator');    out.projection = 'Transverse_Mercator';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'CT_TransverseMercator_Modified_Alaska');    out.projection = 'CT_TransverseMercator_Modified_Alaska';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'Transverse_Mercator_South_Orientated ');    out.projection = 'Transverse_Mercator_South_Orientated ';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'Tunisia_Mining_Grid');    out.projection = 'Tunisia_Mining_Grid';    end
        if strcmp(file{i}(xx(1)+1:xx(2)-1),'VanDerGrinten');    out.projection = 'VanDerGrinten';    end
    elseif (~isempty(findstr(file{i},'Band')) & ~gotBand)        % Get colormap and Type
        gotBand = 1;                                % The 'Band' word may come out more than once
        xx = findstr(file{i},'ColorInterp=');
        switch file{i}(xx+12:end)                   % Get colormap
            case 'Gray'
                out.Cmap = 'gray';
            case 'Palette'
                out.Cmap = 'Palette';
                out.CTable = parse_color_table(file,i+1);
            case {'Red','Green','Blue'}              % Some Tiffs 'Red', 'Green', 'Blue'. Just save the 'Red'
                out.Cmap = 'Red';
            otherwise
                out.Cmap = [];                
        end
        xx = findstr(file{i},'Type=');
        tok = strtok(file{i}(xx+5:end));
        switch tok                   % Get Type
            case {'Byte','Byte,'}
                out.Type = 'Byte';
            case {'Int16','Int16,'}
                out.Type = 'Int16';
            case {'UInt16','UInt16,'}
                out.Type = 'UInt16';
            case {'Int32','Int32,'}
                out.Type = 'Int32';
            case {'UInt32','UInt32,'}
                out.Type = 'UInt32';
            case {'Float','Float,','Float32','Float32,'}
                out.Type = 'Float';
            otherwise
                out.Type = [];                
        end
    elseif ~isempty(findstr(file{i},'NoData Value'))        % Get empty grid value
        xx = findstr(file{i},'NoData Value=');
        out.NoData = file{i}(xx+13:end);                    % It goes out as a text string
    elseif ~isempty(findstr(file{i},'Computed Min/Max'))    % Get Band min/max (more than 1 Band and ??)
        xx = findstr(file{i},'Min/Max=');
        yy = findstr(file{i}(xx+8:end),',');
        out.Zmin = file{i}(xx+8:xx+8+yy(1)-1);              % They go out as a text string
        out.Zmax = file{i}(xx+8+yy(1):end);
    end
end

% Save the lines with the Coordinate System for use in the construction of a .prj file
i_start = 0;    i_stop = 0;
for i=2:length(file)
    if ~isempty(findstr(file{i},'Coordinate System'))
        i_start = i + 1;
    end
    if ~isempty(findstr(file{i},'Origin ='))
        i_stop = i - 1;     break
    end
end
try
    if (i_stop - i_start > 0),  out.Coord_System = file(i_start:i_stop);
    else                        out.Coord_System = [];      end
catch
    out.Coord_System = [];
end

if isempty(out.projection);     out.projection = 'Unknown';     end
if isempty(out.ellipsoid);      out.ellipsoid = 'Unknown';     end
if isempty(out.datum);          out.datum = 'Unknown';     end

% --------------------------------------------------------------------
function deg_dec = strGeo2numGeo(str)
% Translate the string with geo coordinates in the form of (e.g.) '117d38'28.22"W,'
% to the corresponding decimal form
    
% First get read of eventual trailing "N, or "W, or "S or ... BUT get the sign
sign = 1;
if findstr(str,'W');    sign = -1;  end
if findstr(str,'S');    sign = -1;  end
str = str(1:findstr(str,'"')-1);

deg = str(1:findstr(str,'d')-1);
dec_sex = str(findstr(str,'d')+1:end);
minut = dec_sex(1:findstr(dec_sex,'''')-1);
sec = str(findstr(str,'''')+1:end);
deg_dec = str2double(deg)*sign + (str2double(minut)/60)*sign + str2double(sec)/3600*sign;

% --------------------------------------------------------------------
function ctable = parse_color_table(file,i)
% We are searching for the number of entries in a string like "Color Table (RGB with 256 entries)"
j = 1;  n = i + 1;

try
    while isempty(findstr(file{i},'Color Table'))
        n = n + 1;      i = i + 1;
    end
catch   % Didn't find the 'Color Table' string. 
    ctable = [];    return
end

xx = findstr(file{i},'with');
n_colors = str2double( strtok(file{i}(xx+5:end)) );      % Now we know the number color entries
ctable = zeros(n_colors,3);

for k = n:n+n_colors-1
    [tok,rem] = strtok(file{k});
    if strcmp(tok(end),':')
        xx = findstr(rem,',');
        if (length(xx) ~= 3)
            msgbox('Different than assumed format in Color Table. Error on the way','Warning')
        end
        r = rem(1:xx(1)-1);     g = rem(xx(1)+1:xx(2)-1);    b = rem(xx(2)+1:xx(3)-1);
        ctable(j,1) = str2double(r);
        ctable(j,2) = str2double(g);    ctable(j,3) = str2double(b);
        j = j + 1;        n = n + 1;
    else
        % ??
    end
end
