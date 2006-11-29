function write_flederFiles(opt,varargin)
% 'flederize' Matlab graphical objects. 
% NOTE: For every line/point object FM_CMAP and a GEOREF blocks are writen.
% Though it doesn't hurt much, it is an idiot thing

%	Copyright (c) 2004-2006 by J. Luis
%
%	This program is free software; you can redistribute it and/or modify
%	it under the terms of the GNU General Public License as published by
%	the Free Software Foundation; version 2 of the License.
%
%	This program is distributed in the hope that it will be useful,
%	but WITHOUT ANY WARRANTY; without even the implied warranty of
%	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%	GNU General Public License for more details.
%
%	Contact info: w3.ualg.pt/~jluis/mirone
% --------------------------------------------------------------------

switch opt
    case 'geo'
        write_geo(varargin{:});     % Write a .geo block
    case 'dtm'
        write_dtm(varargin{:});     % Write a .dtm object
    case 'shade'
        write_shade(varargin{:});   % Write a .shade object
    case 'main_SD'
        write_main(varargin{:});    % Write a .sd object made of .geo .dtm & .shade blocks
    case 'line_or_points'
        write_lines_or_points(varargin{:});    % Search for and write lines and/or points objects
    case 'line'
        write_line(varargin{:});    % Write line objects
    case 'points'
        write_pts(varargin{:});     % Write point objects
    case 'cmap'
        write_cmap(varargin{:});    % Write a FM_CMAP block
    case 'eof'
        write_eof(varargin{1})      % Write EOF block and close file
end

%----------------------------------------------------------------------------------
function write_geo(fid,mode,limits)
    % Write a GEOREF block
	if (strcmp(mode,'first'))       % The TDR object starts here
        fprintf(fid,'%s\n%s\f\n','%% TDR 2.0 Binary','%%');
	end
	fwrite(fid,[15000 100],'integer*4');     % Tag ID, Data Length
	fwrite(fid,[0 0 1 1 1 1 (1:30)*0],'integer*1');
	fwrite(fid,limits,'real*8');
	fwrite(fid,(1:10)*0,'integer*4');

%----------------------------------------------------------------------------------
function write_dtm(fid,mode,Z,limits)
    % Write a .dtm block
    if (strcmp(mode,'first'))       % The TDR object starts here
	    fprintf(fid,'%s\n%s\f\n','%% TDR 2.0 Binary','%%');
    end
    [m,n] = size(Z);
	fwrite(fid,[1000 m*n*2+30],'integer*4');     % Tag ID, Data Length 
	fwrite(fid,[0 0 1 1 1 2],'integer*1');
	fwrite(fid,[(1:9)*0 2 n m 3 16],'integer*2');      % ?? nDim(?) nCols nRows ?? BitWidth
	fwrite(fid,[limits(5) limits(6)],'real*8');
	fwrite(fid,[0 65535],'uint16');
    Z = flipud( uint16(rot90(scaleto8(Z,16),1)) );
	fwrite(fid,Z,'uint16');

%----------------------------------------------------------------------------------
function write_shade(fid,mode,hFig,hAxes,burnCoasts)
    % Write a .shade block
    if (strcmp(mode,'first'))       % The TDR object starts here
	    fprintf(fid,'%s\n%s\f\n','%% TDR 2.0 Binary','%%');
    end
	img = get(findobj(hAxes,'Type','image'),'CData');
    [m,n,k] = size(img);
	fwrite(fid,[1000 m*n*4+34],'integer*4');     % Tag ID, Data Length
	fwrite(fid,[0 0 1 1 1 2],'integer*1');
	fwrite(fid,[(1:9)*0 2 n m 7 32],'integer*2');      %?? nDim(?) nCols nRows ?? BitWidth
	fwrite(fid,[(1:10)*0 65535 65535],'uint16');
    inplace = false;            % Eventual line burning is not done inplace to not change the Mirone image as well
	if (ndims(img) == 2)
        img = ind2rgb8(img,get(hFig,'Colormap'));
        inplace = true;         % Save to do line burning inplace because img is already a copy
	end
    
    if (burnCoasts)             % Burn coastlines into img (in case they exist)
        img = burnLines(hFig,hAxes,img,inplace);
    end
	img(:,:,4) = uint8(255);                % Tranparency layer
	img = permute(img,[4 3 2 1]);
	fwrite(fid,img,'uint8');

%----------------------------------------------------------------------------------
function write_cmap(fid)
% Write the FM_CMAP block
    fwrite(fid,20010,'integer*4');              % ID of FM_CMAP block
    fwrite(fid,768,'integer*4');                % n of bytes in this block
    fwrite(fid,[0 0 1 1 1 1 (1:18)*0],'integer*1');     % 24 bytes (seams to be an offset)
    pal = uint8(round(jet(256) * 255));         % Make a colormap
    fwrite(fid,pal','uchar');    

%----------------------------------------------------------------------------------
function write_eof(fid)
% Write the EOF block
 	fwrite(fid,999999999,'integer*4');
 	fwrite(fid,[(1:6)*0 1 1 1 1 (1:18)*0],'uchar');
	fclose(fid);

%----------------------------------------------------------------------------------
function write_main(fid,type,hFig,hAxes,Z,limits,burnCoasts)
% Write a basic .sd file. That is with a DTM, a SHADE & a GEO blocks

	fprintf(fid,'%s\n%s\n%s\f\n','%% TDR 2.0 Binary','Created by:     Mirone   ','%%');
    if (strcmp(type,'writePlanarSD') || strcmp(type,'runPlanarSD'))
	    fwrite(fid,[21 39 0 0 2 0 0 0 0 0 1 1 1 1 (1:20)*0],'integer*1');
    else
	    fwrite(fid,[31 39 0 0 0 0 0 0 0 0 1 1 1 1 (1:18)*0],'integer*1');
    end
    
    write_flederFiles('dtm',fid,'add',Z,limits)
    write_flederFiles('shade',fid,'add',hFig,hAxes,burnCoasts)
    write_flederFiles('geo',fid,'add',limits)

%----------------------------------------------------------------------------------
function write_lines_or_points(fid,hFig,hAxes,Z,limits,burnCoasts)
% Look for line and/or points Matlab objects and write them as Fledermaus
% objects. It does so (when it finds them) by calling the corresponding
% function that writes either lines or points objects in Fleder format

if (nargin < 6),    burnCoasts = 1;     end     % Default is to burn coast lines into image
ALLlineHand = findobj(hAxes,'Type','line');
if (~isempty(ALLlineHand))
    h = findobj(ALLlineHand,'Tag','Earthquakes');       % Search first for earthquakes because they have depths
    if (~isempty(h))
        write_flederFiles('points',fid,h,'add',limits,'Earthquakes')   % earthquakes
        ALLlineHand = setxor(ALLlineHand, h);           % h is processed, so remove it from handles list
    end

    if (burnCoasts)     % If we had already burned the coastlines remove them from the ALLlineHand list
        h = findobj(ALLlineHand,'Tag','CoastLineNetCDF');
        if (~isempty(h)),       ALLlineHand = setxor(ALLlineHand, h);     end
        h = findobj(ALLlineHand,'Tag','PoliticalBoundaries');
        if (~isempty(h)),       ALLlineHand = setxor(ALLlineHand, h);     end
        h = findobj(ALLlineHand,'Tag','Rivers');
        if (~isempty(h)),       ALLlineHand = setxor(ALLlineHand, h);     end
    end
    
    % See if we have COASTLINES. If yes they are treated separatly (mainly because od Z and also to create an separate object)
    h = findobj(ALLlineHand,'Tag','CoastLineNetCDF');
    if (~isempty(h))
        z_level = 0;        % It will be a nonsense if the underlying grid is not topographic
        [x,y,z,count] = lines2multiseg(h,z_level);
        line_thick = get(h(1),'LineWidth');       % Line thickness
        line_color = get(h(1),'color');           % Line color
        line_props = [line_thick line_color 0 [1 1 1]];     % 0 means is not a patch and the [1 1 1] is not used
        write_flederFiles('line',fid,'add',x,y,z,count,limits,line_props)
        ALLlineHand = setxor(ALLlineHand, h);     % h are processed, so remove them from handles list
    end
    
    % See if we have national borders. If yes they are treated separatly
    h = findobj(ALLlineHand,'Tag','PoliticalBoundaries');
    if (~isempty(h))
        z_level = 0;        % It will be a nonsense if the underying grid is not topographic
        [x,y,z,count] = lines2multiseg(h,z_level);
        line_thick = get(h(1),'LineWidth');       % Line thickness
        line_color = get(h(1),'color');           % Line color
        line_props = [line_thick line_color 0 [1 1 1]];     % 0 means is not a patch and the [1 1 1] is not used
        write_flederFiles('line',fid,'add',x,y,z,count,limits,line_props)
        ALLlineHand = setxor(ALLlineHand, h);     % h are processed, so remove them from handles list
    end
    
    % See if we have RIVERS. If yes they are treated separatly
    h = findobj(ALLlineHand,'Tag','Rivers');
    if (~isempty(h))
        z_level = 0;        % It will be a nonsense if the underying grid is not topographic
        [x,y,z,count] = lines2multiseg(h,z_level);
        line_thick = get(h(1),'LineWidth');       % Line thickness
        line_color = get(h(1),'color');           % Line color
        line_props = [line_thick line_color 0 [1 1 1]];     % 0 means is not a patch and the [1 1 1] is not used
        write_flederFiles('line',fid,'add',x,y,z,count,limits,line_props)
        ALLlineHand = setxor(ALLlineHand, h);     % h are processed, so remove them from handles list
    end
    
    % See if we have CONTOUR lines. If yes we fish their depths
    h = findobj(ALLlineHand,'Tag','contour');
    if (~isempty(h))
        dz = abs(limits(6) - limits(5)) * 0.01;     % I smell a fleder bug here, so add a small cte to z level
        for (i = 1:length(h))
            z_level = get(h(i),'UserData') + dz;
            [x,y,z,count] = lines2multiseg(h(i),z_level);
            line_thick = get(h(i),'LineWidth');     % Line thickness
            line_color = get(h(i),'color');         % Line color
            line_props = [line_thick line_color 0 [1 1 1]];     % 0 means is not a patch and the [1 1 1] is not used
            write_flederFiles('line',fid,'add',x,y,z,count,limits,line_props)
        end
        ALLlineHand = setxor(ALLlineHand, h);       % h are processed, so remove them from handles list
    end
end

if (~isempty(ALLlineHand))
    % This section deals with repeated line types. The point is that having many objects makes the
    % rendering slow. So I assimilate all lines of the same type (same line thickness and color) into
    % a single mulrisegment line. This makes the rendering much faster.
    LineWidth = get(ALLlineHand,'LineWidth');
    if (iscell(LineWidth)),     LineWidth = cell2mat(LineWidth);    end    
    LineColor = get(ALLlineHand,'Color');
    if (iscell(LineColor)),     LineColor = cell2mat(LineColor);    end
    [tmp,ind] = sortrows([LineWidth LineColor]);
    hands_sort = ALLlineHand(ind);          % Sort also the handles according to the previous sorting cretirea
    difs = diff([tmp(1,:); tmp]);           % Repeat first row to account for the decrease 1 resulting from diff
    [id_row,j] = find(difs ~= 0);           % Find the lines of different type
    id_row = unique(id_row);                % Get rid of repeated values
    id_row = [1; id_row];                   % Make id_row start at one
    id_row(end+1) = size(tmp,1);            % Add the last row as well (for the algo)
    for (i = 1:length(id_row)-1)
        hands{i} = hands_sort(id_row(i):id_row(i+1)-1);
        if (length(hands{i}) > 1)           % Make sure that the following procedure applyies only to repeated line types
            line_thick = get(hands{i}(1),'LineWidth');          % Line thickness
            line_color = get(hands{i}(1),'color');              % Line color
            line_props = [line_thick line_color 0 [1 1 1]];     % 0 means is not a patch and the [1 1 1] is not used
            [x,y,z,count] = lines2multiseg([hands{i}],limits(end));
            write_flederFiles('line',fid,'add',x,y,z,count,limits,line_props)
            ALLlineHand = setxor(ALLlineHand, [hands{i}]);   % hands{i} are processed, so remove them from handles list
        end
    end
        
    % OK, now if we still have lines, they must be of different line type
    for (i = 1:length(ALLlineHand))
        z_level = limits(end);          % Default to z_max (but I have to do something clever)
        [x,y,z,count] = lines2multiseg(ALLlineHand(i),z_level);
        line_thick = get(ALLlineHand(i),'LineWidth');   % Line thickness
        line_color = get(ALLlineHand(i),'color');       % Line color
        line_props = [line_thick line_color 0 [1 1 1]];     % 0 means is not a patch and the [1 1 1] is not used
        write_flederFiles('line',fid,'add',x,y,z,count,limits,line_props)        
    end
end     % end  -> if (~isempty(ALLlineHand))<-

ALLpatchHand = findobj(hAxes,'Type','patch');

if (~isempty(ALLpatchHand))
    z_level = limits(end);          % Default to z_max (but I have to do something clever)

    telhasHand_d = findobj(ALLpatchHand,'Tag','tapete');    % Fish the direct telhas patches
    telhasHand_r = findobj(ALLpatchHand,'Tag','tapete_R');  % And the reverse telhas patches
    line_thick = 2;             % Line thickness
    line_color = [0 0 0];       % Line color
    if (~isempty(telhasHand_d))                     % First the direct telhas
        patch_color = get(telhasHand_d(1),'FaceColor');
        n_vert = 5000;                              % For pre-allocation
        xx = cell(5000,1);  yy = xx;    zz = xx;    % Pre-allocate memory in excess
        count = 1;
        for (i = 1:length(telhasHand_d))            % Loop over direct polarity tapetes
            x = get(telhasHand_d(i),'XData');       y = get(telhasHand_d(i),'YData');
            for (k = 1:size(x,2))                   % Loop over individual telhas in each tapete
                xx1 = x(:,k)';      xx1(end+1) = xx1(1);  
                yy1 = y(:,k)';      yy1(end+1) = yy1(1);
                xx{count} = xx1;    yy{count} = yy1;
                zz{count} = repmat(z_level,1,5);
                count = count + 1;
            end
        end
        count = count - 1;                           % - 1 because the counter incremented one too much
        line_props = [line_thick line_color 1 patch_color];     % The 1 indicates this is a patch
        % Remove eventually unused
        if (count < n_vert)
            xx(count+1:end) = [];    yy(count+1:end) = [];    zz(count+1:end) = [];
        end
        write_flederFiles('line',fid,'add',xx,yy,zz,count,limits,line_props)
        ALLpatchHand = setxor(ALLpatchHand, telhasHand_d);  % telhasHand_d is processed, so remove it from handles list
    end
    
    if (~isempty(telhasHand_r))                     % Now the reverse telhas
        patch_color = get(telhasHand_r(1),'FaceColor');
        xx = cell(5000,1);  yy = xx;    zz = xx;    % Pre-allocate memory in excess
        count = 1;
        for (i = 1:length(telhasHand_r))            % Loop over direct polarity tapetes
            x = get(telhasHand_r(i),'XData');       y = get(telhasHand_r(i),'YData');
            for (k = 1:size(x,2))                   % Loop over individual telhas in each tapete
                xx1 = x(:,k)';      xx1(end+1) = xx1(1);  
                yy1 = y(:,k)';      yy1(end+1) = yy1(1);
                xx{count} = xx1;    yy{count} = yy1;
                zz{count} = repmat(z_level,1,5);
                count = count + 1;
            end
        end
        count = count - 1;                           % - 1 because the counter incremented one too much
        line_props = [line_thick line_color 1 patch_color];     % The 1 indicates this is a patch
        % Remove eventually unused
        if (count < n_vert)
            xx(count+1:end) = [];    yy(count+1:end) = [];    zz(count+1:end) = [];
        end
        write_flederFiles('line',fid,'add',xx,yy,zz,count,limits,line_props)
        ALLpatchHand = setxor(ALLpatchHand, telhasHand_r);  % telhasHand_r is processed, so remove it from handles list
        clear xx yy zz;
    end
end
       
if (~isempty(ALLpatchHand))                 % Now see if we still have more patches
    n_lines = length(ALLpatchHand);
    for (i = 1:n_lines)                     % Do one by one. If they are many, the rendering is slow
        x = get(ALLpatchHand(i),'XData');        y = get(ALLpatchHand(i),'YData');
        line_thick = get(ALLpatchHand(i),'LineWidth');   % Line thickness
        line_color = get(ALLpatchHand(i),'EdgeColor');       % Line color
        
        patch = 1;
        patch_color = get(ALLpatchHand(i),'FaceColor');
        if (strcmp(patch_color,'none'))     % OK, if we have no color we treat it just like an ordinary polyline
            patch = 0;
            patch_color = [1 1 1];          % Not used anyway
        end
        line_props = [line_thick line_color patch patch_color];
        count = length(x);
        z = repmat(z_level,1,count);
        write_flederFiles('line',fid,'add',x,y,z,count,limits,line_props)
    end
end

%----------------------------------------------------------------------------------
function write_line(fid,mode,x,y,z,np,lim_reg,line_props)
% Build an line TDR object
% Se tiver multisegs basta voltar a escrever o np do segmento seguinte + um 0 int*4 e continuar (e o n_byte?)
% Nota o np quase de certeza que e int*8
% LIM_REG   -> is the -R of the map (not of the line, which may be > or <)
% NP        -> Total number of points in this line (if multisegment NP is still the total number of pts)
% LINE_PROPS -> vector with line properties [thickness [color] patch [faceColor]], where color = [r g b]
    ColorBy = 0;                        % 0 -> Solid; 1 -> Line Height (Z); 2 -> Attribute
    code1 = 1;                          % Still don't know what this codes (number of ??)
    lim_line = lim_reg;                 % MERDOSO (nao e assim se lim_reg < lim_line)
    l_thick = line_props(1);
    patch = line_props(5);              % See if we have a polyline or a patch
    if (patch)
        l_color = round(line_props(6:8) * 255);
    else
        l_color = round(line_props(2:4) * 255);
    end

    if (strcmp(mode,'first'))       % The TDR object starts here
	    fprintf(fid,'%s\n%s\f\n','%% TDR 2.0 Binary','%%');
    end
    
    if (iscell(x)),     n_segments = length(x);
    else                n_segments = 1;
    end
    
    % Make sure x & y are row vectors
    if (size(x,1) > 1), x = x';     end
    if (size(y,1) > 1), y = y';     end
    
    % In next line the 2 * 6 term is due to the fact that 'lim' is written twice
    n_byte = 2*6*8 + 3*np*8 + (n_segments-1)*8; % (n_segments-1)*8 accounts for the 'np + 1 zero int*4' for each seg
    n_byte = n_byte + 5*4 + 8+4 + 2*4;          % 
    fwrite(fid,[10525 n_byte],'integer*4');     % ID of block SD_LINES3D and n of bytes in this block
  	fwrite(fid,[0 0 1 1 1 3 (1:18)*0],'integer*1'); % 28 bytes (not counted in n_byte)
    fwrite(fid,[patch n_segments np 0 0],'integer*4');     % n points and some code
  	fwrite(fid,ones(1,8)*205,'uchar');          % ??
   	fwrite(fid,code1,'integer*4');              % ??
	fwrite(fid,[lim_reg lim_line],'real*8');

    if (iscell(x))
		for (i = 1:n_segments)
            np_s = length(x{i});                % Number of points in this segment
            fwrite(fid,[np_s 0],'integer*4');
			fwrite(fid,[x{i}; y{i}; z{i}],'real*8');
		end
    else
        fwrite(fid,[np 0],'integer*4');   % n points
		fwrite(fid,[x; y; z],'real*8');
    end
    
    write_flederFiles('geo',fid,'add',lim_reg)  % Write a GEOREF block
    write_flederFiles('cmap',fid)               % Write a FM_CMAP block

    fwrite(fid,[10526 32],'integer*4');         % ID of block SD_LINES3D_ATB and n bytes in this block
 	fwrite(fid,[0 0 1 1 1 6 (1:18)*0],'integer*1');     % The 6 is a number of version
 	
 	fwrite(fid,[l_color(1:3) 255],'uchar');     % Don't know what is the last 255
    fwrite(fid,[1 1 ColorBy 0 0 1 l_thick],'integer*4'); % First 1 is 'gap', but don't know what are the others

%----------------------------------------------------------------------------------
function write_pts(fid,hand,mode,limits,opt)
% HAND -> handles of the line (points) object
% MODE = FIRST or ADD
% OPT, when it exists, is = 'Earthquakes'
%   TAMBEM NAO SEI O QUE ISTO FAZ SE HOUVER MAIS DE UM CONJUNTO DE PONTOS
    if (nargin == 4),   opt = [];   end
    
    symb = 3;               % 0 -> circle; 1 -> square; 2 -> cross hair; 3 -> cube; 4 -> cylinder; 5 -> sphere; 6 -> point
    PointRad = 0.02;                        % Symbol radius
    LabelSize = 0.502;
    ColorBy = 0;                            % 0 -> Solid; 1 -> Line Height (Z); 2 -> Attribute
    n_col = 3;                              % N of columns
    n_groups = length(hand);                % N of different point ensembles

    if (strcmp(mode,'first'))       % The TDR object starts here
	    fprintf(fid,'%s\n%s\f\n','%% TDR 2.0 Binary','%%');
    end
    
    for (i = 1:n_groups)
        xx = get(hand(i),'XData');      yy = get(hand(i),'YData');
        np = length(xx);                % Number of points in this segment
        PointRad = get(hand(i),'MarkerSize') / 72 * 2.54 / 7;   % Symbol size. The 7 is an ad-hoc corr factor
        if (strcmp(opt,'Earthquakes'))  % Test those first because they have a z
            zz = -double(getappdata(hand(i),'SeismicityDepth')) * 100;    % zz is now in meters positive up
            if (size(zz,1) > 1)         % We need them as a row vector for fwrite
                zz = zz';
            end
        else                            % No Mirone code has yet a ZData property, but who knows
            zz = get(hand(i),'ZData');
        end
        
        if (isempty(zz))                % We need to have a z. Obviously it has to be 0
            zz = repmat(limits(6),1,np);
        else                            % Since we have a z, we need to update limits
            limits(5) = min(zz);        limits(6) = max(zz);
        end
        
	    fwrite(fid,10520,'integer*4');      % 10520 is the SD_POINT3D code
	    fwrite(fid,np*3*8+64,'integer*4');  % data length
	    fwrite(fid,[0 0 1 1 1 1 (1:18)*0],'integer*1');     % ?? o ultimo 1 parece dizer pontos
	    fwrite(fid,[np n_col 0 0],'integer*4');         % Number of points
	    fwrite(fid,limits,'real*8');
	    fwrite(fid,[xx; yy; zz;],'real*8');
    
        write_flederFiles('geo',fid,'add',limits)   % Write a GEOREF block
        write_flederFiles('cmap',fid)               % Write a FM_CMAP block
	
		fwrite(fid,[10521 35],'integer*4');         % 10520 is the SD_POINT3D_ATB code, 35 -> Data Length
		fwrite(fid,[0 0 1 1 1 5 (1:18)*0],'integer*1');
		fwrite(fid,[1 0 symb],'integer*4');
		fwrite(fid,[PointRad LabelSize],'real*4');
		fwrite(fid,[0 1 1],'integer*4');
        cor = uint8(get(hand(i),'MarkerFaceColor')*255);
		fwrite(fid,cor,'uint8');                    % symbol's color
    end
    
%----------------------------------------------------------------------------------
function [x,y,z,count] = lines2multiseg(hands,z_level)
% Convert a collection of lines whose handles are HANDS into a single multiline (cell array) array.
% If HANDS is a scalar, it will check if that line is broken with NaNs. If yes, it will be
% converted into a multiseg cell array. So we can call this function safely even when we don't
% know exactly what is contained in HANDS. The test of if x,... is a single or multisegment line
% are carried out inside the write_line function
% Z_LEVEL, if transmited will be used to set the line height, otherwise ZERO will be used.

if (nargin == 1),   z_level = 0;    end

x = get(hands,'XData');         y = get(hands,'YData');
if (iscell(x))
    n_lines = length(x);        % Number of lines
    z = x;                      % Create a cell array of the same size as x
    for (i = 1:length(x))
        z{i} = ones(1,length(x{i})) * z_level;
    end
else
    n_lines = 1;                % When it's not a cell array n_lines must be one
    z = zeros(1,length(x));
end

id_with_nan = repmat(logical(0),1,n_lines);
count = 0;                          % Counter of the total number of points
for (i = 1:n_lines)
    if (iscell(x))                  % Multi lines case 
        if (any(isnan(x{i})))       % See if we have NaNs in this line. If yes we must treat it as a multiseg
            id_with_nan(i) = 1;     % Yes we have. Mark this handle line to be processed later
            continue                % Jump this line for the time beeing
        end
        count = count + length(x{i});
    else                            % Single line, but with possible NaNs
        if (any(isnan(x)))          % See if we have NaNs in this line. If yes we must treat it as a multiseg
            id_with_nan = 1;        % Yes we have. Mark this handle line to be processed later
            continue                % Jump this line for the time beeing
        end
        count = count + length(x);
    end
end

if (any(id_with_nan))
    if (iscell(x))
        % Remove the cell fields corresponding to lines with NaNs (and add them later)
        x(id_with_nan) = [];        y(id_with_nan) = [];        z(id_with_nan) = [];
    end
    id_where = find(id_with_nan == 1);  % Allways == 1 for single line
    for (i = length(id_where))
        xt = get(hands(id_where(i)),'XData');
        yt = get(hands(id_where(i)),'YData');
        id_nan = find(xt ~= xt);    % Find the NaNs
        id = find(diff(id_nan) == 1) + 1;   % Account for contiguous NaNs
        if (~isempty(id))           % Found contiguous NaNs
            xt(id_nan(id)) = [];
            yt(id_nan(id)) = [];    % Remove them
            id_nan = find(xt ~= xt);% Find the new position of the now non-contiguous NaNs
        end
        id_nan = [0 id_nan];        % Used to make it start at one        
        
        n_segments = length(id_nan)-1;
        xx = cell(n_segments,1);    yy = xx;    zz = xx;
        for (k = 1:n_segments)
           xx{k} = x(id_nan(k)+1:id_nan(k+1)-1);
           yy{k} = y(id_nan(k)+1:id_nan(k+1)-1);
           zz{k} = repmat(z_level,1,length(xx{k}));
           count = count + length(xx{k});
        end
        if (iscell(x))              % Apend to eventual existing lines that had no NaNs
            x = [x; xx];        y = [y; yy];    z = [z; zz];
        else                        % Here we have a single line with NaNs
            x = xx;             y = yy;         z = zz;
        end
    end
end

%----------------------------------------------------------------------------------
function img = burnLines(hFig,hAxes,img,inplace)
% Burn the coastlines directly into de IMG image. We do it because also the Fleder
% doesn't work as advertized. Lines are awfully draped on surfaces

    head = getappdata(hFig,'GMThead');
    head(5:6) = [size(img,2) size(img,1)];
    % See if we have COASTLINES, POLITICALBOUND or RIVERS.
    h{1} = findobj(hAxes,'Type','line','Tag','CoastLineNetCDF');
    h{2} = findobj(hAxes,'Type','line','Tag','PoliticalBoundaries');
    h{3} = findobj(hAxes,'Type','line','Tag','Rivers');
    for (i=1:3)
        if (~isempty(h{i}))        
            xy = coast2pix(h{i}, head);
            line_thick = get(h{i},'LineWidth');       % Line thickness
            line_color = get(h{i},'color') * 255;     % Line color
            lt = 8;     % LINE_TYPE -> 8 connectivity (default)
            if (line_thick <= 1)
                lt = 16;        % antialiased line
            end
            if (inplace)
                cvlib_mex('poly',img,xy,[],line_thick,lt)
            else
                img = cvlib_mex('poly',img,xy,line_color,line_thick,lt);
            end
        end
    end

%----------------------------------------------------------------------------------
function [xy] = coast2pix(hand, lims)
% HAND is a handle of line wich can be broken with NaNs. If yes, it will be converted
% into a multiseg cell array. In the output we have an array of pixel coords (zero based)

x = get(hand,'XData');      y = get(hand,'YData');
x = x(:);                   y = y(:);   % Make sure they are column vectors

if (any(isnan(x)))          % We have NaNs in this line. Treat it as a multiseg
    xt = get(hand,'XData');
    yt = get(hand,'YData');
    id_nan = find(xt ~= xt);        % Find the NaNs
    id = find(diff(id_nan) == 1) + 1;   % Account for contiguous NaNs
    if (~isempty(id))               % Found contiguous NaNs
        xt(id_nan(id)) = [];
        yt(id_nan(id)) = [];        % Remove them
        id_nan = find(xt ~= xt);    % Find the new position of the now non-contiguous NaNs
    end
    id_nan = [0 id_nan];            % Used to make it start at one
    
    n_segments = length(id_nan)-1;
    xy = cell(n_segments,1);
    for (k = 1:n_segments)
        xx = round( localAxes2pix(lims(5),lims(1:2),x(id_nan(k)+1:id_nan(k+1)-1)) -1);  % -1 because we need
        yy = round( localAxes2pix(lims(6),lims(3:4),y(id_nan(k)+1:id_nan(k+1)-1)) -1);  % zero based indexes
        xy{k} = [xx yy];
    end
else
    xy = [xx yy];
end

% -------------------------------------------------------------------------------------
function pixelx = localAxes2pix(dim, x, axesx)
%   Convert axes coordinates to pixel coordinates.
%   PIXELX = AXES2PIX(DIM, X, AXESX) converts axes coordinates
%   (as returned by get(gca, 'CurrentPoint'), for example) into
%   pixel coordinates.  X should be the vector returned by
%   X = get(image_handle, 'XData') (or 'YData').  DIM is the
%   number of image columns for the x coordinate, or the number
%   of image rows for the y coordinate.

	xfirst = x(1);      xlast = x(max(size(x)));	
	if (dim == 1)
        pixelx = axesx - xfirst + 1;        return;
	end
	xslope = (dim - 1) / (xlast - xfirst);
	if ((xslope == 1) & (xfirst == 1))
        pixelx = axesx;
	else
        pixelx = xslope * (axesx - xfirst) + 1;
	end
