function writekml(handles, Z, fname)
% Write a GoogleEarth kml file that will allow displaing the current image in GE
%
% HANDLES -> Can be a Mirone handles or a struct with data as prepared by
%            choosebox.m prepare_poles_struct('kml') function.
% Z       -> Optional (otherwise fished from the HANDLES struct)
% FNAME   -> Optional. Apparently only used when HANDLES is a choosebox struct

%	Copyright (c) 2004-2016 by J. Luis
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

% $Id: writekml.m 7857 2016-04-11 17:29:20Z j $

	n_argin = nargin;
	noFig = true;
	to_delete = [];			% To hold a Fig handle that, if we have to reproject, will be deleted at the end

	% Check for the case where the input is NOT a handles structure but one with the data (choosebox.m uses this)
	if ((n_argin == 1) && isfield(handles, 'nofig'))
		handles.no_file = false;	handles.have_nans = false;
		handles.axes1 = [];			handles.hImg = [];
		handles.image_type = 20;
		handles.is_projected = false;		% Actually this should be made possible one day.
		if (isfield(handles, 'fname'))		% Write the kml file but NOT call GE
			fname = handles.fname;
			n_argin = 3;					% Dirty trick to force the flow into "write file and stop"
		elseif (isfield(handles, 'tmpdir'))
			handles.path_tmp = handles.tmpdir;
		else
			mir_dirs = getappdata(0,'MIRONE_DIRS');
			if (~isempty(mir_dirs))
				handles.path_tmp = [mir_dirs.home_dir filesep 'tmp' filesep];
			else
				error('writekml:directinput', 'Don''t know the adress of Mirone tmp directory')
			end
		end
		noFig = true;
	end

	if (handles.no_file),	return,		end		% This is convoluted. It allows calls by choosebox, but too convoluted.

	if (n_argin <= 2)
		if (n_argin == 1)					% Only HANDLES was transmited (called via clicked callback)
			Z = [];							% Z will be of use to know if need transparency
			if (handles.have_nans),		[X,Y,Z] = load_grd(handles,'silent');	end
		end
		if (handles.is_projected)				% SHIT, we need to convert to geogs under the hood
			handles = guidata(handles.figure1);
			prjSrc  = aux_funs('get_proj_string', handles);			% Get the proj4 string of this projection
			prjDest = '+proj=latlong +datum=WGS84';
			illumComm = getappdata(handles.figure1,'illumComm');	% Need to save this before it's erased below
			Illumin_type = handles.Illumin_type;					%		""
			hMirFigNew = gdal_project(handles, 'hiden', prjSrc, prjDest);
			handles = guidata(hMirFigNew);		% From now on this is the handles we will use
			[X,Y,Z] = load_grd(handles,'silent');
			handles.have_nans = grdutils(Z,'-N');
			if (Illumin_type >= 1 && Illumin_type <= 4)
				if (Illumin_type == 1)
					R = grdgradient_m(Z, handles.head,'-M', illumComm,'-Nt');	% This is the Geog case
				else
					R = grdgradient_m(Z, handles.head, illumComm, '-a1');
				end
				zz = ind2rgb8(get(handles.hImg,'CData'),get(handles.figure1,'ColorMap'));	% zz is now RGB
				zz = shading_mat(zz,R,'no_scale');		set(handles.hImg,'CData',zz)		% and now it is illuminated			
			end
			to_delete = hMirFigNew;			% And now store this guy to be leted at the end
		end
		fname_kml = [handles.path_tmp 'MironeTec.kml'];
		fname_img = [handles.path_tmp 'MironeTec'];
	else
		[PATH,FNAME] = fileparts(fname);
		fname_kml = [PATH filesep FNAME '.kml'];% Make sure it has the .kml extension
		fname_img = [PATH filesep FNAME];       % Extension will be added later
	end

	% Get the image at its own resolution
	img = get(handles.hImg,'CData');
	flipa = false;
	if (strcmp(get(handles.axes1,'YDir'),'normal'))
		img = flipdim(img,1);				% Ghrrrrrrr
		flipa = true;
	end

	% Control transparency
	if (handles.have_nans)					% We need transparency here. Note that NaNs imply that image derives from grid
		fname_img = [fname_img '.png'];		% And we'll use png
		cmap = get(handles.figure1,'Colormap');
		if (ndims(img) == 2)				% Indexed image
			imwrite(img,cmap,fname_img,'Transparency',0)
		else								% Truecolor
			m = size(Z,1);		n = size(Z,2);
			ind = isnan(Z);
			alfa = alloc_mex(m,n,'uint8');
			alfa(~ind) = 255;       clear ind;
			if (flipa),			alfa = flipdim(alfa,1);     end
			imwrite(img,fname_img,'Alpha',alfa);
		end
	elseif (handles.image_type ~= 20)       % Eventual original image transparency was already lost
		if (ndims(img) == 2)                % Indexed image
			fname_img = [fname_img '.png']; % And we still use png because img is indexed
			imwrite(img,get(handles.figure1,'Colormap'),fname_img);
		else
			alfa = get(handles.hImg,'AlphaData');
			if (size(alfa,1) == size(img,1))
				fname_img = [fname_img '.png'];
				if (flipa),		alfa = flipdim(alfa,1);     end
				imwrite(img,fname_img,'Alpha',alfa);
			else
				fname_img = [fname_img '.jpg']; % Now use jpeg because it allows compression
				imwrite(img,fname_img,'Quality',100);
			end
		end
	end
  
	% OPEN kml file
	fid = fopen(fname_kml,'w');
	ploted_PBs = false;         % To know when we need to plot plate boundaries

	fprintf(fid,'%s\n','<?xml version="1.0" encoding="UTF-8"?>');
	fprintf(fid,'%s\n','<kml xmlns="http://earth.google.com/kml/2.1">');
	fprintf(fid,'%s\n','<Document>');
	fprintf(fid,'\t%s%s%s\n','<name>','Mirone Tec','</name>');
	if (handles.image_type ~= 20)
		fprintf(fid,'\t%s\n','<GroundOverlay>');
		fprintf(fid,'\t\t%s\n','<name>Mirone Base Image</name>');
		fprintf(fid,'\t\t%s\n','<Icon>');
		[pato,nome,ext]=fileparts(fname_img);
		fprintf(fid,'\t\t\t%s%s%s\n','<href>',[nome ext],'</href>');
		fprintf(fid,'\t\t%s\n','</Icon>');
		fprintf(fid,'\t\t%s\n','<altitudeMode>clampToGround</altitudeMode>');
		fprintf(fid,'\t\t%s\n','<LatLonBox>');
		fprintf(fid,'\t\t\t%s%s%s\n','<north>',sprintf('%.6f',handles.head(4)),'</north>');
		fprintf(fid,'\t\t\t%s%s%s\n','<south>',sprintf('%.6f',handles.head(3)),'</south>');
		fprintf(fid,'\t\t\t%s%s%s\n','<east>',sprintf('%.6f',handles.head(2)),'</east>');
		fprintf(fid,'\t\t\t%s%s%s\n','<west>',sprintf('%.6f',handles.head(1)),'</west>');
		fprintf(fid,'\t\t%s\n','</LatLonBox>');
		fprintf(fid,'\t%s\n','</GroundOverlay>');
	end

	% Start fishing other candidate elements
	ALLpatchHand = findobj(handles.axes1,'Type','patch');

	% PATCHES
	if (~isempty(ALLpatchHand))
		cores = cell(1,3);
		nPatch = numel(ALLpatchHand);
		for (i = 1:nPatch)                     % Do one by one.
			x = get(ALLpatchHand(i),'XData');
			if (isempty(x)),    continue;       end
			y = get(ALLpatchHand(i),'YData');       z = getappdata(ALLpatchHand(i),'ZData');

			line_thick = get(ALLpatchHand(i),'LineWidth');   % Line thickness

			corA = get(ALLpatchHand(i),'FaceAlpha');
			if (ischar(corA)),  corA = 0;
			else                corA = round(corA*255);
			end

			corF = get(ALLpatchHand(i),'FaceColor');
			if (ischar(corF)),  corF = [255 255 255];   corA = 0;   % To account for 'none' or 'flat'
			else                corF = round(corF(end:-1:1)*255);
			end

			corE = get(ALLpatchHand(i),'EdgeColor');
			if (ischar(corE)),  corE = [255 255 255 255];   % To account for 'none' or 'flat'
			else                corE = [255 round(corE(end:-1:1)*255)];   % first is transparency (none here)
			end
			cores{1} = corA;    cores{2} = corF;    cores{3} = corE;

			writePolygon(fid, x, y, z, cores, 0, 'Patch', line_thick)
		end
	end

	% FISH THE LINES FAMILY (MANY P. DE CRIANCINHAS)
	ALLlineHand = findobj(handles.axes1,'Type','line');
	if (~isempty(ALLlineHand))

		% EARTHQUAKES
		h = findobj(ALLlineHand,'Tag','Earthquakes');       % Search first for earthquakes because ...
		if (~isempty(h))
			if (~ploted_PBs)        % We (well, I) want plate boundaries here
				fprintf(fid,'\t%s\n','<NetworkLink>');
				fprintf(fid,'\t\t%s%s%s\n','<Url><href>',[handles.path_data 'plateBoundaries_PB.kmz'],'</href></Url>');
				fprintf(fid,'\t%s\n','</NetworkLink>');
				ploted_PBs = true;
			end
			%load (['data' filesep 'mirone_icons.mat'],'circle_logic');
			load ([handles.path_data filesep 'mirone_icons.mat'],'circle_logic');
			img = uint8(circle_logic);
			n_groups = length(h);                           % N of different point ensembles
			for (i = 1:n_groups)
				xx = get(h(i),'XData');      yy = get(h(i),'YData');
				PointRad = get(h(i),'MarkerSize') / 72 * 2.54 * 1.5;	% Symbol size. The 1.5 factor is arbitrary
				zz = double(getappdata(h(i),'SeismicityDepth')) / 10;	% zz is now in km
				cor = get(h(i),'MarkerFaceColor');
				cmap = [1 1 1; cor];
				fname_img = [handles.path_tmp 'symb_' sprintf('%d',i) '.png'];
				imwrite(img,cmap,fname_img,'Transparency',0);
				writeStyle(fid,1,fname_img,'Mirone Tec')    % Write style and symbol
				mag = double(getappdata(h(i),'SeismicityMag')) / 10;
				name = sprintf('Mags %.1f - %.1f',min(mag),max(mag));
				writeQuakes(fid,1,xx,yy,zz,mag,PointRad,name)
			end
			ALLlineHand = setxor(ALLlineHand, h);           % h is processed, remove it from the handles list
		end

		% SCATTERS
		h = findobj(ALLlineHand,'Tag','scatter_symbs');
		if (~isempty(h))
			dpis = get(0,'ScreenPixelsPerInch') ;		% screen DPI
			pos = get(handles.axes1,'Position');		ylim = get(handles.axes1,'Ylim');
			escala = diff(ylim)/(pos(4)*2.54/dpis);		% Image units / cm

			fprintf(fid,'\t%s\n','<Folder>');
			for (i = 1:numel(h))
				x = get(h(i),'XData');       y = get(h(i),'YData');     z = getappdata(h(i),'ZData');
				ss = get(h(i),'MarkerSize');
				symb_size = ss / 72 * 2.54;				% Symbol size in cm
				dy = symb_size * escala;
				%rad = geo2dist([pt(1,1) center(1)],[pt(1,2) center(2)],'deg');
				rad = 0.002;
				[latc,lonc] = circ_geo(y,x,rad,[],18);
				z = repmat(z*200,1,numel(latc));
				c = get(h(i),'MarkerFaceColor');
				cores{1} = 255;
				cores{2} = round(c(end:-1:1)*255);
				%cores{2} = [0 0 255];
				cores{3} = [255 round(c(end:-1:1)*255)];

				writePolygon(fid, lonc, latc, z, cores, 1, 'SSs')
			end
			fprintf(fid,'\t%s\n','</Folder>');
			ALLlineHand = setxor(ALLlineHand, h);			% h is processed, remove it from handles list
		end

		% Search for points as Markers but first look for 'DATABASE' symbols
		hDB = findobj(ALLlineHand,'Tag','City_major');
		hDB = [hDB; findobj(ALLlineHand,'Tag','City_other')];
		if (~isempty(hDB))
			ALLlineHand = setxor(ALLlineHand, hDB);			% I'm removing those for now, but I might change my mind
		end
		hDB = findobj(ALLlineHand,'Tag','DSDP');
		hDB = [hDB; findobj(ALLlineHand,'Tag','ODP')];
		if (~isempty(hDB))									% ODP & DSDP & IODP & proposals
			fprintf(fid,'\t%s\n','<NetworkLink>');
			fprintf(fid,'\t\t%s%s%s\n','<Url><href>',[handles.path_data 'iodp_doc.kmz'],'</href></Url>');        
			fprintf(fid,'\t%s\n','</NetworkLink>');
			ALLlineHand = setxor(ALLlineHand, hDB);
		end
		hDB = findobj(ALLlineHand,'Tag','volcano');
		if (~isempty(hDB))									% GOTO SMITHSONIAN FILE
			fprintf(fid,'\t%s\n','<NetworkLink>');
			fprintf(fid,'\t\t%s%s%s\n','<Url><href>',[handles.path_data 'gvp_world.kmz'],'</href></Url>');
			fprintf(fid,'\t%s\n','</NetworkLink>');
			ALLlineHand = setxor(ALLlineHand, hDB);
		end
		hDB = findobj(ALLlineHand,'Tag','hydro');
		if (~isempty(hDB))									% GOTO INTERRIDGE FILE
			fprintf(fid,'\t%s\n','<NetworkLink>');
			fprintf(fid,'\t\t%s%s%s\n','<Url><href>',[handles.path_data 'vents_InterRidge.kmz'],'</href></Url>');
			fprintf(fid,'\t%s\n','</NetworkLink>');
			ALLlineHand = setxor(ALLlineHand, hDB);
		end
		hDB = findobj(ALLlineHand,'Tag','meteor');
		if (~isempty(hDB))									% GOTO http://impacts.rajmon.cz/IDdata.html FILE
			fprintf(fid,'\t%s\n','<NetworkLink>');
			fprintf(fid,'\t\t%s%s%s\n','<Url><href>',[handles.path_data 'Impact_database_2010.kmz'],'</href></Url>');
			fprintf(fid,'\t%s\n','</NetworkLink>');
			ALLlineHand = setxor(ALLlineHand, hDB);
		end
		hDB = findobj(ALLlineHand,'Tag','hotspot');			% FOGSPOTS
		if (~isempty(hDB))
			% Read the whole list. No point in ploting only the visible ones in Mirone window
			fp = fopen([handles.path_data 'hotspots.dat'],'r');
			fgetl(fp);						% Jump the header line
			todos = fread(fp,'*char');     fclose(fp);
			[xx yy fogName fogAge] = strread(todos,'%f %f %s %s');
			clear todos;
			fogNames = cell(numel(fogName),1);
			for (m=1:numel(fogName))
				fogNames{m} = [fogName{m} ', Age = ' fogAge{m} ' Ma'];
			end
			fogNames{1} = [fogNames{1} ' (Yes, you are allowed to lough ... and Very Loud)'];

			writeStyle(fid,1,[handles.path_data 'ghost.gif'],'GostSpots')    % Write style and symbol
			fprintf(fid,'\t%s\n','<Folder>');
			writeSymbols(fid,1,xx,yy,1,'The Ghosts',fogNames)
			fprintf(fid,'\t%s\n','</Folder>');
			ALLlineHand = setxor(ALLlineHand, hDB);
		end

		% Search for plain points as Markers (that is, line with no line - just symbols on vertices)
		h = [findobj(ALLlineHand,'Tag','Symbol') findobj(ALLlineHand,'Tag','Pointpolyline')];
		if (~isempty(h))
			symbol.Marker = get(h,'Marker');
			zz = get(h,'MarkerSize');
			if (~iscell(zz)),	symbol.Size = num2cell(zz,1);
			else				symbol.Size = zz;
			end
			zz = get(h,'MarkerFaceColor');
			if (~iscell(zz)),	symbol.FillColor = {zz};
			else				symbol.FillColor = zz;
			end
			symbol.Marker = char(symbol.Marker);
			symbol.Marker = symbol.Marker(:,1);
			load ([handles.path_data filesep 'mirone_icons.mat'],'star_logic','triangle_logic','losangle_logic',...
				'pentagon_logic','hexagon_logic','circle_logic');
			fprintf(fid,'\t%s\n','<Folder>');
			for (k = 1:numel(symbol.Marker))
				cmap = [0 0 0; symbol.FillColor{k}];
				fname_img = [handles.path_tmp sprintf('ico_%d.png',k)];
				switch symbol.Marker(k)
					case '*',   [img,cmap] = imread([handles.path_data 'AsteriskIcon.gif']);
					case 'o',   img = uint8(circle_logic);
					case 'p',   img = uint8(star_logic);
					case '^',   img = uint8(triangle_logic);
					case '<',   img = uint8(triangle_logic)';
					case '>',   img = fliplr(uint8(triangle_logic)');
					case 'v',   img = flipud(uint8(triangle_logic));
					case 'h',   img = uint8(hexagon_logic);
					case 's',   img = uint8(true(35));
				end
				imwrite(img,cmap,fname_img,'Transparency',0);
				writeStyle(fid,1,fname_img,sprintf('symb_%d',k))		% Write style and symbol
				xx = get(h(k),'XData');      yy = get(h(k),'YData');
				PointRad = get(h(k),'MarkerSize') / 72 * 2.54 * 1.5;	% Symbol size. The 1.5 factor is arbitrary
				writeSymbols(fid,1,xx,yy,PointRad,'o que')
			end
			fprintf(fid,'\t%s\n','</Folder>');
			ALLlineHand = setxor(ALLlineHand, h);			% h is processed, remove it from handles list
		end

		h = findobj(ALLlineHand,'Tag','PB_All');
		if (~isempty(h) && ~ploted_PBs)						% The guy wants the PB and they are not yet ploted
			fprintf(fid,'\t%s\n','<NetworkLink>');
			fprintf(fid,'\t\t%s%s%s\n','<Url><href>',[handles.path_data 'plateBoundaries_PB.kmz'],'</href></Url>');
			fprintf(fid,'\t%s\n','</NetworkLink>');
			ALLlineHand = setxor(ALLlineHand, h);			% h is processed, remove it from handles list
		end

		% If we still have 'line' elements, ... proceed
		for (i = 1:numel(ALLlineHand))
			x = get(ALLlineHand(i),'XData');        y = get(ALLlineHand(i),'YData');
			z = getappdata(ALLlineHand(i),'ZData');
			line_thick = get(ALLlineHand(i),'LineWidth');   % Line thickness
			line_color = get(ALLlineHand(i),'color');       % Line color
			line_color = [255 round(line_color(end:-1:1)*255)];   % first is transparency (none here)
			haveNaNs = isnan(x);
			if (any(haveNaNs))			% Contours do, so here we go again with this NaNs story (...getting tired)
				if (~isnan(x(1))),		x(1) = nan;		end
				if (~isnan(x(end))),	x(end) = nan;	end
				ind = find(isnan(x));
				for (k = 1:numel(ind)-1)
					xx = x( ind(k)+1:ind(k+1)-1 );
					yy = y( ind(k)+1:ind(k+1)-1 );
					if (~isempty(z)),	zz = z( ind(k)+1:ind(k+1)-1 );
					else				zz = z;
					end
					writePolyLine(fid,xx,yy,zz,line_color,line_thick)
				end
			else
				writePolyLine(fid,x,y,z,line_color,line_thick)
			end
		end
	end

	% FISH THE TEXT FAMILY
	ALLtextHand = findobj(handles.axes1,'Type','text');
	if (~isempty(ALLtextHand))
		fprintf(fid,'\t%s\n','<Folder>');
		writeText(fid,1,ALLtextHand,'Texts')
		fprintf(fid,'\t%s\n','</Folder>');
	end

	% OK, Now check if we got any direct data sent as function arguments
	if (noFig)
		try
			% Do we have lines?
			if (isfield(handles, 'line'))
				x = handles.line.x;					y = handles.line.y;
				if (isfield(handles.line, 'z')),	z = handles.line.z;
				else		z = zeros(size(x));
				end
				line_width = 1;
				if (isfield(handles.line, 'width')),	line_width = handles.line.width;	end
				line_color = [255 0 0 0];
				if (isfield(handles.line, 'color')),	line_color = handles.line.color;	end
				writePolyLine(fid,x,y,z,line_color,line_width)
			end

			% Do we have points?
			if (isfield(handles, 'pt'))
				x = handles.pt.x;				y = handles.pt.y;
				symbSize = 0.5;		% Default symbol size in cm (probably to big)
				if (isfield(handles.pt, 'size')),	symbSize = handles.pt.size;		end
				names = '';
				if (isfield(handles.pt, 'str')),	names = handles.pt.str;		end
				writeSymbols(fid,1,x,y,symbSize,'Tralala',names)
			end
		catch
			fclose(fid);
			error('writekml:directinput', lasterr)
		end
	end

	fprintf(fid,'%s\n','</Document>');
	fprintf(fid,'%s','</kml>');
	fclose(fid);

	if (~isempty(to_delete)),	delete(to_delete);	end		% The hidden Mirone Fig holding the reprojected window

	if (n_argin == 1)
		try
			if (ispc),		dos([fname_kml ' &']);
			else			unix(fname_kml);
			end
		catch
			errordlg(lasterr,'Error')
		end
	end

% ------------------------------------------------------------------------------------
function writeStyle(fid,nTab,fname_img,name)
% Write a style block

	sTab = repmat('\t',1,nTab);             sTab_p1 = repmat('\t',1,(nTab+1));
	sTab_p2 = repmat('\t',1,(nTab+2));      sTab_p3 = repmat('\t',1,(nTab+3));

	fprintf(fid,[sTab '%s%s%s\n'],'<name>',name,'</name>');
	fprintf(fid,[sTab '%s\n'],'<Style id="active-cor">');
	fprintf(fid,[sTab_p1 '%s\n'],'<IconStyle>');
	fprintf(fid,[sTab_p2 '%s\n'],'<Icon>');
	fprintf(fid,[sTab_p3 '%s%s%s\n'],'<href>',fname_img,'</href>');
	fprintf(fid,[sTab_p2 '%s\n'],'</Icon>');
	fprintf(fid,[sTab_p1 '%s\n'],'</IconStyle>');
	fprintf(fid,[sTab '%s\n'],'</Style>');

	fprintf(fid,[sTab '%s\n'],'<Style id="inactive-cor">');
	fprintf(fid,[sTab_p1 '%s\n'],'<IconStyle>');
	fprintf(fid,[sTab_p2 '%s\n'],'<Icon>');
	fprintf(fid,[sTab_p3 '%s%s%s\n'],'<href>',fname_img,'</href>');
	fprintf(fid,[sTab_p2 '%s\n'],'</Icon>');
	fprintf(fid,[sTab_p1 '%s\n'],'</IconStyle>');
	fprintf(fid,[sTab_p1 '%s\n'],'<LabelStyle><scale>0</scale></LabelStyle>');
	fprintf(fid,[sTab '%s\n'],'</Style>');

	fprintf(fid,[sTab '%s\n'],'<StyleMap id="cor">');
	fprintf(fid,[sTab_p1 '%s\n'],'<Pair>');
	fprintf(fid,[sTab_p2 '%s\n'],'<key>normal</key>');
	fprintf(fid,[sTab_p2 '%s\n'],'<styleUrl>#inactive-cor</styleUrl>');
	fprintf(fid,[sTab_p1 '%s\n'],'</Pair>');
	fprintf(fid,[sTab_p1 '%s\n'],'<Pair>');
	fprintf(fid,[sTab_p2 '%s\n'],'<key>highlight</key>');
	fprintf(fid,[sTab_p2 '%s\n'],'<styleUrl>#active-cor</styleUrl>');
	fprintf(fid,[sTab_p1 '%s\n'],'</Pair>');
	fprintf(fid,[sTab '%s\n'],'</StyleMap>');

% ------------------------------------------------------------------------------------
function writePolygon(fid, x, y, z, cores, extrude, name, line_thick)
% Write a polygon (patch)

	if (extrude),       outline = 0;    % Do not draw oulines with extrusions
	else                outline = 1;
	end

	x = x(:)';    y = y(:)';        % Make sure they are row vectors
	if (nargin == 3 || isempty(z))
		z = zeros(1,numel(x));
		clampMode = 'clampToGround';
	else 
		z = z(:)';
		clampMode = 'relativeToGround';
	end

	nTab = 1;
	sTab = repmat('\t',1,nTab);             sTab_p1 = repmat('\t',1,(nTab+1));
	sTab_p2 = repmat('\t',1,(nTab+2));      sTab_p3 = repmat('\t',1,(nTab+3));
	sTab_p4 = repmat('\t',1,(nTab+4));

	fprintf(fid,[sTab '%s\n'],'<Placemark>');
	fprintf(fid,[sTab_p1 '%s%s%s\n'],'<name>',name,'</name>');
	fprintf(fid,[sTab_p1 '%s\n'],'<Style>');

	if (nargin == 8)
		fprintf(fid,[sTab_p2 '%s\n'],'<LineStyle>');           % Print polygon outline
		fprintf(fid,[sTab_p2 '\t%s%.2x%.2x%.2x%.2x%s\n'],'<color>',cores{3}(1:4),'</color>');
		fprintf(fid,[sTab_p2 '\t%s%g%s\n'],'<width>',line_thick,'</width>');
		fprintf(fid,[sTab_p2 '%s\n'],'</LineStyle>');
	end

	fprintf(fid,[sTab_p2 '%s\n'],'<PolyStyle>');
	fprintf(fid,[sTab_p3 '%s%.2x%.2x%.2x%.2x%s\n'],'<color>',cores{1},cores{2}(1:3),'</color>');
	fprintf(fid,[sTab_p3 '%s%d%s\n'],'<outline>',outline,'</outline>');
	fprintf(fid,[sTab_p2 '%s\n'],'</PolyStyle>');

	fprintf(fid,[sTab_p1 '%s\n'],'</Style>');
	fprintf(fid,[sTab_p1 '%s\n'],'<Polygon>');
	fprintf(fid,[sTab_p2 '%s%d%s\n'],'<extrude>',extrude,'</extrude>');
	fprintf(fid,[sTab_p2 '%s%s%s\n'],'<altitudeMode>',clampMode,'</altitudeMode>');

	fprintf(fid,[sTab_p2 '%s\n'],'<outerBoundaryIs>');
	fprintf(fid,[sTab_p3 '%s\n'],'<LinearRing>');
	fprintf(fid,[sTab_p4 '%s\n'],'<coordinates>');
	fmt = repmat('%.6f,%.6f,%.0f ',1,4);
	fprintf(fid,[fmt '\n'], [x; y; z]);        % Write the data
	fprintf(fid,[sTab_p4 '%s\n'],'</coordinates>');
	fprintf(fid,[sTab_p3 '%s\n'],'</LinearRing>');
	fprintf(fid,[sTab_p2 '%s\n'],'</outerBoundaryIs>');
	fprintf(fid,[sTab_p1 '%s\n'],'</Polygon>');
	fprintf(fid,[sTab '%s\n'],'</Placemark>');    

% ------------------------------------------------------------------------------------
function writePolyLine(fid,x,y,z,line_color,line_thick,name)
% Write a polyline

	if (nargin == 6),   name = 'Polyline';  end
	x = x(:)';    y = y(:)';        % Make sure they are row vectors
	if (isempty(z)),    z = zeros(1,numel(x));
	else                z = z(:)';
	end

	nTab = 1;
	sTab = repmat('\t',1,nTab);				sTab_p1 = repmat('\t',1,(nTab+1));
	sTab_p2 = repmat('\t',1,(nTab+2));

	fprintf(fid,[sTab '%s\n'],'<Placemark>');
	fprintf(fid,[sTab_p1 '%s%s%s%\n'],'<name>',name,'</name>');
	fprintf(fid,[sTab_p1 '%s\n'],'<Style>');
	fprintf(fid,[sTab_p2 '%s\n'],'<LineStyle>');           % Print polygon outline
	fprintf(fid,[sTab_p2 '\t%s%.2x%.2x%.2x%.2x%s\n'],'<color>',line_color(1:4),'</color>');
	fprintf(fid,[sTab_p2 '\t%s%g%s\n'],'<width>',line_thick,'</width>');
	fprintf(fid,[sTab_p2 '%s\n'],'</LineStyle>');
	fprintf(fid,[sTab_p1 '%s\n'],'</Style>');

	fprintf(fid,[sTab_p1 '%s\n'],'<LineString>');
	% 	fprintf(fid,[sTab_p1 '%s\n'],'<tessellate>1</tessellate>');
	fprintf(fid,[sTab_p2 '%s\n'],'<coordinates>');
	fmt = repmat('%.6f,%.6f,%.0f ',1,4);
	fprintf(fid,[fmt '\n'], [x; y; z]);        % Write the data
	fprintf(fid,[sTab_p2 '%s\n'],'</coordinates>');
	fprintf(fid,[sTab_p1 '%s\n'],'</LineString>');
	fprintf(fid,[sTab '%s\n'],'</Placemark>');

% ------------------------------------------------------------------------------------
function writeQuakes(fid,nTab,x,y,z,mag,scale,name)
% Write earthquakes

	sTab = repmat('\t',1,nTab);             sTab_p1 = repmat('\t',1,(nTab+1));
	sTab_p2 = repmat('\t',1,(nTab+2));      sTab_p3 = repmat('\t',1,(nTab+3));
	sTab_p4 = repmat('\t',1,(nTab+4));

	fprintf(fid,[sTab '%s\n'],'<Folder>');
	fprintf(fid,[sTab_p1 '%s%s%s\n'],'<name>',name,'</name>');

	% OK, now loop ever points (events)
	for (k=1:numel(x))
		fprintf(fid,[sTab_p1 '%s\n'],'<Placemark>');
		fprintf(fid,[sTab_p2 '%s%.1f, Depth %.1f km%s\n'],'<name>M ',mag(k),z(k),'</name>');
		fprintf(fid,[sTab_p2 '%s\n'],'<Snippet maxLines="0"></Snippet>');
		fprintf(fid,[sTab_p2 '%s%s%s\n'],'<styleUrl>','#cor','</styleUrl>');
		fprintf(fid,[sTab_p2 '%s\n'],'<Style>');
		fprintf(fid,[sTab_p3 '%s\n'],'<IconStyle>');
		fprintf(fid,[sTab_p4 '%s%.3f%s\n'],'<scale>',scale,'</scale>');
		fprintf(fid,[sTab_p3 '%s\n'],'</IconStyle>');
		fprintf(fid,[sTab_p2 '%s\n'],'</Style>');

		fprintf(fid,[sTab_p2 '%s\n'],'<Point>');
		fprintf(fid,[sTab_p3 '%s%.4f,%.4f%s\n'],'<coordinates>',x(k),y(k),',0</coordinates>');
		fprintf(fid,[sTab_p2 '%s\n'],'</Point>');
		fprintf(fid,[sTab_p1 '%s\n'],'</Placemark>');        
	end

	fprintf(fid,[sTab '%s\n'],'</Folder>');

% ------------------------------------------------------------------------------------
function writeSymbols(fid,nTab,x,y,scale,nameGroup,names)
% Write symbols coordinates. The symbol icon was already saved with 'Style'

	sTab_p1 = repmat('\t',1,(nTab+1));		sTab_p2 = repmat('\t',1,(nTab+2));
	sTab_p3 = repmat('\t',1,(nTab+3));		sTab_p4 = repmat('\t',1,(nTab+4));
	Pname = false;      % To use when we have Placemark names
	if (nargin == 7 && iscell(names)),      Pname = true;   end
	fprintf(fid,[sTab_p1 '%s%s%s\n'],'<name>',nameGroup,'</name>');
	% OK, now loop over number of points ('grou symbols' have several points)
	for (k = 1:numel(x))
		fprintf(fid,[sTab_p1 '%s\n'],'<Placemark>');
		if (Pname)
			fprintf(fid,[sTab_p1 '%s%s%s\n'],'<name>',names{k},'</name>');
		end
		fprintf(fid,[sTab_p2 '%s%s%s\n'],'<styleUrl>','#cor','</styleUrl>');
		fprintf(fid,[sTab_p2 '%s\n'],'<Style>');
		fprintf(fid,[sTab_p3 '%s\n'],'<IconStyle>');
		fprintf(fid,[sTab_p4 '%s%.3f%s\n'],'<scale>',scale,'</scale>');
		fprintf(fid,[sTab_p3 '%s\n'],'</IconStyle>');
		fprintf(fid,[sTab_p2 '%s\n'],'</Style>');
		fprintf(fid,[sTab_p2 '%s\n'],'<Point>');
		fprintf(fid,[sTab_p3 '%s%.4f,%.4f%s\n'],'<coordinates>',x(k),y(k),',0</coordinates>');
		fprintf(fid,[sTab_p2 '%s\n'],'</Point>');
		fprintf(fid,[sTab_p1 '%s\n'],'</Placemark>');        
	end

% ------------------------------------------------------------------------------------
function writeText(fid,nTab,h,nameGroup)
% Write text strings. Here we better extract position and String info directly from the text handles H
% This function needs improvement as the color and scale are hardwired instead of variable controlled

	sTab_p1 = repmat('\t',1,(nTab+1));		sTab_p2 = repmat('\t',1,(nTab+2));
	sTab_p3 = repmat('\t',1,(nTab+3));

	% Create a style block that doesn't print any marker (but needs improvement)
	fprintf(fid,[sTab_p1 '%s\n'],'<Style id="Mir-1">');
	fprintf(fid,[sTab_p2 '<IconStyle><Icon></Icon></IconStyle>\n']);
	fprintf(fid,[sTab_p2 '%s\n'],'<LineStyle>');
	fprintf(fid,[sTab_p3 '<color>ff000000</color>\n']);
	fprintf(fid,[sTab_p3 '<width>1</width>\n']);
	fprintf(fid,[sTab_p2 '%s\n'],'</LineStyle>');
	fprintf(fid,[sTab_p2 '%s\n'],'<PolyStyle>');
	fprintf(fid,[sTab_p3 '<color>ff80c0ff</color>\n']);
	fprintf(fid,[sTab_p3 '<fill>1</fill>\n']);
	fprintf(fid,[sTab_p3 '<outline>1</outline>\n']);
	fprintf(fid,[sTab_p2 '%s\n'],'</PolyStyle>');
	fprintf(fid,[sTab_p2 '%s\n'],'<LabelStyle>');
	fprintf(fid,[sTab_p3 '%s\n'],'<scale>1</scale>');
	fprintf(fid,[sTab_p3 '<color>bfffffff</color>\n']);
	fprintf(fid,[sTab_p2 '%s\n'],'</LabelStyle>');
	fprintf(fid,[sTab_p1 '%s\n'],'</Style>');

	str = get(h,'String');		pos = get(h,'Position');
	if (~isa(str,'cell'))		% We need to access them as cell arrays
		str = {str};			pos = {pos};
	end

	fprintf(fid,[sTab_p1 '%s%s%s\n'],'<name>',nameGroup,'</name>');
	% OK, now loop ever number of points ('grou symbols' have several points)
	for (k = 1:numel(h))
		fprintf(fid,[sTab_p1 '%s\n'],'<Placemark>');
		fprintf(fid,[sTab_p1 '%s%s%s\n'],'<name>',str{k},'</name>');
		fprintf(fid,[sTab_p2 '%s%s%s\n'],'<styleUrl>','#Mir-1','</styleUrl>');
		fprintf(fid,[sTab_p2 '%s\n'],'<Point>');
		fprintf(fid,[sTab_p3 '%s%.4f,%.4f%s\n'],'<coordinates>',pos{k}(1),pos{k}(2),',0</coordinates>');
		fprintf(fid,[sTab_p2 '%s\n'],'</Point>');
		fprintf(fid,[sTab_p1 '%s\n'],'</Placemark>');        
	end
