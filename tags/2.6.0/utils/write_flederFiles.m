function filename = write_flederFiles(opt,varargin)
% 'flederize' Matlab graphical objects.
%
% Depending of the OPT value, this function builds either:
% OPT = 'writeAll3' a set of three files: .geo, .dtm, .shade as DMagic would do
% OPT = 'writePlanar' directly build a planar Sonar SD file to be used by Fledermaus
% OPT = 'writeSpherical' directly build a spherical Sonar SD file to be used by Fledermaus
% OPT = 'runPlanar' build a planar .sd file 			named "lixoSD.sd" at ...\tmp
% OPT = 'runSpherical' build a spherical .sd file					"
%
% It is also possible to use this as a gateway to some subfunctions. In that case 
% OPT = subfunction, and VARARGIN varies with the function called
% In main_SD must be: VARARGIN = (fid|fname, 'Planar|other', Z, img, limits)
% In write_all3 must be: VARARGIN = (name, flederPlanar, Z, img, limits)

% NOTE: For every line/point object FM_CMAP and a GEOREF blocks are writen.
% Though it doesn't hurt much, it is an idiot thing

%	Copyright (c) 2004-2014 by J. Luis
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

% $Id$

	if (opt(1) == 'w' || opt(1) == 'r')			% Here opt is a "write..." or "run..."
		handles = varargin{1};
		if (handles.validGrid)
			[X,Y,Z,head] = load_grd(handles);
			if (isempty(Z)),	return,		end
		else
			% A bit of cheating. Modify the OPT string to build a geoimage obj 
			if (opt(1) == 'r'),		opt = 'runGeoimg';
			else					opt = 'writeGeoimg';
			end
			head = handles.head;
		end

		% Fish the image
		[img, img2] = getImg(handles);		% If IMG2 is not empty we need to make a "scene"

		if (strcmp(opt,'writeAll3'))		% Write the 3 files and return
			str1 = {'*.*', 'All Files (*.*)'};
			[FileName,PathName] = put_or_get_file(handles,str1,'Name stem for the 3 (.geo, .dtm, .shade) IVS files','put');
			if isequal(FileName,0),		return,		end
			[PATH,FNAME] = fileparts(FileName);    fname = [PathName filesep FNAME];
			set(handles.figure1,'pointer','watch')
			write_all3(fname,handles.flederPlanar,Z,img,head(1:6))
			set(handles.figure1,'pointer','arrow')
			return
		end

		str1 = {'*.sd;*.SD', 'Fledermaus object (*.sd,*.SD)'};
		if (strncmp(opt,'write',5))
			[FileName,PathName] = put_or_get_file(handles,str1,'Name of Fledermaus object','put');
			if isequal(FileName,0),		return,		end
			fname = [PathName FileName];
		else
			fname = [handles.path_tmp 'lixoSD.sd'];
		end

		set(handles.figure1,'pointer','watch')
		[PATH,FNAME,EXT] = fileparts(fname);
		if isempty(EXT),    fname = [fname '.sd'];  end
		if (strcmp(opt(end-2:end),'img'))		% either writeGeoimg or runGeoimg
			fid = fopen(fname,'wb');
			if (~isempty(img))
				write_geoimg(fid, 'first', img)
				write_geo(fid,'add',head(1:6))
			else								% Initiate the TDR
				fprintf(fid,'%s\n%s\f\n','%% TDR 2.0 Binary','%%');
			end
		else
			if (strcmp(opt,'writePlanar') || handles.flederPlanar)
				tipo = 'Planar';
			else
				tipo = 'Spherical';
			end
			vimage = getappdata(handles.axes1,'VIMAGE');
			if (isempty(vimage))
				fid = fopen(fname,'wb');
				write_main(fid, tipo, Z, img, img2, head(1:6));
			else
				fname = [fname '.scene'];		% Will have .sd.scene extension but no big deal
				fid = fopen(fname,'wb');
				write_scene(fid, tipo, 'vimage', Z, img, img2, head(1:6), vimage)
			end
		end
		if (handles.flederBurn ~= 2)		% If not screen capture, see if there are lines & pts to flederize
			write_lines_or_points(fid, handles.axes1, head(1:6), handles.flederBurn);
		end
		write_eof(fid)						% Write EOF block and close the file
		set(handles.figure1,'pointer','arrow')
	end

	% Go trough here when using specific external calls (some are probably not possible to call from outside)
	switch opt
		case 'geo'
			write_geo(varargin{:});     % Write a .geo block
		case 'dtm'
			write_dtm(varargin{:});     % Write a .dtm object
		case 'shade'
			write_shade(varargin{:});   % Write a .shade object
		case 'main_SD'		% The call must be: write_main(fid, tipo, Z, img, img2, limits)
			fid = write_main(varargin{:});		% Write a .sd object made of .geo .dtm & .shade blocks
			write_eof(fid)						% Write EOF block and close the file
		case 'all3'			% The call must be: write_all3(name, flederPlanar, Z, img, limits)
			write_all3(varargin{:})     % Write three files: .geo, .dtm, .shade as DMagic would do
		case 'line_or_points'
			write_lines_or_points(varargin{:});    % Search for and write lines and/or points objects
		case 'line'
			write_line(varargin{:});    % Write line objects
		case 'points'
			if (ischar(varargin{1}))	varargin{1} = fopen(varargin{1},'wb');	end		% Was file name
			write_pts(varargin{:});     % Write point objects
			write_eof(varargin{1})		% Write EOF block and close the file
	end

	if (nargout)	filename = fname;	end		% Otherwise the compiled version would silently error

%----------------------------------------------------------------------------------
function [img, img2] = getImg(handles)
% Fish the image in a Mirone handles and flip it if necessary
% This function is used when only HANDLES was transmited

	if (handles.image_type == 20 && handles.flederBurn ~= 2)		% We don't have any bg image
		img = [];		img2 = [];		return
	end

	% The UD flip horror.
	flipa = false;
	if (strcmp(get(handles.axes1,'Ydir'),'reverse') && handles.flederBurn ~= 2),    	flipa = true;
	elseif (strcmp(get(handles.axes1,'Ydir'),'normal') && handles.flederBurn == 2)		flipa = true;
	end

	img2 = [];		% To hold original image in case of drapping when displayed image was resized
	if (handles.flederBurn == 2)						% "Burn them all"
		hL = findobj(handles.axes1,'Type','line');		hP = findobj(handles.axes1,'Type','patch');
		hT = findobj(handles.axes1,'Type','text');
		if (isempty(hL) && isempty(hP) && isempty(hT))	% No need to SC because image is empty
			img = get(handles.hImg,'CData');
			% Since we didn't do the SC we need to recheck the flip condition
			if ( flipa && strcmp(get(handles.axes1,'Ydir'),'normal') ),		flipa = false;		end
		else
			img = imcapture(handles.axes1,'img',0);		% Do a screen capture
		end
		if (flipa),		img = flipdim(img,1);	end
		if (ndims(img) == 2)
			img = ind2rgb8(img,get(handles.figure1,'Colormap'));
		end
	else
		img = get(handles.hImg,'CData');
		if (flipa),		img = flipdim(img,1);	end
		inplace = false;            % Eventual line burning is not done inplace to not change the Mirone image as well
		if (ndims(img) == 2)
			img = ind2rgb8(img,get(handles.figure1,'Colormap'));
			inplace = true;         % Save to do line burning inplace because img is already a copy
		end
		if (handles.flederBurn == 1)					% Burn coastlines into img (in case they exist, else do nothing)
			img = burnLines(handles.figure1, handles.axes1, img, inplace);
		end
	end

	% See if we are in a "drapped" case with parent image resized. If yes we'll have to build a textureDTM SD
	if ( handles.is_draped && ~isequal([size(img,1) size(img,2)],[size(handles.origFig,1) size(handles.origFig,2)]) )
		img2 = handles.origFig;
		if (flipa),		img2 = flipdim(img2,1);		end
		if (ndims(img2) == 2)
			img2 = ind2rgb8(img2,get(handles.figure1,'Colormap'));
		end
	end

%----------------------------------------------------------------------------------
function write_geo(fid,mode,limits)
% Write a GEOREF block
	if (strcmp(mode,'first'))       % The TDR object starts here
		fprintf(fid,'%s\n%s\n%s\f\n','%% TDR 2.0 Binary','Created by:    Mirone Tech!','%%');
	end
	fwrite(fid,[15000 100],'integer*4');     % Tag ID, Data Length
	fwrite(fid,[0 0 1 1 1 1 (1:30)*0],'uchar');
	fwrite(fid,limits,'real*8');
	fwrite(fid,(1:10)*0,'integer*4');

%----------------------------------------------------------------------------------
function write_dtm(fid,mode,Z,limits)
% Write a .dtm block
	if (strcmp(mode,'first'))       % The TDR object starts here
		fprintf(fid,'%s\n%s\n%s\f\n','%% TDR 2.0 Binary','Created by:    Mirone Tech!','%%');
	end
	[m,n] = size(Z);
	fwrite(fid,[1000 m*n*2+30],'integer*4');		% Tag ID, Data Length 
	fwrite(fid,[0 0 1 1 1 2 (1:18)*0],'uchar');
	fwrite(fid,[2 n m 3 16],'integer*2');			% nDim(?) nCols nRows ?? BitWidth
	fwrite(fid,[limits(5) limits(6)],'real*8');
	fwrite(fid,[0 65535],'uint16');
	Z = flipud( uint16(rot90(scaleto8(Z,16),1)) );
	fwrite(fid,Z,'uint16');

%----------------------------------------------------------------------------------
function write_shade(fid, mode, img, geoimg)
% Write a .shade or a geoimage block
    if (strcmp(mode,'first'))       % The TDR object starts here
		fprintf(fid,'%s\n%s\n%s\f\n','%% TDR 2.0 Binary','Created by:    Mirone Tech!','%%');
    end

	if (nargin == 3)		% Non geoimage
		soma = 34;		nZeros = 10;	bW = [7 32];	nB = 65535;		% Again this is from pure trial
	else
		soma = 28;		nZeros = 7;		bW = [12 8];	nB = 255;		% geoimage object
	end
	
	m = size(img,1);		n = size(img,2);
	fwrite(fid,[1000 m*n*4+soma],'integer*4');		% Tag ID, Data Length
	fwrite(fid,[0 0 1 1 1 2 (1:18)*0],'uchar');
	fwrite(fid,[2 n m bW],'integer*2');				% nDim(?) nCols nRows ?? BitWidth
	fwrite(fid,[(1:nZeros)*0 nB nB],'integer*2');
	img(:,:,4) = uint8(255);						% Tranparency layer
	img = permute(img,[4 3 2 1]);
	fwrite(fid,img,'uint8');

%----------------------------------------------------------------------------------
function write_all3(name_stem, flederPlanar, Z, img, limits)
% Write three files: .geo, .dtm, .shade

	exts = {'geo' 'dtm' 'shade'};
	for (k = 1:3)
		fid = fopen([name_stem '.' exts{k}],'wb');
		fprintf(fid,'%s\n%s\n%s\f\n','%% TDR 2.0 Binary','Created by:    Mirone Tech!','%%');

		if (k == 1)			% .geo
			tipo = '';
			if (flederPlanar),	tipo = 'Planar';	end
			write_sonardtm_block(fid, tipo)
			write_geo(fid,'add',limits);		% Write a .geo block
		elseif (k == 2)		% .dtm
			write_dtm(fid,'add',Z,limits);		% Write a .dtm object
		else				%.shade
			write_shade(fid,'add', img);		% Write a .shade object
		end
		write_eof(fid)							% Write EOF block and close the file
	end

%----------------------------------------------------------------------------------
function fid = write_main(fid, tipo, Z, img, img2, limits)
% Write a basic .sd file. That is with a DTM, a SHADE & a GEO blocks

	if (ischar(fid))		% FID is in fact the file name 
		fid = fopen(fid,'wb');
	end

	fprintf(fid,'%s\n%s\n%s\f\n','%% TDR 2.0 Binary','Created by:    Mirone Tech!','%%');

	if (nargin == 5 || isempty(img2))		% A SD object
		write_sonardtm_block(fid, tipo)
		write_dtm(fid,'add',Z,limits)
		write_shade(fid, 'add', img)
		write_geo(fid,'add',limits)
	elseif (1)								% A textureDTM SD object
		write_texturedtm_block(fid)
		write_shade(fid, 'add', img, 'geoimg')
		write_geo(fid,'add',limits)
		write_dtm(fid,'add',Z,limits)
		write_shade(fid, 'add', img2)
		write_geo(fid,'add',limits)
	else									% A SCENE object
		write_scene_block(fid)
		write_node_block(fid, 'root', 'Root Node', 'Unknown')
		write_geo(fid,'add',limits)
		write_alignparent_block(fid)
		write_geo(fid,'add',limits)
 		write_node_block(fid, 'dtm', 'Surface', 'test1.sd')
		write_geo(fid,'add',limits)
		write_sonardtm_block(fid, tipo)
		write_dtm(fid,'add',Z,limits)
		write_shade(fid, 'add', img2)
		write_geo(fid,'add',limits)
		write_sonardtm_atb_block(fid)
 		write_node_block(fid, 'geoimg', 'Drapped', 'test2.sd')
		limits(5:6) = [0 1];
		write_geo(fid,'add',limits)
		write_geoimg(fid, 'add', img)
		write_geo(fid,'add',limits)
		write_geoimg_atb_block(fid)
	end

%----------------------------------------------------------------------------------
function write_scene(fid, tipo, mode, Z, img, img2, limits, struct_vimage)
% Write a scene (very crude)

%	fprintf(fid,'%s\n%s\n%s\f\n','%% TDR 2.0 Binary','Created by:    Mirone Tech!','%%');
	fprintf(fid,'%s\n%s\n%s\f\n','%% TDR 2.0 Binary','Fledermaus wrote this scene file!','%%');

	write_scene_block(fid)
	write_node_block(fid, 'root', 'Root Node', 'Unknown')
	%write_geo(fid,'add',limits)
	write_geo(fid,'add',[limits(1:4) -10000 limits(end)])
	write_alignparent_block(fid)

	if (~isempty(Z))
		write_geo(fid,'add',[limits(1:4) -10000 limits(end)])	% <
%		write_node_block(fid, 'dtm', 'Surface', 'unknown')
  		write_node_block(fid, 'dtm', 'test1.sd', 'C:|programs|IVS|bin|test1.sd')
		write_geo(fid,'add',limits)		% <
		write_sonardtm_block(fid, tipo)
		write_dtm(fid,'add',Z,limits)
		write_shade(fid, 'add', img)
		write_geo(fid,'add',limits)
		write_sonardtm_atb_block(fid)
	end

	if (strcmp(mode,'vimage'))
		for (k = 1:numel(struct_vimage))
			hLine = struct_vimage(k).hLine;
			if (~ishandle(hLine)),		continue,		end			% Line was deleted
			
			try
				info_img = imfinfo(struct_vimage(k).vimage);
				I = gdalread(struct_vimage(k).vimage,'-U');
				if (strcmp(info_img(1).ColorType,'grayscale') || (strcmp(info_img(1).ColorType,'truecolor') && (ndims(I) ~= 3)) )
					I = ind2rgb8(img2,'Colormap',gray(256));
				elseif (isfield(info_img(1),'ColorTable'))			% Gif images call it 'ColorTable'
					I = ind2rgb8(img2,'Colormap',info_img(1).ColorTable);
				elseif (isfield(info_img(1),'Colormap') && ~isempty(info_img(1).Colormap))
					I = ind2rgb8(img2,'Colormap',info_img(1).Colormap);
				end
			catch
				errordlg(['WriteFleder: Error reading image: ' struct_vimage(k).vimage],'ERROR')
				continue
			end

			[pato, fname] = fileparts(struct_vimage(k).vimage);
 			write_node_block(fid, 'vimage', fname, 'unknown')
			X = get(hLine, 'XData');		Y = get(hLine, 'YData');
			Vlimits = [X(1) Y(1) struct_vimage(k).z_min X(end) Y(end) struct_vimage(k).z_max];
			write_geo(fid,'add',limits)	% <
			write_vimage_block(fid, Vlimits)
			Vlimits = [X(1) X(end) Y(1) Y(end) struct_vimage(k).z_min struct_vimage(k).z_max];
			write_shade(fid, 'add', I, 'geoimg')
			write_geo(fid,'add',Vlimits)
			write_geo(fid,'add',Vlimits)
			write_vimage_atb_block(fid)
		end
	else
		write_node_block(fid, 'geoimg', 'Drapped', 'unknown')
		limits(5:6) = [0 1];
		write_geo(fid,'add',limits)
		write_geoimg(fid, 'add', img2)
		write_geo(fid,'add',limits)
		write_geoimg_atb_block(fid)
	end

%----------------------------------------------------------------------------------
function write_geoimg(fid, mode, img)
% Write a basic image .sd file. That is one with image & a GEO blocks
	if (strcmp(mode,'first'))       % The TDR object starts here
		fprintf(fid,'%s\n%s\n%s\f\n','%% TDR 2.0 Binary','Created by:    Mirone Tech!','%%');
	end
	write_geoimg_block(fid)
	write_shade(fid, 'add', img, 'geoimg')

%----------------------------------------------------------------------------------
function write_sonardtm_block(fid, tipo)
% Write the SD_SONARDTM block
	if (strcmp(tipo,'Planar'))
		fwrite(fid,[10005 2],'integer*4');			% Tag ID, Data Length
	    fwrite(fid,[0 0 1 1 1 1 (1:18)*0 0 0],'uchar');
	else
		fwrite(fid,[10015 0],'integer*4');			% Tag ID, Data Length
		fwrite(fid,[0 0 1 1 1 1 (1:18)*0],'uchar');
	end

%----------------------------------------------------------------------------------
function write_scene_block(fid)
% Write the FM_SCENE block

	fwrite(fid,[20002 116],'integer*4');		% Tag ID, Data Length
	fwrite(fid,[0 0 1 1 1 3 0 0],'uchar');
	fwrite(fid,(1:9)*0,'integer*4');
% 	fwrite(fid,[205 204 204 61 0 64 28 70 0 0 0 0 10 215 35 60 0 36 116 72],'uchar');
% 	fwrite(fid,[0 64 156 69 0 64 156 69 0 0 128 63 (1:19)*0 63 215 179 93 191],'uchar');
	fwrite(fid,[205 204 204 61 0 64 28 70 1 0 1 0 10 215 35 60 0 36 116 72],'uchar');
	fwrite(fid,[72 127 248 66 0 192 90 69 0 0 128 63 (1:16)*0 202 232 227 62 21 61 101 191],'uchar');
	fwrite(fid,[0 0],'integer*4');

	%  	fwrite(fid,[215 179 93 63 0 0 0 63 (1:13)*0 64 156 197 0 0 128 63 0 0 72 66],'uchar');
 	fwrite(fid,[21 61 101 63 202 232 227 62 (1:13)*0 64 156 197 0 0 128 63 0 0 72 66],'uchar');
	fwrite(fid,[9900 0],'integer*4');				% Tag ID, Data Length		-- SD_HIERARCHY
	fwrite(fid,[0 0 1 1 1 1 (1:18)*0],'uchar');

%----------------------------------------------------------------------------------
function write_node_block(fid, tipo, str1, str2)
% Write the FM_NODE block

	dl = 180 + numel(str1) + numel(str2);
	fwrite(fid,[9910 dl],'integer*4');		% Tag ID, Data Length
	fwrite(fid,[0 0 1 1 1 1 (1:18)*0],'uchar');
	fwrite(fid,[-1 0 0 0 0],'integer*4');
	if (strcmp(tipo,'root'))
		%fwrite(fid,[0 0 128 63 0 0 128 63 0 0 128 63],'uchar');		% Non VIMAGE ???
		fwrite(fid,[24 4 41 66 203 90 225 65 0 0 128 63],'uchar');
		fwrite(fid,[0 0 0 3 1 1 1],'integer*4');
	elseif (strcmp(tipo, 'dtm'))
		%fwrite(fid,[0 0 128 63 0 0 128 63 205 204 76 62],'uchar');		% Non VIMAGE ???
		%fwrite(fid,[0 0 0 12 0 1 1],'integer*4');
		fwrite(fid,[0 0 128 63 0 0 128 63],'uchar');
		fwrite(fid,[11 234 229 62 0 0 0 128 0 0 0 128 251 10 141 62],'uchar');
		fwrite(fid,[12 0 1 1],'integer*4');
	elseif (strcmp(tipo, 'geoimg'))
		fwrite(fid,[0 0 128 63 0 0 128 63 0 0 128 63],'uchar');
		fwrite(fid,[0 0 0 36 0 1 1],'integer*4');
	elseif (strcmp(tipo, 'vimage'))
		fwrite(fid,[171 170 42 63 0 205],'uchar');
		fwrite(fid,[204 61 244 106 34 63 0 0],'uchar');		% Somehow it depends on the inserting coods
		fwrite(fid,[0 128 0 205 76 61 25 42 59 190 41 0 0 0],'uchar');		%	"
		fwrite(fid,[0 1 1],'integer*4');
	else
		error('Case unpredicted in write_node_block')
	end
	fwrite(fid,numel(str1),'integer*4');	fprintf(fid,str1);
	fwrite(fid,numel(str2),'integer*4');	fprintf(fid,str2);
	fwrite(fid,(1:8)*0,'integer*4');
	fwrite(fid,[0 16256],'integer*2');
	fwrite(fid,(1:12)*0,'integer*4');
	fwrite(fid,[188 251 29 66 68 4 38 66 17 17 88 193 239 238 39 193 45 208 176 197 0 0 128 63 0 0 0 0],'uchar');

%----------------------------------------------------------------------------------
function write_alignparent_block(fid)
% Write the SD_ALIGNPARENT block
	fwrite(fid,[9955 4],'integer*4');			% Tag ID, Data Length
	fwrite(fid,[0 0 1 1 1 1 (1:18)*0],'uchar');
	fwrite(fid,0,'integer*4');

%----------------------------------------------------------------------------------
function write_geoimg_block(fid)
% Write the SD_GEOIMAGE block
	fwrite(fid,[10530 0],'integer*4');			% Tag ID, Data Length
	fwrite(fid,[0 0 1 1 1 1 (1:18)*0],'uchar');    

%----------------------------------------------------------------------------------
function write_texturedtm_block(fid)
% Write the SD_TEXTUREDTM block
	fwrite(fid,[10050 0],'integer*4');			% Tag ID, Data Length
	fwrite(fid,[0 0 1 1 1 1 (1:18)*0],'uchar');

%----------------------------------------------------------------------------------
function write_vimage_block(fid, limits)
% Write the SD_VIMAGE block
	fwrite(fid,[10555 48],'integer*4');			% Tag ID, Data Length
	fwrite(fid,[0 0 1 1 1 1 (1:18)*0],'uchar');
	fwrite(fid,limits,'real*8');

%----------------------------------------------------------------------------------
function write_sonardtm_atb_block(fid)
% Write the SD_SONARDTM_ATB block
	fwrite(fid,[10007 18],'integer*4');			% Tag ID, Data Length
	fwrite(fid,[0 0 1 1 1 2 (1:18)*0],'uchar');
	fwrite(fid,[1 0 255 0 0 0 10 0 0 0 1 (1:7)*0],'uchar');	% First elem solid (1) mesh (0). Third, transparency  

%----------------------------------------------------------------------------------
function write_geoimg_atb_block(fid)
% Write the SD_GEOIMAGE_ATB block
	fwrite(fid,[10531 16],'integer*4');			% Tag ID, Data Length
	fwrite(fid,[0 0 1 1 1 1 (1:18)*0],'uchar');
	fwrite(fid,[0 0 0 0],'integer*4');

%----------------------------------------------------------------------------------
function write_vimage_atb_block(fid)
% Write the SD_VIMAGE_ATB block
	fwrite(fid,[10556 8],'integer*4');			% Tag ID, Data Length
	fwrite(fid,[0 0 1 1 1 2 (1:18)*0],'uchar');
	fwrite(fid,[0 0],'integer*4');

%----------------------------------------------------------------------------------
function write_cmap(fid, pal)
% Write the FM_CMAP block
	if (nargin == 1),	pal = uint8(round(jet(256) * 255));		end		% default colormap
	fwrite(fid,[20010 768],'integer*4');		% ID of FM_CMAP block & n of bytes in this block
	fwrite(fid,[0 0 1 1 1 1 (1:18)*0],'uchar');	% 24 bytes (seams to be an offset)
	fwrite(fid,pal','uchar');    

%----------------------------------------------------------------------------------
function write_eof(fid)
% Write the EOF block
	fwrite(fid,[999999999 0],'integer*4');
	fwrite(fid,[0 0 1 1 1 1 (1:18)*0],'uchar');
	fclose(fid);

%----------------------------------------------------------------------------------
function write_lines_or_points(fid, hAxes, limits, burnCoasts)
% Look for line and/or points Matlab objects and write them as Fledermaus
% objects. It does so (when it finds them) by calling the corresponding
% function that writes either lines or points objects in Fleder format

if (nargin < 4),    burnCoasts = 1;     end				% Default is to burn coast lines into image

ALLlineHand = findobj(hAxes,'Type','line');
ALLpatchHand = findobj(hAxes,'Type','patch');
if (isempty(ALLlineHand) && isempty(ALLpatchHand))
	return
end

if (~isempty(ALLlineHand))
	h = findobj(ALLlineHand,'Tag','Earthquakes');		% Search first for earthquakes because they have depths
	if (~isempty(h))
		write_pts(fid,h,'add',limits,'Earthquakes')		% earthquakes
		ALLlineHand = setxor(ALLlineHand, h);			% h is processed, so remove it from handles list
	end

	if (burnCoasts)     % If we had already burned the coastlines remove them from the ALLlineHand list
		h = findobj(ALLlineHand,'Tag','CoastLineNetCDF');
		if (~isempty(h)),		ALLlineHand = setxor(ALLlineHand, h);	end
		h = findobj(ALLlineHand,'Tag','PoliticalBoundaries');
		if (~isempty(h)),		ALLlineHand = setxor(ALLlineHand, h);	end
		h = findobj(ALLlineHand,'Tag','Rivers');
		if (~isempty(h)),		ALLlineHand = setxor(ALLlineHand, h);	end
	end
    
	% See if we have COASTLINES. If yes they are treated separatly (mainly because of Z and also to create an separate object)
	h = findobj(ALLlineHand,'Tag','CoastLineNetCDF');
	if (~isempty(h))
		z_level = 0;        % It will be a nonsense if the underlying grid is not topographic
		[x,y,z,count] = lines2multiseg(h,z_level);
		line_thick = get(h(1),'LineWidth');       % Line thickness
		line_color = get(h(1),'color');           % Line color
		line_props = [line_thick line_color 0 [1 1 1]];     % 0 means is not a patch and the [1 1 1] is not used
		write_line(fid,'add',x,y,z,count,limits,line_props)
		ALLlineHand = setxor(ALLlineHand, h);     % h are processed, so remove them from handles list
	end
    
	% See if we have national borders. If yes they are treated separatly
	h = findobj(ALLlineHand,'Tag','PoliticalBoundaries');
	if (~isempty(h))
		z_level = 0;		% It will be a nonsense if the underying grid is not topographic
		[x,y,z,count] = lines2multiseg(h,z_level);
		line_thick = get(h(1),'LineWidth');       % Line thickness
		line_color = get(h(1),'color');           % Line color
		line_props = [line_thick line_color 0 [1 1 1]];     % 0 means is not a patch and the [1 1 1] is not used
		write_line(fid,'add',x,y,z,count,limits,line_props)
		ALLlineHand = setxor(ALLlineHand, h);     % h are processed, so remove them from handles list
	end

	% See if we have RIVERS. If yes they are treated separatly
	h = findobj(ALLlineHand,'Tag','Rivers');
	if (~isempty(h))
		z_level = 0;		% It will be a nonsense if the underying grid is not topographic
		[x,y,z,count] = lines2multiseg(h,z_level);
		line_thick = get(h(1),'LineWidth');       % Line thickness
		line_color = get(h(1),'color');           % Line color
		line_props = [line_thick line_color 0 [1 1 1]];     % 0 means is not a patch and the [1 1 1] is not used
		write_line(fid,'add',x,y,z,count,limits,line_props)
		ALLlineHand = setxor(ALLlineHand, h);     % h are processed, so remove them from handles list
	end
    
	% See if we have CONTOUR lines. If yes we fish their depths
	h = findobj(ALLlineHand,'Tag','contour');
	if (~isempty(h))
		dz = abs(limits(6) - limits(5)) * 0.01;		% I smell a fleder bug here, so add a small cte to z level
		for (i = 1:length(h))
			z_level = get(h(i),'UserData');
			set(h(i),'UserData',z_level + dz)		% <== trick to elevate line by dz, because UD prevails over z_level 
			[x,y,z,count] = lines2multiseg(h(i),z_level);
			set(h(i),'UserData',z_level)			% <== reset original level            
			line_thick = get(h(i),'LineWidth');		% Line thickness
			line_color = get(h(i),'color');			% Line color
			line_props = [line_thick line_color 0 [1 1 1]];     % 0 means is not a patch and the [1 1 1] is not used
			write_line(fid,'add',x,y,z,count,limits,line_props)
		end
		ALLlineHand = setxor(ALLlineHand, h);       % h are processed, so remove them from handles list
	end
end

if (~isempty(ALLlineHand))
	% This section deals with repeated line types. The point is that having many objects makes the
	% rendering slow. So I assimilate all lines of the same type (same line thickness and color) into
	% a single multisegment line. This makes the rendering much faster.
	LineWidth = get(ALLlineHand,'LineWidth');
	if (iscell(LineWidth)),     LineWidth = cell2mat(LineWidth);    end    
	LineColor = get(ALLlineHand,'Color');
	if (iscell(LineColor)),     LineColor = cell2mat(LineColor);    end
	[tmp,ind] = sortrows([LineWidth LineColor]);
	hands_sort = ALLlineHand(ind);			% Sort also the handles according to the previous sorting cretirea
	difs = diff([tmp(1,:); tmp]);			% Repeat first row to account for the decrease 1 resulting from diff
	id_row = find(difs ~= 0);				% Find the lines of different type
	id_row = unique(id_row);				% Get rid of repeated values
	id_row = [1; id_row];					% Make id_row start at one
	id_row(end+1) = size(tmp,1);			% Add the last row as well (for the algo)
	hands = cell(1,numel(id_row)-1);
	for (i = 1:numel(id_row)-1)
		hands{i} = hands_sort(id_row(i):id_row(i+1)-1);
		if (numel(hands{i}) > 1)			% Make sure that the following procedure applyies only to repeated line types
			line_thick = get(hands{i}(1),'LineWidth');			% Line thickness
			line_color = get(hands{i}(1),'color');				% Line color
			line_props = [line_thick line_color 0 1 1 1];		% 0 means is not a patch and the [1 1 1] is not used
			[x,y,z,count] = lines2multiseg(hands{i},limits(end));
			write_line(fid,'add',x,y,z,count,limits,line_props)
			ALLlineHand = setxor(ALLlineHand, hands{i});		% hands{i} are processed, so remove them from handles list
		end
	end

    % OK, now if we still have lines, they must be of different line type
    for (i = 1:numel(ALLlineHand))
        z_level = limits(end);          % Default to z_max (but I have to do something clever)
        [x,y,z,count] = lines2multiseg(ALLlineHand(i),z_level);
        line_thick = get(ALLlineHand(i),'LineWidth');		% Line thickness
        line_color = get(ALLlineHand(i),'color');			% Line color
        line_props = [line_thick line_color 0 1 1 1];		% 0 means is not a patch and the [1 1 1] is not used
		if ( ~strcmp(get(ALLlineHand(i), 'LineStyle'), 'none') )
			write_line(fid,'add',x,y,z,count,limits,line_props)
		else
			write_pts(fid,ALLlineHand(i),'add',limits)
		end
    end
end     % end  -> if (~isempty(ALLlineHand))<-

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
        write_line(fid,'add',xx,yy,zz,count,limits,line_props)
        ALLpatchHand = setxor(ALLpatchHand, telhasHand_d);  % telhasHand_d is processed, so remove it from handles list
    end
    
    if (~isempty(telhasHand_r))						% Now the reverse telhas
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
        write_line(fid,'add',xx,yy,zz,count,limits,line_props)
        ALLpatchHand = setxor(ALLpatchHand, telhasHand_r);  % telhasHand_r is processed, so remove it from handles list
        clear xx yy zz;
    end
end
       
if (~isempty(ALLpatchHand))                 % Now see if we still have more patches
    n_lines = length(ALLpatchHand);
    for (i = 1:n_lines)                     % Do one by one. If they are many, the rendering is slow
        x = get(ALLpatchHand(i),'XData');        y = get(ALLpatchHand(i),'YData');
        if (isempty(x)),    continue;       end
        line_thick = get(ALLpatchHand(i),'LineWidth');   % Line thickness
        line_color = get(ALLpatchHand(i),'EdgeColor');       % Line color
        
        patch = 1;
        patch_color = get(ALLpatchHand(i),'FaceColor');
        if (strcmp(patch_color,'none'))     % OK, if we have no color we treat it just like an ordinary polyline
            patch = 0;
            patch_color = [1 1 1];          % Not used anyway
        end
        line_props = [line_thick line_color patch patch_color];
        count = numel(x);
		z = get(ALLpatchHand(i),'UserData');	% See if we have depth info
		if (isempty(z)),        z = repmat(z_level,1,count);	end
        write_line(fid,'add',x,y,z,count,limits,line_props)
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
		fprintf(fid,'%s\n%s\n%s\f\n','%% TDR 2.0 Binary','Created by:    Mirone Tech!','%%');
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
	fwrite(fid,ones(1,8)*205,'uchar');			% ??
	fwrite(fid,code1,'integer*4');				% ??
	fwrite(fid,[lim_reg lim_line],'real*8');

	if (iscell(x))
		for (i = 1:n_segments)
			np_s = numel(x{i});                % Number of points in this segment
			fwrite(fid,[np_s 0],'integer*4');
			fwrite(fid,[x{i}; y{i}; z{i}],'real*8');
		end
	else
		fwrite(fid,[np 0],'integer*4');			% n points
		fwrite(fid,[x; y; z],'real*8');
	end

	write_geo(fid,'add',lim_reg)				% Write a GEOREF block
	write_cmap(fid)								% Write a FM_CMAP block

	fwrite(fid,[10526 32],'integer*4');         % ID of block SD_LINES3D_ATB and n bytes in this block
	fwrite(fid,[0 0 1 1 1 6 (1:18)*0],'integer*1');     % The 6 is a number of version

	fwrite(fid,[l_color(1:3) 255],'uchar');     % Don't know what is the last 255
	fwrite(fid,[1 1 ColorBy 0 0 1 l_thick],'integer*4'); % First 1 is 'gap', but don't know what are the others

%----------------------------------------------------------------------------------
function write_pts(fid,hand,mode,limits,opt)
% HAND -> handles of the line (points) object
% MODE = FIRST or ADD. Where FIRST indicates that the TDR object starts to created here.
% OPT, when it exists, is = 'Earthquakes'

    if (nargin == 4),   opt = [];   end
    
    symb = 6;	% 0 -> circle; 1 -> square; 2 -> cross hair; 3 -> cube; 4 -> dyamond; 5 -> cylinder; 6 -> sphere; 7 -> point
    PointRad = 0.02;						% Symbol radius
    ColorBy = 1;							% 0 -> Solid; 1 -> Line Height (Z); 2 -> Attribute
    LabelSize = 0.502;
    n_col = 3;								% N of columns
	if (~isa(hand,'cell') && ~ishandle(hand(1)))	% We need cells in this case
		hand = {hand};
	end
	n_groups = numel(hand);					% N of different point ensembles

	if (strcmpi(mode,'first'))				% The TDR object starts here
		fprintf(fid,'%s\n%s\f\n','%% TDR 2.0 Binary','%%');
	end

	for (i = 1:n_groups)
		this_CB = ColorBy;					% Earthquakes are colored by symbol color, others by color scale
		if (ishandle(hand(i)))
			xx = get(hand(i),'XData');	yy = get(hand(i),'YData');
			PointRad = get(hand(i),'MarkerSize') / 72 * 2.54 / 7;   % Symbol size. The 7 is an ad-hoc corr factor
		else
			[m n] = size(hand{i});		% We must allow for column or row vectors in input
			if (m > n)		xx = hand{i}(:,1)';		yy = hand{i}(:,2)';
			else			xx = hand{i}(1,:);		yy = hand{i}(2,:);
			end
		end
		if (strcmp(opt,'Earthquakes'))  % Test those first because they have a z
			zz = -double(getappdata(hand(i),'SeismicityDepth')) * 100;    % zz is now in meters positive up
			if (size(zz,1) > 1)         % We need them as a row vector for fwrite
				zz = zz';
			end
			this_CB = 0;	symb = 6;	% Use spheres for epicenters and color by symbol (solid)
		else
			if (ishandle(hand(i)))
				zz = get(hand(i),'ZData');
				if (isempty(zz))	zz = getappdata(hand(i),'ZData');	end		% Try this too
				if (isempty(zz))	zz = get(hand(i),'UserData');	end		% Try this too
			else
				if (min(m,n) == 3)	% Have Z
					if (m > n)		zz = hand{i}(:,3)';
					else			zz = hand{i}(3,:);
					end
				end
			end
		end

		np = numel(xx);					% Number of points in this segment
		if (isempty(zz))				% We need to have a z. Obviously it has to be 0
			zz = repmat(limits(6),1,np);
		else							% Since we have a z, we need to update limits
			limits(5) = min(zz);        limits(6) = max(zz);
		end

		fwrite(fid,10520,'integer*4');		% 10520 is the SD_POINT3D code
		%fwrite(fid,np*3*8+64,'integer*4');	% data length
		fwrite(fid,np*3*8+20,'integer*4');	% data length
		fwrite(fid,[0 0 1 1 1 2 (1:18)*0],'integer*1');
		fwrite(fid,[np n_col 0 0 1],'integer*4');	% Number of points
		fwrite(fid,[xx; yy; zz;],'real*8');

		limits(1:4) = [min(xx) max(xx) min(yy) max(yy)];	% Required by Fleder7 (???)
		write_geo(fid,'add',limits)					% Write a GEOREF block
		if (this_CB > 0),	write_cmap(fid),	end	% Write a FM_CMAP block (no need if color is solid)

		fwrite(fid,[10521 35],'integer*4');			% 10520 is the SD_POINT3D_ATB code, 35 -> Data Length
		fwrite(fid,[0 0 1 1 1 5 (1:18)*0],'integer*1');
		fwrite(fid,[1 this_CB symb],'integer*4');
		fwrite(fid,[PointRad LabelSize],'real*4');
		fwrite(fid,[0 1 1],'integer*4');
		if (ishandle(hand(i)))
			cor = uint8(get(hand(i),'MarkerFaceColor')*255);
		else
			cor = [255 255 255];	% ??
		end
		fwrite(fid,cor,'uint8');					% symbol's color
	end

%----------------------------------------------------------------------------------
function [x,y,z,count] = lines2multiseg(hands,z_level)
% Convert a collection of lines whose handles are HANDS into a single multiline (cell array) array.
% If HANDS is a scalar, it will check if that line is broken with NaNs. If yes, it will be
% converted into a multiseg cell array. So we can call this function safely even when we don't
% know exactly what is contained in HANDS. The test of if x,... is a single or multisegment line
% are carried out inside the write_line function
% Z_LEVEL, if transmited will be used to set the line height, otherwise ZERO will be used.
% If lines have Z in UserData AND take precedence over transmited z_level

	if (nargin == 1),   z_level = 0;    end

	x = get(hands,'XData');		y = get(hands,'YData');		z = get(hands,'UserData');
	got_Z = false;
	
	if (iscell(x))
		n_lines = length(x);			% Number of lines
		if (~isempty(z{1}) && isnumeric(z{1})),		got_Z = true;	end		% Line(s) have Z info
		if (~got_Z)						% Common line type (2D)
			z = x;						% Create a cell array of the same size as x
			for (i = 1:n_lines)
				z{i} = ones(1,numel(x{i})) * z_level;
			end
		end
	else
		n_lines = 1;					% When it's not a cell array n_lines must be one
		if (~isempty(z) && isnumeric(z)),	got_Z = true;	end		% Line(s) have Z info
		if ( got_Z && (numel(x) ~= numel(z)) ),		z = repmat(z,size(x));		end
		if (~got_Z),		z = zeros(1,numel(x));	end
	end

	id_with_nan = false(1,n_lines);
	count = 0;							% Counter of the total number of points
	for (i = 1:n_lines)
		if (iscell(x))					% Multi lines case 
			if (any(isnan(x{i})))		% See if we have NaNs in this line. If yes we must treat it as a multiseg
				id_with_nan(i) = 1;		% Yes we have. Mark this handle line to be processed later
				continue				% Jump this line for the time beeing
			end
			count = count + length(x{i});
		else							% Single line, but with possible NaNs
			if (any(isnan(x)))			% See if we have NaNs in this line. If yes we must treat it as a multiseg
				id_with_nan = 1;		% Yes we have. Mark this handle line to be processed later
				continue				% Jump this line for the time beeing
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
		for (i = 1:numel(id_where))
			xt = get(hands(id_where(i)),'XData');
			yt = get(hands(id_where(i)),'YData');
			zt = z;
			id_nan = find(xt ~= xt);    % Find the NaNs
			id = find(diff(id_nan) == 1) + 1;   % Account for contiguous NaNs
			if (~isempty(id))           % Found contiguous NaNs
				xt(id_nan(id)) = [];
				yt(id_nan(id)) = [];    % Remove them
				zt(id_nan(id)) = [];
				id_nan = find(xt ~= xt);% Find the new position of the now non-contiguous NaNs
			end
			id_nan = [0 id_nan];        % Used to make it start at one        
			
			n_segments = length(id_nan)-1;
			xx = cell(n_segments,1);    yy = xx;    zz = xx;
			for (k = 1:n_segments)
				xx{k} = xt(id_nan(k)+1:id_nan(k+1)-1);
				yy{k} = yt(id_nan(k)+1:id_nan(k+1)-1);
				%zz{k} = repmat(z_level,1,numel(xx{k}));
				zz{k} = zt(id_nan(k)+1:id_nan(k+1)-1);					% <== Idem
				count = count + numel(xx{k});
			end

			if (iscell(x))              % Apend to eventual existing lines that had no NaNs
				x = [x; xx];		y = [y; yy];	z = [z; zz];
			else                        % Here we have a single line with NaNs
				x = xx;				y = yy;			z = zz;
			end
		end
	end

%----------------------------------------------------------------------------------
function img = burnLines(hFig, hAxes, img, inplace)
% Burn the coastlines directly into de IMG image. We do it because also the Fleder
% doesn't work as advertized. Lines are awfully drapped on surfaces

	handMir = guidata(hFig);
	head = handMir.head;
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
			lt = 8;				% LINE_TYPE -> 8 connectivity (default)
			if (line_thick <= 1)
				lt = 16;		% antialiased line
			end
			if (inplace)
				cvlib_mex('poly',img,xy,line_color,line_thick,lt)
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

if (any(isnan(x)))				% We have NaNs in this line. Treat it as a multiseg
	id_nan = find(x ~= x);				% Find the NaNs
	id = find(diff(id_nan) == 1) + 1;	% Account for contiguous NaNs
	if (~isempty(id))					% Found contiguous NaNs
		x(id_nan(id)) = [];
		y(id_nan(id)) = [];				% Remove them
		id_nan = find(x ~= x);			% Find the new position of the now non-contiguous NaNs
	end	
	id_nan = [0; id_nan];				% Used to make it start at one

    n_segments = length(id_nan)-1;
    xy = cell(n_segments,1);
    for (k = 1:n_segments)
        xx = round( getPixel_coords(lims(5),lims(1:2),x(id_nan(k)+1:id_nan(k+1)-1)) -1);  % -1 because we need
        yy = round( getPixel_coords(lims(6),lims(3:4),y(id_nan(k)+1:id_nan(k+1)-1)) -1);  % zero based indexes
        xy{k} = [xx yy];
    end
else
	xy = [x y];
end

% -------------------------------------------------------------------------------------
function pix_coords = getPixel_coords(img_length, XData, axes_coord)
% Convert coordinates from axes (real coords) to image (pixel) coordinates.
% IMG_LENGTH is the image width (n_columns)
% XDATA is the image's [x_min x_max] in axes coordinates
% AXES_COORD is the (x,y) coordinate of the point(s) to be converted

	slope = (img_length - 1) / (XData(end) - XData(1));
	if ((XData(1) == 1) && (slope == 1))
		pix_coords = axes_coord;
	else
		pix_coords = slope * (axes_coord - XData(1)) + 1;
	end
