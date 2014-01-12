function out = structuring_elem(x,y,p)
% structuring_elem  editor for creating structuring elements for morphology ops.
% 
% usage:
% structuring_elem(x,y,pixel)
%     x - horizontal size (7 as deault)
%     y - horizontal size (7 as deault)
%     pixel - size of each dot (15 as deafult)
% 
% structuring_elem
%     same as structuring_elem(7,7,15)
%
%	OUT is a structure with these fields
%		operation  - One of: 'dilate' 'erode' 'open' 'close' 'gradient' 'tophat' 'blackhat'
%		rows	- number of rows
%		cols	- number of columns
%		iterations - number of iterations
%		values  - the structuring element (a double with 1s & 0s)
%		shape   - One of: 'custom' 'rectangle' 'ellipse' 'diamond' 'cross'
%		anchorX - Relative horizontal offset of the anchor point
%		anchorY - Relative vertical offset of the anchor point
%
% Inspired on iconeditor of  Elmar Tarajan [MCommander@gmx.de]

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

	switch nargin
		case 0
			H.x = 7;
			H.y = 7;
			H.p = 15;
		case 3
			H = struct('x',x, 'y',y, 'p',p);
	end

	operations = {'dilate' 'erode' 'open' 'close' 'gradient' 'tophat' 'blackhat'};

	% main figure
	H.nan = [.75 .74 .70 ; .83 .81 .78 ; .90 .89 .85 ; .83 .81 .78];
	ssz = get(0,'ScreenSize');
	DX = 90;
	figSizeW = max(max(H.x*H.p+14, 222), H.x*5+60) + 17 + DX;
	figSizeH = max((H.y + 1)*H.p + 85, 190);
	location = (ssz(3:4) - [figSizeW figSizeH]) / 2;

	gotPrevFig = findobj('Type','figure','Name','Structuring element', 'Tag','StructElem');
	if (isempty(gotPrevFig))	gotPrevFig = false;		end

	if (~gotPrevFig)
		H.gcf = figure('menubar','none', 'resize','off', 'NumberTitle','off', ...
						'units','pixel', 'pos',[location figSizeW figSizeH], ...
						'Name','Structuring element', 'Tag','StructElem', ...
						'visible','off', 'WindowStyle','modal');
	else
		x = H.x;	y = H.y;	p = H.p;	% Backup those
		H = getappdata(gotPrevFig,'H');
		H.gcf = gotPrevFig;
		H.x = x;	H.y = y;	H.p = p;
		pos = get(H.gcf, 'Pos');
		pos(3) = figSizeW;		pos(4) = figSizeH;
		if ((pos(1) + pos(3)) > ssz(3))		% Check that it doesn't screen overflows
			pos(1) = pos(1) - ( ssz(3) - (pos(1) + pos(3)) );
		end
		if ((pos(2) + pos(4)) > ssz(4))
			pos(2) = pos(2) - ( ssz(4) - (pos(2) + pos(4)) );
		end
		set(H.gcf, 'pos',pos, 'UserData', []);	% Resize to fit with the new dims
	end

	% paint area
	pos = [(figSizeW-H.x*H.p-11+DX)/2 22 H.x*H.p+14 H.y*H.p+37];
	if (~gotPrevFig)
		H.bgAx = axes('parent',H.gcf,'units','pixel','pos',pos','XTick',[],'YTick',[],'Color',[.45 .45 .45]);
		H.gca = axes('parent',H.gcf,'units','pixel','pos',[pos(1)+6 pos(2)+22 H.x*H.p H.y*H.p], ...
					 'Xlim',[0 H.x*H.p],'Ylim',[0 H.y*H.p],'vis','off');
	else
		set(H.gca, 'pos',[pos(1)+6 pos(2)+22 H.x*H.p H.y*H.p],'Xlim',[0 H.x*H.p],'Ylim',[0 H.y*H.p],'vis','off');
		set(H.bgAx, 'pos',pos);
	end

	% paint patch
	xx = H.p*reshape(repmat([0:H.x-1;1:H.x;1:H.x;0:H.x-1],H.y,1),4,H.x*H.y);
	yy = H.p*repmat([(H.y-1):-1:0;(H.y-1):-1:0;H.y:-1:1;H.y:-1:1],1,H.x);
	clr = repmat(reshape(H.nan,[4,1,3]), 1, H.x * H.y);
	if (~gotPrevFig)
		H.draw = patch(xx,yy,clr,'EdgeColor',[.45 .45 .45],'Parent',H.gca);
	else
		set(H.draw, 'XData',xx, 'YData',yy, 'CData',clr)
	end

	% icon's
	img = icondata;
	if ((pos(1)+28) < (floor(pos(1)+pos(3)/2))),		xs = pos(1)+7;		d = 6;
	else												xs = pos(1)-25;		d = 2;
	end

	if (~gotPrevFig)
		for (i = 1:2)
			H.iconaxes(i) = axes('parent',H.gcf,'units','pixel','pos',[xs+(i-1)*(15+d) pos(2)+4 15 15]);
			H.icon(i) = image(uint8(repmat(img(:,:,i),[1 1 3])));
			set(H.iconaxes(i),'vis','off', 'XTick',[],'XColor',[.9 .9 .9],'YTick',[],'YColor',[.9 .9 .9])
		end
		set(H.iconaxes(1),'Vis','on')
		set(H.icon,'buttondownfcn',{@iconclick,H})
		xs = max( ceil(pos(1)+pos(3)/2), pos(1)+H.x*H.p-89+3*20 );
		for (i = 3:4)
			H.iconaxes(i) = axes('parent',H.gcf,'units','pixel','pos',[xs+(i-3)*(15+d) pos(2)+4 15 15]);
			H.icon(i)= image(uint8(228-repmat(img(:,:,i),[1 1 3])),'UserData',uint8(repmat(img(:,:,i),[1 1 3])));
			set(H.iconaxes(i),'vis','off')
		end

		cfg = struct('parent', H.gcf, 'style', 'text', 'enable','inactive', 'fontsize', 8, ...
			'backgroundcolor',[.45 .45 .45], 'foregroundcolor', [.9 .9 .9], 'horizontalalignment','left');
		H.pos = uicontrol(cfg,'Pos',[pos(1)+6 pos(2)+H.y*H.p+22 35 12]);
		cfg.foregroundcolor = [0 0 0];
		cfg.backgroundcolor = get(H.gcf,'Color');
		cfg.fontsize = 10;
		cfg.fontweight = 'bold';
		cfg.string = 'X';
		cfg.FontName = 'Helvetica';
		H.txt = uicontrol(cfg,'Pos',[33 figSizeH-22 12 15]);
	else
		for (i = 1:2),		set(H.iconaxes(i),'pos',[xs+(i-1)*(15+d) pos(2)+4 15 15]),	end
		xs = max( ceil(pos(1)+pos(3)/2), pos(1)+H.x*H.p-89+3*20 );
		for (i = 3:4),		set(H.iconaxes(i),'pos',[xs+(i-3)*(15+d) pos(2)+4 15 15]),	end
		set(H.pos, 'Pos',[pos(1)+6 pos(2)+H.y*H.p+22 40 12])
		set(H.txt, 'Pos',[33 figSizeH-22 12 15])
	end

	if (~gotPrevFig)
		H.editNrow = uicontrol('parent',H.gcf,'style','edit','units','pixel', ...
			'pos',[6 figSizeH-25 25 21], 'String',H.y, 'Tooltip','Number of rows (Odd numbers please)');
		H.editNcol = uicontrol('parent',H.gcf,'style','edit','units','pixel', ...
			'pos',[45 figSizeH-25 25 21], 'String',H.x, 'Tooltip','Number of columns (Odd numbers please)');
		H.resize = uicontrol('parent',H.gcf,'style','pushbutton','units','pixel', 'string','Apply', ...
			'pos',[6 figSizeH-40 64 15], 'Call',{@resize,H},'Tooltip','Resize figure by new dims');
		H.OK = uicontrol('parent',H.gcf,'style','pushbutton','units','pixel', 'string','OK', ...
			'pos',[figSizeW-41 1 40 21], 'FontWeight','bold', 'Call',@ok_CB);

		cfg.string = 'Iterations';
		uicontrol(cfg,'Pos',[6 figSizeH-68 70 18]);
		H.editIter = uicontrol('parent',H.gcf,'style','edit','units','pixel', ...
			'pos',[76 figSizeH-68 20 21], 'String','1', 'Tooltip','Number of iterations');

		dz = 17;	rad_y0 = figSizeH-90;	bg_color = get(H.gcf,'Color');
		props = {'parent',H.gcf,'style','radio','units','pixel','backgroundcolor',bg_color};
		H.rad(1) = uicontrol(props{:}, 'pos',[6 rad_y0 80 dz], 'str','dilate','UserData',1);
		H.rad(2) = uicontrol(props{:}, 'pos',[6 rad_y0-dz 80 dz], 'str','erode','UserData',2);
		H.rad(3) = uicontrol(props{:}, 'pos',[6 rad_y0-2*dz 80 dz], 'str','open','Value',1,'UserData',3,'Tooltip','Erosion + Dilation');
		H.rad(4) = uicontrol(props{:}, 'pos',[6 rad_y0-3*dz 80 dz], 'str','close','UserData',4,'Tooltip','Dilation + Erosion');
		H.rad(5) = uicontrol(props{:}, 'pos',[6 rad_y0-4*dz 80 dz], 'str','gradient','UserData',5,'Tooltip', 'Dilate - Erode');
		H.rad(6) = uicontrol(props{:}, 'pos',[6 rad_y0-5*dz 80 dz], 'str','tophat','UserData',6,'Tooltip','Original - open');
		H.rad(7) = uicontrol(props{:}, 'pos',[6 rad_y0-6*dz 80 dz], 'str','blackhat','UserData',7,'Tooltip','close - Original');
		set(H.rad,'Call',{@radios_CB,H})

		H.push = zeros(1,5);
		for (i = 5:-1:1)
			H.push(i) = uicontrol('parent',H.gcf,'style','pushbutton', ...
				'cdata',ones(H.y,H.x,3)*NaN,'backgroundcolor',[.8 .8 .8], ...
				'units','pixel','pos',[figSizeW-(6-i)*(21+5)-5 figSizeH-25 21 21], ...
				'Call',{@setpushbutton,H,i});
		end
		ico = uint8(ones(19)*255);	ico(10,:) = 0;	ico(:,10) = 0;
		set(H.push(5), 'CData', repmat(ico,[1 1 3]), 'Tooltip', 'Cross')
		set(H.push(4), 'CData', uint8(zeros(19,19,3)), 'Tooltip', 'Rectangle')
		[x,y] = meshgrid(-9:9);		ico = uint8( double((abs(x) + abs(y)) > 9) * 255);
		set(H.push(3), 'CData', repmat(ico,[1 1 3]), 'Tooltip', 'Diamond')
		ico = uint8( double(((x.^2) + (y.^2)) > 70) * 255);
		set(H.push(2), 'CData', repmat(ico,[1 1 3]), 'Tooltip', 'Disc')
		ico = uint8( round(rand(17))*255 );
		set(H.push(1), 'CData', repmat(ico,[1 1 3]), 'Tooltip', 'Custom')
	else
		set(H.editNrow, 'pos',[6 figSizeH-25 25 21], 'String',H.y)
		set(H.editNcol, 'pos',[45 figSizeH-25 25 21], 'String',H.x)
		set(H.resize, 'pos',[6 figSizeH-40 64 15], 'Call',{@resize,H})
		set(H.OK, 'pos',[figSizeW-41 1 40 21])
		for (i = 5:-1:1)
			set(H.push(i), 'pos',[figSizeW-(6-i)*(21+5)-5 figSizeH-25 21 21], 'Call',{@setpushbutton,H,i})
		end
	end

	% -------------- prepare "undo"-feature ----------------------------------------
	setappdata(H.gcf,'undoid',1)
	setappdata(H.gcf,'undo',{get(H.draw,'FaceVertexCData')});

	% --------------- prepare main callbacks ---------------------------------------
	set(H.draw,'buttondownfcn',{@mouse_down,H})
	set(H.gcf,'WindowButtonMotionFcn',{@mouse_move,H,[],0},'vis','on','HandleVisibility','callback')
	setappdata(H.gcf,'H',H)

	% Set a 3x3 rectangle as default
	setpushbutton(H.push(1),[],H,6)

	if (gotPrevFig),	return,		end

	uiwait(H.gcf);			% UIWAIT makes it wait for user response
	if (nargout)
		try
			values = get(H.gcf,'UserData');
			r = get(H.rad, 'Val');
			r = logical([r{:}]);
			[ny, nx] = size(values);
			out = struct('operation',operations{r}, 'values', values, 'shape', 'custom', 'cols', nx, 'rows', ny, 'iterations',1);
			out.anchorX = floor((nx-1) / 2);
			out.anchorY = floor((ny-1) / 2);

			% Test if we have any of the default OpenCV structuring element types
			test = true(ny,nx);						donne = false;
			if (isequal(values,test)),				out.shape = 'rectangle';	donne = true;
			else									test = makeStrel(nx,ny,'cross');
			end
			if (~donne && isequal(values,test)),	out.shape = 'cross';		donne = true;
			else									test = makeStrel(nx,ny,'ellipse');
			end
			if (~donne && isequal(values,test)),	out.shape = 'ellipse';		end
			iter = round(str2double(get(H.editIter,'string')));
			if (~isnan(iter) && iter > 1)			out.iterations = iter;		end
		catch
			out = [];
		end
	end
	if (ishandle(H.gcf)),	delete(H.gcf),		end

%-------------------------------------------------------------------------------
function  radios_CB(obj,event,H)
% Set all the other radios to 0
	if (~get(obj,'Value')),		set(obj,'Val',1),		return,		end
	ind = (((1:7) - get(obj, 'UserData')) == 0);
	set(H.rad(~ind), 'Val', 0)

%-------------------------------------------------------------------------------
function cdata = ok_CB(obj,event)
	H = getappdata(get(obj,'Parent'),'H');
	cdata = get(H.draw,'CData');
	cdata(1,any(cdata(:,:,1)-repmat(cdata(1,:,1),4,1)),:) = NaN;
	cdata = reshape(cdata(1,:,:),[H.y H.x 3]);
	cdata = cdata(:,:,1);
	ind = isnan(cdata);
	cdata(ind) = 0;		cdata(~ind) = 1;

	% Wipe out full zero row/columns arround the core element (inefficient, but the array is small)
	while (~any(cdata(:,end))),		cdata(:,end) = [];		end		% Right
	while (~any(cdata(:,1))),		cdata(:,1) = [];		end		% Left
	while (~any(cdata(1,:))),		cdata(1,:) = [];		end		% Top
	while (~any(cdata(end,:))),		cdata(end,:) = [];		end		% Bot
	set(H.gcf,'UserData',cdata)
	uiresume(H.gcf);

%-------------------------------------------------------------------------------
function mouse_down(obj,cnc,H)
	a = get(H.gcf,'SelectionType');
	b = find(strcmp(get(H.iconaxes,'visible'),'on'));

	if (strcmp(a,'normal') && b == 1)		% PAINT
	   color = zeros(4,3);

	elseif (strcmp(a,'normal') && b == 2)	% CLEAR
	   color = H.nan;

	elseif (strcmp(a,'alt'))				% CLEAR (right mouse button)
	   if (find(strcmp(get(H.iconaxes,'visible'),'on')) == 4)
		  pos = floor(get(H.gca,'CurrentPoint')/H.p);
		  pos = pos(1)*H.y+(H.y+1)-pos(3)-1;
		  img = get(H.draw,'cdata');
		  clr1 = img(:,pos,:);
		  clr2 = reshape(H.nan,[4 1 3]);
		  if ~isequal(clr1,clr2)
			 set(gco,'Cdata',recfill(img,pos,H.y,clr1,clr2))
			 undo_feature('add',[],H)
		  end
		  return
	   end
	   color = H.nan;
	   set(H.iconaxes,'visible','off')
	   set(H.iconaxes(1),'visible','on')

	else
	   return
	end

	set(gcf,'WindowButtonMotionFcn',{@mouse_move,H,color,1},'WindowButtonUpFcn',{@mouse_up,H})
	mouse_move([],[],H,color,1)

%-------------------------------------------------------------------------------
function mouse_move(fig,cnc,H,color,flag)
	pos = floor(get(H.gca,'CurrentPoint')/H.p);
	id = (pos(1)*H.y+(H.y-pos(3)-1))*4+1;
	clr = get(H.draw,'FaceVertexCData');
	if ~(H.x-pos(1)<1 || H.y-pos(3)<1 || pos(1)<0 || pos(3)<0)
		if flag
			clr(id:id+3,:) = color;
			set(H.draw,'FaceVertexCData',clr);
		end
		set(H.pos,'string',sprintf('%02d,%02d',pos(1)+1,pos(3)+1))
	else
		set(H.pos,'string','')
	end

%-------------------------------------------------------------------------------
function mouse_up(fig,cnc,H) 
	set(fig,'WindowButtonMotionFcn',{@mouse_move,H,[],0},'WindowButtonUpFcn','')
	undo_feature('add',[],H)

%-------------------------------------------------------------------------------
function setpushbutton(obj,cnc,H,n)
	%H = getappdata(H.gcf,'H');
	setappdata(H.gcf,'undoid',1)
	undo_feature(obj,[],H,-1)

	switch n
		case 1		% Custom
			return
		case 2		% Ellipse
			nhood = makeStrel(H.x,H.y,'ellipse');
		case 3		% Diamond
			nhood = makeStrel(H.x,H.y,'diamond');
		case 4		% Rectangle
			nhood = true(H.y, H.x);
		case 5		% Cross
			nhood = makeStrel(H.x,H.y,'cross');
		case 6		% This case is not "public" but used to set a default rect of 3x3
			nhood = false(H.y, H.x);
			row_c = floor((H.y+1) / 2);
			col_c = floor((H.x+1) / 2);
			nhood(row_c-1:row_c+1, col_c-1:col_c+1) = true(3, 3);
	end

	clr = get(H.draw,'FaceVertexCData');
	ind = 1:4:4*numel(nhood);
	ind(~nhood) = [];
	ind = [ind; ind+1; ind+2; ind+3];
	clr(ind(:),:) = 0;

	set(H.draw,'FaceVertexCData',clr);
	undo_feature('add',[],H)

%-------------------------------------------------------------------------------
function img = makeStrel(nx,ny,type)
%
	r = floor((nx-1) / 2);
	A = floor((nx-1) / 2);
	B = floor((ny-1) / 2);
	if (strncmp(type, 'ell', 3))		% Ellipse
		[x,y] = meshgrid(-A:A,-B:B);
		img = ((x.^2) / A^2 + (y.^2) / B^2) <= 1;
	elseif (strncmp(type, 'rec', 3))	% Rectangle
		img = true(ny, nx);
	elseif (strncmp(type, 'dia', 3))	% Diamond
		[x,y] = meshgrid(-A:A,-B:B);
		img = (abs(x) + abs(y)) <= r;
	elseif (strncmp(type, 'cro', 3))	% Cross
		img = false(ny, nx);
		row_c = floor((ny+1) / 2);
		col_c = floor((nx+1) / 2);
		img(row_c,:) = true;
		img(:,col_c) = true;
	end

%-------------------------------------------------------------------------------
function img = recfill(img,pos,step,clr1,clr2)
	id = floor((pos-1)/step)*step+1:pos;
	tmp1 = max(id(find(~all(all(img(:,id,:) == repmat(clr1,[1 length(id) 1]),3)))+1));
	if isempty(tmp1)
		tmp1 = id(1);
	end

	id = pos:ceil(pos/step)*step;
	tmp2 = min(id(find(~all(all(img(:,id,:) == repmat(clr1,[1 length(id) 1]),3)))-1));
	if isempty(tmp2)
		tmp2 = id(end);
	end

	img(:,tmp1:tmp2,:) = repmat(clr2,[1 length(tmp1:tmp2)]);

	for (n = tmp1:tmp2)
	   if n-step>0 && all(all(img(:,n-step,:)==clr1,3))
		  img = recfill(img,n-step,step,clr1,clr2);
	   end

	   if ( n + step <= size(img,2) && all(all(img(:,n+step,:) == clr1,3)) )
		  img = recfill(img,n+step,step,clr1,clr2);
	   end
	end

%-------------------------------------------------------------------------------
function undo_feature(obj,cnc,H,step)
	undo = getappdata(H.gcf,'undo');
	undoid = getappdata(H.gcf,'undoid');
	if ~ishandle(obj)
	   switch obj
		  case 'reset'
			 setappdata(H.gcf,'undoid',1)
			 setappdata(H.gcf,'undo',{get(H.draw,'FaceVertexCData')});
			 set(H.icon(3),'CData',uint8(228-double(get(H.icon(3),'UserData'))),'buttondownfcn','')
			 set(H.icon(4),'CData',uint8(228-double(get(H.icon(4),'UserData'))),'buttondownfcn','')

		  case 'add'
			 undo = undo(1:undoid);
			 undo{undoid+1} = get(H.draw,'FaceVertexCData');
			 setappdata(H.gcf,'undo',undo);
			 setappdata(H.gcf,'undoid',undoid+1);

			 set(H.icon(3),'CData',get(H.icon(3),'UserData'),'buttondownfcn',{@undo_feature,H,-1});
			 set(H.icon(4),'CData',uint8(228-double(get(H.icon(4),'UserData'))),'buttondownfcn','')
	   end
	else
	   undoid = min(max(undoid+step,1),length(undo));
	   set(H.draw,'FaceVertexCData',undo{undoid})
	   setappdata(H.gcf,'undoid',undoid);

	   if (undoid < length(undo))
		  set(H.icon(4),'CData',get(H.icon(4),'UserData'),'buttondownfcn',{@undo_feature,H,1});
	   else
		  set(H.icon(4),'CData',uint8(228-double(get(H.icon(4),'UserData'))),'buttondownfcn','')
	   end

	   if (undoid > 1)
		  set(H.icon(3),'CData',get(H.icon(3),'UserData'),'buttondownfcn',{@undo_feature,H,-1});
	   else
		  set(H.icon(3),'CData',uint8(228-double(get(H.icon(3),'UserData'))),'buttondownfcn','')
	   end

	end

%-------------------------------------------------------------------------------
function resize(obj,event,H)
% Resize the figure according to the new dimension request.
	nRows = round(str2double(get(H.editNrow,'String')));
	nCols = round(str2double(get(H.editNcol,'String')));
	if (isnan(nRows) || nRows < 1)
		set(H.editNrow,'String', H.y),		return
	end
	if (isnan(nCols) || nCols < 1)
		set(H.editNrow,'String', H.x),		return
	end
	if (nRows == H.y && nCols == H.x),		return,		end

	% Make sure they are odd numbers
	if ((fix(nRows / 2) - nRows / 2) == 0),		nRows = nRows + 1;	end
	if ((fix(nCols / 2) - nCols / 2) == 0),		nCols = nCols + 1;	end
	set(H.editNrow,'String',nRows)
	set(H.editNrow,'String',nCols)

	menor = min(nRows,nCols);
	if (menor <= 3),		p = 32;
	elseif (menor <= 5)		p = 28;
	elseif (menor <= 7)		p = 24;
	elseif (menor <= 9)		p = 20;
	elseif (menor <= 11)	p = 16;
	else					p = 14;
	end
	structuring_elem(nCols, nRows, p)

%-------------------------------------------------------------------------------
function iconclick(obj,cnc,H)
	set(H.iconaxes,'visible','off')
	set(get(obj,'parent'),'visible','on')

%-------------------------------------------------------------------------------
function img = icondata
%
% PENCIL
img(:,:,1) = ...
 [ 114 114 114 114 114 114 114 114 114 114 114 114 114 114 114 114
   114 114 114 114 114 114 114 114 114 114 118 165 147 144 152 120
   114 114 114 114 114 114 114 114 114 114 168 165 225 227 149 125
   114 114 114 114 114 114 114 114 114 186 212 225 227 227 219 144
   114 114 114 114 114 114 114 114 186 237 226 174 207 207 219 146
   114 114 114 114 114 114 114 186 237 244 217 207 157 171 145 114
   114 114 114 114 114 114 186 237 244 220 214 238 180 143 114 114
   114 114 114 114 114 186 237 244 220 214 238 180 143 114 114 114
   114 114 114 114 144 237 244 220 214 238 180 143 114 114 114 114
   114 114 114 137 230 244 220 214 238 180 143 114 114 114 114 114
   114 113 137 203 169 220 214 238 180 143 114 114 114 114 114 114
   114 109 225 231 125 153 196 156 131 114 114 114 114 114 114 114
   114 109 227 170 152 103  96 118 114 114 114 114 114 114 114 114
   114  66 156 207 207 139 118 114 114 114 114 114 114 114 114 114
   114  64  66 105 118 131 114 114 114 114 114 114 114 114 114 114
   114 114 114 114 114 114 114 114 114 114 114 114 114 114 114 114 ];
%
% CLEAR
img(:,:,2) = ...
 [ 114 114 114 114 114 114 114 114 110  97  86 107 114 114 114 114
   114 114 114 114 114 114 114  95  81  74  86  86  80 109 114 114
   114 114 114 114 114 114  98  80  86 102 100 128 129  92 109 114
   114 114 114 114 114 103  93  90  91  89 100 109 151 148  82 114
   114 114 114 114 121 144  95  95  91 102 111 120 129 182 125 107
   114 114 114 129 194 240 170  95 104 113 123 132 140 196 148 108
   114 114 138 208 240 234 230 171 116 125 134 143 165 206 143 114
   114 124 206 240 235 231 233 235 191 135 144 168 209 184 120 114
   114 188 235 242 231 233 236 238 241 211 170 211 195 137 114 114
   114 209 243 232 234 237 239 242 243 246 242 208 154 114 114 114
   114 206 242 237 237 240 242 244 247 252 232 165 114 114 114 114
   114 161 232 247 240 242 244 248 252 236 168 114 114 114 114 114
   114 118 215 244 249 247 250 252 239 171 114 114 114 114 114 114
   114 114 123 217 237 244 244 240 171 114 114 114 114 114 114 114
   114 114 114 116 152 188 172 132 114 114 114 114 114 114 114 114
   114 114 114 114 114 114 114 114 114 114 114 114 114 114 114 114];
%
% UNDO
img(:,:,3) = ...
 [ 114 114 114 114 114 114 114 114 114 114 114 114 114 114 114 114
   114 114 114 114 114 114 114 114 114 114 114 114 114 114 114 114
   114 114 114 114 114 114 114 114 114 114 114 114 114 114 114 114
   111 103 114 114 114 114 104 112 112 105 104 101 101 114 114 114
   131 114 103 114 103 117 142 204 204 208 189 149  97  98 114 114
   131 208 127 107 122 207 237 223 219 215 210 202 162  93  99 114
   131 255 214 123 219 238 229 230 236 222 210 198 187 134  90 114
   131 255 254 250 242 234 239 235 129 108 104 114 114 158  94  97
   131 255 249 246 240 244 228 114 103 114 114 114  99 127 115  93
   131 255 249 245 245 252 115 104 114 114 114 114 100 124 137  90
   131 255 251 246 242 249 226 113 103 114 114 114 114  96 166  88
   131 255 254 252 250 248 255 219 110 114 114 114 114  97 162  87
   132 160 156 150 145 140 134 114 107 104 114 114 114  99  91  95
   114 103 103 114 114 114 114 114 114 114 114 114 114 114 114 114
   114 114 114 114 114 114 114 114 114 114 114 114 114 114 114 114
   114 114 114 114 114 114 114 114 114 114 114 114 114 114 114 114 ];

% REDO
img(:,:,4) = fliplr(img(:,:,3));
