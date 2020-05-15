function varargout = masker(varargin)
% Helper Window to run masking ops using second grid/images

%	Copyright (c) 2004-2020 by J. Luis
%
%             DO WHAT THE FUCK YOU WANT TO PUBLIC LICENSE
%                     Version 2, December 2004
% 
%  Everyone is permitted to copy and distribute verbatim or modified
%  copies of this license document, and changing it is allowed as long
%  as the name is changed.
% 
%             DO WHAT THE FUCK YOU WANT TO PUBLIC LICENSE
%    TERMS AND CONDITIONS FOR COPYING, DISTRIBUTION AND MODIFICATION
% 
%   0. You just DO WHAT THE FUCK YOU WANT TO.
%
%	Contact info: w3.ualg.pt/~jluis/mirone
% --------------------------------------------------------------------

	if (nargin > 1 && ischar(varargin{1}))
		gui_CB = str2func(varargin{1});
		[varargout{1:nargout}] = feval(gui_CB,varargin{2:end});
	else
		h = masker_OF(varargin{:});
		if (nargout),	varargout{1} = h;   end
	end

% ---------------------------------------------------------------------------------
function hObject = masker_OF(varargin)
	if (isempty(varargin)),		hObject = [];	return,		end
	if (varargin{1}.no_file),	hObject = [];	return,		end
 
	hObject = figure('Vis','off');
	masker_LayoutFcn(hObject);
	handles = guihandles(hObject);
	move2side(hObject,'right')

	handMir = varargin{1};
	handles.hMirFig = handMir.figure1;
	handles.hMirImg = handMir.hImg;
	handles.last_dir = handMir.last_dir;
	handles.mask = [];

	[masks, handles.hMaskFigs] = fish_MirMasks(handles);
	if (~isempty(masks))
		set(handles.listbox, 'String', masks)
	else								% Resize fig since no floating masks are around
		pos = get(handles.listbox, 'Pos');		H = pos(2);
		pos = get(handles.figure1, 'Pos');		pos(4) = H - 1;
		set(handles.figure1, 'Pos', pos)
	end

	if (handMir.validGrid)
		set(handles.check_alpha, 'Vis', 'off')
	else
		set(handles.edit_value, 'String', 0)
		if (size(get(handMir.hImg, 'CData'),3) == 1)
			set(handles.edit_value, 'Enable', 'off')	% Indexed imges can only use the bg color 
		else
			s = sprintf('Enter a scalar in the [0 255] range\nor a color like 100/20/250\nThe scalar 0 means use bg color');
			set(handles.edit_value, 'Tooltip', s)
		end
	end
	guidata(hObject, handles);
	set(hObject,'Visible','on');
	if (nargin > 1),	external_drive(handles, 'masker', varargin{2:end}),	end

% -----------------------------------------------------------------------------------------
function listbox_CB(hObject, handles)
% ...
	val = get(hObject, 'Val');
	if (val == 1),	return,		end
	thisHand = guidata(handles.hMaskFigs(val-1));
	if (thisHand.validGrid)
		handles.mask = logical(getappdata(thisHand.figure1,'dem_z'));
	else
		handles.mask = logical(get(thisHand.hImg,'CData'));
	end
	guidata(handles.figure1, handles)

% -----------------------------------------------------------------------------------------
function edit_fname_CB(hObject, handles)
	fname = get(hObject,'String');
	if (isempty(fname)),	return,		end
	push_fname_CB(handles.push_fname, handles, fname)

% -----------------------------------------------------------------------------------------
function push_fname_CB(hObject, handles, opt)
	if (nargin == 3)
		[img, handles_out, handles.att] = mirone('FileOpenNewImage_CB', handles, opt);
	else
		[img, handles_out, handles.att] = mirone('FileOpenNewImage_CB', handles);
		set(handles.edit_fname,'Str',handles_out.fileName)		
	end
	handles.mask = logical(img);
	handles.head = handles_out.head;
	guidata(handles.figure1, handles)

% -----------------------------------------------------------------------------------------
function push_OK_CB(hObject, handles)
% ...
	if (isempty(handles.mask)),		errordlg('Yes, yes. And what about giving me a mask?','Error'),		end
	if (get(handles.check_revert, 'Val'))
		mask = ~handles.mask;
	else
		mask = handles.mask;
	end
	handMir = guidata(handles.hMirFig);
	if (handMir.validGrid)
		val = single(str2double(get(handles.edit_value, 'String')));
		Z = getappdata(handles.hMirFig,'dem_z');
		if (size(Z,1) ~= size(mask,1) || size(Z,2) ~= size(mask,2))
			errordlg('Maks and target do not have the same size. Impossible to continue', 'Error')
			return
		end
		Z(mask) = val;
		tmp.head = handMir.head;
		if (val < tmp.head(5)),	tmp.head(5) = val;	end
		if (val > tmp.head(6)),	tmp.head(6) = val;	end
		tmp.geog = handMir.geog;
		name = get(handles.hMirFig, 'Name');
 		ind = strfind(name, '@');
		if (~isempty(ind)),	name = ddewhite(name(1:ind(end)-1));	end
		[pato, fname] = fileparts(name);
		tmp.name = [fname 'Masked'];
		mirone(Z, tmp, handles.hMirFig)
	else
		img = get(handMir.hImg, 'CData');
		if (size(img,3) == 1)
			img(mask) = 0;			% 0 means bg color
		else
			val_str = ddewhite(get(handles.edit_value, 'String'));
			ind = strfind(val_str, '/');
			if (numel(ind) == 2 && (ind(1) ~= 1 && ind(2) < numel(val_str)))	% Accept entries like 12/123/234
				val = [str2double(val_str(1:ind(1)-1)) str2double(val_str(ind(1)+1:ind(2)-1)) str2double(val_str(ind(2)+1:end))];
				img = maskRGB(img, mask, uint8(val));
			else
				if (strcmp(val_str, '0'))		% 0 means we'll use the bg_color
					img = maskRGB(img, mask, uint8(handMir.bg_color(1,:) * 255));
				else
					img = maskRGB(img, mask, uint8(str2double(val_str)));
				end
			end
		end
		if (get(handles.check_alpha, 'Val'))
			new_mask = alloc_mex(size(mask,1), size(mask,2), 'uint8', 255);
			new_mask(mask) = 0;
			set(handles.hMirImg, 'AlphaData', new_mask)
		else
			set(handles.hMirImg, 'AlphaData', 1)
		end
		set(handMir.hImg, 'CData', img)
	end

% -----------------------------------------------------------------------------------------
function img = maskRGB(img, mask, val)
% Mask the IMG RGB array at positions where MASK == 1 to VAL
	if (numel(val) == 3)
		t = img(:,:,1);		t(mask) = val(1);	img(:,:,1) = t;
		t = img(:,:,2);		t(mask) = val(2);	img(:,:,2) = t;
		t = img(:,:,3);		t(mask) = val(3);	img(:,:,3) = t;
	else
		t = img(:,:,1);		t(mask) = val;		img(:,:,1) = t;
		t = img(:,:,2);		t(mask) = val;		img(:,:,2) = t;
		t = img(:,:,3);		t(mask) = val;		img(:,:,3) = t;
	end

% -----------------------------------------------------------------------------------------
function [masks, hFigs] = fish_MirMasks(handles)
% Fish Mirone figs with masks and return their names, or empty if not found
	masks = [];
	hFigs = findobj(0,'type','figure');						% Fish all figures
	hFigs = aux_funs('figs_XOR', handles.hMirFig, hFigs);	% Get all unique Mirone Figs
	if (isempty(hFigs)),	return,		end
	masks = false(1, numel(hFigs));
	[H,W,layers] = size(get(handles.hMirImg, 'CData'));
	for (k = 1:numel(hFigs))
		thisHand = guidata(hFigs(k));
		this_im = get(thisHand.hImg,'CData');
		[H_,W_,layer] = size(this_im);
		if (layer == 1 && W == W_ && H == H_),	masks(k) = true;	end
	end
	hFigs = hFigs(masks);
	if (isempty(hFigs)),	masks = [];	return,		end		% None of them is a mask with the same size as IMG
	nomes = get(hFigs,'name');
	if (~isa(nomes,'cell')),	nomes = {nomes};	end
	for (k = 1:numel(hFigs))
		ind = strfind(nomes{k}, '@');
		if (~isempty(ind)),	nomes{k} = ddewhite(nomes{k}(1:ind(end)-1));	end
		[pato, nomes{k}] = fileparts(nomes{k});
	end
	masks = [{' '}; nomes];

% ---------------------------------------------------
function masker_LayoutFcn(h1)

set(h1, 'Position',[520 605 321 211],...
	'Color',get(0,'factoryUicontrolBackgroundColor'),...
	'MenuBar','none',...
	'Name','Masker',...
	'NumberTitle','off',...
	'Resize','off',...
	'HandleVisibility','callback',...
	'Tag','figure1');

uicontrol('Parent',h1, 'Position',[11 183 190 15],...
	'HorizontalAlignment','left',...
	'String','Use this mask file (click to select)',...
	'Style','text');

uicontrol('Parent',h1, 'Position',[10 131 301 51],...
	'Callback',@masker_uiCB,...
	'String','',...
	'Style','listbox',...
	'Value',1,...
	'Tag','listbox');

uicontrol('Parent',h1, 'Position',[10 81 281 21],...
	'BackgroundColor',[1 1 1],...
	'Callback',@masker_uiCB,...
	'HorizontalAlignment','left',...
	'String','',...
	'Style','edit',...
	'Tag','edit_fname');

uicontrol('Parent',h1, 'Position',[10 102 100 15],...
	'HorizontalAlignment','left',...
	'String','Load mask from file',...
	'Style','text',...
	'Tag','text_mask_file');

uicontrol('Parent',h1, 'Position',[10 53 140 15],...
	'HorizontalAlignment','left',...
	'String','Replace masked values by:',...
	'Style','text',...
	'Tag','text_mask_val');

uicontrol('Parent',h1, 'Position',[287 81 23 21],...
	'Callback',@masker_uiCB,...
	'String','...',...
	'Tag','push_fname');

uicontrol('Parent',h1, 'Position',[147 50 90 21],...
	'BackgroundColor',[1 1 1],...
	'String','NaN',...
	'Style','edit',...
	'Tooltip','All masked values in grid will be replaced by this',...
	'Tag','edit_value');

uicontrol('Parent',h1, 'Position',[10 28 78 15],...
	'String','Revert mask',...
	'Style','checkbox',...
	'Tooltip','Negate mask. That is, swapp zeros and ones',...
	'Tag','check_revert');

uicontrol('Parent',h1, 'Position',[10 9 100 15],...
	'String','Add alpha layer',...
	'Style','checkbox',...
	'Tooltip','Add mask to image alpha layer',...
	'Tag','check_alpha');

uicontrol('Parent',h1, 'Position',[250 10 60 23],...
	'Callback',@masker_uiCB,...
	'FontSize',10,...
	'FontWeight','bold',...
	'String','OK',...
	'Tag','push_OK');

function masker_uiCB(hObject, evt)
% This function is executed by the callback and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'],hObject, guidata(hObject));
