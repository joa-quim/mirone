function varargout = show_palette(varargin)
% Plot a color palette either as an idependent figure or inside the main window

%	Copyright (c) 2004-2015 by J. Luis
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

% $Id: show_palette.m 4678 2015-02-04 13:52:12Z j $

    handMir = varargin{1};
    tipo = varargin{2};
    
    % Some tests
	if (handMir.no_file),		return,		end
	if (~handMir.validGrid && ndims(get(handMir.hImg,'CData')) == 3)
        msgbox('True color images do not use color palettes. Bye','Warning');    return
	end
    
    cmap = get(handMir.figure1,'Colormap');
	axUnits = get(handMir.axes1, 'Units');			set(handMir.axes1, 'Units', 'pixels');
	posAxParent = get(handMir.axes1,'Pos');
	posFigParent = get(handMir.figure1,'Pos');
	screen = get(0,'ScreenSize');				% In some cases we need this
	
	min_max = handMir.head(5:6);
	if (min_max(1) == 0 && min_max(2) == 0)		% Happens for example with the stacks
		Z = getappdata(handMir.figure1,'dem_z');
		if (isempty(Z))
			errordlg('Something wrong. I have no info to create a color palette.','WarnError'),		return
		end
		min_max = grdutils(Z,'-L')';
	end
    
    if (strcmp(tipo,'Float'))            % Create a new Figure with the colorbar
        
        figPos = [100 100 100 min(420,posAxParent(4)+10)];
		hObject = figure('Units','Pixels',...
			'Position',figPos, 'PaperUnits','centimeters',...
			'Color',get(0,'factoryUicontrolBackgroundColor'),...
			'KeyPressFcn',@figure1_KeyPressFcn,...
			'MenuBar','none','Resize','off',...
			'Name','','NumberTitle','off',...
            'Colormap',cmap,'Visible','off');
    
		% Try to position this figure glued to the right of calling figure
		ecran = get(0,'ScreenSize');
		xLL = posFigParent(1) + posFigParent(3) + 6;
		xLR = xLL + figPos(3);
		if (xLR > ecran(3))         % If figure is partially out, bring totally into screen
			xLL = ecran(3) - figPos(3);
		end
		yLL = (posFigParent(2) + posFigParent(4)/2+12) - figPos(4) / 2;
		set(hObject,'Pos',[xLL yLL figPos(3:4)])

		% Create the color bar axes and image object
        axPos = [10 5 40 figPos(4)-10];
		hAx = axes('Units','Pixels','pos',axPos,'Parent',hObject);
		image([1 10], min_max, (1:size(cmap,1))');
		set(hAx,'XTick',[],'YAxisLocation','right','YDir','normal')
		set(hObject,'Visible','on');
		if (nargout),   varargout{1} = hObject;     end

    elseif (strcmp(tipo,'At'))				% Create a colorbar right next to image's right side

		origFigWidth = posFigParent(3);		% To store original figure's width before palette insertion
		hUict = handMir.PalAt;
		% If we already have a colorbar, remove it
		if (strcmp(get(hUict,'Check'),'on'))
			ud = get(hUict,'Userdata');
			posFigParent(3) = ud(2);		% Figure's width before palette insertion 
			if (posFigParent(3) == screen(3)),      posFigParent(3) = posFigParent(3)-1;    end     % The elastic
			set(handMir.figure1,'pos',posFigParent)
			delete(ud(1))					% Delete the color palette
			set(hUict,'Checked','off')
			return
		end

		% Create the color bar axes at the right side of the main image axes
        barW = 20;			% Colorbar width
        marg = 8;			% Margin between Image and colorbar (must be enoug for slider width)
        axPos = [posAxParent(1)+posAxParent(3)+marg posAxParent(2) barW min(420,posAxParent(4))];
		hAx = axes('Units','Pixels','pos',axPos,'Parent',handMir.figure1,'vis','off');
		image([1 10], min_max, (1:size(cmap,1))');
		set(hAx,'XTick',[],'YAxisLocation','right','YDir','normal','HandleVisibility','off','Tag','MIR_CBat')

		h_Ylabel = get(hAx,'Ylabel');	set(h_Ylabel,'units','pixels')
		Ylabel_pos = get(h_Ylabel,'Extent');

        gutterRight = posFigParent(3) - (posAxParent(1) + posAxParent(3));
		if (gutterRight < Ylabel_pos(1)+marg)			% Grow Figure's width to acomodate the colorbar
			posFigParent(3) = posFigParent(3) + (Ylabel_pos(1)+marg - gutterRight - 3); % -3 is empiric
			if (posFigParent(3) == screen(3)),      posFigParent(3) = posFigParent(3)-1;    end     % The elastic
			set(handMir.figure1, 'Position', posFigParent);
		end
		set(hUict,'Userdata',[hAx origFigWidth])		% Save it so that we can restore upon colorbar deletion
		set(hAx,'Units','normalized','Vis','on')
		set(hUict,'Checked','on')

    elseif (strcmp(tipo,'In'))					% Create a colorbar inside the image, at it's right side

		hUict = handMir.PalIn;
		% If we already have a colorbar, remove it
		if (strcmp(get(hUict,'Check'),'on'))
			ud = get(hUict,'Userdata');
			delete(ud)
			set(hUict,'Checked','off')
			return
		end

		% Create the color bar axes at the right (in)side of the main image axes
        barW = 20;      % Colorbar width
        axPos = [posAxParent(1)+posAxParent(3)-barW posAxParent(2) barW posAxParent(4)];
		hAx = axes('Units','Pixels','pos',axPos,'Parent',handMir.figure1,'Vis','off');
		image([1 10], min_max, (1:size(cmap,1))');
		set(hAx,'XTick',[],'YDir','normal','HandleVisibility','off','Units','normalized','Vis','on','Tag','MIR_CBin')
		set(hUict,'Userdata',hAx)            % Save it so that we can restore upon colorbar deletion
		set(hUict,'Checked','on')

    end
    
	set(handMir.axes1, 'Units', axUnits);   % Do it here because the "At side" case needs it in Pixels

% --- Executes on key press over figure1 with no controls selected.
function figure1_KeyPressFcn(hObject, eventdata)
	if isequal(get(hObject,'CurrentKey'),'escape')      % Check for "escape"
        delete(hObject);
	end
