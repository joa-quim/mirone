function fig_hdl = message_win(option,in2)
%MESSAGE_WIN  Message function.
%
%   Usage:  message_win('create','Some text')   for the first call
%           message_win('add','More text')      for subsequent calls. When the window message
%                                               is full of text a slider will be added.
%   NOTE:   This apparently doesn't work well when the "More text" is long enough to fill the
%           edit box in one call.

%   Based on the dmsgfun Matlab function, but with a slider and some other modifs
%   Joaquim Luis

%------------------------
tag_fig = 'Wdmsgfig';
tag_txt = 'Axe_info';

fig_hdl = wfindobj('figure','tag',tag_fig);
switch option
    case 'create'
        if isempty(fig_hdl)
            win_units  = 'pixels';
            win_height = 500;
            Screen_Size = get(0,'ScreenSize');
			defFigPos = get(0,'DefaultfigurePosition');
            win_width = defFigPos(3);
            pos_win = [Screen_Size(3)-5-win_width 40 win_width win_height];
			if Screen_Size(4)<800 , pos_win(2) = 20; end
			axe_col = 'w'; 
            fig_hdl = colordef('new','none');
            set(fig_hdl, 'MenuBar','none', 'Name','Information window',...
                    'Visible','off', 'Unit',win_units, 'Position',pos_win,...
                    'Color',axe_col, 'NumberTitle','off', 'Tag',tag_fig);
			bord  = 10;
            p_text = [bord bord/5 win_width-2*bord win_height-2*bord/5];
            txt_hdl = uicontrol('Parent',fig_hdl,...
                    'Style','text',      ...
                    'Visible','off',     ...
                    'Units',win_units,   ...
                    'Position',p_text,   ...
					'FontWeight','bold', ...
					'FontSize',8,       ...
					'Max',40,            ...
					'HorizontalAlignment','left',...
                    'BackgroundColor',[1 1 1],   ...
					'ForegroundColor',[0 0 0],   ...
                    'Tag',tag_txt);
            set([fig_hdl,txt_hdl],'units','normalized','Visible','on');
            set(fig_hdl,'Name','Message window');
        else
            set(fig_hdl, 'HandleVisibility', 'on');
            txt_hdl = findobj(fig_hdl,'Tag',tag_txt);
            figure(fig_hdl);
        end
		set(txt_hdl,'String',in2);
        set(fig_hdl, 'HandleVisibility', 'off');
    case 'add'
        txt_hdl = findobj(fig_hdl,'Style','Text');
        h_slide = findobj(fig_hdl,'Style','slider');
        old_text = get(txt_hdl,'String');
        post = get(txt_hdl,'Position');
        extent = get(txt_hdl,'Extent');
        scal = 5;
        if (extent(4) > 0.999 && isempty(h_slide))   % Text to big to fit. Add a slider to the figure
            pos=[0.97 0 .03 1];
            cb_slide_step = {@slide_step,txt_hdl,post,scal};
            uicontrol(fig_hdl,'style','slider','units','normalized','position',pos,...
                'callback',cb_slide_step,'min',0,'max',scal,'Value',scal);
        elseif (extent(4) > 1)          % Text is biger than figure. Scroll it up and move slider down
            slid_val = scal - (extent(4) - 1);
            set(h_slide,'Value',slid_val)
            new_text = char(old_text,in2);
            new_pos = [post(1) post(2) post(3) scal - slid_val + 1];
            set(txt_hdl,'String',new_text,'Position',new_pos)
        else
            new_text = char(old_text,in2);
            set(txt_hdl,'String',new_text);
        end
    case 'close'
        delete(fig_hdl);
end

function slide_step(obj,eventdata,h,pos,scal)
new_pos = [pos(1) pos(2) pos(3) scal-get(gcbo,'value')+1];
set(h,'Position',new_pos)
