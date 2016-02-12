function varargout = uisetdate(varargin)
% uisetdate is designed to select any date among the past current and future years.
% uisetdate by itself uses the current date and year as a starting point and returns a string
%           containing the selected date in 'dd-mmm-yyyy' format.
% uisetdate(date) uses date as a starting point. date must be in 'dd-mmm-yyyy' format.
%
% [datet,Y,M,D,nD]=uisetdate returns a string containing the date plus the year, the month, the day
%                            and the number of days between the selected date and the 1st January.
%
%   example:
%      if you select the 5th of August of year 2004, the returned values would be
%
%     '05-Aug-2004'
%              2004
%                 8
%                 5
%               218  (218 days between 5th of August and 1st January)
%
% To change the year, just type + or - key while pointer is above the figure (unit step change) or
% type y to select a given year. If you close the figure, all the outputs will be empty. Figure appearance
% may be changed in the "init" function. 
%
%  Luc Masset (2004)  e-mail: luc.masset@ulg.ac.be

% Rewritten to use handles and be more robust.
%	Coffeeright (c) 2004-2016 by J. Luis
%
%	Contact info: w3.ualg.pt/~jluis/mirone
% --------------------------------------------------------------------

% $Id: uisetdate.m 7790 2016-02-12 01:33:03Z j $

	[datet,Y,M,D,nD] = init;
	varargout{1} = datet;
	varargout{2} = Y;   varargout{3} = M;
	varargout{4} = D;   varargout{5} = nD;

%------------------------------------------------------------------------------
function [datet,Y,M,D,nD] = init(datet)
	if (nargin == 0),	datet = date;     end

	listJ = {'Su','Mo','Tu','We','Th','Fr','Sa'};
	listM = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};    % month list

	%year, month and day
	Y = str2double(datet(end-3:end));		D = str2double(datet(1:2));
	M = datet(end-7:end-5);					M = strcmp(M,listM);
	M = find(M);

	hFig = figure('units','pixels','position',[0 0 290 330],'menubar','none','numbertitle','off', ...
				'name','Calendar','resize','off','KeyPressFcn',@changeyear, ...
				'Color',get(0,'factoryUicontrolBackgroundColor'),'Tag','figure1', ...
				'Visible','off','DefaultUIControlFontName','arial','DefaultUIControlFontSize',10);
	move2side(hFig,'center');

	%frame buttons
	uicontrol('style','frame','units','pixels','position',[2 2 286 215]);
	uicontrol('style','frame','units','pixels','position',[2 2 286 185]);
	uicontrol('style','frame','units','pixels','position',[2 222 286 66]);
	uicontrol('style','frame','units','pixels','position',[2 292 286 36]);

	%current date button
	tts = 'Use +/- keys to change year by unit step. Use y key to set the year';
	uicontrol('style','text','units','pixels','position',[10 298 270 20],'string',datet, ...
				'horizontalalignment','center','fontsize',12,'Tag','date', 'tooltipstring',tts);

	%validate button
	uicontrol('style','pushbutton','units','pixels','position',[245 300 30 20], ...
				'string','OK','tooltipstring','Validate current date', ...
				'Tag', 'OK', ...
				'Call',@uisetdate_uiCB)

	%static text buttons for day name
	for (i = 1:7)
		pos=[10+40*(i-1) 190 30 20];
		uicontrol('style','text','units','pixels','position',pos,'string',listJ{i}, 'horizontalalignment','center');
	end

	%------------ Give a Pro look (3D) to the frame boxes  --------
	new_frame3D(hFig, NaN)
	set(hFig,'Visible','on')

	%figure appdata
	setappdata(hFig,'year',Y);  setappdata(hFig,'month',M);
	setappdata(hFig,'day',D);   setappdata(hFig,'SelectColor','b')

	handles = guihandles(hFig);
	guidata(hFig, handles);

	update(hFig, handles)		% update buttons and text

	uiwait(hFig);

	%wait for temp button to be deleted
	if (~ishandle(hFig))
		datet = [];   Y = [];   M = [];   D = [];   nD = [];  return
	end

	%compute outputs
	Y = getappdata(hFig,'year');    M = getappdata(hFig,'month');   D = getappdata(hFig,'day');
	datet = datestr([Y M D 0 0 0],'dd-mmm-yyyy');
	nD = 0;  %indice du jour
	for (i = 1:M-1),	nD = nD + eomday(Y,i);  end
	nD = nD+D;
	delete(hFig)

%------------------------------------------------------------------------------
function OK_CB(hObject, handles)
% ...
	uiresume(handles.figure1)

%------------------------------------------------------------------------------
function update(hObject, handles)
	% Update buttons and text when changing year, month or day

	hFig = handles.figure1;

	%delete old buttons
	delete(findobj('tag','changeday','type','uicontrol','parent',hFig))

	%year, month, day
	Y = getappdata(hFig,'year');       M = getappdata(hFig,'month');
	D = getappdata(hFig,'day');        Dmax = eomday(Y,M);
	D = min([D Dmax]);
	setappdata(hFig,'day',D)

	%current month calendar
	C = calendar(Y,M);

 	listM = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
	if (isempty(findobj('tag','changemonth','type','uicontrol','parent',hFig)))		% First time. Create the month buttons
		for (i = 1:2)
			for (j = 1:6)
				pos = [15+45*(j-1) 260-(i-1)*30 35 20];
				st  = listM{6*(i-1)+j};
				uicontrol('style','togglebutton','units','pixels','position',pos,'string',st, ...
						'Tag','changemonth', 'Call',@uisetdate_uiCB)
			end
		end
	end

	for (i = 1:size(C,1))		% day buttons
		for (j = 1:7)
			if C(i,j)
				pos = [10+40*(j-1) 160-(i-1)*30 30 20];
				st  = sprintf('%d', C(i,j));
				uicontrol('style','togglebutton','units','pixels','position',pos,'string',st, ...
					'Tag','changeday', 'Call',@uisetdate_uiCB)
			end
		end
	end

	%selected month
	scolor = getappdata(hFig,'SelectColor');
	set(findobj('tag','changemonth','type','uicontrol','parent',hFig),'value',0,'foregroundcolor','k')
	h = findobj('tag','changemonth','string',listM{M},'type','uicontrol','parent',hFig);
	set(h,'value',1,'foregroundcolor',scolor)

	%selected day
	h = findobj('tag','changeday','string',num2str(D),'type','uicontrol','parent',hFig);
	set(h,'value',1,'foregroundcolor',scolor)

	%update current date text
	h = findobj('tag','date','type','uicontrol','parent',hFig);
	set(h,'string',datestr([Y M D 0 0 0],'dd-mmm-yyyy'))

%------------------------------------------------------------------------------
function changeday_CB(hObject, handles)
	D = str2num(get(hObject,'string'));
	setappdata(handles.figure1,'day',D)
	update(hObject, handles)

%------------------------------------------------------------------------------
function changemonth_CB(hObject, handles)
	listM = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
	M = get(hObject, 'string');
	M = find(strcmp(M,listM));
	setappdata(handles.figure1,'month',M)
	update(hObject, handles)

%------------------------------------------------------------------------------
function changeyear(hObject, evt)
	handles = guidata(hObject);
	Y = getappdata(handles.figure1,'year');
	cc = get(handles.figure1,'currentcharacter');
	switch cc,
		case '+',   Y=Y+1;
		case '-',   Y=Y-1;
		case 'y',
			def = {sprintf('%i',Y)};
			answer = inputdlg({'Year:'},'Set current year',1,def);
			if isempty(answer),		return,		end
			Y = str2num(answer{1});
			if isempty(Y),	return,		end
			Y = round(Y);
		otherwise
			return
	end
	setappdata(handles.figure1,'year',Y)
	update(hObject, handles)

%------------------------------------------------------------------------------
function uisetdate_uiCB(hObject, eventdata)
% This function is executed by the callback and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'], hObject, guidata(hObject));
