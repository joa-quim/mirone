function [leghandle,labelhandles,outH,outM]=legend_(varargin)
%LEGEND Graph legend.
%   LEGEND(string1,string2,string3, ...) puts a legend_ on the current plot
%   using the specified strings as labels. LEGEND works on line graphs,
%   bar graphs, pie graphs, ribbon plots, etc.  You can label any
%   solid-colored patch or surface object.  The fontsize and fontname for
%   the legend_ strings matches the axes fontsize and fontname.
%
%   LEGEND(H,string1,string2,string3, ...) puts a legend_ on the plot
%   containing the handles in the vector H using the specified strings as
%   labels for the corresponding handles.
%
%   LEGEND(M), where M is a string matrix or cell array of strings, and
%   LEGEND(H,M) where H is a vector of handles to lines and patches also
%   works.
%
%   LEGEND(AX,...) puts a legend_ on the axes with handle AX.
%
%   LEGEND OFF removes the legend_ from the current axes.
%   LEGEND(AX,'off') removes the legend_ from the axis AX.
%
%   LEGEND HIDE makes legend_ invisible.
%   LEGEND(AX,'hide') makes legend_ on axis AX invisible.
%   LEGEND SHOW makes legend_ visible.
%   LEGEND(AX,'show') makes legend_ on axis AX visible.
%
%   LEGEND BOXOFF sets appdata property legendboxon to 'off' making
%   legend_ background box invisible when the legend_ is visible.
%   LEGEND(AX,'boxoff') sets appdata property legendboxon to 'off for axis AX
%   making legend_ background box invisible when the legend_ is visible.
%   LEGEND BOXON sets appdata property legendboxon to 'on' making
%   legend_ background box visible when the legend_ is visible.
%   LEGEND(AX,'boxon') sets appdata property legendboxon to 'off for axis AX
%   making legend_ background box visible when the legend_ is visible.
%
%   LEGH = LEGEND returns the handle to legend_ on the current axes or
%   empty if none exists.
%
%   LEGEND with no arguments refreshes all the legends in the current
%   figure (if any).  LEGEND(LEGH) refreshes the specified legend_.
%
%   LEGEND(...,Pos) places the legend_ in the specified
%   location:
%       0 = Automatic "best" placement (least conflict with data)
%       1 = Upper right-hand corner (default)
%       2 = Upper left-hand corner
%       3 = Lower left-hand corner
%       4 = Lower right-hand corner
%      -1 = To the right of the plot
%
%   To move the legend_, press the left mouse button on the legend_ and drag
%   to the desired location. Double clicking on a label allows you to edit
%   the label.
%
%   [LEGH,OBJH,OUTH,OUTM] = LEGEND(...) returns a handle LEGH to the
%   legend_ axes; a vector OBJH containing handles for the text, lines,
%   and patches in the legend_; a vector OUTH of handles to the
%   lines and patches in the plot; and a cell array OUTM containing
%   the text in the legend_.
%
%   LEGEND will try to install a ResizeFcn on the figure if it hasn't been
%   defined before.  This resize function will try to keep the legend_ the
%   same size.
%
%   Examples:
%       x = 0:.2:12;
%       plot(x,bessel(1,x),x,bessel(2,x),x,bessel(3,x));
%       legend_('First','Second','Third');
%       legend_('First','Second','Third',-1)
%
%       b = bar(rand(10,5),'stacked'); colormap(summer); hold on
%       x = plot(1:10,5*rand(10,1),'marker','square','markersize',12,...
%                'markeredgecolor','y','markerfacecolor',[.6 0 .6],...
%                'linestyle','-','color','r','linewidth',2); hold off
%       legend_([b,x],'Carrots','Peas','Peppers','Green Beans',...
%                 'Cucumbers','Eggplant')
%
%   See also PLOT.

%   D. Thomas 5/6/93
%             9/6/95
%   Rich Radke 6/17/96
%   Copyright 1984-2002 The MathWorks, Inc.
%   $Revision: 5.90 $  $Date: 2002/05/20 20:35:57 $

%   Private syntax:
%
%     LEGEND('ShowLegendPlot',h) is called from MOVEAXIS to set the gco to
%        the plot the legend goes with.

narg = nargin;
isresize(0);

%--------------------------
% Parse inputs
%--------------------------

% Determine the legend_ parent axes (ha) is specified
if narg > 0 && ~isempty(varargin{1}) && ~isa(varargin{1}, 'char') &&ishandle(varargin{1}) && strcmp(get(varargin{1}(1),'type'),'axes') % legend_(ax,...)
	ha = varargin{1}(1);
	varargin(1) = []; % Remove from list
	narg = narg - 1;
	if islegend(ha) % Use the parent
		ud = get(ha,'userdata');
		if isfield(ud,'PlotHandle')
			ha = ud.PlotHandle;
		else
			warning('Can''t put a legend on a legend.')
			if nargout>0, leghandle=[]; labelhandles=[]; outH=[]; outM=[]; end
			return
		end
	end
else
	ha = [];
end

if (numel(varargin) == 1 && isa(varargin{1}, 'char'))
	varargin = {varargin};				% Just to simplify with the if tests bellow
end

if (narg == 0)					% legend_
	if nargout==1 % h = legend_
		if isempty(ha)
			leghandle = find_legend(find_gca);
		else
			leghandle = find_legend(ha);
		end
	elseif nargout==0 % legend_
		if isempty(ha)
			update_all_legends
		else
			update_legend(find_legend(ha));
		end
	else % [h,objh,...] = legend_
		if isempty(ha)
			[leghandle,labelhandles,outH,outM] = find_legend(find_gca);
		else
			[leghandle,labelhandles,outH,outM] = find_legend(ha);
		end
		if (nargout>3 && ischar(outM)), outM = cellstr(outM); end
	end
	return
	
elseif narg==1 && strcmp(varargin{1}(1),'off')		% legend_('off') or legend_(AX,'off')
	if isempty(ha)
		delete_legend(find_legend(find_gca))
	else
		delete_legend(find_legend(ha))
	end
	if nargout>0, error('Too many outputs.'); end
	return
	
	% Legend hide or legend_(ax,'hide')
elseif narg==1 && strcmp(varargin{1}(1),'hide')
	if isempty(ha)
		leg = find_legend(find_gca);
	else
		leg = find_legend(ha);
	end
	
	% if legend_ axes are already invisible but
	% some of it's children are visible
	% then make children invisible
	if strcmp(get(leg,'visible'),'off')
		legch = get(leg,'children');
		if any(strcmp(get(legch,'visible'),'on'))
			set(legch,'visible','off');
		end
	else
		set(leg,'visible','off')
	end
	
	return
	
	% Legend show or legend_(ax,'show')
elseif narg==1 && strcmp(varargin{1}(1),'show')
	if isempty(ha)
		leg = find_legend(find_gca);
	else
		leg = find_legend(ha);
	end
	
	legch = get(leg,'children');
	
	% get legendboxon from appdata if there, otherwise it's on.
	if isappdata(leg,'legendboxon')
		legboxon = getappdata(leg,'legendboxon');
	else
		legboxon = 'on';
	end
	
	% set legend_ axes visibility
	if strcmp(legboxon,'on')
		set(leg,'visible','on');
	else
		set(leg,'visible','off');
	end
	% make sure children are visible
	set(legch,'visible','on');
	return
	
	% Legend boxoff or legend_(ax,'boxoff')
elseif narg==1 && strcmp(varargin{1}(1),'boxoff')
	if isempty(ha)
		leg = find_legend(find_gca);
	else
		leg = find_legend(ha);
	end
	% set legendboxon appdata to off
	setappdata(leg,'legendboxon','off');
	legch = get(leg,'children');
	
	% check for legend_ hidden
	if strcmp(get(leg,'visible'),'off') && ~any(strcmp(get(legch,'visible'),'on'))
		return
	else
		set(leg,'visible','off');		% make legend_ axis invisible
		set(legch,'visible','on');		% make children visible
	end
	return
	
	% Legend boxon or legend_(ax,'boxon')
elseif narg==1 && strcmp(varargin{1}(1),'boxon')
	if isempty(ha)
		leg = find_legend(find_gca);
	else
		leg = find_legend(ha);
	end
	
	% set legendboxon appdata to off
	setappdata(leg,'legendboxon','on');
	legch = get(leg,'children');
	
	% check for legend_ hidden
	if strcmp(get(leg,'visible'),'off') && ~any(strcmp(get(legch,'visible'),'on'))
		return
	else
		% make legend_ axis visible
		set(leg,'visible','on');
		% make sure children are visible
		set(legch,'visible','on');
	end
	return
	
elseif narg==1 && islegend(varargin{1}(1)) % legend_(legh)
	[hl,labelhandles,outH,outM] = update_legend(varargin{1});
	if nargout>0, leghandle = hl; end
	if (nargout>3 && ischar(outM)), outM = cellstr(outM); end
	return
	
end


% Look for legendpos code
if isa(varargin{end},'double')
	legendpos = varargin{end};
	varargin(end) = [];
else
	legendpos = [];
end

% Determine the active children (kids) and the strings (lstrings)
if narg < 1
	error('Not enough input arguments.');
elseif ishandle(varargin{1}) % legend_(h,strings,...)
	kids = varargin{1};
	if isempty(ha)
		ha=get(varargin{1}(1),'Parent');
		if ~strcmp(get(ha,'type'),'axes')
			error('Handle must be an axes or child of an axes.');
		end
	end
	if narg==1, error('A string must be supplied for each handle.'); end
	lstrings = getstrings(varargin(2:end));
else % legend_(strings,...) or legend_(linetype,string,...)
	if isempty(ha), ha=find_gca; end
	kids = getchildren(ha);
	lstrings = getstrings(varargin);
end

% Set default legendpos
if isempty(legendpos)
	if ~isequal(get(ha,'view'),[0 90])
		legendpos = -1;  % To the right of axis is default for 3-D
	else
		legendpos = 1;   % Upper right is default for 2-D
	end
end

% Remove any existing legend_ on this plot
if isempty(ha)
	hl = find_legend;
else
	hl = find_legend(ha);
end
if ~isempty(hl)
	ud = get(hl,{'userdata'});
	for i=1:length(ud)
		if isfield(ud{i},'PlotHandle') && ud{i}.PlotHandle == ha
			delete_legend(hl)
		end
	end
end

if isempty(kids)
	warning('Plot empty.')
	if nargout>0
		leghandle = []; labelhandles = []; outH = []; outM = [];
	end
	return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% make_legend
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[hl,labelhandles] = make_legend(ha,kids,lstrings,legendpos);

if nargout > 0, leghandle=hl; end
if nargout > 2, outH = kids; end
if nargout > 3
	if ischar(lstrings), lstrings = cellstr(lstrings); end
	outM = lstrings;
end

%--------------------------------
function [hl,hobjs,outH,outM] = find_legend(ha)
%FIND_LEGEND Return current legend_ handle or error out if none.
if nargin==0
	hFig=find_gcf;
else
	hFig=find_gcf(ha);
end

hAx = findobj(hFig,'type','axes');
hl=[];
for i=1:length(hAx)
	if islegend(hAx(i))
		hl(end+1)=hAx(i);
	end
end
hobjs = [];

if nargin>0 && strcmp(get(ha,'type'),'axes')
	if length(ha)~=1
		error('Requires a single axis handle.');
	end
	ud = get(hl,{'userdata'});
	for i=1:length(ud)
		if isfield(ud{i},'PlotHandle') && ud{i}.PlotHandle == ha
			hl = hl(i);
			udi = ud{i};
			hobjs = udi.LabelHandles;
			outH  = udi.handles;
			outM  = udi.lstrings;
			return
		end
	end
	hl = []; % None found
	hobjs = [];
	outH = [];
	outM = [];
end

%-------------------------------
function tf = isresize(setting)
persistent s
if nargin==1
	s = setting;
end
if nargout==1
	tf = s;
end

%--------------------------------
function hf = find_gcf(ha)
%FIND_GCF Find gcf.
%   FIND_GCF Returns the callback figure if there is one otherwise
%   it returns the current figure.
if nargin==1 && strcmp(get(ha,'type'),'axes')
	hf = get(ha,'parent');
else
	if isresize
		hf = gcbf;
		if isempty(hf)
			hf = gcf;
		end
	else
		hf = gcf;
	end
end

%---------------------------------
function ha = find_gca(ha)
%FIND_GCA Find gca (skipping legend_)
if nargin==0
	fig = find_gcf;
else
	fig = find_gcf(ha);
end
ha = get(fig,'currentaxes');
if isempty(ha), ha = gca; end
if islegend(ha)
	ud = get(ha,'userdata');
	if isfield(ud,'PlotHandle')
		ha = ud.PlotHandle;
		% Make sure legend isn't the gca
		set(fig,'currentaxes',ud.PlotHandle)
	end
end

%--------------------------------------------------------------------------------------------------------
function PlaceLegendOnTop(hf,hl,ha)
%PlaceLengendOpTop  Make sure the legend_ is on top of its axes.
ord = findobj(allchild(hf),'flat','type','axes');
axpos = find(ord==ha);
legpos = find(ord==hl);
axlegdist = axpos - legpos;

% legend_ needs to be the next axes type child above its plotaxes
% if it's higher than that, move it back down where it belongs.
% I'm not sure if this is really needed, but it IS the old
% (expected) behavior.  Just above should be good enough
% For now I'm commenting out the stacking down.

if axlegdist>1
	% this may be needed if legend_ being more than one child above its
	% axes causes problems (overlapping axes?), but hopefully it wont
	axes% since using it will cause the flashing toolbar bug to return.
	% uistack(hl,'down',axlegdist-1);
	
	% need to move legend_ up if its stack order number is larger
	% than (or same as - impossible?) that of the plot axes.  Infact
	% this may not even be needed, because legend_ is never called to
	% do anything to an old legend_ - it always creates a new one, which
	% will always be on top.  So this whole thing may be unnecessary.
elseif axlegdist<1
	uistack(hl,'up',1-axlegdist);
end

%------------------------------
function info = legend_info(ha,hf,Kids,lstrings)
%LEGEND_INFO Get legend_ info from parent axes, Kids, and strings
%   INFO = LEGEND_INFO(HA,KIDS,STRINGS) returns a structure array containing
%      objtype  -- Type of object 'line', 'patch', or 'surface'
%      label    -- label string
%      linetype -- linetype;
%      edgecol  -- edge color
%      facecol  -- face color
%      lnwidth  -- line width
%      marker   -- marker
%      marksize -- markersize
%      markedge -- marker edge color
%      markface -- marker face color (not used for 'line')

defaultnverts = 4;

linetype = {};
edgecol = {};
facecol = {};
lnwidth = {};
marker = {};
marksize = {};
markedge = {};
markface = {};
nverts = {};

% These 8 variables are the important ones.  The only ambiguity is
% edgecol/facecol.  For lines, edgecol is the line color and facecol
% is unused.  For patches, edgecol/facecol mean the logical thing.

Kids = Kids(:);  %  Reshape so that we have a column vector of handles.
lstrings = lstrings(:);

% Check for valid handles
nonhandles = ~ishandle(Kids);
if any(nonhandles)
	%  warning('Some invalid handles were ignored.')
	Kids(nonhandles) = [];
end
if ~isempty(Kids)
	badhandles = ~(strcmp(get(Kids,'type'),'patch') | strcmp(get(Kids,'type'),'line') | strcmp(get(Kids,'type'),'surface'));
	if any(badhandles)
		warning('Some handles to non-lines and/or non-solid color objects were ignored.')
		Kids(badhandles) = [];
	end
end

objtype = get(Kids,{'type'});
nk = length(Kids);
nstack = length(lstrings);
n = min(nstack,nk);

% Treat empty strings as a single space
for i=1:nstack
	if isempty(lstrings{i})
		lstrings{i} = ' ';
	end
end

% Truncate kids if necessary to match the number of strings
objtype = objtype(1:n);
Kids = Kids(1:n);

for i=1:n
	linetype = [linetype,get(Kids(i),{'LineStyle'})];
	if strcmp(objtype{i},'line')
		edgecol = [edgecol,get(Kids(i),{'Color'})];
		facecol = [facecol,{'none'}];
		nverts = [nverts,{defaultnverts}];
	elseif strcmp(objtype{i},'patch') || strcmp(objtype{i},'surface')
		[e,f] = patchcol(Kids(i));
		edgecol = [edgecol,{e}];
		facecol = [facecol,{f}];
		nverts = [nverts,{length(get(Kids(i),'xdata'))}];
	end
	lnwidth = [lnwidth,get(Kids(i),{'LineWidth'})];
	marker = [marker,get(Kids(i),{'Marker'})];
	marksize = [marksize,get(Kids(i),{'MarkerSize'})];
	markedge = [markedge,get(Kids(i),{'MarkerEdgeColor'})];
	markface = [markface,get(Kids(i),{'MarkerFaceColor'})];
end

if n < nstack     % More strings than handles
	objtype(end+1:nstack) = {'none'};
	linetype(end+1:nstack) = {'none'};
	edgecol(end+1:nstack) = {'none'};
	facecol(end+1:nstack) = {'none'};
	lnwidth(end+1:nstack) = {get(hf,'defaultlinelinewidth')};
	marker(end+1:nstack) = {'none'};
	marksize(end+1:nstack) = {get(hf,'defaultlinemarkersize')};
	markedge(end+1:nstack) = {'auto'};
	markface(end+1:nstack) = {'auto'};
	nverts(end+1:nstack) = {defaultnverts};
end

% Limit markersize to axes fontsize
fonts = get(ha,'fontsize');
marksize([marksize{:}]' > fonts & strcmp(objtype(:),'line')) = {fonts};
marksize([marksize{:}]' > fonts/2 & strcmp(objtype(:),'patch')) = {fonts/2};

% Package everything into the info structure
info = struct('objtype',objtype(:),'label',lstrings(:),...
	'linetype',linetype(:),'edgecol',edgecol(:),...
	'facecol',facecol(:),'lnwidth',lnwidth(:),'marker',marker(:),...
	'marksize',marksize(:),'markedge',markedge(:),'markface',markface(:),'nverts',nverts(:));


%-----------------------------
function update_all_legends
%UPDATE_ALL_LEGENDS Update all legends on this figure
legh = find_legend;
for i=1:length(legh)
	update_legend(legh(i));
end

%-------------------------------
function [hl,objh,outH,outM] = update_legend(legh)
%UPDATE_LEGEND Update an existing legend
if isempty(legh)
	hl = [];
	objh = [];
	return
end
if length(legh)~=1
	error('Can only update one legend_ at a time.')
end

ud = get(legh,'userdata');
if ~isfield(ud,'LegendPosition')
	warning('No legend to update.')
	hl = []; objh = [];
	return
end

moved = DidLegendMove(legh);

units = get(legh,'units');
set(legh,'units','points')
oldpos = get(legh,'position');

% Delete old legend_
delete_legend(legh)

% Make a new one
if moved || length(ud.legendpos)==4
	[hl,objh] = make_legend(ud.PlotHandle,ud.handles,ud.lstrings,oldpos);
else
	[hl,objh] = make_legend(ud.PlotHandle,ud.handles,ud.lstrings,ud.legendpos);
end
set(hl,'units',units)
outH = ud.handles;
outM = ud.lstrings;

%----------------------------------------------
function moved = DidLegendMove(legh)
% Check to see if the legend_ has been moved
ud = get(legh,'userdata');
if isstruct(ud) && isfield(ud,'LegendPosition')
	units = get(legh,'units');
	set(legh,'units','normalized')
	pos = get(legh,'position');
	set(legh,'units','pixels')
	tol = pos ./ get(legh,'position')/2;
	if any(abs(ud.LegendPosition - pos) > max(tol(3:4)))
		moved = 1;
	else
		moved = 0;
	end
	set(legh,'units',units)
else
	moved = 1;
end

%----------------------------------------------
function moved = DidAxesMove(legh)
% Check to see if the axes has been moved
ud = get(legh,'userdata');
ax = ud.PlotHandle;
if isfield(ud,'PlotPosition')
	units = get(ax,'units');
	set(ax,'units','normalized')
	pos = get(ax,'position');
	set(ax,'units','pixels')
	tol = pos ./ get(ax,'position')/2;
	if any(abs(ud.PlotPosition - pos) > max(tol(3:4)))
		moved = 1;
	else
		moved = 0;
	end
	set(ax,'units',units)
else
	moved = 0;
end

%----------------------------
function delete_legend(ax)
%DELETE_LEGEND Remove legend_ from plot
if isempty(ax) || ~ishandle(ax), return, end
ax = ax(1);

ud=get(ax,'UserData');
if isfield(ud,'DeleteProxy') && ishandle(ud.DeleteProxy)
	delete(ud.DeleteProxy)
end

%-------------------------------
function [lpos,axpos] = getposition(ha,legendpos,llen,lhgt)
%GETPOS Get position vector from legendpos code
stickytol=1;
cap=get(ha,'position');
edge = 5; % 5 Point edge -- this number also used in make_legend

if length(legendpos)==4
	% Keep the top at the same place
	Pos = [legendpos(1) legendpos(2)+legendpos(4)-lhgt];
else
	switch legendpos
		case 0
			Pos = lscan(ha,llen,lhgt,0,stickytol);
		case 1
			Pos = [cap(1)+cap(3)-llen-edge cap(2)+cap(4)-lhgt-edge];
		case 2
			Pos = [cap(1)+edge cap(2)+cap(4)-lhgt-edge];
		case 3
			Pos = [cap(1)+edge cap(2)+edge];
		case 4
			Pos = [cap(1)+cap(3)-llen-edge cap(2)+edge];
		otherwise
			Pos = -1;
	end
end

if isequal(Pos,-1)
	axpos=[cap(1) cap(2) cap(3)-llen-.03 cap(4)];
	lpos=[cap(1)+cap(3)-llen+edge cap(4)+cap(2)-lhgt llen lhgt];
	if any(axpos<0) || any(lpos<0)
		warning('Insufficient space to draw legend_.')
		if any(axpos<0), axpos = []; end
	end
else
	axpos=[];
	lpos=[Pos(1) Pos(2) llen lhgt];
end

%--------------------------------------------
function Pos = lscan(ha,wdt,hgt,tol,stickytol)
%LSCAN  Scan for good legend_ location.

% Calculate tile size
cap=get(ha,'Position'); % In Point coordinates
xlim=get(ha,'Xlim');
ylim=get(ha,'Ylim');
H=ylim(2)-ylim(1);
W=xlim(2)-xlim(1);

dh = 0.03*H;
dw = 0.03*W;
Hgt = hgt*H/cap(4);
Wdt = wdt*W/cap(3);
Thgt = H/max(1,floor(H/(Hgt+dh)));
Twdt = W/max(1,floor(W/(Wdt+dw)));
dh = (Thgt - Hgt)/2;
dw = (Twdt - Wdt)/2;

% Get data, points and text

Kids=get(ha,'children');
Xdata=[];Ydata=[];
for i=1:length(Kids)
	type = get(Kids(i),'type');
	if strcmp(type,'line')
		xk = get(Kids(i),'Xdata');
		yk = get(Kids(i),'Ydata');
		n = length(xk);
		if n < 100 && n > 1
			xk = interp1(xk,linspace(1,n,200));
			yk = interp1(yk,linspace(1,n,200));
		end
		Xdata=[Xdata,xk];
		Ydata=[Ydata,yk];
	elseif strcmp(type,'patch') || strcmp(type,'surface')
		xk = get(Kids(i),'Xdata');
		yk = get(Kids(i),'Ydata');
		Xdata=[Xdata,xk(:)'];
		Ydata=[Ydata,yk(:)'];
	elseif strcmp(get(Kids(i),'type'),'text')
		tmpunits = get(Kids(i),'units');
		set(Kids(i),'units','data')
		tmp=get(Kids(i),'Position');
		ext=get(Kids(i),'Extent');
		set(Kids(i),'units',tmpunits);
		Xdata=[Xdata,[tmp(1) tmp(1)+ext(3)]];
		Ydata=[Ydata,[tmp(2) tmp(2)+ext(4)]];
	end
end
in = finite(Xdata) & finite(Ydata);
Xdata = Xdata(in);
Ydata = Ydata(in);

% Determine # of data points under each "tile"
xp = (0:Twdt:W-Twdt) + xlim(1);
yp = (0:Thgt:H-Thgt) + ylim(1);
wtol = Twdt / 100;
htol = Thgt / 100;
pop = zeros(length(yp), length(xp));
for j=1:length(yp)
	for i=1:length(xp)
		pop(j,i) = sum(sum((Xdata > xp(i)-wtol) & (Xdata < xp(i)+Twdt+wtol) & (Ydata > yp(j)-htol) & (Ydata < yp(j)+Thgt+htol)));
	end
end

if all(pop(:) == 0), pop(1) = 1; end

% Cover up fewest points.  After this while loop, pop will
% be lowest furthest away from the data
while any(pop(:) == 0)
	newpop = filter2(ones(3),pop);
	if all(newpop(:) ~= 0)
		break;
	end
	pop = newpop;
end
[j,i] = find(pop == min(pop(:)));
xp =  xp - xlim(1) + dw;
yp =  yp - ylim(1) + dh;
Pos = [cap(1)+xp(i(end))*cap(3)/W
	cap(2)+yp(j(end))*cap(4)/H];

%--------------------------------
function Kids = getchildren(ha)
%GETCHILDREN Get children that can have legends
%   Note: by default, lines get labeled before patches;
%   patches get labeled before surfaces.

linekids = findobj(ha,'type','line');
surfkids = findobj(ha,'type','surface');
patchkids = findobj(ha,'type','patch');

if ~isempty(linekids)
	goodlk = ones(1,length(linekids));
	for i=1:length(linekids)
		if ((isempty(get(linekids(i),'xdata')) || isallnan(get(linekids(i),'xdata'))) && ...
			(isempty(get(linekids(i),'ydata')) || isallnan(get(linekids(i),'ydata'))) && ...
			(isempty(get(linekids(i),'zdata')) || isallnan(get(linekids(i),'zdata'))) )
			goodlk(i) = 0;
		end
	end
	linekids = linekids(logical(goodlk));
end

if ~isempty(surfkids)
	goodsk = ones(1,length(surfkids));
	for i=1:length(surfkids)
		if ((isempty(get(surfkids(i),'xdata')) || isallnan(get(surfkids(i),'xdata'))) && ...
			(isempty(get(surfkids(i),'ydata')) || isallnan(get(surfkids(i),'ydata'))) && ...
			(isempty(get(surfkids(i),'zdata')) || isallnan(get(surfkids(i),'zdata'))) )
			goodsk(i) = 0;
		end
	end
	surfkids = surfkids(logical(goodsk));
end

if ~isempty(patchkids)
	goodpk = ones(1,length(patchkids));
	for i=1:length(patchkids)
		if ((isempty(get(patchkids(i),'xdata')) || isallnan(get(patchkids(i),'xdata'))) && ...
			(isempty(get(patchkids(i),'ydata')) || isallnan(get(patchkids(i),'ydata'))) && ...
			(isempty(get(patchkids(i),'zdata')) || isallnan(get(patchkids(i),'zdata'))) )
			goodpk(i) = 0;
		end
	end
	patchkids = patchkids(logical(goodpk));
end

Kids = flipud([surfkids ; patchkids ; linekids]);

%----------------------------
function allnan = isallnan(d)
nans = isnan(d);
allnan = all(nans(:));

%----------------------------
function s = getstrings(c)
%GETSTRINGS Get strings from legend_ input
%   S = GETSTRINGS(C) where C is a cell array containing the legend_
%   input arguments.  Handles three cases:
%      (1) legend_(M) -- string matrix
%      (2) legend_(C) -- cell array of strings
%      (3) legend_(string1,string2,string3,...)
%   Returns a cell array of strings
if length(c)==1 % legend_(M) or legend_(C)
	s = cellstr(c{1});
elseif iscellstr(c)
	s = c;
else
	error('Legend labels must be strings.');
end

%-----------------------------------
function  out=ctorgb(arg)
%CTORGB Convert color string to rgb value
switch arg
	case 'y', out=[1 1 0];
	case 'm', out=[1 0 1];
	case 'c', out=[0 1 1];
	case 'r', out=[1 0 0];
	case 'g', out=[0 1 0];
	case 'b', out=[0 0 1];
	case 'w', out=[1 1 1];
	otherwise, out=[0 0 0];
end


%----------------------------------
function  [edgecol,facecol] = patchcol(h)
%PATCHCOL Return edge and facecolor from patch handle
cdat = get(h,'Cdata');
facecol = get(h,'FaceColor');
if strcmp(facecol,'interp') || strcmp(facecol,'texturemap')
	if ~all(cdat == cdat(1))
		warning(['Legend not supported for patches with FaceColor = ''',facecol,''''])
	end
	facecol = 'flat';
end
if strcmp(facecol,'flat')
	if size(cdat,3) == 1       % Indexed Color
		k = find(finite(cdat));
		if isempty(k)
			facecol = 'none';
		end
	else                       % RGB values
		facecol = reshape(cdat(1,1,:),1,3);
	end
end

edgecol = get(h,'EdgeColor');
if strcmp(edgecol,'interp')
	if ~all(cdat == cdat(1))
		warning('Legend not supported for patches with EdgeColor = ''interp''.')
	end
	edgecol = 'flat';
end
if strcmp(edgecol,'flat')
	if size(cdat,3) == 1      % Indexed Color
		k = find(finite(cdat));
		if isempty(k)
			edgecol = 'none';
		end
	else                      % RGB values
		edgecol = reshape(cdat(1,1,:),1,3);
	end
end

%------------------------------------------------------------------------------------------------------
function show_plot(legend_ax)
%Set the axes this legend_ goes with to the current axes

if islegend(legend_ax)
	ud = get(legend_ax,'userdata');
	if isfield(ud,'PlotHandle')
		set(find_gcf(legend_ax),'currentaxes',ud.PlotHandle)
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tf=islegend(ax)
if length(ax) ~= 1 || ~ishandle(ax)
	tf = false;
else
	tf=strcmp(get(ax,'tag'),'legend_');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tOut=singleline(tIn)
%converts cellstrs and 2-d char arrays to
%\n-delimited single-line text

if ischar(tIn)
	if size(tIn,1)>1
		nRows=size(tIn,1);
		cr=char(10);
		cr=cr(ones(nRows,1));
		tIn=[tIn,cr]';
		tOut=tIn(:)';
		tOut=tOut(1:end-1); %remove trailing \n
	else
		tOut=tIn;
	end
elseif iscellstr(tIn)
	tOut=singleline(char(tIn));
else
	tOut='';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function expandplot(ha,ud,legendpos)
% Expand plot to cover rectangle that includes legend_
% but when a colorbar is present only expand up to
% the left edge of the colorbar minus a margin.

edge=5;

% get all the figure colorbars
cbars = findobj(get(ha,'parent'),'tag','Colorbar');
cbarfound = 0;
i=1;
% find vertical colorbar with plothandle ha
% and get its position
while ~cbarfound && i<=length(cbars)
	cbud = get(cbars(i),'userdata');
	if isfield(cbud,'PlotHandle')
		if isequal(cbud.PlotHandle,ha)
			cbpos = get(cbars(i),'Position');
			if cbpos(3)<cbpos(4) && strcmpi(get(cbars(i),'beingdeleted'),'off')  % vertical
				cbarfound=1;
				% get its position in point units
				oldunits = get(cbars(i),'units');
				set(cbars(i),'units','points');
				cbppos = get(cbars(i),'Position');
				set(cbars(i),'units',oldunits);
			end
		end
	end
	i = i+1;
end

% get axes position
oldunits = get(ha,'units');
set(ha,'units','points');
axespos = get(ha,'Position');

% get legend_ position
set(ud.LegendHandle,'units','points');
legpos  = get(ud.LegendHandle,'Position');

startpt = max(min(axespos(1:2),legpos(1:2)),[0 0]);
endpt   = max(axespos(1:2)+axespos(3:4), legpos(1:2)+legpos(3:4));
% adjust endpoint for vertical colorbar if one was found
if cbarfound
	endpt(1) = min(endpt(1),cbppos(1)-edge);
end
newpos = [startpt (endpt-startpt)];


legright = legpos(1) + legpos(3);
axright = axespos(1) + axespos(3);
if legendpos<0
	if cbarfound || legright>axright
		% subtract edge size
		newpos(3)=newpos(3)-edge;
	elseif  axright>legright && (axright-legright)<edge
		% subtract edge less distance from right of legend_ to right of axes
		newpos(3)=newpos(3) - (edge-(axright-legright));
	end
end

set(ha,'position',newpos);
set(ha,'units',oldunits);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout=record_size(hPlot,plotSize)

if nargin<2
	oldUnits=get(hPlot,'Units');
	set(hPlot,'Units','normalized');
	plotSize=get(hPlot,'Position');
	set(hPlot,'Units',oldUnits);
end

if nargout==1
	varargout{1}=plotSize;
else
	hLegend=find_legend(hPlot);
	ud = get(hLegend,'UserData');
	if (isfield(ud,'legendpos') && length(ud.legendpos)==1 && ud.legendpos<1) || ...
		(isfield(ud,'NegativePositionTripped') && ud.NegativePositionTripped)
		ud.PlotHandle=hPlot;
		ud.PlotPosition=plotSize;
		set(hLegend,'UserData',ud);
	end
end

% ----------------------------------------------------------------------------------------------------
function moveaxis_(hObj, evt, arg)
%MOVEAXIS Used by LEGEND to enable dragging of legend_.
%   This function is for private use of LEGEND, and is
%   not intended for use by other programs.  It is
%   automatically installed as the ButtonDownFcn of the
%   legend_ axes and its children, enabling dragging of
%   legends with the mouse.

% get the figure handle in a way that's robust to hidden handles
fig = gcbf;
% don't do anything if we're not called from a callback:
if isempty(fig)
	return
end

% don't do anything if we're supposed to be zooming
tmp = zoom_j(fig,'getmode');
if (isequal(tmp,'in') || isequal(tmp,'on')),		return,		end

if (nargin == 2)		% we're getting called as a result of the ButtonDownFcn:
	lax = find_legend_axes(gcbo);
	% if we can't identify a legend_ as the ancestor of the object that was clicked on, bail out.
	if isempty(lax)
		return
	end
	st=get(fig,'SelectionType');
	if (strcmp(st,'normal'))
		moveaxis_([], [], 1);
	elseif strcmp(st,'open') && strcmp(get(gcbo,'Type'),'text')
		edit_legend(gcbo)
	end
	
elseif (arg == 1)		% set up for a drag
	mad.oldFigUnits = get(fig,'units');             set(fig,'units','pixels');
	mad.oldFigPointer = get(fig,'pointer');         set(fig,'pointer','fleur');
	mad.legendAxes = find_legend_axes(gcbo);
	mad.oldAxesUnits = get(mad.legendAxes,'units'); set(mad.legendAxes,'units','pixels');
	pnt=get(fig,'currentpoint');
	pos = get(mad.legendAxes,'position');
	mad.DELTA = pos(1:2) - pnt;
	mad.oldWindowButtonMotionFcn = get(fig,'windowbuttonmotionfcn');
	set(fig,'windowbuttonmotionfcn',{@moveaxis_,2});
	mad.oldWindowButtonUpFcn = get(fig,'windowbuttonupfcn');
	set(fig,'windowbuttonupfcn',{@moveaxis_,3});
	setappdata(fig,'moveaxisData',mad);
	
elseif (arg == 2)		% mouse motion during a drag
	mad = getappdata(fig,'moveaxisData');
	pos=get(mad.legendAxes,'position');
	set(mad.legendAxes,'position',[get(fig,'currentpoint')+mad.DELTA pos(3:4)]);
	
elseif (arg == 3)		% mouse release - done dragging
	mad = getappdata(fig,'moveaxisData');
	set(fig,'WindowButtonMotionfcn',mad.oldWindowButtonMotionFcn,...
		'WindowButtonUpFcn',mad.oldWindowButtonUpFcn,...
		'pointer', mad.oldFigPointer,...
		'units', mad.oldFigUnits);
	set(mad.legendAxes, 'units', mad.oldAxesUnits);
	show_plot(mad.legendAxes)
	rmappdata(fig,'moveaxisData');
end

function lax = find_legend_axes(obj)
lax = [];
while ~isempty(obj)
	if strcmp(get(obj,'type'),'axes') && strcmp(get(obj,'tag'),'legend_')
		lax = obj;
		break
	end
	obj = get(obj,'parent');
end

%------------------------
function edit_legend(gco)
%Edit a legend

if ~strcmp(get(gco,'type'),'text'), return, end
legh = get(gco,'parent');

% Determine which string was clicked on
units = get(gco,'units');
set(gco,'units','data')
cp = get(legh,'currentpoint');
ext = get(gco,'extent');
nstrings = size(get(gco,'string'),1);

% The k-th string (from the top) was clicked on
k = floor((ext(4) - cp(1,2))/ext(4)*nstrings) + 1;
ud = get(legh,'userdata');
nrows = cellfun('size',ud.lstrings,1);
crows = cumsum(nrows);

% Determine which string in the cell array was clicked on
active_string = floor(interp1([0 cumsum(nrows)+1],[1:length(nrows) length(nrows)],k));
if isnan(active_string), return, end

% Disable legend_ buttondownfcn's
%savehandle = findobj(legh,'buttondownfcn','moveaxis');
savehandle = legh;
set(savehandle,'buttondownfcn','')

% Make a editable string on top of the legend_ string
pos = get(gco,'position');
y = ext(4) - (crows(active_string)-nrows(active_string)/2)*ext(4)/nstrings;
pos(2) = ext(2) + y;

% Coverup text
CoverHandle = copyobj(gco,legh);
color = get(legh,'color');
if ischar(color), color = get(get(legh,'parent'),'color'); end
set(CoverHandle,'color',color);

% Make editable case
TextHandle = copyobj(gco,legh);
set(TextHandle,'string',char(ud.lstrings{active_string}),'position',pos, 'Editing','on')
waitfor(TextHandle,'Editing');

% Protect against the handles being destroyed during the waitfor
if ishandle(CoverHandle)
	delete(CoverHandle)
end
if ishandle(TextHandle) && ishandle(legh) && ishandle(savehandle)
	newLine = get(TextHandle,'String');
	delete(TextHandle)
	
	ud.lstrings{active_string}=newLine;
	set(legh,'UserData',ud)
	set(gco,'units',units)
	
	% Enable legend_ buttondfcn's
	set(savehandle,'buttondownfcn',@moveaxis_)
	resize_legend(legh);
	%update_legend(legh);
end

%----------------------------
function resize_legend(legh)
%RESIZE_LEGEND Resize all legend_ in this figure

ud = get(legh,'userdata');
units = get(legh,'units');
set(legh,'units','normalized')

if ~isfield(ud,'LegendPosition') || ~ishandle(ud.LegendHandle)
	warning('No legend_ to update.')
	return
end

moved = DidLegendMove(legh);

set(legh,'units','points')
oldpos = get(legh,'position');

% Update the legend_
if moved || length(ud.legendpos)==4
	[hl,objh] = make_legend(ud.PlotHandle,ud.handles,ud.lstrings,oldpos,ud,'resize');
else
	[hl,objh] = make_legend(ud.PlotHandle,ud.handles,ud.lstrings,ud.legendpos,ud,'resize');
end
% make_legend returns empty on error
if ~isempty(hl)
	set(hl,'units',units)
end

%-------------------------------
function [hl,labelhandles] = make_legend(ha,Kids,lstrings,legendpos,ud,resize)
%MAKE_LEGEND Make legend given parent axes, kids, and strings
%
%   MAKE_LEGEND(...,ud) is called from the resizeFcn. In this case
%   just update the position of the legend_ pieces instead of recreating it from scratch.

ud.PlotHandle = ha;
hf = get(ha,'parent'); % Parent figure
doresize = 0;
if (nargin>=6 && isequal(resize,'resize')), doresize = 1; end

% Get the legend_ info structure from the inputs
info = legend_info(ha,hf,Kids,lstrings);

% Remember current state
hfold = find_gcf(ha);
haold = find_gca(ha);
punits=get(hf,'units');
aunits=get(ha,'units');
% Remember Figure Default Text Font Units and Size
oldFigDefaultTextFontUnits = get(hf,'DefaultTextFontUnits');
oldFigDefaultTextFontSize = get(hf,'DefaultTextFontSize');

if strncmp(get(hf,'NextPlot'),'replace',7)
	set(hf,'NextPlot','add')
	oldNextPlot = get(hf,'NextPlot');
else
	oldNextPlot = '';
end
set(ha,'units','points');
set(hf,'units','points');

if ~doresize
	textStyleSource=ha;
	tInterp = 'tex';
else
	textStyleSource=ud.LabelHandles(1);
	tInterp = get(textStyleSource,'interpreter');
end

% Determine size of legend_ in figure points
oldUnits= get(textStyleSource,'FontUnits');
set(textStyleSource,'FontUnits','points');
fontn = get(textStyleSource,'fontname');
fonts = get(textStyleSource,'fontsize');
fonta = get(textStyleSource,'fontangle');
fontw = get(textStyleSource,'fontweight');
set(textStyleSource,'FontUnits',oldUnits);
% Set figure Default Text Font Size and Units
% Otherwise these can interact badly with legend_ fontsize
set(hf,'DefaultTextFontUnits','points');
set(hf,'DefaultTextFontSize',fonts);

% Symbols are the size of 3 numbers, the temphackytext tag is used so that
% scribefiglisten will know not to turn zoom or rotate3d off when it is added or deleted
h = text(0,0,'123',...
	'fontname',fontn,...
	'fontsize',fonts,...
	'fontangle',fonta,...
	'fontweight',fontw,...
	'visible','off',...
	'tag','temphackytext',...
	'units','points',...
	'parent',ha,...
	'interpreter',tInterp);
ext = get(h,'extent');
lsym = ext(3);
loffset = lsym/3;
delete(h);

% Make box big enough to handle longest string
% the temphackytext tag is used so
% that scribefiglisten will know not
% to turn zoom or rotate3d off when it is added
% or deleted
h=text(0,0,{info.label},...
	'fontname',fontn,...
	'fontsize',fonts,...
	'fontangle',fonta,...
	'fontweight',fontw,...
	'units','points',...
	'visible','off',...
	'tag','temphackytext',...
	'interpreter',tInterp,...
	'parent',ha);

ext = get(h,'extent');
width = ext(3);
height = ext(4)/size(get(h,'string'),1);
margin = height*0.075;
delete(h);

llen = width + loffset*3 + lsym;
lhgt = ext(4) + 2*margin;

%We reposition an axes if its position is becoming -1 or if setting
%a position of 1,2,3,4 when the position has been -1 in the past.
repositionAxes=false;
if length(legendpos)==1
	legendpos=round(legendpos); %deal with non-integer numbers
	if legendpos<0 || legendpos>4
		legendpos=-1; %do this in case someone passes in something not in [-1,4]
		ud.NegativePositionTripped = true;
		repositionAxes=true;
		
		% Remember old axes position if resizing to -1
		if ~doresize
			ud.PlotPosition = record_size(ha);
		end
	else
		if isfield(ud,'NegativePositionTripped') && ud.NegativePositionTripped
			repositionAxes=true;
			if isfield(ud,'PlotPosition')
				ud=rmfield(ud,'PlotPosition');
			end
		end
		ud.NegativePositionTripped=false;
	end
end

% If resizing a plot, temporarily set the
% axes position to cover a rectangle that also encompasses the old legend_.
if doresize && repositionAxes
	expandplot(ha,ud,legendpos)
end

[lpos,axpos] = getposition(ha,legendpos,llen,lhgt);

ud.legendpos = legendpos;

% Shrink axes if necessary
if ~isempty(axpos)
	set(ha,'Position',axpos)
end

%
% Create legend_ object
%
if strcmp(get(ha,'color'),'none')
	acolor = get(hf,'color');
else
	acolor = get(ha,'color');
end

if ~doresize
	% Create legend_ axes and LegendDeleteProxy object (an
	% invisible text object in target axes) so that the
	% legend_ will get deleted correctly.
	ud.DeleteProxy = text('parent',ha, 'visible','off', 'tag','LegendDeleteProxy', 'handlevisibility','off');
	
	hl=axes('units','points',...
		'position',lpos,...
		'box','on',...
		'nextplot','add',...
		'xtick',[-1],...
		'ytick',[-1],...
		'xticklabel','',...
		'yticklabel','',...
		'xlim',[0 1],...
		'ylim',[0 1], ...
		'clipping','on',...
		'color',acolor,...
		'tag','legend_',...
		'view',[0 90],...
		'climmode',get(ha,'climmode'),...
		'clim',get(ha,'clim'),...
		'parent',hf);
	
	set(hl,'units','normalized')
	setappdata(hl,'NonDataObject',[]); % Used by DATACHILDREN.M
	ud.LegendPosition = get(hl,'position');
	set(ud.DeleteProxy,'userdata',hl);
else
	hl = ud.LegendHandle;
	labelhandles = ud.LabelHandles;
	set(hl,'units','points','position',lpos);
	set(hl,'units','normalized')
	ud.LegendPosition = get(hl,'position');
end

%
% Draw text description above legend_
%
nrows = size(char(info.label),1);

% draw text one on chunk so that the text spacing is good
label = char(info.label);
top = (1-max(1,size(label,1)))/2;
if ~doresize
	texthandles = text('parent',hl,...
		'units','data',...
		'position',[1-(width+loffset)/llen,1-(1-top)/(nrows+1)],...
		'string',char(info.label),...
		'fontname',fontn,...
		'fontweight',fontw,...
		'fontsize',fonts,...
		'fontangle',fonta,...
		'ButtonDownFcn',@moveaxis_);
else
	texthandles = ud.LabelHandles(1);
	
	oldFontUnits=get(texthandles,'FontUnits');
	oldErr=lasterr;
	try
		set(texthandles,...
			'String',{info.label},...
			'units','data',...
			'position',[1-(width+loffset)/llen,1-(1-top)/(nrows+1)],...
			'fontname',fontn,...
			'fontsize',fonts,...
			'FontUnits','points');
		set(texthandles,'FontUnits',oldFontUnits);
	catch
		lasterr(oldErr);
		% make HL empty so caller knows we got toasted
		hl = [];
		return
	end
end

% adjust text position
ext = get(texthandles,'extent');
centers = linspace(ext(4)-ext(4)/nrows,0,nrows)+ext(4)/nrows/2 + 0.4*(1-ext(4));
edges = linspace(ext(4),0,nrows+1) + 0.4*(1-ext(4));
indent = [1 1 -1 -1 1] * ext(4)/nrows/7.5;

%
% Draw lines and / or styles and labels for each legend_ item
%

% start handleIndex at 2 cuz labels and line styles follow the text object from above in the handle list
handleIndex = 2;
r = 1;
objhandles = [];
nstack = length(info);

for i=1:nstack
	p = [];
	if strcmp(info(i).objtype,'line')
		
		% draw lines with markers like this: --*--
		
		% draw line
		if ~doresize
			p = line('parent',hl,...
				'xdata',(loffset+[0 lsym])/llen,...
				'ydata',[centers(r) centers(r)],...
				'linestyle',info(i).linetype,...
				'marker','none',...
				'tag',singleline(info(i).label),...
				'color',info(i).edgecol, ...
				'linewidth',info(i).lnwidth,...
				'ButtonDownFcn',@moveaxis_,...
				'SelectionHighlight','off');
			
			markerHandle=line('parent',hl,...
				'xdata',(loffset+lsym/2)/llen,...
				'ydata',centers(r),...
				'color',info(i).edgecol, ...
				'HitTest','off',...
				'linestyle','none',...
				'linewidth',info(i).lnwidth,...
				'marker',info(i).marker,...
				'markeredgecolor',info(i).markedge,...
				'markerfacecolor',info(i).markface,...
				'markersize',info(i).marksize,...
				'ButtonDownFcn',@moveaxis_);
			
			% set line handle after legendmarker handle so setting
			% marker handle won't cause change to kid line markers
			p=[p;markerHandle];
			
		else
			
			% For legends created pre R12, don't try to draw lines where linestyle is none.
			p1 = ud.LabelHandles(handleIndex);
			if strcmp(get(p1, 'tag'), singleline(info(i).label)) || ~strcmp(info(i).linetype,'none')
				set(p1, 'xdata',(loffset+[0 lsym])/llen, 'ydata',[centers(r) centers(r)]);
				handleIndex = handleIndex+1;
				p = p1;
			end
			
			% For legends created pre R12, don't try to draw markers where markerstyle is none.
			if strcmp(get(p, 'tag'), singleline(info(i).label)) || ~strcmp(info(i).linetype,'none')
				p2 = ud.LabelHandles(handleIndex);
				set(p2, 'xdata',(loffset+lsym/2)/llen, 'ydata',centers(r));
				handleIndex = handleIndex+1;
				p = [p;p2];
			end
			
		end
		
	elseif strcmp(info(i).objtype,'patch') || strcmp(info(i).objtype,'surface')
		% draw patches
		
		% Adjusting ydata to make a thinner box will produce nicer
		% results if you use patches with markers.
		
		% set patch xdata depending on n vertices in axes patch object
		
		if info(i).nverts == 1
			pxdata = (loffset + (lsym/2))/llen;
			pydata = (((2*edges(r)) + (1*edges(r+1)))/3) - indent;
			mksize = lsym/2.3;
		else
			pxdata = (loffset+[0 lsym lsym 0 0])/llen;
			pydata = [edges(r) edges(r) edges(r+1) edges(r+1) edges(r)]-indent;
			mksize = info(i).marksize;
		end
		
		if ~doresize
			p = patch('parent',hl,...
				'xdata',pxdata,...
				'ydata',pydata,...
				'linestyle',info(i).linetype,...
				'edgecolor',info(i).edgecol, ...
				'facecolor',info(i).facecol,...
				'linewidth',info(i).lnwidth,...
				'tag',singleline(info(i).label),...
				'marker',info(i).marker,...
				'markeredgecolor',info(i).markedge,...
				'markerfacecolor',info(i).markface,...
				'markersize',mksize,...
				'SelectionHighlight','off',...
				'ButtonDownFcn',@moveaxis_);
		else
			p = ud.LabelHandles(handleIndex);
			set(p,'xdata',(loffset+[0 lsym lsym 0 0])/llen, 'ydata',[edges(r) edges(r) edges(r+1) edges(r+1) edges(r)]-indent);
			handleIndex = handleIndex+1;
		end
		
		if strcmp(info(i).facecol,'flat') || strcmp(info(i).edgecol,'flat')
			c = get(Kids(i),'cdata');
			k = min(find(finite(c)));
			if ~isempty(k)
				set(p,'cdata',c(k)*ones(1,5),'cdatamapping',get(Kids(i),'cdatamapping'));
			end
		end
	end
	
	% p will be empty when resizing a pre R12 legend_ with lines using linestyle none.
	if ~isempty(p)
		objhandles = [objhandles;p];
	end
	
	r = r + max(1,size(info(i).label,1));
end

% set both of these cuz label handles is an output argument to this function
labelhandles = [texthandles;objhandles];
ud.LabelHandles = labelhandles;

% Clean up a bit
set(hf,'DefaultTextFontUnits',oldFigDefaultTextFontUnits);
set(hf,'DefaultTextFontSize',oldFigDefaultTextFontSize);
set(hf,'units',punits)
set(ha,'units',aunits)
if (hfold ~= hf) && ~doresize, figure(hfold); end
if ~isempty(oldNextPlot)
	set(hf,'nextplot',oldNextPlot)
end
ud.handles = Kids;
ud.lstrings = {info.label};

ud.LegendHandle = hl;
set(hl, 'ButtonDownFcn',@moveaxis_, 'interruptible','on', 'busyaction','queue', 'userdata',ud);

if ~doresize
	PlaceLegendOnTop(hf,hl,ha)
end

set(hf,'currentaxes',haold);  % this should be last
