function hh = clabel_j(cs,varargin)
%CLABEL Contour plot elevation labels.
%   CLABEL(CS,H) adds height labels to the current contour
%   plot.  The labels are rotated and inserted within the contour
%   lines.  CS and H are the contour matrix output and object handle
%   outputs from CONTOUR, CONTOUR3, or CONTOURF.
%
%   CLABEL(CS,H,V) labels just those contour levels given in
%   vector V.  The default action is to label all known contours.
%   The label positions are selected randomly.
%
%   CLABEL(CS) or CLABEL(CS,V) places
%   contour labels as above, except that the labels are drawn as
%   plus signs on the contour with a nearby height value.
%
%   H = CLABEL(...) returns handles to the TEXT (and possibly LINE)
%   objects created.  The UserData property of the TEXT objects contain
%   the height value for each label.
%
%   CLABEL(...,'text property',property_value,...) allows arbitrary 
%   TEXT property/value pairs to specified for the label strings.
%
%   One special property ('LabelSpacing') is also available to specify 
%   the spacing between labels (in points). This defaults to 144, or
%   2 inches. 
%
%   Uses code by R. Pawlowicz to handle inline contour labels.
%
%   Example
%      subplot(1,3,1), [cs,h] = contour(peaks); clabel(cs,h,'labelspacing',72)
%      subplot(1,3,2), cs = contour(peaks); clabel(cs)
%      subplot(1,3,3), [cs,h] = contour(peaks); 
%      clabel(cs,h,'fontsize',15,'color','r','rotation',0)
%
%   See also CONTOUR, CONTOUR3, CONTOURF.

%   Thanks to R. Pawlowicz (IOS) rich@ios.bc.ca for the algorithm used
%   in 'inline_labels' so that clabel can produce inline labeling.

%   Copyright 1984-2002 The MathWorks, Inc. 
%   $Revision: 5.39 $  $Date: 2002/06/05 17:52:51 $

% Modified by R Pawlowicz to allow for text properties as in 
% extcontour code 14/5/97
% 28/10/97 - modified to work in map contouring
%  9/01/98 - improved calculation of gaps for line labels
%   Fix by Eric Firing, efiring@soest.hawaii.edu, 4/97, to
%   make the rotation angles correct when XDir and/or YDir are
%   reverse.

% Modified by J. Luis to add the text handles to appdata 'LabelHands' and removed 'manual' case

if nargin == 0
    error('Not enough input arguments.')
end
if min(size(cs)) > 2
    error('First input must be a valid contour description matrix.')
end
threeD = IsThreeD(gca);

if nargin == 1,
  h = plus_labels(threeD,cs);
else
  if ~isempty(varargin{1}(1)) && ishandle(varargin{1}(1)) && ...
    (strcmp(get(varargin{1}(1),'type'),'line') || strcmp(get(varargin{1}(1),'type'),'patch')),
    h = inline_labels(cs,varargin{:});
  else
    h = plus_labels(threeD,cs,varargin{:});
  end
end

if nargout>0, hh = h; end
if ~ishold, 
  if threeD, view(3), else view(2), end
end

%--------------------------------------------------------------
function H = inline_labels(CS,h,varargin)
%
% Draw the labels along the contours and rotated to match the local slope.
%

% To open up space in the contours, we rely on the order in which
% the handles h are created in CONTOUR3.  If CONTOUR3 changes you
% might need to change the algorithm below.

% Author: R. Pawlowicz IOS rich@ios.bc.ca
%         12/12/94
%         changes - R. Pawlowicz 14/5/97 - small bug in "that ole' 
%         matlab magic" fixed, also another in manual selection
%         of locations.

v=[]; inargs=zeros(1,length(varargin));
axHand = gca;   figHand = get(axHand,'Parent');

if nargin>=3 && ~ischar(varargin{1}),
  v = varargin{1};
  inargs(1)=1;
end;

lab_int=72*2;  % label interval (points)

for k=find(inargs==0),
 if strncmpi(varargin{k},'lab',3),
   inargs([k k+1])=1;
   lab_int=varargin{k+1};
 end;
end;

varargin(find(inargs))=[]; 

if strcmp(get(h(1),'type'),'patch') && ~strcmp(get(h(1),'facecolor'),'none'),
  isfilled = 1;
else
  isfilled = 0;
end

%% EF 4/97
if (strcmp(get(axHand, 'XDir'), 'reverse')), XDir = -1; else XDir = 1; end
if (strcmp(get(axHand, 'YDir'), 'reverse')), YDir = -1; else YDir = 1; end
%%

% Compute scaling to make sure printed output looks OK. We have to go via
% the figure's 'paperposition', rather than the the absolute units of the
% axes 'position' since those would be absolute only if we kept the 'units'
% property in some absolute units (like 'points') rather than the default
% 'normalized'.

UN=get(axHand,'units');
if (UN(1:3)=='nor'),
  UN=get(figHand,'paperunits');
  set(figHand,'paperunits','points');
  PA=get(figHand,'paperposition');
  set(figHand,'paperunits',UN);
  PA=PA.*[get(axHand,'position')];
else
  set(axHand,'units','points');
  PA=get(axHand,'pos');
  set(axHand,'units',UN); 
end

% Find beginning of all lines

lCS=size(CS,2);

if ~isempty(get(axHand,'children')),
  XL=get(axHand,'xlim');    YL=get(axHand,'ylim');
else
  iL=[];
  k=1;
  XL=[Inf -Inf];
  YL=[Inf -Inf];
  while (k<lCS),
    x=CS(1,k+(1:CS(2,k)));
    y=CS(2,k+(1:CS(2,k)));
    XL=[ min([XL(1),x]) max([XL(2),x]) ];
    YL=[ min([YL(1),y]) max([YL(2),y]) ]; 
    iL=[iL k];
    k=k+CS(2,k)+1;
  end;
  set(gca,'xlim',XL,'ylim',YL);
end;

Aspx=PA(3)/diff(XL);  % To convert data coordinates to paper (we need
                      % to do this
Aspy=PA(4)/diff(YL);  % to get the gaps for text the correct size)

H=[];

% Set up a dummy text object from which you can get text extent info
H1=text(XL(1),YL(1),'dummyarg','units','points','visible','off',varargin{:});

% Get labels all at once to get the length of the longest string.  
% This allows us to call extent only once, thus speeding up this routine
labels = getlabels(CS);
% Get the size of the label
set(H1,'string',repmat('9',1,size(labels,2)),'visible','on',varargin{:})
EX=get(H1,'extent'); set(H1,'visible','off')
len_lab=EX(3)/2;

ii=1; k = 0;
while (ii<lCS),
  k = k+1;

  if ~isfilled && k>length(h), error('Not enough contour handles.'); end

  l=CS(2,ii);
  x=CS(1,ii+(1:l));
  y=CS(2,ii+(1:l));

  lvl=CS(1,ii);

  %RP - get rid of all blanks in label
  lab=labels(k,labels(k,:)~=' ');
  %RP - scale label length by string size instead of a fixed length
  len_lab=EX(3)/2*length(lab)/size(labels,2);
  
  % RP28/10/97 - Contouring sometimes returns x vectors with 
  % NaN in them - we want to handle this case!
  sx=x*Aspx;
  sy=y*Aspy;
  d=[0 sqrt(diff(sx).^2 +diff(sy).^2)]; 
  % Determine the location of the NaN separated sections
  section = cumsum(isnan(d));
  d(isnan(d))=0;
  d=cumsum(d);
  
  len_contour = max(0,d(l)-3*len_lab);
  slop = (len_contour - floor(len_contour/lab_int)*lab_int);
  start = 1.5*len_lab + max(len_lab,slop)*rands(1); % Randomize start
  psn=[start:lab_int:d(l)-1.5*len_lab];
  lp=size(psn,2);
  
  if (lp>0) & isfinite(lvl)  && (isempty(v) | any(abs(lvl-v)/max(eps+abs(v)) < .00001)),
 
    Ic=sum( d(ones(1,lp),:)' < psn(ones(1,l),:),1 );
    Il=sum( d(ones(1,lp),:)' <= psn(ones(1,l),:)-len_lab,1 );
    Ir=sum( d(ones(1,lp),:)' < psn(ones(1,l),:)+len_lab,1 );
 
    % Check for and handle out of range values
    out = (Ir < 1 | Ir > length(d)-1) | ...
          (Il < 1 | Il > length(d)-1) | ...
          (Ic < 1 | Ic > length(d)-1);
    Ir = max(1,min(Ir,length(d)-1));
    Il = max(1,min(Il,length(d)-1));
    Ic = max(1,min(Ic,length(d)-1));

    % For out of range values, don't remove datapoints under label
    Il(out) = Ic(out);
    Ir(out) = Ic(out);
    
    % Remove label if it isn't in the same section
    bad = (section(Il) ~= section(Ir));
    Il(bad) = [];
    Ir(bad) = [];
    Ic(bad) = [];
    psn(:,bad) = [];
    out(bad) = [];
    lp = length(Il);
    in = ~out;
    
    if ~isempty(Il)
      % Endpoints of text in data coordinates
      wl=(d(Il+1)-psn+len_lab.*in)./(d(Il+1)-d(Il));
      wr=(psn-len_lab.*in-d(Il)  )./(d(Il+1)-d(Il));
      xl=x(Il).*wl+x(Il+1).*wr;
      yl=y(Il).*wl+y(Il+1).*wr;
   
      wl=(d(Ir+1)-psn-len_lab.*in)./(d(Ir+1)-d(Ir));
      wr=(psn+len_lab.*in-d(Ir)  )./(d(Ir+1)-d(Ir));
      xr=x(Ir).*wl+x(Ir+1).*wr;
      yr=y(Ir).*wl+y(Ir+1).*wr;
     
      trot=atan2( (yr-yl)*YDir*Aspy, (xr-xl)*XDir*Aspx )*180/pi; %% EF 4/97
      backang=abs(trot)>90;
      trot(backang)=trot(backang)+180;
      
      % Text location in data coordinates 
      wl=(d(Ic+1)-psn)./(d(Ic+1)-d(Ic));
      wr=(psn-d(Ic)  )./(d(Ic+1)-d(Ic));    
      xc=x(Ic).*wl+x(Ic+1).*wr;
      yc=y(Ic).*wl+y(Ic+1).*wr;
  
      % Shift label over a little if in a curvy area
      shiftfrac=.5;
      
      xc=xc*(1-shiftfrac)+(xr+xl)/2*shiftfrac;
      yc=yc*(1-shiftfrac)+(yr+yl)/2*shiftfrac;
      
      % Remove data points under the label...
      % First, find endpoint locations as distances along lines
    
      dr=d(Ir)+sqrt( ((xr-x(Ir))*Aspx).^2 + ((yr-y(Ir))*Aspy).^2 );
      dl=d(Il)+sqrt( ((xl-x(Il))*Aspx).^2 + ((yl-y(Il))*Aspy).^2 );
    
      % Now, remove the data points in those gaps using that
      % ole' MATLAB magic. We use the sparse array stuff instead of
      % something like:  
      %        f1=zeros(1,l); f1(Il)=ones(1,lp);
      % because the sparse functions will sum into repeated indices,
      % rather than just take the last accessed element - compare
      %   x=[0 0 0]; x([2 2])=[1 1]
      % with
      %   x=full(sparse([1 1],[2 2],[1 1],1,3))
      % (bug fix in original code 18/7/95 - RP)
  
      f1=full(sparse(ones(1,lp),Il,ones(1,lp),1,l));
      f2=full(sparse(ones(1,lp),Ir,ones(1,lp),1,l));
      irem=find(cumsum(f1)-cumsum(f2))+1;
      x(irem)=[];
      y(irem)=[];
      d(irem)=[];
      l=l-size(irem,2);
      
      % Put the points in the correct order...
      
      xf=[x(1:l),xl,repmat(NaN,size(xc)),xr];
      yf=[y(1:l),yl,yc,yr];
  
      [df,If]=sort([d(1:l),dl,psn,dr]);
    
      % ...and draw.
      % Here's where we assume the order of the h(k).  
      %
  
      z = get(h(k),'zdata');
      if ~isfilled, % Only modify lines or patches if unfilled
        set(h(k),'xdata',[xf(If) NaN],'ydata',[yf(If) NaN])
  
        % Handle contour3 case (z won't be empty).
        if ~isempty(z), 
          set(h(k),'zdata',[]) % Work around for bug in face generation
          % Set z to a constant while preserving the location of NaN's
          set(h(k),'zdata',z(1)+0*get(h(k),'xdata'))
        end
  
        if strcmp(get(h(k),'type'),'patch')
          set(h(k),'cdata',lvl+[0*xf(If) nan])
        end
      end
  
      H_tmp = [];
      for jj=1:lp,
        % Handle contour3 case (z won't be empty).
        if ~isempty(z),
          H_tmp = [H_tmp; text(xc(jj),yc(jj),z(1),lab,'rotation',trot(jj), ...
               'verticalAlignment','middle','horizontalAlignment','center',...
               'clipping','on','userdata',lvl,varargin{:})];
        else
          H_tmp = [H_tmp; text(xc(jj),yc(jj),lab,'rotation',trot(jj), ...
               'verticalAlignment','middle','horizontalAlignment','center',...
               'clipping','on','userdata',lvl,varargin{:})];       
        end
      end
      H = [H; H_tmp];

      if (~isempty(H_tmp))
        setappdata(h(k),'LabelHands',H_tmp)         % JL
      end
      
    end % ~isempty(Il)
  else
    if ~isfilled, % Only modify lines or patches if unfilled
      %
      % Here's another place where we assume the order of the h(k)
      %
      set(h(k),'xdata',[x NaN],'ydata',[y NaN])

      % Handle contour3 case (z won't be empty).
      z = get(h(k),'zdata');
      if ~isempty(z), 
        set(h(k),'zdata',[]) % Work around for bug in face generation
        % Set z to a constant while preserving the location of NaN's
        set(h(k),'zdata',z(1)+0*get(h(k),'xdata'))
      end
      if strcmp(get(h(k),'type'),'patch')
         set(h(k),'cdata',lvl+[0*x nan])
      end
    end
  end;
  
  ii=ii+1+CS(2,ii);
end;
  
% delete dummy string
delete(H1);
%-------------------------------------------------------

%-------------------------------------------------------
function h = plus_labels(threeD,cs,varargin)
%
% Draw the labels as plus symbols next to text (v4 compatible)
%

%    RP - 14/5/97
%    Clay M. Thompson 6-7-96
%    Charles R. Denham, MathWorks, 1988, 1989, 1990.
cax = gca;
choice = 0;

if nargin > 2
    choice = 1;
    v = sort(varargin{1}(:));
    varargin(1)=[];
end

[mcs, ncs] = size(cs);

% Find range of levels.
k = 1; i = 1;
while k <= ncs
   levels(i) = cs(1,k);
   i = i + 1;
   k = k + cs(2,k) + 1;
end
crange = max(abs(levels));
cdelta = abs(diff(levels)); 
cdelta = min(cdelta(cdelta > eps))/max(eps,crange); % Minimum significant change
if isempty(cdelta), cdelta = 0; end

k = 0; n = 0; flip = 0; h = [];

while (1)       % Select a labeling point randomly.

	k = k + n + 1; if k > ncs, break, end
	c = cs(1, k); n = cs(2, k);
	if choice
        f = find(abs(c-v)/max(eps+abs(v)) < .00001);
        okay = length(f) > 0;
	else
        okay = 1;
	end
	if okay
        r = rands(1);
        j = fix(r.* (n - 1)) + 1;
        if flip, j = n - j; end
        flip = ~flip;
        if n == 1    % if there is only one point
            xx = cs(1, j+k); yy = cs(2, j+k);
        else
            x1 = cs(1, j+k); y1 = cs(2, j+k);
            x2 = cs(1, j+k+1); y2 = cs(2, j+k+1);
            xx = (x1 + x2) ./ 2; yy = (y1 + y2) ./ 2;  % Test was here; removed.
        end
	end

% Label the point.

   if okay
      % Set tiny labels to zero.
      if abs(c) <= 10*eps*crange, c = 0; end
      % Determine format string number of digits
      if cdelta > 0, 
        ndigits = max(3,ceil(-log10(cdelta)));
      else
        ndigits = 3;
      end
      s = num2str(c,ndigits);
      hl = line('xdata',xx,'ydata',yy,'marker','+','erasemode','none');
      ht = text(xx, yy, s, 'verticalalignment', 'bottom', ...
             'horizontalalignment', 'left','erasemode','none', ...
             'clipping','on','userdata',c,varargin{:});
      if threeD, 
        set(hl,'zdata',c);
        set(ht,'position',[xx yy c]);
      end
      h = [h;hl];
      h = [h;ht];
   end
end
%-------------------------------------------------------

%-------------------------------------------------------
function labels = getlabels(CS)
%GETLABELS Get contour labels
v = []; i =1;
while i < size(CS,2), 
  v = [v,CS(1,i)];
  i = i+CS(2,i)+1;
end
labels = num2str(v');

%---------------------------------------------------
function threeD = IsThreeD(cax)
%ISTHREED  True for a contour3 plot
hp = findobj(cax,'type','patch');
if isempty(hp), hp = findobj(gca,'type','line'); end
if ~isempty(hp),
  % Assume a contour3 plot if z data not empty
  threeD = ~isempty(get(hp(1),'zdata'));
else
  threeD = 0;
end

%----------------------------------------
function r = rands(n)
%RANDS Stateless rand
%   R = RANDS(N) is the same as calling RAND(INP) except
%   that the state of the random number generator is unaffected.

Currstate=rand('state');
rand('state',sum(100*clock));
r = rand(n);
rand('state',Currstate); % resets to original state.
