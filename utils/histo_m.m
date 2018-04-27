function  varargout = histo_m(opt,varargin)
% This hacked hist functions work also with single y
	switch opt
		case 'hist'
			[varargout{1:nargout}] = hist(varargin{:});
		case 'bar'
			[varargout{1}] = bar(varargin{:});
	end

% -------------------------------------------------------------------------
function [no,xo] = hist(y,x,min_max, opt)
%   N = HIST(Y) bins the elements of Y into 10 equally spaced containers
%   and returns the number of elements in each container.  If Y is a
%   matrix, HIST works down the columns.
%
%   N = HIST(Y,M), where M is a scalar, uses M bins.
%
%   N = HIST(Y,X), where X is a vector, returns the distribution of Y
%   among bins with centers specified by X.  Note: Use HISTC if it is
%   more natural to specify bin edges instead.
%
%   [N,X] = HIST(...) also returns the position of the bin centers in X.
%
%   HIST(...) without output arguments produces a histogram bar plot of the results.
%   H = HIST(..., 'whatever') produces a histogram bar plot and returns a vector of patch handles

%   J.N. Little 2-06-86
%   Copyright 1984-2002 The MathWorks, Inc. 
%   J. Luis
%   Revised 13-6-05. Make it work with singles, added the min_max optional argument and some cleaning
%	        30-3-17. Add the 'OPT' option to return the bar handles for the histogram case.

	n_in = nargin;		ret_hand = false;
	if (n_in == 0)
		error('Requires one or two input arguments.')
	end
	if (n_in == 1)
		x = 10;
	elseif (n_in == 2 && isa(y, 'char'))
		x = 10;		ret_hand = true;	n_in = n_in - 1;
	elseif (n_in == 2 && numel(x) == 1)
		n_in = n_in - 1;
	elseif (n_in == 3 && isa(min_max, 'char'))
		ret_hand = true;	n_in = n_in - 1;
	else
		ret_hand = true;	n_in = n_in - 1;
	end
	
	if (n_in == 3)            % Min & Max of y was transmited in argument
		miny = min_max(1);
		maxy = min_max(2);
		got_min_max = 1;
		if (isempty(x)),    x = 10;     end
	else
		got_min_max = 0;
	end
	if min(size(y))==1, y = y(:); end

	if (~got_min_max)
		miny = min(y(:));
		maxy = max(y(:));
		if (~isa(miny,'double')),   miny = double(miny);    maxy = double(maxy);    end
	end
	if (length(x) == 1)
		if (miny == maxy)
			miny = miny - floor(x/2) - 0.5; 
			maxy = maxy + ceil(x/2) - 0.5;
		end
		binwidth = (maxy - miny) ./ x;
		xx = miny + binwidth*(0:x);
		xx(length(xx)) = maxy;
		x = xx(1:length(xx)-1) + binwidth/2;
	else
		xx = x(:)';
		binwidth = [diff(xx) 0];
		xx = [xx(1)-binwidth(1)/2 xx+binwidth/2];
		xx(1) = min(xx(1),miny);
		xx(end) = max(xx(end),maxy);
	end
	% Shift bins so the internal is ( ] instead of [ ).
	%xx = full(real(xx)); y = full(real(y)); % For compatibility
	bins = xx + max(eps,eps*abs(xx));
	nn = histc(y,[-inf bins],1);

	% Combine first bin with 2nd bin and last bin with next to last bin
	nn(2,:) = nn(2,:)+nn(1,:);
	nn(end-1,:) = nn(end-1,:)+nn(end,:);
	nn = nn(2:end-1,:);

	if (nargout == 0 || ret_hand)
		no = bar(x,nn,'hist');			% If bar handles required, return them in 'no'
		set(gca, 'xlim',[miny maxy])
	else
		if (min(size(y)) == 1) % Return row vectors if possible.
			no = nn';   xo = x;
		else
			no = nn;    xo = x';
		end
	end

%--------------------------------------------------------------------------------
function xo = bar(varargin)
%BAR Bar graph.
%    BAR(X,Y) draws the columns of the M-by-N matrix Y as M groups of N
%    vertical bars.  The vector X must be monotonically increasing or
%    decreasing.
%
%    BAR(Y) uses the default value of X=1:M.  For vector inputs, BAR(X,Y)
%    or BAR(Y) draws LENGTH(Y) bars.  The colors are set by the colormap.
%
%    BAR(X,Y,WIDTH) or BAR(Y,WIDTH) specifies the width of the bars. Values
%    of WIDTH > 1, produce overlapped bars.  The default value is WIDTH=0.8
%
%    BAR(...,'grouped') produces the default vertical grouped bar chart.
%    BAR(...,'stacked') produces a vertical stacked bar chart.
%    BAR(...,LINESPEC) uses the line color specified (one of 'rgbymckw').
%
%    H = BAR(...) returns a vector of patch handles.
%
%    Use SHADING FACETED to put edges on the bars.  Use SHADING FLAT to
%    turn them off.
%
%    Examples: subplot(3,1,1), bar(rand(10,5),'stacked'), colormap(cool)
%              subplot(3,1,2), bar(0:.25:1,rand(5),1)
%              subplot(3,1,3), bar(rand(2,3),.75,'grouped')

%    C.B Moler 2-06-86
%    Copyright 1984-2002 The MathWorks, Inc. 

	[msg,x,y,xx,yy,linetype,plottype,barwidth,equal] = makebars(varargin{:});
	if (~isempty(msg)),		error(msg),		end

	state = uisuspend_j(gcf);				% Remember initial figure state
	hImg = findobj(gcf, 'type', 'image');	% To be used if a Mirone figure is the gcf
	if (~isempty(hImg))
		cax = newplot(hImg);
	else									% Used with a new or still empty figure 
		cax = gca;
	end
	next = lower(get(cax,'NextPlot'));
	hold_state = ishold;
	edgec = get(gcf,'defaultaxesxcolor');
	facec = 'flat';
	h = zeros(1, size(xx,2)); 
	cc = ones(size(xx,1),1);
	if ~isempty(linetype), facec = linetype; end
	for i=1:size(xx,2)
		numBars = (size(xx,1)-1)/5;
		f = 1:(numBars*5);
		f(1:5:(numBars*5)) = [];
		f = reshape(f, 4, numBars);
		f = f';
		v = [xx(:,i) yy(:,i)];
		h(i) = patch('faces', f, 'vertices', v, 'cdata', i*cc, 'FaceColor',facec,'EdgeColor',edgec);
	end
	%if length(h)==1, set(cax,'clim',[1 2]), end
	set(h,'FaceColor',[.8 .8 .8])
	if (~equal)
		hold on;
		plot(x(:,1),zeros(size(x,1),1),'*')
	end
	if (~hold_state) 
		% Set ticks if less than 16 integers
		if (all(all(floor(x)==x)) && (size(x,1) < 16))
			set(cax,'xtick',x(:,1))
		end
		hold off,  set(cax,'NextPlot',next); %view(2);
		set(cax,'Layer','Bottom','box','on')
		% Turn off edges when they start to overwhelm the colors
		if (size(xx,2)*numBars > 150)
			set(h,{'edgecolor'},get(h,{'facecolor'}));
		end
	end
	uirestore_j(state, 'nochildren');         % Restore the figure's initial state
	if (nargout==1),   xo = h;     end

%--------------------------------------------------------------------------------
function [msg,x,y,xx,yy,linetype,plottype,barwidth,arg8] = makebars(varargin)
%MAKEBARS Make data for bar charts.
%   [MSG,X,Y,XX,YY,LINETYPE,PLOTTYPE,BARWIDTH,EQUAL] = MAKEBARS(X,Y) 
%   returns X and Y, the original data values, XX and YY, the data
%   values formatted to be plotted by one of the BARxx m-files (BAR,
%   BARH, BAR3, BAR3H).
%   
%   LINETYPE returns the color desired for the plot.
%   PLOTTYPE determines whether the plot will be grouped (PLOTTYPE=0),
%   stacked (PLOTTYPE=1), or detached (PLOTTYPE=2--only for 3-D plots).
%   BARWIDTH has the bar width (normalized to one).
%
%   [MSG,X,Y,XX,YY,LINETYPE,PLOTTYPE,BARWIDTH,ZZ] = MAKEBARS(X,Y)
%   does the same as above, except for the final parameter.
%   ZZ is the z-axis data for 3-D plots; used in BAR3 and BAR3H.
%
%   [...] = MAKEBARS(X,Y,WIDTH) or MAKEBARS(Y,WIDTH) returns the
%   specified width given in WIDTH.  The default is 0.8.
%   [...] = MAKEBARS(...,'grouped') returns the data in the form
%   so that the information will be plotted in groups.
%   [...] = MAKEBARS(...,'detached') {3-D only} returns the data
%   such that the information will be plotted detached.
%   [...] = MAKEBARS(...,'stacked') returns the data such that the
%   information will be plotted in stacked form.
%   [...] = MAKEBARS(...,'hist') creates centered bars that touch.
%   [...] = MAKEBARS(...,'histc') creates bars that touch edges.
%   EQUAL is true if spacing of the data is equal; false otherwise.
%   EQUAL is always true except for 'hist' and 'histc' plottypes.

	% Initialize everything
	x = []; y=[]; xx=[]; yy=[]; arg8 = [];

	barwidth = 0.8;		% Normalized width of bar.
	groupwidth = 0.8;	% Normalized width of groups.
	linetype = [];		% Assume linetype is not specified
	ishist = 0;			% Assume this isn't a histogram

	nin = nargin;

	if strcmp(varargin{nin},'3')
		threeD = 1;
		nin = nin - 1; 
		plottype = 2; % Detached plot default
	else
		threeD = 0;
		plottype = 0; % Grouped plot default
	end

	if ischar(varargin{nin})	% Try to parse this string as a color
		[ls,co,mark,msg] = colstyle(varargin{nin});
		if isempty(msg), linetype = co; nin = nin - 1; end
	end

	if ischar(varargin{nin})	% Process 'grouped','stacked' or 'detached' string.
		kind = [lower(varargin{nin}) 'g'];
		if kind(1)=='g'			% grouped
			plottype = 0;
		elseif kind(1)=='s'		% stacked
			plottype = 1;
		elseif threeD && kind(1) == 'd'	% detached
			plottype = 2;
		elseif kind(1)=='h'		% histogram
			if strcmpi(varargin{nin},'histc')
				ishist = -1;	barwidth = 1;
			else
				ishist = 1;		barwidth = 1;
			end
		else
			msg = sprintf('Unrecognized option "%s".',varargin{nin});
			return
		end
		nin = nin-1;
	end

	% Parse input arguments.
	if (nin>1) && (length(varargin{nin})==1) && (length(varargin{nin-1})>1)
		% If last argument is a scalar and next to last isn't then last argument must be the barwidth.
		barwidth = varargin{nin};
		[msg,x,y] = xychk(varargin{1:nin-1});
	else
		[msg,x,y] = xychk(varargin{1:nin});
	end

	% Make sure x is monotonic
	[x,ndx] = sort(x);
	if (min(size(y)) == 1),		y = y(ndx);
	else,						y = y(ndx,:);
	end

	if ~isempty(msg), return, end

	% Expand x to be the same size as y.
	if min(size(y))>1, x = x(:,ones(size(y,2),1)); end

	% Make sure vector arguments are columns
	if min(size(y))==1, x = x(:); y = y(:); end

	[n,m] = size(y);

	% Remove y values to 0 where y is NaN;
	k = find(isnan(y));
	if ~isempty(k), y(k) = 0; end

	if threeD
		z = y; % Make variables consistent with 3-D bar graph
		y = x;

		if (m == 1 || plottype ~= 0)
			groupwidth = 1;
		else
			groupwidth = min(groupwidth,m/(m+1.5));
		end

		nn = 6*n; mm = 4*m;
		zz = zeros(nn,mm);
		yy = zz;
		xx = zeros(nn,4);

		% Define xx
		xx(:,3:4) = 1;
		if (plottype == 0)
			xx = (xx-0.5)*barwidth*groupwidth/m;
		else
			xx = (xx-0.5)*barwidth*groupwidth;
		end

		% Define heights
		zz(2:6:nn,2:4:mm) = z;
		zz(2:6:nn,3:4:mm) = z;
		zz(3:6:nn,2:4:mm) = z;
		zz(3:6:nn,3:4:mm) = z;

		if (plottype == 1 && m > 1)		% Stacked
			z = cumsum(z.').';
			zz = zz + [zeros(nn,4) z(ones(6,1)*(1:n),ones(4,1)*(1:m-1))];
		end     

		if length(y)==1 || max(diff(y))==0
			equal = []; % Special case
		else
			equal = max(abs(diff(diff(y)))) <= max(max(abs(y)))*sqrt(eps);
		end

		%       
		% Determine beginning of bars (t) and bar spacing (delta)
		%       
		if isempty(equal)		% Scalar case and special case
			delta = ones(size(y));
			t = y - 0.5*delta;
			x = 1;
		elseif equal
			if plottype ~= 0		% Stacked or detached
				delta = ones(n,1) * (max(y) - min(y)) * groupwidth / (n-1);
				t = y - 0.5*delta;
			else % grouped
				delta = ones(n,1) * (max(y) - min(y)) * groupwidth / (n-1) / m ;
				t = y - 0.5*delta*m + (ones(n,1)*(0:m-1)).*delta;
			end 
			if (plottype == 2),		x = 1:m;
			else,					x = ones(1,m);
			end
		else % Spacing is unequal.
			if ishist==1 % Width of bin is average of difference between points.
				dy = diff(y); dy = (dy([1 1:end],:)+dy([1:end end],:))/2;
				t = [(y(1:end-1,:)+y(2:end,:))/2;y(end,:)+(y(end,:)-y(end-1,:))/2]-dy/2;
			elseif ishist==-1 % Width of bin is difference between edges
				dy = [diff(y,1);zeros(1,size(y,2))];
				t = y + dy*groupwidth/2;
			else % Width of bin is minimum of difference between points.
				dy = ones(n,1)*min(diff(y));
				t = y;
			end 
			if plottype ~= 0		% Stacked or detached
				delta = dy * groupwidth;
				t = t - delta/2;
			else % Grouped
				delta = dy * groupwidth / m ;
				t = t - delta/2*m + (ones(n,1)*(0:m-1)).*delta;
			end
			if (plottype == 2)
				x = 1:m;
			else
				x = ones(1,m);
			end
		end     
		t = t(:,ones(4,1)*(1:m));
		delta = delta(:,ones(4,1)*(1:m));
		yy(1:6:nn,:) = t + (1-barwidth)/2.*delta;
		yy(2:6:nn,:) = t + (1-barwidth)/2.*delta;
		yy(3:6:nn,:) = t + (1+barwidth)/2.*delta;
		yy(4:6:nn,:) = t + (1+barwidth)/2.*delta;
		yy(5:6:nn,:) = t + (1-barwidth)/2.*delta;

		% Insert NaN's so distinct bars are drawn
		ii1 = [(1:6:nn) (4:6:nn) (5:6:nn) (6:6:nn)];
		ii2 = 6:6:nn;
		xx(ii1,[1 4]) = nan;
		xx(ii2,[2 3]) = nan;
		yy(ii1,[(1:4:mm) (4:4:mm)]) = nan;
		yy(ii2,[(2:4:mm) (3:4:mm)]) = nan;
		zz(ii1,[(1:4:mm) (4:4:mm)]) = nan;
		zz(ii2,[(2:4:mm) (3:4:mm)]) = nan;
		arg8 = zz;
	else
		nn = 5*n;
		yy = zeros(nn+1,m);
		xx = yy;
		if plottype && (m>1)
			yc = cumsum(y')';
			ys = [zeros(n,1),yc(:,1:end-1)];     

			yy(2:5:nn,:) = ys;
			yy(3:5:nn,:) = yc;
			yy(4:5:nn,:) = yc;
			yy(5:5:nn,:) = ys;
		else
			yy(3:5:nn,:) = y;
			yy(4:5:nn,:) = y;
		end

		if (m==1 || plottype),	groupwidth = 1;
		else,					groupwidth = min(groupwidth,m/(m+1.5));
		end

		equal = max(abs(diff(diff(x)))) <= max(max(abs(x)))*sqrt(eps);
		if length(x)==1 || max(diff(x))==0, equal=[]; end % Special case 

		%
		% Determine beginning of bars (t) and bar spacing (delta)
		%
		dt = [];
		t = x;
		if isempty(equal)				% Scalar case and special case
			delta = 1;
			if (ishist~=-1), dt = 0.5*delta; end
		elseif equal
			if plottype
				delta = ones(n,1) * (max(x) - min(x)) * groupwidth / (n-1);
				if (ishist~=-1), dt = 0.5*delta; end
			else
				delta = ones(n,1) * (max(x) - min(x)) * groupwidth / (n-1) / m ;
				t = x + (ones(n,1)*(0:m-1)).*delta;
				if (ishist~=-1), dt = 0.5*delta*m; end
			end
		else % Spacing is unequal.
			if (ishist == 1)			% Width of bin is average of difference between points.
				dx = diff(x); dx = (dx([1 1:end],:)+dx([1:end end],:))/2;
				t = [(x(1:end-1,:)+x(2:end,:))/2;x(end,:)+(x(end,:)-x(end-1,:))/2]-dx/2;
			elseif (ishist == -1)		% Width of bin is difference between edges
				dx = diff(x,1); dx = dx([(1:end) end],:);
				%t = x + dx*groupwidth/2;
			else % Width of bin is minimum of difference between points.
				dx = ones(n,1)*min(diff(x));
			end
			if plottype
				delta = dx * groupwidth;
				if (ishist~=-1), dt = delta/2; end
			else
				delta = dx * groupwidth / m ;
				t = t + (ones(n,1)*(0:m-1)).*delta;
				if (ishist~=-1), dt = delta/2*m; end
			end
		end
		if (~isempty(dt)), t = t - dt; end

		xx(1:5:nn,:) = t;
		xx(2:5:nn,:) = t + (1-barwidth)/2.*delta;
		xx(3:5:nn,:) = t + (1-barwidth)/2.*delta;
		xx(4:5:nn,:) = t + (1+barwidth)/2.*delta;
		xx(5:5:nn,:) = t + (1+barwidth)/2.*delta;
		xx(nn+1,:) = xx(nn,:);

		if (ishist~=1) && (ishist~=-1), equal = 1; end
		arg8 = equal;
	end
