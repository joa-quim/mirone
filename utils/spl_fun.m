%#mex
function    varargout = spl_fun(opt,varargin)

switch opt
    case 'csaps'
        [output,p] = csaps(varargin{:});
        varargout{1} = output;  varargout{2} = p;
    case 'fnval'
        varargout{1} = fnval(varargin{:});
    case 'fnder'
        varargout{1} = fnder(varargin{:});
    case 'ppual'
        varargout{1} = ppual(varargin{:});
end

%----------------------------------------------------------------------------------
function t = brk2knt(breaks,mults)
%BRK2KNT Breaks with multiplicities into knots.
%   Copyright 1987-2003 C. de Boor and The MathWorks, Inc. 
%   $Revision: 1.12 $

s = sum(mults);
if s==0
   t = [];
else
   li = length(breaks);
      % make sure there is a multiplicity assigned to each break,
      % and drop any break whose assigned multiplicity is not positive.
   if length(mults)~=li, mults = repmat(mults(1),1,li); s = mults(1)*li;
   else
      fm = find(mults<=0);
      if ~isempty(fm), breaks(fm)=[]; mults(fm)=[]; li = length(breaks); end
   end
   mm = zeros(1,s);
   mm(cumsum([1 reshape(mults(1:li-1),1,li-1)])) = ones(1,li);
   t = breaks(cumsum(mm));
end

%----------------------------------------------------------------------------------
function [x,y,sizeval,w,origint,p,tolred] = chckxywp(x,y,nmin,w,p,adjtol)
%CHCKXYWP Check and adjust input for *AP*1 commands.
%   Copyright 1984-2003 The MathWorks, Inc.
%   $Revision: 1.3 $  $Date: 2003/02/13 22:22:18 $

% make sure X is a vector:
if iscell(x)||length(find(size(x)>1))>1
   error('SPLINES:CHCKXYWP:Xnotvec','X must be a vector.'), end

% make sure X is real:
if ~all(isreal(x))
   x = real(x);
   warning('SPLINES:CHCKXYWP:Xnotreal', 'Imaginary part of complex data sites ignored.')
end

% deal with NaN's and Inf's among the sites:
nanx = find(~isfinite(x));
if ~isempty(nanx)
   x(nanx) = [];
   warning('SPLINES:CHCKXYWP:NaNs', ...
           'All data points with NaN or Inf as their site will be ignored.')
end

n = length(x);
if nargin>2&&nmin>0, minn = nmin; else minn = 2; end
if n<minn
   error('SPLINES:CHCKXYWP:toofewpoints', ...
   'There should be at least %g data sites.',minn)
end

% re-sort, if needed, to ensure nondecreasing site sequence:
tosort = false;
if any(diff(x)<0), tosort = true; [x,ind] = sort(x); end

nstart = n+length(nanx);
% if Y is ND, reshape it to a matrix by combining all dimensions but the last:
sizeval = size(y);
yn = sizeval(end); sizeval(end) = []; yd = prod(sizeval);
if length(sizeval)>1
   y = reshape(y,yd,yn);
else
   % if Y happens to be a column matrix, of the same length as the original X,
   % then change Y to a row matrix
   if yn==1&&yd==nstart
      yn = yd; y = reshape(y,1,yn); yd = 1; sizeval = yd;
   end
end
y = y.'; x = reshape(x,n,1);

% make sure that sites, values and weights match in number:

if nargin>2&&~nmin % in this case we accept two more data values than
                   % sites, stripping off the first and last, and returning
		   % them separately, in W, for use in CSAPE1.
   switch yn
   case nstart+2, w = y([1 end],:); y([1 end],:) = [];
      if ~all(isfinite(w)),
         error('SPLINES:CHCKXYWP:InfY', 'Some of the end condition values fail to be finite.')
      end
   case nstart, w = [];
   otherwise
      error('SPLINES:CHCKXYWP:XdontmatchY', ...
           ['The number of sites, %g, does not match the number of', ' values, %g.'], nstart, yn)
   end
else
   if yn~=nstart
      error('SPLINES:CHCKXYWP:XdontmatchY', ...
           ['The number of sites, %g, does not match the number of', ' values, %g.'], nstart, yn)
   end
end

nonemptyw = nargin>3&&~isempty(w);
if nonemptyw
   if length(w)~=nstart
      error('SPLINES:CHCKXYWP:weightsdontmatchX', ...
       ['The number of weights, %g, does not match the number of', ' sites, %g.'], length(w), nstart)
   else
      w = reshape(w,1,nstart);
   end
end

roughnessw = exist('p','var')&&length(p)>1;
if roughnessw
   if tosort
      warning('SPLINES:CHCKXYWP:cantreorderrough', ...
           'Since data sites are not ordered, roughness weights are ignored.')
      p = p(1);
   else
      if length(p)~=nstart
         error('SPLINES:CHCKXYWP:rweightsdontmatchX', ...
	 ['The number of roughness weights is incompatible with the', ' number of sites, %g.'], nstart)
      end
   end
end

%%% remove values and error weights corresponding to nonfinite sites:
if ~isempty(nanx), y(nanx,:) = []; if nonemptyw, w(nanx) = []; end
   if roughnessw  % as a first approximation, simply ignore the
                  % specified weight to the left of any ignored point.
      p(max(nanx,2)) = [];
   end
end
if tosort, y = y(ind,:); if nonemptyw, w = w(ind); end, end

% deal with nonfinites among the values:
nany = find(sum(~isfinite(y),2));
if ~isempty(nany)
   y(nany,:) = []; x(nany) = []; if nonemptyw, w(nany) = []; end
   warning('SPLINES:CHCKXYWP:NaNs', ...
           'All data points with NaNs or Infs in their value will be ignored.')
   n = length(x);
   if n<minn
      error('SPLINES:CHCKXYWP:toofewX', 'There should be at least %g data sites.',minn)
   end
   if roughnessw  % as a first approximation, simply ignore the
                  % specified weight to the left of any ignored point.
      p(max(nany,2)) = [];
   end
end

if nargin==3&&nmin, return, end % for SPAPI, skip the averaging

if nargin>3&&isempty(w) %  use the trapezoidal rule weights:
   dx = diff(x);
   if any(dx), w = ([dx;0]+[0;dx]).'/2;
   else,       w = ones(1,n);
   end
   nonemptyw = ~nonemptyw;
end

tolred = 0;
if ~all(diff(x)) % conflate repeat sites, averaging the corresponding values
                 % and summing the corresponding weights
   mults = knt2mlt(x);
   for j=find(diff([mults;0])<0).'
		if nonemptyw
         temp = sum(w(j-mults(j):j));
			if nargin>5
				tolred = tolred + w(j-mults(j):j)*sum(y(j-mults(j):j,:).^2,2); 
			end
			y(j-mults(j),:) = (w(j-mults(j):j)*y(j-mults(j):j,:))/temp;
			w(j-mults(j)) = temp;
			if nargin>5
				tolred = tolred - temp*sum(y(j-mults(j),:).^2);
			end
		else
			y(j-mults(j),:) = mean(y(j-mults(j):j,:),1);
		end
   end
      
   repeats = find(mults);
   x(repeats) = []; y(repeats,:) = []; if nonemptyw, w(repeats) = []; end
   if roughnessw  % as a first approximation, simply ignore the
                  % specified weight to the left of any ignored point.
      p(max(repeats,2)) = [];
   end
   n = length(x);
   if n<minn, error('SPLINES:CHCKXYWP:toofewX', 'There should be at least %g data sites.',minn), end
end

if nargin<4, return, end


% remove all points corresponding to relatively small weights (since a
% (near-)zero weight in effect asks for the corresponding datum to be dis-
% regarded while, at the same time, leading to bad condition and even
% division by zero).
origint = []; % this will be set to x([1 end]).' in case the weight for an end
             % data point is near zero, hence the approximation is computed
             % without that endpoint.
if nonemptyw
   ignorep = find( w <= (1e-13)*max(abs(w)) );
   if ~isempty(ignorep)
      if ignorep(1)==1||ignorep(end)==n, origint = x([1 end]).'; end
      x(ignorep) = []; y(ignorep,:) = []; w(ignorep) = []; 
      if roughnessw
                     % as a first approximation, simply ignore the
                     % specified weight to the left of any ignored point.
         p(max(ignorep,2)) = [];
      end
      n = length(x);
      if n<minn
        error('SPLINES:CHCKXYWP:toofewX', ...
	     ['There should be at least %g data points with positive',...
	       ' weights.'],minn)
      end
   end
end

%----------------------------------------------------------------------------------
function [output,p] = csaps(x,y,p,xx,w)
%CSAPS Cubic smoothing spline.
%   Copyright 1987-2003 C. de Boor and The MathWorks, Inc.
%   $Revision: 1.24 $

if nargin<3||isempty(p), p = -1; end
if nargin<4, xx = []; end
if nargin<5, w = []; end

if iscell(x)     % we are to handle gridded data
   m = length(x);
   sizey = size(y);
   if length(sizey)<m
     error('SPLINES:CSAPS:toofewdims',...
          ['If X is a cell-array of length m, then Y must have', ...
            ' at least m dimensions.'])
   end

   if length(sizey)==m,  % grid values of a scalar-valued function
     if issparse(y), y = full(y); end 
     sizey = [1 sizey]; 
   end

   sizeval = sizey(1:end-m); sizey = [prod(sizeval), sizey(end-m+(1:m))];
   y = reshape(y, sizey); 
   
   if ~iscell(p)  % because of the possibility of weighted roughness measures
                  % must have P be a cell array in the multivariate case.
      if length(p)~=m, p = repmat(p(1),1,m); end
      p = num2cell(p);
   end
   if isempty(w), w = cell(1,m); end

   v = y; sizev = sizey;
   for i=m:-1:1   % carry out coordinatewise smoothing
      [cs,p{i}] = csaps1(x{i}, reshape(v,prod(sizev(1:m)),sizev(m+1)), ...
                  p{i}, [], w{i});
      [b,v,l,k] = ppbrk(cs);
      breaks{i} = b;
      sizev(m+1) = l*k; v = reshape(v,sizev);
      if m>1
         v = permute(v,[1,m+1,2:m]); sizev(2:m+1) = sizev([m+1,2:m]);
      end
   end
   % At this point, V contains the tensor-product pp coefficients;
   % It remains to make up the formal description:
   output = ppmak(breaks, v);
   if length(sizeval)>1, output = fnchg(output,'dz',sizeval); end
   if ~isempty(xx)
      output = fnval(output,xx);
   end
else             % we have univariate data
   [output,p] = csaps1(x,y,p,xx,w);
end
%varargout{1} = output;      varargout{2} = p;

%----------------------------------------------------------------------------
function [output,p] = csaps1(x,y,p,xx,w)
%CSAPS1 univariate cubic smoothing spline

n=length(x); if isempty(w), w = ones(1,n); end
[xi,yi,sizeval,w,origint,p] = chckxywp(x,y,2,w,p);
n = size(xi,1); yd = size(yi,2); dd = ones(1,yd);

dx = diff(xi); divdif = diff(yi)./dx(:,dd);
if n==2 % the smoothing spline is the straight line interpolant
   pp=ppmak(xi.',[divdif.' yi(1,:).'],yd); p = 1;
else % set up the linear system for solving for the 2nd derivatives at  xi .
     % this is taken from (XIV.6)ff of the `Practical Guide to Splines'
     % with the diagonal matrix D^2 there equal to diag(1/w) here.
     % Make use of sparsity of the system.

   dxol = dx;
   if length(p)>1 
      lam = p(2:end).'; p = p(1);
      dxol = dx./lam;
   end

   R = spdiags([dxol(2:n-1), 2*(dxol(2:n-1)+dxol(1:n-2)), dxol(1:n-2)],...
                                         -1:1, n-2,n-2);
   odx=1./dx;
   Qt = spdiags([odx(1:n-2), -(odx(2:n-1)+odx(1:n-2)), odx(2:n-1)], ...
                                                0:2, n-2,n);
   % solve for the 2nd derivatives
   W = spdiags(1./w(:),0,n,n);
   Qtw = Qt*spdiags(1./sqrt(w(:)),0,n,n);
   if p<0 % we are to determine an appropriate P
      QtWQ = Qtw*Qtw.'; p = 1/(1+trace(R)/(6*trace(QtWQ)));
          % note that the resulting  p  behaves like
          %   1/(1 + w_unit*x_unit^3/lambda_unit)
          % as a function of the various units chosen
      u=((6*(1-p))*QtWQ+p*R)\diff(divdif);
   else
      u=((6*(1-p))*(Qtw*Qtw.')+p*R)\diff(divdif);
   end
   clear divdif
   % ... and convert to pp form
   % Qt.'*u=diff([0;diff([0;u;0])./dx;0])
   yi = yi - ...
    (6*(1-p))*W*diff([zeros(1,yd)
                 diff([zeros(1,yd);u;zeros(1,yd)])./dx(:,dd)
                 zeros(1,yd)]);
   c3 = [zeros(1,yd);p*u;zeros(1,yd)];
   clear u
   c2=diff(yi)./dx(:,dd)-dxol(:,dd).*(2*c3(1:n-1,:)+c3(2:n,:));
   if exist('lam','var')
      dxtl = dx.*lam;
      pp=ppmak(xi.',...
      reshape([(diff(c3)./dxtl(:,dd)).',3*(c3(1:n-1,:)./lam(:,dd)).', ...
                                    c2.',yi(1:n-1,:).'], (n-1)*yd,4),yd);
   else
      pp=ppmak(xi.',...
      reshape([(diff(c3)./dx(:,dd)).',3*c3(1:n-1,:).',c2.',yi(1:n-1,:).'],...
                                                            (n-1)*yd,4),yd);
   end
end

if ~isempty(origint), pp = fnchg(pp,'int',origint); end
if length(sizeval)>1, pp = fnchg(pp,'dz',sizeval); end

if isempty(xx)
   output = pp;
else
   output = fnval(pp,xx);
end

%----------------------------------------------------------------------------------
function g = fn2fm(f,form,sconds)
%FN2FM Convert to specified form.
%   Copyright 1987-2003 C. de Boor and The MathWorks, Inc.
%   $Revision: 1.18 $

%%%%%%%%%%%%%%%%%%%%  array versions of univariate forms (from FNBRK)
%   case 10, out1 = 'ppform, univariate';
%   case 11, out1 = 'B-form, univariate';
%   case 12, out1 = 'BBform, univariate';
%   case 15, out1 = 'polynomial in Newton form';
%   case 25, out1 = 'thinplate spline';

if nargin==1 % set FORM to convert the form in F into its current version
   if isstruct(f)
      form = f.form;
   else
      switch f(1)
      case 10,  form = 'pp';
      case 11,  form = 'B-';
      case 12,  form = 'BB';
      case 15,  form = 'NP';
      case 25,  form = 'st-tp00';
      otherwise
         error('SPLINES:FN2FM:unknownform','Input is of unknown (function) form.')
      end
   end
else
   if length(form)<2
      error('SPLINES:FN2FM:unknownform',...
           ['The specified form, ''',form,''', is not (yet) recognized.']),
   end 
end

switch form(1:2)
case 'pp'
   if isstruct(f)
      switch f.form(1:2)
      case 'pp',          g = f;
      case 'rp',          g = f; g.form = 'pp'; g.dim = g.dim+1;
      case 'rB',          f.form = 'B-'; f.dim = f.dim+1; g = sp2pp(f);
      case {'B-','BB'},   g = sp2pp(f);
      end
   else
      switch f(1)
      case 10,      [b,c,l,k,d] = ppbrk(f); g = ppmak(b,c,d);
      case {11,12}, g = sp2pp(f);
      case 15
       error('SPLINES:FN2FM:nonnewtonhere',...
             'Use NP2PP directly to convert from Newton form to ppform.')
      end
   end
case {'sp','B-'}
   form = 'B-';
   if isstruct(f)
      switch f.form(1:2)
      case 'pp'
         if nargin<3, g = pp2sp(f);
         else         g = pp2sp(f,sconds);
         end
      case 'rp',      f.form = 'pp'; f.dim = f.dim+1; g = pp2sp(f);
      case 'rB',      g = f; g.form = 'B-'; g.dim = g.dim+1;
      case 'B-',      g = f;
      case 'BB',      g = pp2sp(sp2pp(f));
      end
   else
      switch f(1)
      case 10, g = pp2sp(f);
      case 11, [knots,coefs] = spbrk(f); g = spmak(knots,coefs);
      case 12, g = pp2sp(sp2pp(f));
      end
   end
case {'bb','BB'}
   form = 'BB';
   if isstruct(f)
      switch f.form(1:2)
      case 'pp', g = sp2bb(pp2sp(f));          
      case 'rp', f.form = 'pp'; f.dim = f.dim+1; g = sp2bb(pp2sp(f));
      case 'rB', f.form = 'B-'; f.dim = f.dim+1; g = sp2bb(f);
      case 'B-', g = sp2bb(f);
      case 'BB',  g = f;
      end
   else
      switch f(1)
      case 10, g = sp2bb(pp2sp(f));
      case 11, g = sp2bb(f);
      case 12, [knots,coefs] = spbrk(f); g = spmak(knots,coefs); g.form = 'BB';
      end
   end
case 'st' 
   if isstruct(f)
      switch f.form(1:2)
      case 'st', g = f;
      end
   else
      switch f(1)
      case 25, [ce,co] = stbrk(f); g = stmak(ce,co,'tp00');
      end
   end
case 'rp'  % switch first into ppform if need be
   if isstruct(f)
      switch f.form(1:2)
      case 'pp',      g = f;
      case 'rp',      g = f; g.dim = g.dim+1;
      case 'rB',      g = sp2pp(fn2fm(f,'B-'));
      case {'B-','BB','sp','bb'}, g = sp2pp(f);
      end
   else
      switch f(1)
      case 10,      [b,c,l,k,d] = ppbrk(f); g = ppmak(b,c,d);
      case {11,12}, g = sp2pp(f);
      end
   end
   if exist('g','var')
      g.form = 'rp'; 
      if length(g.dim)>1
         warning('SPLINES:FN2FM:noNDrats', ...
	 ['While the given function has values of size [', num2str(g.dim), ...
	  '], the function returned is ', num2str(prod(g.dim)-1), ...
	  '-vector valued.'])
	 g.dim = prod(g.dim);
      end
      g.dim = g.dim-1;
      if g.dim<1
         error('SPLINES:FN2FM:ratneedsmoredim',...
               'A rational spline must have more than one component.')
      end
   end
case 'rB'  % switch first into B-form if need be
   if isstruct(f)
      switch f.form(1:2)
      case 'pp'
         if nargin<3, g = pp2sp(f);
         else         g = pp2sp(f,sconds);
         end
      case 'rp',      g = fn2fm(fn2fm(f,'pp'),'B-');
      case 'rB',      g = f; g.dim = g.dim+1;
      case {'B-','BB','sp','bb'}, g = f;
      end   
   else 
      switch f(1)
      case 10, g = pp2sp(f);
      case 11, g = spmak(fnbrk(f,'knots'),fnbrk(f,'coefs'));
      case 12, g = pp2sp(sp2pp(f));
      end
   end
   if exist('g','var')
      g.form = 'rB';
      if length(g.dim)>1
         warning('SPLINES:FN2FM:noNDrats', ...
	 ['While the given function has values of size [', num2str(g.dim), ...
	  '], the function returned is ', num2str(prod(g.dim)-1), ...
	  '-vector valued.'])
	 g.dim = prod(g.dim);
      end
      g.dim = g.dim-1;
      if g.dim<1
         error('SPLINES:FN2FM:ratneedsmoredim',...
               'A rational spline must have more than one component.')
      end
   end
case 'NP' % at present, this form only exists in the old version
   if isstruct(f)
      switch f.form(1:2)
      case 'NP', g = f;
      end
   else
      switch f(1)
      case 15, g = f;
      end
   end
case 'MA'     % convert univariate spline to old, nonstructure ppform.
   if isstruct(f)&&~iscell(f.order)
      switch f.form(1:2)
      case 'pp' 
         g = [10 f.dim f.pieces f.breaks(:).' f.order f.coefs(:).'];
      case {'B-','BB'}
         f = sp2pp(f);
         g = [10 f.dim f.pieces f.breaks(:).' f.order f.coefs(:).'];
      otherwise
         error('SPLINES:FN2FM:notintoMA',...
               'Cannot convert the given function into MATLAB''s ppform.')
      end
   else
      g = f;
   end
case 'ol'    % convert univariate structured form to corresponding 
             % formerly used array form.
   if ~isstruct(f)
      error('SPLINES:FN2FM:olneedsstruct',...
            'When FORM=''ol'', F must be a structure.')
   else
      if length(f.order)>1
         error('SPLINES:FN2FM:nonstructneedsuni',...
              ['Reversion to array is only possible for',...
                ' univariate functions.'])
      else
         switch f.form(1:2)
         case 'pp' 
            g = [10 f.dim f.pieces f.breaks(:).' f.order f.coefs(:).'];
         case {'B-','BB'}
            g = [11 f.dim f.number f.coefs(:).' f.order f.knots(:).'];
            if (strcmp(f.form(1:2),'BB')),		g(1) = 12;	end
         otherwise
            error('SPLINES:FN2FM:unknownfn','Unknown function type encountered.')
         end
      end
   end
otherwise
   error('SPLINES:FN2FM:unknownfn',...
    ['The specified form, ''',form,''', is not (yet) recognized.']),
end

if ~exist('g','var')
   error('SPLINES:FN2FM:notintothatform',...
        ['Cannot convert the given function into ',form,'form.'])
end

%----------------------------------------------------------------------------------
function varargout = fnbrk(fn,varargin)
%FNBRK Name or part(s) of a form.
%   Copyright 1987-2003 C. de Boor and The MathWorks, Inc.
%   $Revision: 1.21 $ 

if nargin>1
   np = max(1,nargout); % FNBRK(FN,PART) may be part of an expression
   if np <= length(varargin)
      varargout = cell(1,np);
   else
      error('SPLINES:FNBRK:moreoutthanin', ...
            'Too many output arguments for the given input.')
   end
end 

if ~isstruct(fn)    % this branch should eventually be abandoned
   switch fn(1)
   % curves:
      case 10, fnform = 'ppform, univariate, array format';
      case 11, fnform = 'B-form, univariate, array format';
      case 12, fnform = 'BBform, univariate, array format';
      case 15, fnform = 'polynomial in Newton form';
   % surfaces:
      case 20, fnform = 'ppform, bivariate tensor product, array format';
      case 21, fnform = 'B-form, bivariate tensor product, array format';
      case 22, fnform = 'BBform, bivariate, array format';
      case 24, fnform = 'polynomial in shifted power form, bivariate';
      case 25, fnform = 'thin-plate spline, bivariate';
   % matrices:
      case 40, fnform = 'almost block diagonal form';
      case 41, fnform = 'spline version of almost block diagonal form';
   % 42 = 'factorization of spline version of almost block diagonal form'
   %      (not yet implemented)
   
   % multivariate:
      case 94, fnform = ...
                  'polynomial in shifted normalized power form, multivariate';
      otherwise
         error('SPLINES:FNBRK:unknownform','Input is of unknown (function) form.')
   end
   
   if nargin>1 %  return some parts if possible
      switch fn(1)
      case 10, [varargout{:}] = ppbrk(fn,varargin{:});
      case {11,12}, [varargout{:}] = spbrk(fn,varargin{:});
      otherwise
         error('SPLINES:FNBRK:unknownpart',...
              ['Parts for ',fnform,' are not (yet) available.'])
      end
   else        % print available information
      if nargout
         error('SPLINES:FNBRK:partneeded','You need to specify a part to be returned.')
      else
         fprintf(['The input describes a ',fnform,'\n\n'])
         switch fn(1)
         case 10, ppbrk(fn);
         case {11,12}, spbrk(fn);
         otherwise
            fprintf('Its parts are not (yet) available.\n')
         end
      end
   end
   return
end   
 
% we reach this point only if FN is a structure.
pre = fn.form(1:2);
switch pre
case {'pp','rp','st'}
case {'B-','BB'},     pre = 'sp';
case 'rB',            pre = 'rs';
otherwise
   error('SPLINES:FNBRK:unknownform','Input is of unknown (function) form.')
end

if nargin>1
   if(strcmp(pre,'pp'))             % This is only for working with the compiler
      varargout{:} = ppbrk(fn,varargin{:});
   else
      eval(['[varargout{:}] = ',pre,'brk(fn,varargin{:});'])
   end
else
   if nargout
      error('SPLINES:FNBRK:partneeded','You need to specify a part to be returned.')
   else
      fprintf(['The input describes a ',fn.form(1:2),'form\n\n'])
      eval([pre,'brk(fn)'])
   end
end

%----------------------------------------------------------------------------------
function fn = fnchg(fn,part,value)
%FNCHG Change part(s) of a form.
%   Copyright 1987-2003 C. de Boor and The MathWorks, Inc.
%   $Revision: 1.3 $  $Date: 2003/02/13 22:22:42 $

if ~ischar(part)
     error('SPLINES:FNCHG:partnotstr','PART must be a string.'), end

switch part(1)
case 'd'
   if part(end)~='z'
      oldvalue = fnbrk(fn,'dim');
      if prod(value)~=prod(oldvalue)
         error('SPLINES:FNCHG:wrongdim', ...
	      ['The specified target dimension, [',num2str(value),...
               '], does not match the present target dimension, [',...
               num2str(oldvalue),'].'])
      end
   end
   fn.dim = value;
case 'i', fn = fnbrk(fn,value);
otherwise
   error('SPLINES:FNCHG:wrongpart','The part ''', part, ''' cannot be reset via FNCHG.')
end

%----------------------------------------------------------------------------------
function fprime = fnder(f,dorder)
%FNDER Differentiate a function.
%   Copyright 1987-2003 C. de Boor and The MathWorks, Inc.
%   $Revision: 1.19 $

if ~isstruct(f), f = fn2fm(f); end

try %treat the function as vector-valued if it is not
   sizeval = fnbrk(f,'dim');
   if length(sizeval)>1, f = fnchg(f,'dz',prod(sizeval)); end
catch
 error('SPLINES:FNDER:unknownfn',['Cannot handle the given fooorm, ',f.form,'.'])
end

if nargin<2, dorder=1; end

switch f.form(1:2)
case 'pp' % the function is in ppform:
	[breaks,coefs,l,k,d]=ppbrk(f);
	if iscell(breaks) % the function is multivariate
		m = length(k);
		if length(dorder)~=m
			error('SPLINES:FNDER:ordermustbevec', ...
				['DORDER should be a ' num2str(m) '-vector.'])
		end
		sizec = [d,l.*k]; %size(coefs);
      for i=m:-1:1
         dd = prod(sizec(1:m));
         dpp = fnderp(ppmak(breaks{i},reshape(coefs,dd*l(i),k(i)),dd), ...
                      dorder(i));
         breaks{i} = dpp.breaks; sizec(m+1) = dpp.pieces*dpp.order;
         coefs = reshape(dpp.coefs,sizec);
         if m>1
             coefs = permute(coefs,[1,m+1,2:m]);
             sizec(2:m+1) = sizec([m+1,2:m]);
         end
      end
      fprime = ppmak(breaks,coefs,sizec);
	else
		fprime = fnderp(f,dorder);
	end

case {'B-','BB'} % the function is in B-form or BB-form;
                 % omit trivial B-spline terms.
   [knots,coefs,n,k,d]=spbrk(f);
   if iscell(knots)       % the function is multivariate
      m = length(knots);
		if length(dorder)~=m
			error('SPLINES:FNDER:ordermustbevec', ...
				['DORDER should be a ' num2str(m) '-vector.'])
		end
      sizec = [d,n];% size(coefs);
      for i=m:-1:1
         dsp = fnderb(spmak(knots{i},...
            reshape(coefs,prod(sizec(1:m)),sizec(m+1))),dorder(i));
         knots{i} = dsp.knots; sizec(m+1) = dsp.number; 
         coefs = reshape(dsp.coefs,sizec); 
         if m>1
            coefs = permute(coefs,[1,m+1,2:m]);
            sizec(2:m+1) = sizec([m+1,2:m]);
         end
      end
      fprime = spmak(knots,coefs,sizec);
   else
      fprime = fnderb(f,dorder);
   end
case {'rp','rB'}
  error('SPLINES:FNDER:notforrat',...
       'FNDER does not work for rational splines. Use FNTLR instead.')
case 'st'
   if strcmp(f.form,'st-tp00')
      if length(dorder)~=2||sum(dorder)>1
         error('SPLINES:FNDER:onlyfirstpartial',...
              ['At present, only the first order partial derivatives of', ...
                ' a bivariate thin-plate spline can be generated.'])
      end        
      if sum(dorder)==0, fprime = f; return, end
      [centers,coefs] = stbrk(f);
      if dorder(1)==1  % we are to differentiate wrto the first argument
         type = 'tp10'; coefs(:,[end-1,end]) = [];     
      else             % we are to differentiate wrto the second argument
         type = 'tp01'; coefs(:,[end-2,end]) = [];     
      end 
      fprime = stmak(centers,coefs,type,fnbrk(f,'interv'));
   else
       error('SPLINES:FNDER:notforst',...
            ['FNDER does not (yet) work for st functions of type ',...
                 f.form(4:end),'.'])
   end
  
otherwise
   error('SPLINES:FNDER:unknownfn','F does not appear to describe a function.')
end

if length(sizeval)>1, fprime = fnchg(fprime,'dz',sizeval); end

function fprime = fnderp(f,dorder)
%FNDERP Differentiate a univariate function in ppform.
[breaks,coefs,l,k,d]=ppbrk(f);
if k<=dorder
   fprime=ppmak([breaks(1) breaks(l+1)],zeros(d,1));
elseif dorder<0    % we are to integrate
   fprime = f;
   for j=1:(-dorder)
      fprime = fnint(fprime);
   end
else
   knew=k-dorder;
   for j=k-1:-1:knew
      coefs=coefs.*repmat([j:-1:j-k+1],d*l,1);
   end
   fprime=ppmak(breaks,coefs(:,1:knew),d);
end

function fprime = fnderb(f,dorder)
%FNDERB Differentiate a univariate function in B-form.

[t,a,n,k,d]=spbrk(f);
if k<=dorder
   fprime=spmak(t,zeros(d,n));
elseif dorder<0    % we are to integrate
   fprime = f;
   for j=1:(-dorder)
      fprime = fnint(fprime);
   end
else
   knew=k-dorder;
   for j=k-1:-1:knew
      tt=t(j+1+[0:n])-t(1:n+1); z=find(tt>0); nn=length(z);
      temp=(diff([zeros(1,d);a.'; zeros(1,d)])).';
      a=temp(:,z)./repmat(tt(z)/j,d,1);
      t=[t(z) t(n+2:n+j+1)]; n=nn;
   end
   fprime=spmak(t,a);
end

%----------------------------------------------------------------------------------
function intgrf = fnint(f,ifa)
%FNINT Integrate a function.
%   Copyright 1987-2003 C. de Boor and The MathWorks, Inc.
%   $Revision: 1.17 $

if ~isstruct(f), f = fn2fm(f); end

try %treat the function as vector-valued if it is not
   sizeval = fnbrk(f,'dim');
   if length(sizeval)>1, f = fnchg(f,'dz',prod(sizeval)); end
catch
   error('SPLINES:FNINT:unknownfn',['Cannot handle the given form, ',f.form,'.'])
end

switch f.form(1:2)
case 'pp'      % the function is in ppform:
   
   if length(f.order)>1
         error('SPLINES:FNINT:useFNDER','Use FNDER with negative derivative order.'), end
      
   [breaks,coefs,l,k,d]=ppbrk(f);
   coefs=coefs./repmat([k:-1:1],d*l,1);
   if nargin==1, ifa = zeros(d,1); end
   if l<2
      intgrf=ppmak(breaks,[coefs ifa],d);
   else
      % evaluate each integrated polynomial at the right endpoint of its
      % interval (this is adapted from PPUAL)
      xs=diff(breaks(1:l));index=[1:l-1];
      if d>1
         xs = reshape(repmat(xs,d,1),1,(l-1)*d);
         index = reshape(1+repmat(d*index,d,1)+repmat([-d:-1].',1,l-1), ...
                        (l-1)*d,1);
      end
      vv=xs.*coefs(index,1).';
      for i=2:k
         vv = xs.*(vv + coefs(index,i).');
      end
      if (d>1)
         junk=zeros(d,l-1);junk(:)=vv;last=(cumsum([ifa junk].')).';
      else
         last=cumsum([ifa,vv]);
      end

      intgrf=ppmak(breaks,[coefs(:,1:k) last(:)],d);
   end
case {'B-','BB'}   % the function is in B-form or BBform.
   if length(f.order)>1
      error('SPLINES:FNINT:useFNDER','Use FNDER with negative derivative order.'), end
   
   % Set it up so that it would be correct on the interval [t(1) .. t(n+k)].
   % There is no way to make it correct everywhere since the integral of
   % a spline over the interval [t(1) .. t(n+k)] need not be zero.
   [t,a,n,k,d]=spbrk(f);
   index = find(diff(t)>0);      % increase multiplicity of last knot to  k
   needed = index(length(index)) - n; % =  k+1 - (n+k - index(length(index));
   if (needed > 0)
      t = [t repmat(t(n+k),1,needed)]; a = [a zeros(d,needed)]; n = n+needed;
   end
   if nargin>1 % if a left-end value is specified, increase left-end knot
               % multiplicity to k+1, making the additional coefficients
               %  0 , then add IFA to all coefficients of the integral.
      needed = k - index(1);
      intgrf = spmak([repmat(t(1),1,needed+1) t t(n+k)], ...
         cumsum([ifa,zeros(d,needed),a.*repmat((t(k+[1:n])-t(1:n))/k,d,1)],2));
   else
      intgrf = spmak([t t(n+k)], ...
                     cumsum(a.*repmat((t(k+[1:n])-t(1:n))/k,d,1),2));
   end
case {'rB','rp'}
   error('SPLINES:FNINT:notforrat','FNINT does not work for rational splines.')
case 'st'
   error('SPLINES:FNINT:notforst','FNINT does not work for st functions.')
otherwise  % unknown representation
   error('SPLINES:FNINT:unknownfn','F does not appear to describe a function.')
end

if length(sizeval)>1, intgrf = fnchg(intgrf,'dz',sizeval); end

%----------------------------------------------------------------------------------
function fnew = fnrfn(f,varargin)
%FNRFN Insert additional sites into the partition for F.
%   Copyright 1987-2003 C. de Boor and The MathWorks, Inc.
%   $Revision: 1.20 $

if ~isstruct(f), f = fn2fm(f); end

switch f.form(1:2)
case {'B-','BB'},  fnew = sprfn(f,varargin{:});
case 'pp',         fnew = pprfn(f,varargin{:});
case 'rB',         fnew = fn2fm(sprfn(fn2fm(f,'B-'),varargin{:}),'rB');
case 'rp',         fnew = fn2fm(pprfn(fn2fm(f,'pp'),varargin{:}),'rp');
case 'st'
   error('SPLINES:FNPLT:notforst','FNRFN does not work for st functions.')
otherwise
   error('SPLINES:FNPLT:unknownfn','F does not appear to describe a function.')
end

%----------------------------------------------------------------------------------
function v = fnval(f,varargin)
%FNVAL Evaluate a function.
%   Copyright 1987-2003 C. de Boor and The MathWorks, Inc.
%   $Revision: 1.20 $

%   Hacked version to work with the compiler (in the 'pp' case)
%   MUST reside in root directory

if ~isstruct(f)
   if isstruct(varargin{1})
      temp = f; f = varargin{1}; varargin{1} = temp;
   else
      f = fn2fm(f);
   end
end

try
   [m, sizeval] = fnbrk(f,'var','dim');     % Fails in compiled version
catch
   if(strcmp(f.form(1:2),'pp'))             % This is only for working with the compiler
      [m, sizeval] = ppbrk(f,'var','dim');
   else
      error('SPLINES:FNVAL:unknownfn',...
      ['Cannot handle the given form, ',f.form,'.'])
   end
end
    % record, then adjust, size of site array and of function values.
if ~iscell(varargin{1})
   sizex = size(varargin{1});
   if m>1
      if sizex(1)~=m
         error('SPLINES:FNVAL:wrongsizex', ...
	      ['Each X(:,j) must be a ',num2str(m),'-vector.'])
      end
   sizex(1) = [];
   end
   if length(sizex)>2
      varargin{1} = reshape(varargin{1},sizex(1),prod(sizex(2:end)));
   elseif length(sizex)==2&&sizex(1)==1, sizex = sizex(2); end
else
   sizex = cellfun('length',varargin{1});
end
if length(sizeval)>1, f = fnchg(f,'dz',prod(sizeval));
else if sizeval==1&&length(sizex)>1; sizeval = []; end
end

switch f.form(1)
case 'B',  ff = @spval;
case 'p',  ff = @ppual;
case 'r',  ff = @rsval;
case 's',  ff = @stval;
otherwise
   error('SPLINES:FNVAL:unknownfn','Unknown function type encountered.')
end
v = reshape(feval(ff,f,varargin{:}),[sizeval,sizex]);

%----------------------------------------------------------------------------------
function [xi,m] = knt2brk(t)
%KNT2BRK From knots to breaks and their multiplicities.
%   Copyright 1987-2003 C. de Boor and The MathWorks, Inc. 
%   $Revision: 1.13 $

difft = diff(t); if any(difft<0) t = sort(t); difft = diff(t); end
[r,c] = size(t);
if r>1 % make sure to return vectors of the same kind
   index = [1;find(difft>0)+1];
   m = diff([index;r*c+1]);
else
   index = [1 find(difft>0)+1];
   m = diff([index r*c+1]);
end
xi = t(index);

%----------------------------------------------------------------------------------
function [m,t] = knt2mlt(t)
%KNT2MLT Knot multiplicities.
%   Copyright 1987-2003 C. de Boor and The MathWorks, Inc. 
%   $Revision: 1.14 $

[r,c] = size(t);
if r*c<2 m=0; return, end

difft = diff(t); if any(difft<0) t = sort(t); difft = diff(t); end
index = zeros([r,c]); index(2:end) = difft==0;
m = cumsum(index);
zz = find(diff(index)<0);
if isempty(zz) return, end

z = zeros([r,c]);       pt = m(zz);
pt(2:end) = diff(pt);   z(zz+1) = pt;
m = m - cumsum(z);

%----------------------------------------------------------------------------------
function sp = pp2sp(pp,sconds)
%PP2SP Convert from ppform to B-form.
%   Copyright 1987-2003 C. de Boor and The MathWorks, Inc.
%   $Revision: 1.17 $

if nargin<2, sconds = []; end

if ~isstruct(pp), pp = fn2fm(pp); end

sizeval = fnbrk(pp,'dim');
if length(sizeval)>1, pp = fnchg(pp,'dz',prod(sizeval)); end

if iscell(pp.breaks) % we are dealing with a multivariate spline

   [breaks,coefs,l,k,d] = ppbrk(pp);
   m=length(k);
   if isempty(sconds), sconds = cell(1,m);
   elseif ~iscell(sconds), sconds = num2cell(repmat(sconds,1,m));
   end
   sizec = [d,l.*k]; %size(coefs);
   for i=m:-1:1
      dd = prod(sizec(1:m));
      spi = pp2sp(ppmak(breaks{i},reshape(coefs,dd*l(i),k(i)),dd),...
                           sconds{i});
      knots{i} = spi.knots;  sizec(m+1) = spi.number;
      coefs = reshape(spi.coefs,sizec);
      if m>1
         coefs = permute(coefs,[1,m+1,2:m]); sizec = sizec([1,m+1,2:m]);
      end
   end
   sp = spmak(knots,coefs,sizec);
else
   sp = pp2sp1(pp,sconds);
end

if length(sizeval)>1, sp = fnchg(sp,'dz',sizeval); end

function sp = pp2sp1(pp,sconds)
%PP2SP1 Convert univariate spline from ppform to B-form.

mustguess = 0;
  % if SCONDS is not specified, or else if its first entry is in (0..1),
  % the proper continuity across each interior knot is to be guessed, using,
  % in the second case, SCONDS(1) as the tolerance for it.
if isempty(sconds), mustguess = 1; tol = 1.e-12;
elseif (0<sconds(1)&&sconds(1)<1), mustguess = 1; tol = sconds(1);
end

[breaks,coefs,l,k,d] = ppbrk(pp);

if mustguess
                                 % guess at smoothness across breaks
   if l==1, sconds = [];
   else % evaluate each piece (but the last) at its right endpoint
      x = breaks(2:l)-breaks(1:l-1);
      if d>1 % repeat each point D times if necessary
         x = repmat(x,d,1);
      end
      x = x(:); a = coefs(1:(d*(l-1)),:);
      for ii=k:-1:2
         for i=2:ii
            a(:,i) = x.*a(:,i-1)+a(:,i);
         end
      end
      % now, at each interior break, look for the smallest i with
      %  |a(:,k-i)-coefs(:+d,k-i)|  >
      %                   >  tol*max(a(:,k-i),coefs(:,k-i),coefs(:+d,k-i))
      % taking i = k  if there is none.
      % if d>1, one would have to take the smallest such i over the d
      % functions involved.

      % first get the sizes
      temp = 1:d*(l-1); temp = [temp;temp+d;temp+l*d];
      tmp = abs([coefs;a]);
      maxes = reshape(max(reshape(tmp(temp,:),3,d*(l-1)*k)),d*(l-1),k);

      % then do the comparison
      tmp = repmat(1:k,d*(l-1),1);
      index = find(abs(coefs(d+1:d*l,:)-a)<=tol*maxes);
      tmp(index) = zeros(length(index),1);
      sconds = k - max(tmp.');
      if d>1, sconds = min(reshape(sconds,d,l-1)); end
   end
end

if (length(sconds)~=l-1)
error('SPLINES:PP2SP:condsdontmatchbreaks', ...
      'Number of smoothness conditions must match number of interior breaks.')
end

mults = k - sconds(:).';
knots = brk2knt(breaks,[k mults k]);
rights = cumsum([1 k mults]);   %     RIGHTS(j) is the first place in the
                                %     knot sequence at which BREAKS(j)
                                %     appears (if it appears at all).

n = length(knots)-k;

% At each break  tau = BREAKS(i) ,i=1,...,l , use the de Boor-Fix formula
%
%       a_j = sum_{r=1:k} (-)^{r-1} D^{r-1} psi_j(tau) D^{k-r-1}f(tau)
%
%                                   for  t_j+ \le tau \le t_{j+k-1}-
%
% with    psi_j(t) := (knots(j+1)-t) ... (knots(j+k-1)-t)/(k-1)!
% to compute the coefficients of the  k  B-splines having the interval
%  [BREAKS(i) .. BREAKS(i+1)]  in their support. Different break intervals may,
% in this way, provide a value for the B-spline coefficient; in that case,
% choose a `most accurate' one.
% Generate the needed derivatives of  psi_j  by differentiating,  k-1  times
% with respect to  t , the program
%         v = 1
%         for i=1:(k-1)
%            v = v*(knots(j+i)-t)
%         end
% for the calculation of  (k-1)! psi(t) .

nx = k*l;
%      Each break is used k times, hence
xindex = reshape(repmat(1:l,k,1),1,nx);
%  TINDEX((j-1)*k + r)  is the index of the B-spline having the interval
%  [BREAKS(j) .. BREAKS(j+1))  as the  (k+1-r)th  knot interval in its support,
%  r=1,...,k ; i.e.,
tindex = reshape(repmat(rights(2:l+1),k,1)-repmat((k:-1:1).',1,l),1,nx);

values = zeros(k,nx);
values(k,:) = ones(1,nx);
for j=1:k-1
   xx = knots(tindex+j) -breaks(xindex);
   for i=(k-j):k-1
      values(i,:) = values(i,:).*xx + ((k-i)/i)*values(i+1,:);
   end
   values(k,:) = values(k,:).*xx;
end
% In the above, a straight-forward inner loop, with the second term being
%  -(k-i)*VALUES() , would generate in  VALUES(r,:)  the  (k-r)th derivative of
%  (k-1)!psi , while  COEFS(:,s)  contains  D^{k-s}f/(k-s)! . We want to sum
%  (-)^{k-1-r} D^{k-1-r}psi D^r f =
%                     = (-)^{k-1-r} (r!/(k-1)!) values(r+1,:) coefs(:,k-r)
% over  r = 0, ..., k-1 . In particular,
%       for  r = k-1 , the weight is 1,
%       for  r = k-2 , the weight is -1/(k-1),
%       for  r = k-3 , the weight is 1/((k-1)(k-2)) = (-1/(k-1))(-1/(k-2)),
%        etc.
%   In other words, we can incorporate these weights into the calculations
% of VALUES by changing the weight in the inner loop, from  -(k-i)  to
% +(k-i)/i , as was done above.
%
%   Here is the summing:

ac = repmat((1:d).',1,nx)+repmat(d*(xindex-1),d,1);
av = repmat(1:nx,d,1);
coefs = reshape(sum(coefs(ac(:),k:-1:1).'.*values(:,av(:)),1),d,nx);

%  Now, choose, for each B-spline coefficient, the best of possibly several
% ways of computing it, namely the one for which the corresponding  psi
% vector is the smallest (in 1-norm).
% Start this by presetting the (k,n)-array NORMS to inf, then setting
%  VALUES(i,j)  equal to the 1-norm of  VALUES(:,r)  if that vector is
% the one that computes coefficient  j  from the i-th knot interval in the
% support of the  j-th B-spline.
norms = repmat(inf,1,k*n);
vindex = ...
reshape(repmat(k*(rights(2:l+1)-k),k,1)+repmat((1-k:k+1:(k*(k-1))).',1,l),1,nx);
norms(vindex) = sum(abs(values),1);
%  ... Then, for each  j , find the number  INDEX(j)  so that
%  NORMS(INDEX(j),j)  is the smallest element in  NORMS(:,j) .
%  (worry about the one-row vs more-than-one-row discontinuity in MATLAB)
if (k==1)
   index = ones(1,n);
else
   [ignore,index] = min(reshape(norms,k,n));
end
%  ... Finally, for each  j=1:n , determine which column of VALUES this
% smallest number came from and choose the corresponding column of COEFS
% as the j-th B-spline coefficient:
norms(vindex) = 1:nx;
sp = spmak(knots, coefs(:,norms((0:n-1)*k+index)));

% Check correctness by comparing input  pp  with  sp2pp(sp) ??

%----------------------------------------------------------------------------------
function varargout = ppbrk(pp,varargin)
%PPBRK Part(s) of a ppform.
%   Copyright 1987-2003 C. de Boor and The MathWorks, Inc.
%   $Revision: 1.21 $

if ~isstruct(pp)
   if pp(1)~=10
      error('SPLINES:PPBRK:unknownfn',...
      'The input array does not seem to describe a function in ppform.')
   else
      ppi = pp;
      di=ppi(2); li=ppi(3); bi=reshape(ppi(3+[1:li+1]),1,li+1); ki=ppi(5+li);
      pp.coefs=reshape(ppi(5+li+[1:di*li*ki]),di*li,ki);
      pp.form = 'pp'; pp.dim = di; pp.breaks = bi; pp.pieces = li; 
      pp.order = ki;
   end
end 

if ~isequal(pp.form,'pp')
   error('SPLINES:PPBRK:notpp',...
   'The input does not seem to describe a function in ppform.')
end
if nargin>1 % we have to hand back one or more parts
   np = max(1,nargout);
   if np>length(varargin)
      error('SPLINES:PPBRK:moreoutthanin', ...
            'Too many output arguments for the given input.')
   end
   varargout = cell(1,np);
   for jp=1:np
      part = varargin{jp};
      if ischar(part)
         if isempty(part)
	    error('SPLINES:PPBRK:partemptystr',...
	    'Part specification should not be an empty string.')
	 end
         switch part(1)
         case 'f',       out1 = [pp.form,'form'];
         case 'd',       out1 = pp.dim;
         case {'l','p'}, out1 = pp.pieces;
         case 'b',       out1 = pp.breaks;
         case {'o','k'}, out1 = pp.order;
         case 'c',       out1 = pp.coefs;
	 case 'v',       out1 = length(pp.order);
         case 'g',       % if the spline is univariate, scalar-valued,
                         % return the coefs in the form needed in the ppform
                         % used in PGS.
            if length(pp.dim)>1||pp.dim>1||iscell(pp.order)
               error('SPLINES:PPBRK:onlyuniscalar', ...
                     '''%s'' is only available for scalar-valued',...
                              ' univariate pp functions.',part)
            else
               k = pp.order;
               out1 = (pp.coefs(:,k:-1:1).').* ...
	                repmat(cumprod([1 1:k-1].'),1,pp.pieces);
            end
         case 'i'
            if iscell(pp.breaks)
               for i=length(pp.order):-1:1
                  out1{i} = pp.breaks{i}([1 end]); end
            else
               out1 = pp.breaks([1 end]);
            end
         otherwise
            error('SPLINES:PPBRK:unknownpart','''%s'' is not part of a ppform.',part)
         end
      elseif isempty(part)
	 out1 = pp;
      else % we are to restrict PP to some interval or piece
	 sizeval = pp.dim; if length(sizeval)>1, pp.dim = prod(sizeval); end
         if iscell(part)  % we are dealing with a tensor-product spline
   
            [breaks,c,l,k,d] = ppbrk(pp); m = length(breaks);
            sizec = [d,l.*k]; %size(c);
            if length(sizec)~=m+1
	       error('SPLINES:PPBRK:inconsistentfn', ...
	       'Information in PP is inconsistent.'),
            end
            for i=m:-1:1
               dd = prod(sizec(1:m));
               ppi = ppbrk1(ppmak(breaks{i},reshape(c,dd*l(i),k(i)),dd),...
                           part{i}) ;
               breaks{i} = ppi.breaks; sizec(m+1) = ppi.pieces*k(i);
               c = reshape(ppi.coefs,sizec);
               if m>1
                  c = permute(c,[1,m+1,2:m]);
                  sizec(2:m+1) = sizec([m+1,2:m]);
               end
            end
            out1 = ppmak(breaks,c, sizec);
   
         else  % we are dealing with a univariate spline
   
            out1 = ppbrk1(pp,part);
         end
         if length(sizeval)>1, out1 = fnchg(out1,'dz',sizeval); end
      end
      varargout{jp} = out1;
   end
else
   if nargout==0
     if iscell(pp.breaks) % we have a multivariate spline and, at present,
                          % I can't think of anything clever to do; so...
       disp(pp)
     else
       disp('breaks(1:l+1)'),        disp(pp.breaks)
       disp('coefficients(d*l,k)'),  disp(pp.coefs)
       disp('pieces number l'),      disp(pp.pieces)
       disp('order k'),              disp(pp.order)
       disp('dimension d of target'),disp(pp.dim)
     end
   else
      varargout = {pp.breaks, pp.coefs, pp.pieces, pp.order, pp.dim};
   end
end

function pppart = ppbrk1(pp,part)
%PPBRK1 restriction of pp to some piece or interval

if isempty(part)||ischar(part), pppart = pp; return, end

if size(part,2) > 1 , % extract the part relevant to the interval 
                      % specified by  part =: [a b]  
   pppart = ppcut(pp,part(1,1:2));
else                  % extract the part(1)-th polynomial piece of pp (if any)
   pppart = pppce(pp,part(1));
end

function ppcut = ppcut(pp,interv)
%PPCUT returns the part of pp  specified by the interval interv =: [a b]  

xl = interv(1); xr = interv(2); if xl>xr, xl = xr; xr = interv(1);  end
if xl==xr
   warning('SPLINES:PPBRK:PPCUT:trivialinterval', ...
           'No changes made since the given end points are equal.')
   ppcut = pp; return
end
 
%  the first pol. piece is  jl ,
% the one responsible for argument  xl
jl=pp.pieces; index=find(pp.breaks(2:jl)>xl); 
                                   % note that the resulting  index  ...
if (~isempty(index)), jl=index(1); % ... is shifted down by one  ...
end                                % ... because of  breaks(2: ...
%  if xl ~= breaks(jl), recenter the pol.coeffs.
x=xl-pp.breaks(jl);
di = pp.dim;
if x ~= 0
   a=pp.coefs(di*jl+[1-di:0],:);
   for ii=pp.order:-1:2
      for i=2:ii
         a(:,i)=x*a(:,i-1)+a(:,i);
      end
   end
   pp.coefs(di*jl+[1-di:0],:)=a;
end
 
%  the last pol. piece is  jr ,
% the one responsible for argument  xr .
jr=pp.pieces;index=find(pp.breaks(2:jr+1)>=xr); 
                                   % note that the resulting ...
if (~isempty(index)), jr=index(1); % index  is shifted down by
end                                % ... one because of  breaks(2: ...
 
%  put together the cut-down  pp
di = pp.dim;
ppcut = ppmak([xl pp.breaks(jl+1:jr) xr], ...
                        pp.coefs(di*(jl-1)+[1:di*(jr-jl+1)],:),di);

function pppce = pppce(pp,j)
%PPPCE returns the j-th polynomial piece of pp  (if any).

%  if  pp  has a  j-th  piece, ...
if (0<j)&&(j<=pp.pieces)  %             ...  extract it
   di = pp.dim;
   pppce = ppmak([pp.breaks(j) pp.breaks(j+1)], ...
              pp.coefs(di*j+[1-di:0],:),di);
else
   error('SPLINES:PPBRK:wrongpieceno', ...
   'The given pp function does not have %g pieces.',j);
end

%----------------------------------------------------------------------------------
function pp = ppmak(breaks,coefs,d)
%PPMAK Put together a spline in ppform.
%   Copyright 1987-2003 C. de Boor and The MathWorks, Inc.
%   $Revision: 1.20 $

if nargin==0
   breaks=input('Give the (l+1)-vector of breaks  >');
   coefs=input('Give the (d by (k*l)) matrix of local pol. coefficients  >');
end

sizec = size(coefs);

if iscell(breaks)  % we are dealing with a tensor-product spline
   if nargin>2
      if prod(sizec)~=prod(d)
        error('SPLINES:PPMAK:coefsdontmatchsize', ...
	      'The coefficient array is not of the explicitly specified size.')
      end, sizec = d;
   end
   m = length(breaks);
   if length(sizec)<m
      error('SPLINES:PPMAK:coefsdontmatchbreaks', ...
           ['If BREAKS is a cell-array of length m, then COEFS must ',...
             'have at least m dimensions.'])
   end
   if length(sizec)==m,  % coefficients of a scalar-valued function
      sizec = [1 sizec];
   end
   sizeval = sizec(1:end-m); sizec = [prod(sizeval), sizec(end-m+(1:m))];
   coefs = reshape(coefs, sizec);

   d = sizec(1);
   for i=m:-1:1
      l(i) = length(breaks{i})-1;
      k(i) = fix(sizec(i+1)/l(i));
      if k(i)<=0||k(i)*l(i)~=sizec(i+1)
         error('SPLINES:PPMAK:piecesdontmatchcoefs', ...
	       ['The specified number %g of polynomial pieces is', ...
                ' incompatible\nwith the total number %g of coefficients', ...
                ' supplied in variable %g.'], l(i),sizec(i+1),i)
      end
      breaks{i} = reshape(breaks{i},1,l(i)+1);
   end
else
  if nargin<3
     if isempty(coefs)
        error('SPLINES:PPMAK:emptycoefs','The coefficient sequence is empty.')
     end
     sizeval = sizec(1:end-1);
     d = prod(sizeval); kl = sizec(end);
     l=length(breaks)-1;k=fix(kl/l);
     if (k<=0)||(k*l~=kl)
       error('SPLINES:PPMAK:piecesdontmatchcoefs', ...
            ['The specified number %g of polynomial pieces is',...
             ' incompatible\nwith the total number %g of coefficients',...
             ' supplied.'],l,kl);
        elseif any(diff(breaks)<0)
        error('SPLINES:PPMAK:decreasingbreaks', ...
	      'The break sequence should be nondecreasing.')
     elseif breaks(1)==breaks(l+1)
        error('SPLINES:PPMAK:extremebreakssame', ...
	      'The extreme breaks should be different.')
     else
        % the ppformat expects coefs in array  (d*l) by k, while the standard
        % input supplies them in an array d by (k*l) . This requires the
        % following shuffling, from  D+d(-1+K + k(-1+L))=D-d +(K-k)d + dkL
        % to  D+d(-1+L + l(-1+K)=D-d +(L-l)d + dlK .
        % This used to be handled by the following:
        % c=coefs(:); temp = ([1-k:0].'*ones(1,l)+k*ones(k,1)*[1:l]).';
        % coefs=[1-d:0].'*ones(1,kl)+d*ones(d,1)*(temp(:).');
        % coefs(:)=c(coefs);
        % Thanks to multidimensional arrays, we can now simply say
        coefs = reshape(permute(reshape(coefs,[d,k,l]),[1,3,2]),d*l,k);
     end
  else % in the univariate case, a scalar D only specifies the dimension of
       % the target and COEFS must be a matrix (though that is not checked for);
       % but if D is a vector, then it is taken to be the intended size of
       % COEFS whatever the actual dimensions of COEFS might be.
     if length(d)==1
        k = sizec(end); l = prod(sizec(1:end-1))/d;
     else
	if prod(d)~=prod(sizec)
	   error('SPLINES:PPMAK:coefsdontmatchsize', ...
	        ['The size of COEFS, [',num2str(sizec), ...
	         '], does not match the specified size, [',num2str(d),'].'])
        end
        k = d(end); l = d(end-1); d(end-1:end) = [];
	if isempty(d), d = 1; end
     end
     if l+1~=length(breaks)
        error('SPLINES:PPMAK:coefsdontmatchbreaks', ...
	      'COEFS indicates %g piece(s) while BREAKS indicates %g.', ...
	l, length(breaks)-1), end
     sizeval = d;
  end
  breaks = reshape(breaks,1,l+1);
end
pp.form = 'pp';
pp.breaks = breaks;
pp.coefs = coefs;
pp.pieces = l;
pp.order = k;
pp.dim = sizeval;

%----------------------------------------------------------------------------------
function ppout = pprfn(pp,varargin)
%PPRFN Insert additional breaks into a ppform.
%   Copyright 1987-2003 C. de Boor and The MathWorks, Inc.
%   $Revision: 1.19 $

[var,sizeval] = fnbrk(pp,'var','dim');
if length(sizeval)>1, pp = fnchg(pp,'dz',prod(sizeval)); end

if var>1   % we are dealing with a multivariate spline

   if nargin>1&&~iscell(varargin{1})
      error('SPLINES:PPRFN:addbreaksnotcell', ...
            'If PP is multivariate, then ADDBREAKS must be a cell-array.')
   end

   [b,c,l,k,d] = ppbrk(pp);
   m = length(l);
   coefs = c; sizec = [d,l.*k]; %size(coefs);
   for i=m:-1:1   % carry out coordinatewise breaks refinement
      dd = prod(sizec(1:m));
      if nargin>1
         ppi = ...
	   pprfn1(ppmak(b{i},reshape(coefs,dd*l(i),k(i)),dd),varargin{1}{i});
      else
         ppi = pprfn1(ppmak(b{i},reshape(coefs,dd*l(i),k(i)),dd));
      end
      b{i} = ppi.breaks; sizec(m+1) = ppi.pieces*ppi.order; 
      coefs = reshape(ppi.coefs,sizec);  
      coefs = permute(coefs,[1,m+1,2:m]); sizec(2:m+1) = sizec([m+1,2:m]);
   end

   % At this point, COEFS contains the tensor-product pp coefficients;
   % also, the various break sequences in B will have been updated. 
   % It remains to return information:
   ppout = ppmak(b, coefs, sizec);

else             % univariate spline refinement
   ppout = pprfn1(pp,varargin{:});
end
if length(sizeval)>1, ppout = fnchg(ppout,'dz',sizeval); end

function ppout = pprfn1(pp,breaks)
%PPRFN1 Insert additional breaks into a univariate ppform.

if nargin<2||(ischar(breaks)&&breaks(1)=='m') % we must supply the midpoints
                                          % of all nontrivial knot intervals
   oldbreaks = ppbrk(pp,'breaks');
   breaks = (oldbreaks(1:end-1)+oldbreaks(2:end))/2;
end

if isempty(breaks), ppout = pp; return, end
breaks = sort(breaks(:).'); lb = length(breaks);

[b,c,l,k,d] = ppbrk(pp);

index0 = find(breaks<b(1)); l0 = length(index0);
% any of these become left-end points for new pieces to the left of the first
% piece, with their coefs all computed from the first piece, i.e., jl(j) = 1.

index2 = find(breaks>b(l+1)); l2 = length(index2);

% now look at the entries of BREAKS in [B(1) .. B(L+1)]. Any of these which
% are not equal to some B(j) become new left-end points, with the coefs
% computed from the relevant piece.
index1 = (l0+1):(lb-l2);
if isempty(index1)
   index = index1;
   jl = [ones(1,l0),repmat(l,1,l2)];
else
   pointer = sorted(b(1:l+1),breaks(index1));
   % find any BREAKS(j) not in B.
   % For them, the relevant left-end point is B(POINTER(INDEX)).
   index = find(b(pointer)~=breaks(index1));
   jl = [ones(1,l0),pointer(index),repmat(l,1,l2)];
end
ljl = length(jl);
     % If all entries of BREAKS are already in B, then just return the input.
if ljl==0, ppout = pp; return, end

% if there are any BREAKS to the right of B(L+1), then B(L+1) and all but the
% rightmost of these must become left-end points, with coefs computed from
% the last piece, i.e., JL(j) = L  for these, and the rightmost BREAKS
% becomes the new right endpoint of the basic interval.
if l2>0
   tmp = breaks(lb);
   breaks(lb:-1:(lb-l2+1)) = [breaks(lb-1:-1:(lb-l2+1)),b(l+1)];
   b(l+1) = tmp;
end

% These are all the additional left-end points:
addbreaks = breaks([index0,index1(index),index2]);
% Now compute the new coefficients in lockstep:
x = addbreaks - b(jl);
if d>1 % repeat each point D times if necessary
   x = repmat(x,d,1);
   omd = (1-d:0).'; jl = repmat(d*jl,d,1)+repmat(omd,1,ljl);
end
a = c(jl,:); x = x(:);
for ii=k:-1:2
   for i=2:ii
      a(:,i) = x.*a(:,i-1)+a(:,i);
   end
end

% Now, all that's left is to insert the coefficients appropriately.
% First, get the enlarged breaks sequence:
newbreaks = sort([b, addbreaks]);
% This should be of length  L + length(JL)  +  1, requiring
newc = zeros(d*(length(newbreaks)-1),k);
if d>1
   newc(repmat(d*sorted(newbreaks,b(1:l)),d,1)+repmat(omd,1,l),:) = c;
   newc(repmat(d*sorted(newbreaks,addbreaks),d,1)+repmat(omd,1,ljl),:) = a;
else
   newc(sorted(newbreaks,b(1:l)),:) = c;
   newc(sorted(newbreaks,addbreaks),:) = a;
end
ppout = ppmak(newbreaks,newc,d);

%----------------------------------------------------------------------------------
function v = ppual(pp,x,left,deriv)
%PPUAL Evaluate function in ppform.
%   Copyright 1987-2002 C. de Boor and The MathWorks, Inc.
%   $Revision: 1.18 $

if ~isstruct(pp), pp = fn2fm(pp); end

if iscell(pp.breaks)   % we are dealing with a multivariate spline
   [breaks,coefs,l,k,d] = ppbrk(pp); m = length(breaks);
   if nargin>2 % set up left appropriately
      if ~iscell(left)
         temp = left; left = cell(1,m); [left{:}] = deal(temp);
      end
   else      left = cell(1,m);
   end

   if iscell(x)  % evaluation on a mesh
      if length(x)~=m, error(['X should specify a(n) ',num2str(m), '-dimensional grid.']), end
      v = coefs; sizev = [d,l.*k]; % size(coefs); 
      nsizev = zeros(1,m);
      for i=m:-1:1
         nsizev(i) = length(x{i}(:)); dd = prod(sizev(1:m));
         v = reshape(ppual1(...
              ppmak(breaks{i},reshape(v,dd*l(i),k(i)),dd), x{i}, left{i}, [] ), [sizev(1:m),nsizev(i)]);
         sizev(m+1) = nsizev(i);
         if m>1
            v = permute(v,[1,m+1,2:m]); sizev(2:m+1) = sizev([m+1,2:m]);
         end
      end
      if (d > 1),        v = reshape(v,[d,nsizev]);
      else               v = reshape(v,nsizev);      end
 
   else          % evaluation at scattered points
      % locate the scattered data in the break sequences:
      [mx,n] = size(x);
      if mx~=m, error(['Each X(:,j) must be a ',num2str(m),'-vector.']), end
      ix = zeros(m,n);
      for i=1:m
         [ox,iindex] = sort(x(i,:));
         ix(i,iindex) = get_index(breaks{i},ox,left{i});
      end
      
      % ... and now set up lockstep polynomial evaluation
      % %%First,  select the relevant portion of the coefficients array. 
      % This has the additional pain that now there are k(i) coefficients
      % for the i-th univariate interval.  
      % The coefficients sit in the (m+1)-dimensional array COEFS, with
      % the (i+1)st dimension containing the coefficients in the i-th
      % dimension, and organized to have first the highest coefficients
      % for each interval, then the next-highest, etc (i.e., as if coming
      % from an array of size [l(i),k(i)]).
      % ix(:,j) is the index vector for the lower corner of j-th point
      % The goal is to extract, for the j-th point, the requisite coefficients
      % from the equivalent one-dimensional array for COEFS, computing a
      % base index from ix(:,j), and adding to this the same set of offsets
      % computed from the l(i) and k(i).
      
      temp = l(1)*[0:k(1)-1]'; 
      for i=2:m
         lt = length(temp(:,1));
         temp = [repmat(temp,k(i),1),reshape(repmat(l(i)*[0:k(i)-1],lt,1),k(i)*lt,1)];
      end
      % also take care of the possibility that the function in PP is vector-valued 
      lt = length(temp(:,1));
      temp = [reshape(repmat([0:d-1].',1,lt),d*lt,1) temp(repmat(1:lt,d,1),:)];

      temp = num2cell(1+temp,1);
      offset = repmat(reshape(sub2ind(size(coefs),temp{:}),d*prod(k),1),1,n);
      
      temp = num2cell([ones(n,1) ix.'],1);
      base = repmat(sub2ind(size(coefs),temp{:}).',d*prod(k),1)-1;
      v = reshape(coefs(base+offset),[d,k,n]);
      
      % ... then do a version of local polynomial evaluation
      for i=m:-1:1
         s = reshape(x(i,:) - breaks{i}(ix(i,:)),[1,1,n]);
         otherk = d*prod(k(1:i-1));
         v = reshape(v,[otherk,k(i),n]);
         for j=2:k(i)
            v(:,1,:) = v(:,1,:).*s(ones(otherk,1),1,:)+v(:,j,:);
         end
         v(:,2:k(i),:) = [];
      end
      v = reshape(v,d,n);
   end
else
   if nargin<3,     left = [];  deriv = []; end
   if nargin<4,     deriv = []; end
   v = ppual1(pp,x,left,deriv);
end

function v = ppual1(pp,x,left,deriv)
%PPUAL1 Evaluate univariate function in ppform.
[mx,nx] = size(x); lx = mx*nx; xs = reshape(x,1,lx);
%  if necessary, sort XS
tosort = 0;
if any(diff(xs)<0),   tosort = 1;   [xs,ix] = sort(xs);   end

%  take apart PP
[breaks,c,l,k,d] = ppbrk(pp);
%  if there are no points to evaluate at, return empty matrix of appropriate size.
if lx==0, v = zeros(d,0); return, end

% for each data site, compute its break interval
index = get_index(breaks,xs,left);

% now go to local coordinates ...
xs = xs-breaks(index);
if d>1 % ... replicate XS and INDEX in case PP is vector-valued ...
   xs = reshape(xs(ones(d,1),:),1,d*lx);
   index = d*index; temp = [-d:-1].';
   index = reshape(1+index(ones(d,1),:)+temp(:,ones(1,lx)), d*lx, 1 );
end
% ... and apply nested multiplication:
if isempty(deriv)                                           % evaluate a X sites
    v = c(index,1).';
    for i=2:k
        v = xs.*v + c(index,i).';
    end
elseif (strcmp(deriv,'first'))                              % Take 1 st derivative
    for j=1:l       % this is j=1:L, not j=1:1
        c(j,1:3)=polyder(c(j,:));
    end
    v = c(index,1).';
    for i=2:k-1
        v = xs.*v + c(index,i).';
    end
elseif (strcmp(deriv,'second'))
    for j=1:l,        c(j,1:3)=polyder(c(j,:));     end     % Take the 1 st derivative
    for j=1:l,        c(j,1:2)=polyder(c(j,1:3));   end     % Now take the 2 nd derivative
    v = c(index,1).';
    for i=2:k-2
        v = xs.*v + c(index,i).';
    end
end

v = reshape(v,d,lx);
if tosort>0,  v(:,ix) = v; end
v = reshape(v,d*mx,nx);

function ind = get_index(mesh,sites,left)
%GET_INDEX appropriate mesh intervals for given ordered data sites
if isempty(left)|left(1)~='l'
   ind = max(sorted(mesh(1:end-1),sites),1);
else
   ind = fliplr(max(length(mesh)-sorted(-fliplr(mesh(2:end)),-sites),1));
end

%----------------------------------------------------------------------------------
function values = rsval(rs,varargin)
%RSVAL Evaluate rational spline.
%   Copyright 1987-2003 C. de Boor and The MathWorks, Inc.
%   $Revision: 1.12 $  $Date: 2003/02/13 23:15:10 $

if ~isstruct(rs)
   error('SPLINES:RSVAL:fnnotstruct','The first argument must be a structure.'), end

% treat rs as a spline in B-form or ppform:
switch rs.form(2)
case 'B', temp = fnval(fn2fm(rs,'B-'),varargin{:});
case 'p', temp = fnval(fn2fm(rs,'pp'),varargin{:});
otherwise
   error('SPLINES:RSVAL:fnnotrat','The function is not a rational spline.')
end

% divide each resulting (d+1)-vector by its last entry (having changed any
% zero entry to 1), then return only the first d entries:
d = fnbrk(rs,'dim'); temp(d+1,find(temp(d+1,:)==0)) = 1;

sv = size(temp);
if length(sv)<3 newsize = [d*sv(1)/(d+1),sv(2)];
else            newsize = [d,sv(2:end)];
end

temp = reshape(temp,d+1,prod(sv)/(d+1));
values = reshape(temp(1:d,:)./temp(repmat(d+1,d,1),:),newsize);

%----------------------------------------------------------------------------------
function pointer = sorted(meshsites, sites)
%SORTED Locate sites with respect to meshsites.
%   Copyright 1987-2003 C. de Boor and The MathWorks, Inc.
%   $Revision: 1.16 $

[ignored,index] = sort([meshsites(:).' sites(:).']);
pointer = find(index>length(meshsites))-(1:length(sites));

%----------------------------------------------------------------------------------
function [v,b] = sprpp(tx,a)
%SPRPP Right Taylor coefficients from local B-coefficients.
%   Copyright 1987-2003 C. de Boor and The MathWorks, Inc. 
%   $Revision: 1.14 $

k = length(a(1,:)); km1 = k-1; b = a;
for r=1:km1
   for i=1:k-r
      b(:,i) =(tx(:,i+km1).*b(:,i)-tx(:,i+r-1).*b(:,i+1))./...
               (tx(:,i+km1)-tx(:,i+r-1));
   end
end

%  Use differentiation at  0  to generate the derivatives
v = b;
for r=2:k
   factor = (k-r+1)/(r-1);
   for i=k:-1:r
      v(:,i) = (v(:,i) - v(:,i-1))*factor./tx(:,i+k-r);
   end
end
v = v(:,k:-1:1);

%----------------------------------------------------------------------------------
function spnew = sprfn(sp,varargin)
%SPRFN Insert additional knots into B-form of a spline.
%   Copyright 1987-2003 C. de Boor and The MathWorks, Inc.
%   $Revision: 1.21 $

sizeval = fnbrk(sp,'dim');
if length(sizeval)>1, sp = fnchg(sp,'dz',prod(sizeval)); end
   
if fnbrk(sp,'var')>1    % we are dealing with a multivariate spline

   if nargin>1&&~iscell(varargin{1})
      error('SPLINES:SPRFN:addknotsnotcell', ...
            'If SP is m-variate, then ADDKNOTS must be an m-cell-array.')
   end

   [t,a,n,k,d] = spbrk(sp);
   m = length(n);
   coefs = a; sizec = [d,n];
   for i=m:-1:1   % carry out coordinatewise knot refinement
      if nargin>1
         spi = sprfn1(spmak(t{i},reshape(coefs,prod(sizec(1:m)),sizec(m+1))),...
                                  varargin{1}{i});
      else
         spi = sprfn1(spmak(t{i},reshape(coefs,prod(sizec(1:m)),sizec(m+1))));
      end 
      t{i} = spi.knots; sizec(m+1) = spi.number; 
      coefs = reshape(spi.coefs, sizec);
      coefs = permute(coefs,[1,m+1,2:m]); sizec(2:m+1) = sizec([m+1,2:m]);
   end
   % At this point, COEFS contains the tensor-product B-spline coefficients;
   % also, the various knot sequences will have been updated. 
   % It remains to return information:
   spnew = spmak(t, coefs,sizec);

else             % univariate spline refinement
   spnew = sprfn1(sp,varargin{:});
end

if length(sizeval)>1, spnew = fnchg(spnew,'dz',sizeval); end

function spnew = sprfn1(sp,addknots)
%SPRFN1 Insert additional knots into B-form of a univariate spline.

if nargin<2||(ischar(addknots)&&addknots(1)=='m') % we must supply the midpoints
                                          % of all nontrivial knot intervals
   breaks = knt2brk(fnbrk(sp,'knots'));
   addknots = (breaks(1:end-1)+breaks(2:end))/2;
end

if isempty(addknots), spnew = sp; return, end
addknots = sort(addknots(:).'); ladd = length(addknots); 

[t,a,n,k,d] = spbrk(sp);

% retain only distinct points, but record their input multiplicity
index = [1 find(diff(addknots)>0)+1];
inmults = diff([index ladd+1]); sortedadds = addknots;
addknots = addknots(index); ladd = length(index);

% compute the current multiplicity of ADDKNOTS in T.
indexr = sorted(t, addknots);
temp = n+k - sorted(-t,-addknots); mults = indexr -  temp(ladd:-1:1);
% ... then reduce INMULTS to make certain that output has knots of multiplicity at most  k .
excess = subplus(inmults+mults-k);
if any(excess>0)
   inmults = inmults - excess;
   index = find(inmults>0);
   if isempty(index)
      warning('SPLINES:SPRFN:alreadyfullmult', ...
              'All additional knots occur already to full multiplicity.')
      spnew = sp; return
   end
   warning('SPLINES:SPRFN:excessmult', ...
           'Insertion of some knot(s) would have caused excess multiplicity.')
   if length(index)<ladd
      addknots = addknots(index); mults = mults(index);
      inmults = inmults(index); ladd = length(addknots);
   end
end

% if the endknot multiplicity is to be increased and/or there are knots
% outside the current basic interval, do it now
lamin=1;
index = find(addknots<t(1));
if ~isempty(index) % there are knots to be put to the left
   totals = sum(inmults(index));
   t = [sortedadds(1:totals) t];
   a = [zeros(d,totals) a]; n = n + totals; indexr = indexr + totals;
   lamin = 1+length(index);
elseif addknots(1)==t(1)
   t = t([ones(1,inmults(1)) 1:(n+k)]);
   a = [zeros(d,inmults(1)) a]; n = n+inmults(1); indexr = indexr + inmults(1);
   lamin = 2;
end
lamax = ladd;
index = find(addknots>t(n+k));
if ~isempty(index) % there are knots to be put to the right
   totals = sum(inmults(index));
   t = [t sortedadds(length(sortedadds)+((1-totals):0))];
   a = [a zeros(d,totals)]; n = n + totals;
   lamax = ladd-length(index);
elseif addknots(ladd)==t(n+k)
   t = t([1:(n+k) repmat(n+k,1,inmults(ladd))]);
   a = [a zeros(d,inmults(ladd))]; n = n+inmults(ladd);
   lamax = ladd-1;
end

% Increase endknot multiplicity to  k , to avoid difficulties.
index = find(diff(t)>0); addl = k-index(1); addr = index(end)-n;
if ( addl>0 || addr>0 )
   t = t([ones(1,addl) 1:(n+k) repmat(n+k,1,addr)]);
   a = [zeros(d,addl) a zeros(d,addr)];
   n = n+addl+addr; indexr = indexr+addl;
end

% Ready for knot insertion, one at a time.
for la=lamin:lamax
   for mm=1:inmults(la)
      newa = a(:,[1:indexr(la)-k+1 indexr(la)-k+1:n]);
      newt = [t((1:indexr(la))) addknots(la) t((indexr(la)+1:n+k))];

      for j=(indexr(la)-k+2):(indexr(la)-mults(la))
         newa(:,j) = ...
             (a(:,j-1)*(t(j+k-1)-addknots(la)) + a(:,j)*(addknots(la)-t(j)))/...
             (          t(j+k-1)                                     -t(j) );
      end
      t = newt; a = newa; n = n+1; indexr = indexr+1; mults(la) = mults(la)+1;
   end
end

if addl>0||addr>0 % remove again those additional end knots and coefficients
   a(:,[1:addl n+((1-addr):0)]) = [];
   t([1:addl n+k+((1-addr):0)]) = []; n = n - addl-addr;
end
spnew = spmak(t,a);

%----------------------------------------------------------------------------------
function spline = spmak(knots,coefs,sizec)
%SPMAK Put together a spline in B-form.
%   Copyright 1987-2003 C. de Boor and The MathWorks, Inc.
%   $Revision: 1.17 $

if nargin==0;
   knots = input('Give the vector of knots  >');
   coefs = input('Give the array of B-spline coefficients  >');
end

if nargin>2
   if numel(coefs)~=prod(sizec)
     error('SPLINES:SPMAK:coefsdontmatchsize', ...
           'The coefficient array is not of the explicitly specified size.')
   end
else
   if isempty(coefs)
      error('SPLINES:SPMAK:emptycoefs','The coefficient array is empty.')
   end
   sizec = size(coefs);
end

m = 1; if iscell(knots), m = length(knots); end
if length(sizec)<m
   error('SPLINES:SPMAK:coefsdontmatchknots', ...
        ['According to KNOTS, the function is %g-dimensional;\n',...
          'hence COEFS must be at least %g-dimensional.'],m,m)
end
if length(sizec)==m,  % coefficients of a scalar-valued function
   sizec = [1 sizec];
end

% convert ND-valued coefficients into vector-valued ones, retaining the
% original size in SIZEVAL, to be stored eventually in SP.DIM .
sizeval = sizec(1:end-m); sizec = [prod(sizeval), sizec(end-m+(1:m))];
coefs = reshape(coefs, sizec);

if iscell(knots), % we are putting together a tensor-product spline
   [knots,coefs,k,sizec] = chckknt(knots,coefs,sizec);
else            % we are putting together a univariate spline
   [knots,coefs,k,sizec] = chckknt({knots},coefs,sizec); knots = knots{1};
end

spline.form = 'B-';
spline.knots = knots;
spline.coefs = coefs;
spline.number = sizec(2:end);
spline.order = k;
spline.dim = sizeval;
% spline = [11 d n coefs(:).' k knots(:).'];

function [knots,coefs,k,sizec] = chckknt(knots,coefs,sizec)
%CHCKKNT check knots, omit trivial B-splines

for j=1:length(sizec)-1
   n = sizec(j+1); k(j) = length(knots{j})-n;
   if k(j)<=0, error('SPLINES:SPMAK:knotsdontmatchcoefs', ...
                     'There should be more knots than coefficients.'), end
   if any(diff(knots{j})<0)
      error('SPLINES:SPMAK:knotdecreasing',...
      'The knot sequence should be nondecreasing.')
   end
   if knots{j}(1)==knots{j}(end)
      error('SPLINES:SPMAK:extremeknotssame',...
      'The extreme knots should be different.')
   end

   % throw out trivial B-splines:
   index = find(knots{j}(k(j)+[1:n])-knots{j}(1:n)>0);
   if length(index)<n
      oldn = n; n = length(index);
      knots{j} = reshape(knots{j}([index oldn+[1:k(j)]]),1,n+k(j));
      coefs = ...
          reshape(coefs, [prod(sizec(1:j)),sizec(j+1),prod(sizec(j+2:end))]);
      sizec(j+1) = n; coefs = reshape(coefs(:,index,:),sizec);
   end
end

%----------------------------------------------------------------------------------
function varargout = spbrk(sp,varargin)
%SPBRK Part(s) of a B-form or a BBform.
%   Copyright 1987-2003 C. de Boor and The MathWorks, Inc.
%   $Revision: 1.24 $

if ~isstruct(sp)
  if sp(1)~=11&&sp(1)~=12
     error('SPLINES:SPBRK:notBform', ...
     'The input array does not seem to describe a function in B-form.')
  else
     di=sp(2);ni=sp(3);
     ci=reshape(sp(3+[1:di*ni]),di,ni);
     kk=sp(4+di*ni);ki=sp(4+di*ni+[1:kk+ni]);
     sp = spmak(ki,ci);
  end
end

if length(sp.form)~=2||sp.form(1)~='B'
   error('SPLINES:SPBRK:notBform','The input does not seem to describe a spline in B-form.')
end
if nargin>1 % we have to hand back one or more parts
   lp = max(1,nargout); % SPBRK(SP,PART) may be part of an expression
   if lp>length(varargin)
      error('SPLINES:SPBRK:moreoutthanin','Too many output arguments for the given input.')
   end
   varargout = cell(1,lp);
   for jp=1:lp
      part = varargin{jp};
      if ischar(part)
         if isempty(part)
	    error('SPLINES:SPBRK:partemptystr',...
	    'Part specification should not be an empty string.')
	 end
         switch part(1)
         case 'f',       out1 = [sp.form,'form'];
         case 'd',       out1 = sp.dim;
         case 'n',       out1 = sp.number;
         case {'k','t'}, out1 = sp.knots;
         case 'o',       out1 = sp.order;
         case 'c',       out1 = sp.coefs;
	 case 'v',       out1 = length(sp.order);
         case 'i', % this must be treated differently in multivariate case
            if iscell(sp.knots)
               for i=length(sp.knots):-1:1  % loop backward to avoid redef.
                  out1{i} = sp.knots{i}([1 end]);
               end
            else
               out1 = sp.knots([1 end]);
            end
	 case 'b', % this must be treated differently in multivariate case
	    if iscell(sp.knots)
               for i=length(sp.knots):-1:1  % loop backward to avoid redef.
                  out1{i} = knt2brk(sp.knots{i});
               end
            else
               out1 = knt2brk(sp.knots);
            end
         otherwise
            error('SPLINES:SPBRK:wrongpart',['''',part,''' is not part of a B-form.'])
         end
      elseif isempty(part)
	 out1 = sp;
      else
         if iscell(part)  % we must be dealing with a tensor-product spline
            c = sp.coefs; knots = sp.knots; m = length(knots);
            sizec = size(c);
            if length(sizec)~=m+1 % trouble because of trailing singleton dims
               sizec = [sp.dim,sp.number]; c = reshape(c,sizec);
            end
            for i=m:-1:1
               dd = prod(sizec(1:m));
               spi = spcut(spmak(knots{i},reshape(c,dd,sp.number(i))), part{i});
               knots{i} = spi.knots; sizec(m+1) = spi.number;
               c = reshape(spi.coefs,sizec);
               if m>1
                  c = permute(c,[1,m+1,2:m]);
                  sizec(2:m+1) = sizec([m+1,2:m]);
               end
            end
            out1 = spmak(knots,c,sizec);
   
         else             % we must be dealing with a univariate spline
            out1 = spcut(sp,part);
         end
      end
      varargout{jp} = out1;
   end
else
   if nargout==0
     if iscell(sp.knots) % we have a multivariate spline and, at present,
                         % I can't think of anything clever to do; so...
       disp(sp)
     else
       disp('knots(1:n+k)'),disp(sp.knots),
       disp('coefficients(d,n)'),disp(sp.coefs),
       disp('number n of coefficients'),disp(sp.number),
       disp('order k'),disp(sp.order),
       disp('dimension d of target'),disp(sp.dim),
     end
   else
    varargout = {sp.knots,sp.coefs, sp.number, sp.order, sp.dim};
   end
end
function out1 = spcut(sp,interv)
%SPCUT change the basic interval

if isempty(interv)||ischar(interv), out1 = sp; return, end

sizei = size(interv);
if sizei(2)>1 % we are to change the basic interval
   tl = interv(1,1); tr = interv(1,2);
   if tl==tr
      warning('SPLINES:SPBRK:SPCUT:trivial_interval', ...
	         'No changes made since the given end points are equal.')
      out1 = sp; return
   end
   if tl>tr, tl = tr; tr = interv(1); end

   index = sorted(sp.knots,[tl,tr]); mults = knt2mlt(sp.knots);
   if tl<sp.knots(1),      m1 = 1;
   elseif tl==sp.knots(1), m1 = 0;
   else,                   m1 = sp.order;                    
      if tl==sp.knots(index(1)), m1 = m1-mults(index(1))-1; end
   end
   if tr>sp.knots(end),      m2 = 1;
   elseif tr==sp.knots(end), m2 = 0;
   else,                     m2 = sp.order;                    
      if tr==sp.knots(index(2)), m2 = m2-mults(index(2))-1; end
   end
   sp = fnrfn(sp, [repmat(tl,1,m1),repmat(tr,1,m2)]);
   index = sorted(sp.knots,[tl tr]);
   if sp.knots(end)>tr
      sp = spmak(sp.knots(1:index(2)),sp.coefs(:,1:(index(2)-sp.order)));
   end
   if sp.knots(1)<tl
      sp = spmak(sp.knots(index(1)-sp.order+1:end), ...
                 sp.coefs(:,index(1)-sp.order+1:end));
   end
   out1 = sp;
else 
   error('SPLINES:SPBRK:partnotinterv',...
   'The given part, %g, does not specify an interval.',interv)
end

%----------------------------------------------------------------------------------
function pp = sp2pp(spline)
%SP2PP Convert from B-form to ppform.
%   Copyright 1987-2003 C. de Boor and The MathWorks, Inc.
%   $Revision: 1.17 $

if ~isstruct(spline), spline = fn2fm(spline); end

sizeval = fnbrk(spline,'dim');
if length(sizeval)>1, spline = fnchg(spline,'dz',prod(sizeval)); end

if iscell(spline.knots)   % we are dealing with a multivariate spline
   [t,a,n,k,d] = spbrk(spline);
   m = length(k);
   coefs = a; sizec = [prod(d),n]; % size(coefs);
   for i=m:-1:1
      ppi = sp2pp1(spmak(t{i},reshape(coefs,prod(sizec(1:m)),n(i))));
      breaks{i} = ppi.breaks;  sizec(m+1) = ppi.pieces*k(i);
      coefs = reshape(ppi.coefs,sizec);
      if m>1
         coefs = permute(coefs,[1,m+1,2:m]); sizec = sizec([1,m+1,2:m]);
      end
   end
   pp = ppmak(breaks,coefs,sizec);
else
   pp = sp2pp1(spline);
end

if length(sizeval)>1, spline = fnchg(spline,'dz',sizeval); end

%----------------------------------------------------------------------------------
function pp = sp2pp1(spline)
%  Take apart the  spline
[t,a,n,k,d] = spbrk(spline);

%  and augment the knot sequence so that first and last knot each have multiplicity  k .
index = find(diff(t)>0); addl = k-index(1); addr = index(end)-n;
if (addl>0||addr>0)
   t = [repmat(t(1),1,addl) t(:).' repmat(t(n+k),1,addr)];
   a = [zeros(d,addl) a zeros(d,addr)];
end

%  From this, generate the pp description.
inter = find( diff(t)>0 ); l = length(inter);
if k>1
   temp = repmat(inter,d,1); dinter = temp(:);
   tx = repmat([2-k:k-1],d*l,1)+repmat(dinter,1,2*(k-1)); tx(:) = t(tx);
   tx = tx-repmat(t(dinter).',1,2*(k-1)); a = a(:);
   temp = repmat(d*inter,d,1)+repmat([1-d:0].',1,l); dinter(:) = temp(:);
   b = repmat(d*[1-k:0],d*l,1)+repmat(dinter,1,k); b(:) = a(b);
   c = sprpp(tx,b);
else temp = a(:,inter); c = temp(:);
end

%   put together the  pp
pp = ppmak([t(inter) t(inter(end)+1)],c,d);

%----------------------------------------------------------------------------------
function spline = sp2bb(spline)
%SP2BB Convert from B-form to BBform.
%   Copyright 1987-2003 C. de Boor and The MathWorks, Inc.
%   $Revision: 1.16 $

if ~isstruct(spline), spline = fn2fm(spline); end

sizeval = fnbrk(spline,'dim');
if length(sizeval)>1, spline = fnchg(spline,'dz',prod(sizeval)); end
if iscell(spline.knots)   % we are dealing with a multivariate spline
   [t,coefs,n,k,d] = spbrk(spline);
   m=length(k);
   sizec = [d,n]; %size(coefs);
   for i=m:-1:1
      spi = sp2bb1(spmak(t{i},reshape(coefs,prod(sizec(1:m)),n(i))));
      knots{i} = spi.knots;  sizec(m+1) = spi.number;
      coefs = reshape(spi.coefs,sizec);
      if m>1
         coefs = permute(coefs,[1,m+1,2:m]); sizec = sizec([1,m+1,2:m]);
      end
   end
   spline = spmak(knots,coefs,sizec); spline.form = 'BB';
      
else
   spline = sp2bb1(spline);
end
if length(sizeval)>1, spline = fnchg(spline,'dz',sizeval); end

%----------------------------------------------------------------------------------
function spline = sp2bb1(spline)
%SP2BB1 Convert univariate spline from B-form to BBform.

[xi,m] = knt2brk(spbrk(spline,'knots'));
spline = sprfn(spline, brk2knt(xi,subplus(spbrk(spline,'order') - m)));
spline.form = 'BB';

%----------------------------------------------------------------------------------
function v = spval(sp,x,left)
%SPVAL Evaluate function in B-form.
%   Copyright 1987-2003 C. de Boor and The MathWorks, Inc.
%   $Revision: 1.20 $

if ~isstruct(sp)
   error('SPLINES:SPVAL:fnnotstruct','SP must be a structure.'), end
if iscell(sp.knots)  % we are dealing with a tensor product spline
   [t,a,n,k,d] = spbrk(sp); m = length(t);
   if nargin>2 % set up left appropriately
      if ~iscell(left)
         temp = left; left = cell(1,m); [left{:}] = deal(temp);
      end
   else
      left = cell(1,m);
   end

   if iscell(x)  % evaluation on a mesh
      v = a; sizev = [d,n]; nsizev = zeros(1,m);
      for i=m:-1:1
         nsizev(i) = length(x{i}(:));
         v = reshape(...
         spval1(spmak(t{i},reshape(v,prod(sizev(1:m)),sizev(m+1))), ...
                 x{i},left{i}),   [sizev(1:m),nsizev(i)]);
         sizev(m+1) = nsizev(i);
         if m>1
            v = permute(v,[1,m+1,2:m]); sizev(2:m+1) = sizev([m+1,2:m]);
         end
      end
      if d>1
         v = reshape(v,[d,nsizev]);
      else
         v = reshape(v,nsizev);
      end
   else          % evaluation at scattered points; this will eventually be done directly here.
      v = ppual(sp2pp(sp),x);
   end
else                 % we are dealing with a univariate spline
   if nargin<3, left = []; end
   v = spval1(sp,x,left);
end

%----------------------------------------------------------------------------------
function v = spval1(sp,x,left)
%SPVAL1 Evaluate univariate function in B-form.

[mx,nx] = size(x); lx = mx*nx; xs = reshape(x,1,lx);
%  If necessary, sort XS:
tosort = 0;
if any(diff(xs)< 0)
   tosort = 1; [xs,ix] = sort(xs);
end

%  Take apart spline:
[t,a,n,k,d] = spbrk(sp);
%  If there are no points to evaluate at, return empty matrix of appropriate
%  size:
if lx==0, v = zeros(d,0); return, end

%  Otherwise, augment the knot sequence so that first and last knot each
%  have multiplicity  >= K . (AUGKNT would not be suitable for this
%  since any change in T must be accompanied by a corresponding change
%  in A.)

index = find(diff(t)>0); addl = k-index(1); addr = index(length(index))-n;
if ( addl>0 || addr>0 )
   npk = n+k; t = t([ones(1,addl) 1:npk npk(ones(1,addr))]);
   a = [zeros(d,addl) a zeros(d,addr)];
   n = n+addl+addr;
end

% For each data point, compute its knot interval:
if isempty(left)||left(1)~='l'
   index = max(sorted(t(1:n),xs),k);
else
   index = fliplr(max(n+k-sorted(-fliplr(t(k+1:n+k)),-xs),k));
end

% Now, all is ready for the evaluation.
if  k>1  % carry out in lockstep the first spline evaluation algorithm
         % (this requires the following initialization):
   dindex = reshape(repmat(index,d,1),d*lx,1);
   tx =reshape(t(repmat([2-k:k-1],d*lx,1)+repmat(dindex,1,2*(k-1))),d*lx,2*(k-1));
   tx = tx - repmat(reshape(repmat(xs,d,1),d*lx,1),1,2*(k-1));
   dindex = reshape(repmat(d*index,d,1)+repmat([1-d:0].',1,lx),d*lx,1);
   b = repmat([d*(1-k):d:0],d*lx,1)+repmat(dindex,1,k);
   a = a(:); b(:) = a(b);

   % (the following loop is taken from SPRPP)

   for r = 1:k-1
      for i = 1:k-r
         b(:,i) = (tx(:,i+k-1).*b(:,i)-tx(:,i+r-1).*b(:,i+1)) ./ ...
                  (tx(:,i+k-1)    -    tx(:,i+r-1));
      end
   end

   v = reshape(b(:,1),d,lx);
else     % the spline is piecewise constant, hence ...
   v = a(:,index);
end

if tosort>0, v(:,ix) = v; end

% Finally, zero out all values for points outside the basic interval:
index = find(x<t(1)|x>t(n+k));
if ~isempty(index)
   v(:,index) = zeros(d,length(index));
end
v = reshape(v,d*mx,nx);

function varargout = stbrk(st, varargin)
%STBRK Part(s) of an stform.
%   Copyright 1987-2003 C. de Boor and The MathWorks, Inc.
%   $Revision: 1.5 $

if ~isstruct(st)
   if st(1)~=25
      error('SPLINES:STBRK:notst',...
      'The input array does not seem to describe a function in stform.')
   end
   dce = st(2); nce = st(3);
   dco = st(4+dce*nce); nco = st(5+dce*nce);
   st = stmak( reshape(st(3+(1:dce*nce)),dce,nce), ...
               reshape(st(5+dce*nce+(1:dco*nco)),dco,nco), 'tp00' );
end

if length(st.form)<2||~isequal(st.form(1:2),'st')
   error('SPLINES:STBRK:unknownfn','The input does not seem to describe a function in stform.')
end

if nargin>1 % we have to hand back one or more parts
   lp = max(1,nargout);
   if lp>length(varargin)
      error('SPLINES:STBRK:moreoutthanin','Too many output arguments for the given input.')
   end
   varargout = cell(1,lp);
   for jp=1:lp
      part = varargin{jp};
      if ischar(part)
         if isempty(part)
	    error('SPLINES:STBRK:partemptystr','Part specification should not be an empty string.')
	 end
         switch part(1)
         case 'f',       out1 = 'stform';
         case 'c'
            if length(part)<2
               error('SPLINES:STBRK:ambiguouspart',['For Part %g, ', ...
	              'did you mean ''centers'' or ''coefficients''?'],jp)
            else
               switch part(2)
               case 'e', out1 = st.centers;
               case 'o', out1 = st.coefs;
               end
            end
         case 'd', out1 = st.dim;
         case 'i', out1 = st.interv;
         case 'n', out1 = st.number;
         case 't', out1 = st.form(4:end);
         case 'v', out1 = size(st.centers,1);
         otherwise
            error('SPLINES:STBRK:wrongpart',['''',part,'''',' is not part of an stform.'])
         end
      else % part is expected to be a basic interval spec.
         if ~iscell(part)
            if size(st.centers,1)==1, part = {part};
            else
               error('SPLINES:STBRK:inarg2notcell',...
	                'The second argument is expected to be a cell array.')
            end
         end
         lpart = length(part);
         if lpart~= size(st.centers,1)
            error('SPLINES:STBRK:toofewinterv',...
	                ['A basic interval spec. must have as many entries', ...
                   ' as the function has variables.'])
         end 
         for j=1:length(part)
	    if isempty(part{j}), part{j} = st.interv{j};
            elseif ~isequal(size(part{j}),[1 2])||part{j}(1)>=part{j}(2)   
               error('SPLINES:STBRK:wronginterv',...
	                ['A basic interval must be a cell array,\n', ...
                      'with each nonempty entry of the form\n', ...
                      '[a b] for some a < b.'],0)
            end
         end
         st.interv = part; out1 = st;
      end
      varargout{jp} = out1;
   end
else
   if nargout==0
      disp(st)
   else
      varargout = {st.centers, st.coefs, st.form(4:end), st.interv};
   end
end

%----------------------------------------------------------------------------------
function colmat = stcol(centers, x, varargin)
%STCOL Scattered translates collocation matrix.
%   Copyright 1987-2003 C. de Boor and The MathWorks, Inc.
%   $Revision: 1.4 $  $Date: 2003/02/13 23:16:05 $

transp = 0; type = 'tp';
for j=1:nargin-2
   narg = varargin{j};
   if ~isempty(narg)&&ischar(narg)
      if isequal(narg(1:2),'tr'), transp = 1;
      else, type = narg;
      end
   end
end

[d,nc] = size(centers);
[dx,nx] = size(x);
if dx~=d
   error('SPLINES:STCOL:wrongsizeX', ...
   [' x  should be of size [',num2str(d),',nx].']), end
if transp
   temp = repmat(1:nx,nc,1);
   points = x(:,temp) - reshape(repmat(centers(:),1,nx),d,nx*nc);
else
   temp = nc; nc = nx; nx = temp;
   temp = repmat(1:nx,nc,1);
   points = reshape(repmat(x(:),1,nx),d,nx*nc) - centers(:,temp);
end
ap2 = sum(points.*points,1);

switch d
case 1
   ap = sqrt(ap2); 
   colmat = reshape(ap2.*ap,nc,nx); 
case 2
   ap2(find(ap2==0)) = 1;
   switch type
   case 'tp'
      colmat = reshape(ap2.*log(ap2),nc,nx); 
   case 'tp00'
      if transp
         colmat = [reshape(ap2.*log(ap2),nc,nx) ; x; ones(1,nx)]; 
      else
         colmat = [reshape(ap2.*log(ap2),nc,nx), x.', ones(nc,1)]; 
      end
   case 'tp10'
      if transp
         colmat = [reshape(2*(points(1,:).*(log(ap2)+1)),nc,nx); ones(1,nx)]; 
      else
         colmat = [reshape(2*(points(1,:).*(log(ap2)+1)),nc,nx), ones(nc,1)]; 
      end
   case 'tp01'
      if transp
         colmat = [reshape(2*(points(2,:).*(log(ap2)+1)),nc,nx); ones(1,nx)]; 
      else
         colmat = [reshape(2*(points(2,:).*(log(ap2)+1)),nc,nx), ones(nc,1)]; 
      end
   otherwise
      error('SPLINES:STCOL:unknowntype',['Cannot handle st functions of type ',type,' (yet).']) 
   end
otherwise
   error('SPLINES:STCOL:atmostbivar',...
   'Cannot handle functions of more than two variables (yet).')
end

%----------------------------------------------------------------------------------
function st = stmak(centers,coefs,type,interv)
%STMAK Put together a function in stform.
%   Copyright 1987-2003 C. de Boor and The MathWorks, Inc.
%   $Revision: 1.4 $  $Date: 2003/02/13 23:16:07 $

if nargin<3||isempty(type)
   type = 'tp';
else
   if ~strcmp(type(1:2),'tp')
      error('SPLINES:STMAK:unknowntype',['Cannot handle the specified type, ''',type,'''.'])
   end
   for (j=length(type):-1:3)      der(j-2) = str2num(type(j));   end
end

[dce, nce] = size(centers); [dco, nco] = size(coefs);

% check that nco and nce are consistent with the specified type.
if ~exist('der','var')
   exces = 0;
else
   switch dce
   case 1
      exces = 2-der;
   case 2
      switch sum(der)
      case 0, exces = 3;
      case 1, exces = 1;
      case 2, exces = 0;
      end   
   otherwise
      error('SPLINES:STMAK:atmostbivar','Cannot handle functions of more than 2 variables.')
   end
end

if nco~=nce+exces
   error('SPLINES:STMAK:centersdontmatchcoefs',['Number of centers and coefficients',...
         ' is inconsistent with the specified type.']) 
end

st.form = ['st-',type];     st.centers = centers;
st.coefs = coefs;           st.ncenters = nce;
st.number = nco;            st.dim = dco;

if nargin<4||isempty(interv)
   % For want of some better idea, define the basic interval as the 
   % bounding box of the centers, i.e., the smallest axiparallel 
   % (hyper-)rectangle that contains all the centers.
   % This can always be altered by STBRK(st,interv).
   if nce==1
      interv = {centers(1,1)+[0 1], centers(2,1)+[0 1]};
   else
      for (j=dce:-1:1)
         interv{j} = [min(centers(j,:)), max(centers(j,:))];
      end
   end
   st.interv = interv;
else
   try
      st = stbrk(st,interv);
   catch
      error('SPLINES:STMAK:wronginterv',['INTERV must be a cell array of size ',...
            '[1,size(CENTERS,1)],\n with each INTERV{i} ','of the form [a,b] with a<=b.'],0)
   end
end

%----------------------------------------------------------------------------------
function values = stval(st,x)
%STVAL Evaluate function in stform.
%   Copyright 1987-2003 C. de Boor and The MathWorks, Inc.
%   $Revision: 1.4 $  $Date: 2003/02/13 23:16:09 $

[centers, coefs] = stbrk(st);

if iscell(x) % we must determine the gridpoints from the given univariate meshes
    [xx,yy] = ndgrid(x{1},x{2});
    nx = [length(x{1}),length(x{2})];
    x = [reshape(xx,1,prod(nx));reshape(yy,1,prod(nx))];
else
    [mx,nx] = size(x);
    if mx~=2
        if nx==2 %switch the two
            x = x.'; nx = mx;
        else
            error('SPLINES:STVAL:wrongsizex','X must be of size [2,...], or else a cell array.')
        end
    end
end

d = size(coefs,1); lx = size(x,2); values = zeros(d,lx);

% avoid use of out-of-core memory by doing calculations in small enough 
% pieces if need be; the upper limit, of 100000, is well below what it
% has to be, but the resulting time penalty for this undershot is
% negligible compared to the cost of this calculation:

lefttodo = lx*size(centers,2);
segments = ceil(lefttodo/100000);
llx = round(linspace(0,lx,segments+1));

for j=1:segments 
   values(:,llx(j)+1:llx(j+1)) = ...
        coefs*stcol(centers, x(:,llx(j)+1:llx(j+1)), 'tr', st.form(4:end)); 
end

if d>1||length(nx)==1, nx = [d nx]; end
values = reshape(values,nx);

%----------------------------------------------------------------------------------
function y = subplus(x)
%SUBPLUS Positive part.
%                                  x , if  x>=0
%   y  = subplus(x) := (x)_{+}  =               ,
y=max(x,zeros(size(x)));
