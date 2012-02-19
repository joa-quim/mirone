function varargout = findcell(k,c,varargin)
% [cixn,cixo] = findcell(key,ca,'-opt')
% cix         = findcell(key,ca,'-opt')
%
%	returns cell number:offset
%	of any occurrence of <key> in the cell array <ca>
%
% key	search pattern
%	format
%		string
%		numeric
%
% ca	a row or col array of cell(s)
%	format
%		string(s)
%		numeric
%	classes must not be mixed
%
% cix	result vector(s)/struct
%	format
%		struct					[def]
%		[cell_nr offset]			[opt:     -p]
%		[cell_nr] [offset]			[nargout:  2]
%
% opt	parameter	action
%  -a :	-		search ACROSS cells
%  -p :	-		output <cix> as vector(s)	[def: struct]
%
% examples
%	ac={'xa';'xxa';'foo';'xaxa';'goo'};
%	findcell('a',ac)
%		% cell nr:offset   1   2
%		% cell nr:offset   2   3
%		% cell nr:offset   4   2
%		% cell nr:offset   4   4
%	nc={[1:10];pi;pi*[1:10]};
%	findcell(pi,nc)
%		% cell nr:offset   2   1
%		% cell nr:offset   3   1
%	findcell([pi pi],nc);
%		% findcell> key not found
%	findcell([pi pi],nc,'-a');
%		% cell nr:offset   2   1
%	zc={['this' 0 'is0'];[0 0 0 'this' 0 'is0']};
%	findcell(['is' 0],zc);
%		% cell nr:offset   1   3
%		% cell nr:offset   2   6
%
% remark
%	program based on a problem by
%	CSSM member <michael robbins>
%	<Vectorizing FINDSTR> (3/11/03)

% search engine
%	since <findstr> does not find <nan>s, we use <nan>s as
%	stop markers after each cell to prevent <across> results,
%	which is the default behavior
%	cellarray
%		cell1 {nan}
%		cell2 {nan}
%		...
%		cellN {nan}
%
%	since <nan>s cannot be converted to <char>s, we must
%	convert cells containing string arrays to double

% created:
%	us	12-Mar-2003
% modified:
%	us	14-Mar-2003 20:40:37	/ TMW

if	nargout
	varargout = cell(nargout,1);
end
if	nargin < 2 || isempty(k) || isempty(c)
    error('Wrong number of input args')
end
[p,k,c] = chk_input(k,c);
if (p.res), return;     end

% get/set options
mat=[];
p = get_opt(p,varargin{:});

% precondition key/cells
c=c(:);
cs=length(c);
cl=length([c{:}]);

if (p.cflg)			% ... do NOT search across <cell>s
	if (ischar(k) && all(p.isc))
		k=double(k);
		inn=cellfun('length',c)';
		c=[c num2cell(0*ones(cs,1))];
		c=reshape(c.',1,2*cs);
		in=cellfun('length',c)';
		c=double([c{:}]);
		c(cumsum(inn+1))=nan;
	elseif	isnumeric(k) && all(p.isn)
		c=[c num2cell(nan*ones(cs,1))];
		c=reshape(c.',1,2*cs);
		in=cellfun('length',c)';
		c=[c{:}];
	else
		fprintf('findcell> unexpected error');
		return
	end
else				% ...	do search across <cell>s
	in=cellfun('length',c)';
	c=[c{:}];
end

% find indices
ix = strfind(k,c);
if	~isempty(ix)
	[mx,mn]=meshgrid(ix,cumsum(in)-in);
	crow=sum(mx>mn,1)';
	if	p.cflg
		crow=.5*(crow-1)+1;
	end
	cix=find(crow==fix(crow));
	ccol=mx-mn;
	ccol(ccol<=0)=nan;
	ccol=min(ccol,[],1)';
	mat=[crow(cix) ccol(cix)];
end

% prepare output
if (isempty(mat))
	%disp('findcell> key not found');
	return
end
if (nargout == 1)
	if	p.pflg
		varargout{1}.par=p;
		varargout{1}.csiz=cl;
		varargout{1}.cn=mat(:,1);
		varargout{1}.co=mat(:,2);
		varargout{1}.cno=mat;
	else
		varargout{1}=mat;
	end
	elseif	nargout == 2
		varargout{1}=mat(:,1);
		varargout{2}=mat(:,2);
	else	% no output requested
		fprintf('cell nr:offset %8d %8d\n',mat.');
end

%--------------------------------------------------------------------------------
function [p,k,c] = chk_input(k,c)
% must do extensive checking ...
p.res=0;
p.key=k;
p.cs=numel(c);
p.cl=length(c);
if	p.cs ~= p.cl
	txt=sprintf('%d/',size(c));
	fprintf('findcell> input must be a ROW or COL array of cells [%s]',txt(1:end-1))
	p.res=1;
	return
end
if	~iscell(c)
	fprintf('findcell> input must be a CELL array [%s]',class(c));
	p.res=2;
	return
end
p.isc=cellfun('isclass',c,'char');
p.isn=cellfun('isclass',c,'double');
if	sum(p.isc) ~= p.cs && sum(p.isn) ~= p.cs
	fprintf('findcell> input contains invalid or mixed classes');
	p.res=3;
	return
end
if	any(isnan(k))
	fprintf('findcell> numeric key must NOT include <NaN>s');
	p.res=4;
	return
end
if	ischar(k) && all(p.isn)
	p.res=5;
	fprintf('findcell> input mismatch: key(char) / cells(double)');
	return
end
if	isnumeric(k) && all(p.isc)
	k=char(k);
end
%--------------------------------------------------------------------------------
function p = get_opt(p,varargin)
p.cflg=1;	% do NOT search across <cell>s
p.pflg=1;	% output is a <struct>

if	nargin
	if	any(strcmp('-a',varargin))
		p.cflg=0;
	end
	if	any(strcmp('-p',varargin))
		p.pflg=0;
	end
end
