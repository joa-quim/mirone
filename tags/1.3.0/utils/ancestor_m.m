function p = ancestor_m(p,type,varargin)
%ANCESTOR  Get object ancestor.
%    P = ANCESTOR(H,TYPE) returns the handle of the closest ancestor
%    of h of one of the types in TYPE, or empty if none exists. TYPE
%    may be a single string (single type) or cell array of strings
%    (types). If H is a vector of handles the P is a cell array of the
%    same length as H and P{n} is the ancestor of H(n). If H has one
%    of the specified types then the ancestor of H is H itself.
%    P = ANCESTOR(H,TYPE,'TOPLEVEL') finds the highest level ancestor of
%    one of the types in TYPE
%
%    If H is not an Handle Graphics object, ANCESTOR returns empty.
%
%  Examples:
%    p = ancestor(gca,'figure');
%    p = ancestor(gco,{'hgtransform','hggroup','axes'},'toplevel');

%   Copyright 1984-2004 The MathWorks, Inc.
%   $Revision: 1.1.6.9 $  $Date: 2004/08/16 01:47:07 $

%nargchk(2,4,nargin);

% if ~ishghandle(p)
%     p = [];
%     return;
% end

if (length(p) > 1)
    n = length(p);
    pv = cell(n,1);
    for k=1:n
        pv{k} = ancestor_m(p(k),type,varargin{:});
    end
    p = pv;
    return
end

if (nargin == 2)        % ancestor(h,type)
    while ~isempty(p) && ~isatype(p,type)
        p = get(handle_j(p),'parent');
    end
elseif nargin==3 % ancestor(h,type,'toplevel')
    P=[];
    if isatype(p,type)
        P = p;
    end
    while ~isempty(p)
        p = get(handle_j(p),'parent');
        if isatype(p,type)
            P = p;
        end
    end
    p=P;
end

%-------------------------------------------------------------%
function istype=isatype(h,type)

istype = false;
if ischar(type) % ancestor(h,'type'..)
    if strcmpi(get(h,'type'),type) | isa(handle_j(h),type)
        istype=true;
    end
else % ancestor(h,{'type1','type2',..}...)
    % % make sure it's a cell array
    % type = cellstr(type);
    if any(strcmpi(get(h,'type'),type))
        istype = true;
    else
        % check each cell
        for k=1:length(type)
            if isa(handle_j(h),type{k})
                istype=true;
            end
        end
    end
end

%-------------------------------------------------------------%
function type = handle_j(h)
% Attempt to short-circuit the fact that the built in function HANDLE
% is not compilable and is not documented. Apparently it returns the TYPE
% of the handles H as an object class. Here I just return the handle in input
% because the ANCESTOR function seams to be called within Mirone only with H
% as a figure handle. In this case, things seams to work
% and I hope to not f... the ancestor functionality
type = h;









