function [yvec,xvec,trimpts] = trimpatch(yvec,ylim,xvec,xlim)

%TRIMPATCH  Trim patch vector data exceeding projection limits
%
%  Purpose
%  TRIMPATCH will identify points in patch vector data
%  which exceed projection limits.  The projection limits
%  are defined by the lower and upper inputs.  Points outside
%  the projection range are replaced with their respective limits,
%  thereby trimming them from the display.  If a patch lies completely
%  outside the trim limits, it is completely replaced with NaNs.
%
%  Synopsis
%           [ymat,xmat,trimpts] = trimpatch(ymat,ylim,xmat,xlim)
%
%  Copyright 1996-2003 The MathWorks, Inc.

%  Argument tests

if nargin ~= 4;  error('Incorrect number of arguments');  end

%  Dimension tests
if ndims(xvec) > 2 || ndims(yvec) > 2
    error('Pages not allowed for xmat or ymat')
elseif ~isequal(size(xvec),size(yvec))
    error('Inconsistent xmat and ymat dimensions')
end

xvec = xvec(:);
yvec = yvec(:);

%  Find the individual patches
indx = find(isnan(xvec) | isnan(yvec));
if isempty(indx);   indx = length(xvec)+1;   end

trimpts = [];

for i = 1:length(indx)
    if (i == 1),    startloc = 1;
    else            startloc = indx(i-1)+1;
    end
    endloc   = indx(i)-1;

    indices = (startloc:endloc)';   %  Indices will be empty if NaNs are
	                                %  neighboring in the vector data.
    if ~isempty(indices)            %  Should not happen, but test just in case

        %  Patches which lie completely outside the trim window.
        %  If at least one point of the patch edge does not lie with the
        %  specified window limits, then the entire patch is trimmed.
	
        if ~any(xvec(indices) >= min(xlim) & xvec(indices) <= max(xlim) & ...
             yvec(indices) >= min(ylim)  & yvec(indices) <= max(ylim))
            trimpts = [trimpts; indices ones(size(indices)) yvec(indices) xvec(indices)];
		    xvec(indices) = NaN;    yvec(indices) = NaN;
        end
	
        %  Need to only test along edge since patch must lie somehow within
        %  the window.  Make sure that the original data is saved before
        %  the points are replaced with the limit data.
	
        %  Points off the bottom
        loctn = find( xvec(indices) < min(xlim) );
        if ~isempty(loctn)
            trimpts = [trimpts; indices(loctn) ones(size(loctn)) yvec(indices(loctn)) xvec(indices(loctn))];
			xvec(indices(loctn)) = min(xlim);
        end
	
        %  Points off the top
        loctn = find( xvec(indices) > max(xlim) );
        if ~isempty(loctn)
            trimpts = [trimpts; indices(loctn) ones(size(loctn)) yvec(indices(loctn)) xvec(indices(loctn))];
			xvec(indices(loctn)) = max(xlim);
        end
	
        %  Points off the left
        loctn = find( yvec(indices) < min(ylim) );
        if ~isempty(loctn)
            trimpts = [trimpts; indices(loctn) ones(size(loctn)) yvec(indices(loctn)) xvec(indices(loctn))];
            yvec(indices(loctn)) = min(ylim);
        end
	
        %  Points off the right
        loctn = find( yvec(indices) > max(ylim) );
        if ~isempty(loctn)
            trimpts = [trimpts; indices(loctn) ones(size(loctn)) yvec(indices(loctn)) xvec(indices(loctn))];
            yvec(indices(loctn)) = max(ylim);
        end
    end
end
