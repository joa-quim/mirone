function [xy_prj, msg, opt_R] = geog2projected_pts(handles,xy,lims,pad)
    % project a 2 column array in geogs to the projection currently stored in 'ProjWKT' or 'ProjGMT'
    %
    % XY , a [lon lat] 2 column array.
    % LIMS, if present & not empty,  must be a [x_min x_max y_min y_max] vector. 
    %       Normally, the image four corners coords
    % PAD is a scalar with its value in projected coords. It is used as a 'skirt' boundary buffer
    % In case we have both a 'ProjGMT' & 'ProjWKT' the later prevails since
    % it is supposed to correctly describe the data projection.
    % If none proj string is found, XY_PRJ = XY since nothing was done
    %
    % On output OPT_R is a -R GMT option string which will be a crude estimation
    % in case of a mapproject conversion (that story of -Rg). We need this OPT_R to
    % use in calls to shoredump_m

    msg = '';       xy_prj = [];    opt_R = [];

    n_arg = nargin;
    if (n_arg < 2)
        msg = 'geog2projected_pts: Bad usage. Need at least 2 arguments in';
        return
    elseif (n_arg < 3)
        lims = [];      pad = 0;
    elseif (n_arg < 3)
        pad = 0;
    end

    % Check that input data is a 2 column vector
    if (size(xy,2) > 2)
        msg = 'geog2projected_pts: ERROR, data input must be a 2 column array';
        return
    end
    
    if (~isempty(lims) && numel(lims) ~= 4)
        msg = 'geog2projected_pts: ERROR, "lims" vector must have 4 elements';
        return
    end
    
    % Fish eventual proj strings
    projGMT = getappdata(handles.figure1,'ProjGMT');
    projWKT = getappdata(handles.axes1,'ProjWKT');

    if (~isempty(projWKT))
        projStruc.SrcProjWKT = projWKT;
        if (~isempty(lims) && pad ~= 0)
        	% Get rid of data that are outside the map limits
            lims = lims + [-pad pad -pad pad];      % Extend limits by PAD value
            ind = (xy(:,1) < lims(1) | xy(:,1) > lims(2));
	        xy(:,ind) = [];
	        ind = (xy(:,2) < lims(3) | xy(:,2) > lims(4));
	        xy(:,ind) = [];
        end
        xy_prj = ogrproj(xy, projStruc);
        if (nargout == 3)
            x_min = min(xy_prj(:,1));        x_max = max(xy_prj(:,1));
            y_min = min(xy_prj(:,2));        y_max = max(xy_prj(:,2));
            opt_R = ['-R' sprintf('%f',x_min) '/' sprintf('%f',x_max) '/' sprintf('%f',y_min) '/' sprintf('%f',y_max)];
        end
    elseif (~isempty(projGMT))
        if (~handles.geog)      % Currently the projGMT has not any inverse projection defined in it
            msg = 'geog2projected_pts: This operation is currently possible only for geographic type data';
            return
        end
        out = mapproject_m(lims,'-R-180/180/0/80','-I','-F',projGMT{:});    % Convert lims back to geogs
        x_min = min(out(:,1));        x_max = max(out(:,1));
        y_min = min(out(:,2));        y_max = max(out(:,2));
      	opt_R = ['-R' sprintf('%f',x_min) '/' sprintf('%f',x_max) '/' sprintf('%f',y_min) '/' sprintf('%f',y_max)];
        if (pad ~= 0)
        	% Get rid of data that are outside the map limits
            lims = lims + [-pad pad -pad pad];      % Extend limits by PAD value
            ind = (xy(:,1) < lims(1) | xy(:,1) > lims(2));
	        xy(:,ind) = [];
	        ind = (xy(:,2) < lims(3) | xy(:,2) > lims(4));
	        xy(:,ind) = [];
        end
        xy_prj = mapproject_m(xy, opt_R, '-F', projGMT{:});
    else
        xy_prj = xy;
        msg = '0';          % Signal that nothing has been done and output = input
        if (nargout == 3)   % Input in geogs, we need opt_R for shoredump_m (presumably)
            opt_R = ['-R' sprintf('%f/%f/%f/%f',lims(1),lims(2),lims(3),lims(4))];
        end
        return
    end
	