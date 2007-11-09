function [Z,Z_rect,handles] = roi_filtering(handles, Z, head, Z_rect, r_c, mask, border)
 
is_rect = 0;        done = 0;
if (ischar(mask)),   is_rect = 1;    end
if (nargin == 6),   border = 'no';  end

prompt = {'Enter number of filter rows' ,'Enter number of filter cols', 'Maximum allowed change (z units)'};
def = {num2str(3) num2str(3) num2str(10)};
resp  = inputdlg(prompt,'Median Filtering',[1 30; 1 30; 1 30],def);    pause(0.01)
if isempty(resp);    set(handles.figure1,'pointer','arrow');    return;     end
Z_rect = double(Z_rect);      % It has to be

have_local_nans = 0;
if (handles.have_nans)              % If the grid has NaNs
    rect_nans = isnan(Z_rect);
    if (any(rect_nans(:)))             % Check if we have NaNs inside the Z_rect rectangle
        have_local_nans = 1;
    else
        clear rect_nans;
    end
end

Z_filt = img_fun('medfilt2',Z_rect,[str2double(resp{1}) str2double(resp{2})]);
max_z_cut = str2double(resp{3});
idx_dif_p = find((Z_rect - Z_filt) > max_z_cut);
idx_dif_n = find((Z_rect - Z_filt) < -max_z_cut);
Z_filt(idx_dif_p) = Z_rect(idx_dif_p) - max_z_cut;
Z_filt(idx_dif_n) = Z_rect(idx_dif_n) + max_z_cut;
clear idx_dif_p idx_dif_n;

if (is_rect)
    if (strcmp(border,'no'))     % We don't need to do border re-interpolation
		hfc = fix(str2double(resp{1})/2);   % Half filter width. Don't use the filtered data on beyond
		hfr = fix(str2double(resp{2})/2);   % the transition zone. It may produce sharp edges.
		if (isa(Z,'single'))
			Z_rect = single(Z_filt);		clear Z_filt;
		elseif (isa(Z,'int16'))
			Z_rect = int16(Z_filt);			clear Z_filt;
		elseif (isa(Z,'uint16'))
			Z_rect = uint16(Z_filt);		clear Z_filt;
		else		% Should never come here
			Z_rect = single(Z_filt);		clear Z_filt;
		end
		Z(r_c(1)+hfr:r_c(2)-hfr,r_c(3)+hfc:r_c(4)-hfc) = Z_rect(hfr+1:end-hfr,hfc+1:end-hfc);
        done = 1;
    else            % User asked for border re-interpolation. We must create a mask grid
        mask = true(size(Z_rect));
    end
end

% handles.Z_back = Z(r_c(1):r_c(2),r_c(3):r_c(4));    % For the undo op
% handles.r_c = r_c;

if (done),  return;     end     % Ractangle without border re-interpolation

hf = 3;

mask_b = mask_border(mask,hf);
mask_b = mask_b(:);

X = (head(1) + (r_c(3)-1)*head(8)):head(8):(head(1) + (r_c(4)-1)*head(8));
Y = (head(3) + (r_c(1)-1)*head(9)):head(9):(head(3) + (r_c(2)-1)*head(9));
opt_R = sprintf('-R%.10f/%.10f/%.10f/%.10f', X(1), X(end), Y(1), Y(end));
opt_I = sprintf('-I%.10f/%.10f',head(8),head(9));
Z_filt(~mask) = Z_rect(~mask);      % Outside mask zone Z_filt is set to original value
Z_filt = Z_filt(:);     Z_filt(mask_b) = [];
[X,Y] = meshgrid(X,Y);
X = X(:);   X(mask_b) = [];
Y = Y(:);   Y(mask_b) = [];
Z_rect = surface_m(X,Y,Z_filt,opt_R,opt_I);    clear X Y Z_filt;
if (have_local_nans)            % If we had NaNs inside original Z_rect, restore them
    Z_rect(rect_nans) = NaN;
end

if (isa(Z,'single'))
	Z(r_c(1):r_c(2),r_c(3):r_c(4)) = single(Z_rect);
elseif (isa(Z,'int16'))
	Z(r_c(1):r_c(2),r_c(3):r_c(4)) = int16(Z_rect);
elseif (isa(Z,'uint16'))
	Z(r_c(1):r_c(2),r_c(3):r_c(4)) = uint16(Z_rect);
else		% Should never come here
	Z(r_c(1):r_c(2),r_c(3):r_c(4)) = single(Z_rect);
end

% -----------------------------------------------------------------------
function mask_b = mask_border(mask,w)
% Takes the MASK matrix which has ones defining the ROI and computes
% MASK_B that is composed of ones at the border of the ROI and W nodes
% to the inside

	[m,n] = size(mask);
	mask_b = [repmat(logical(0),m,w) mask(:,1:end-w)];      % shift mask rigth
	mask_in = mask & mask_b;
	mask_b = [mask(:,w+1:end) repmat(logical(0),m,w)];      % shift mask left
	mask_in = mask_in & mask_b;
	mask_b = [mask(w+1:end,:); repmat(logical(0),w,n)];      % shift mask up
	mask_in = mask_in & mask_b;
	mask_b = [repmat(logical(0),w,n); mask(1:end-w,:)];      % shift mask down
	mask_in = mask_in & mask_b;
	mask_b = xor(mask_in,mask);
