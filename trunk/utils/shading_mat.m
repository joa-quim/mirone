function    img = shading_mat(img,R,scale)
% New version that calls a mex function named mex_illuminate
% With this we gained in performance by a factor > 2

if (nargin == 2),    scale = 'scale';    end

[rows,cols,cores] = size(img);
if (strcmp(scale,'scale'))      % Scale R into [-1 1] interval
    Rmax = max(R(:));   Rmin = min(R(:));
    esc = 1 / (Rmax -Rmin);
    R = (-1 + 2*((R-Rmin)*esc)) * 0.95;
end

[r,g,b] = mex_illuminate(img,R);        % mex_illuminate should output img
img = reshape([r g b],rows,cols,3);
