function new_pal = cdf2pal(z_min,z_max,z_mean,z_std,pal,n)
% Use the normalized CDF (cumulative distribution function) created by cdf_cpt to
% re-compute PAL as a continuous-color-palette, with a non-linear histogram-equalized
% The idea is beased on grd2cpt. If n == 64 the WS "stupid bug fix" may not be enough

if (nargin == 5),   n = 32;     end
cdf = cdf_cpt(z_min,z_max,z_mean,z_std,n);

n_int = length(cdf) - 1;            % Number of intervals
xlin = linspace(0,1,n_int+1)';
dif_lin = diff(xlin);
dif_cdf = diff(cdf);
n_colors = round(256 / n_int);      % Number of colors in each interval

pal_split = cell(1,n_int-1);
for i=1:n_int-1
    pal_split{i} = pal( ((i-1)*n_colors+1:i*n_colors),:);
end
pal_split{n_int} = pal(i*n_colors:end,:);

f = dif_cdf ./ dif_lin;
new_pal = cell(1,n_int);
for i=1:n_int
    n = round(f(i)*length(pal_split{i}));
    x = linspace(0,dif_lin(i), length(pal_split{i}));
    xi = linspace(0,dif_lin(i),n);
    new_pal{i} = interp1(x,pal_split{i},xi,'linear','extrap');
end

new_pal = cat(1,new_pal{:});
if (length(new_pal) ~= 256)
    x = linspace(1,length(new_pal),256)';
    new_pal = interp1(new_pal,x);
end

% Insure that all colormap intensities are be between 0 and 1.
new_pal(new_pal < 0) = 0;
new_pal(new_pal > 1) = 1;

%-------------------------------------------------------------------------------
function cdf = cdf_cpt(z_min,z_max,z_mean,z_std,n)
% Compute the normalized CDF (cumulative distribution function)
step = 1 / n;
x = linspace(-3,3,2000)';
p = 0.5 * erfcore(-x ./ sqrt(2),1);
pp = (step:step:1-step)';
p_int = interp1(p,x,pp);

% This is what WS called "stupid bug fix" in grd2cpt (may be not enough if n == 64)
if ((z_mean + p_int(1)*z_std) <= z_min || (z_mean + p_int(end)*z_std) >= z_max)
    z_mean = 0.5 * (z_min + z_max);
    z_std = (z_max - z_mean) / 1.5;
end

cdf = z_mean + p_int * z_std;
cdf = [z_min; cdf; z_max];

% Now normalize the cdf distribution 
b = 1 / (cdf(end) - cdf(1));
a = -cdf(1) * b;
cdf = a + b * cdf;
cdf(1) = 0;     cdf(end) = 1;   % Prevent roundoff errors
