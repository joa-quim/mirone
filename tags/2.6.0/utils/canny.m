function [eout,thresh] = canny(varargin)
%CANNY Find edges in intensity image (cvlib_mex powered version).
%
% Re-written version to work in singles by the help of cvlib_mex
%
%   BW = CANNY(I,THRESH) specifies sensitivity thresholds for the
%   Canny method. THRESH is a two-element vector in which the first element
%   is the low threshold, and the second element is the high threshold. If
%   you specify a scalar for THRESH, this value is used for the high
%   threshold and 0.4*THRESH is used for the low threshold. If you do not
%   specify THRESH, or if THRESH is empty ([]), EDGE chooses low and high 
%   values automatically.
%
%   BW = CANNY(I,THRESH,SIGMA) specifies the Canny method, using
%   SIGMA as the standard deviation of the Gaussian filter. The default
%   SIGMA is 1; the size of the filter is chosen automatically, based on SIGMA. 
%
%   [BW,thresh] = CANNY(I,...) returns the threshold values as a two-element vector.


[a,thresh,sigma] = parse_inputs(varargin{:});
is_indexed = false;

% Transform to a double precision intensity image if necessary
if (isa(a,'double') || isa(a,'logical'))
	a = single(a);
elseif isa(a,'uint8')
	if (~is_indexed)
		a = single(a);
		cvlib_mex('CvtScale', a, 1. / 255);
	else
		cvlib_mex('addS',single(a), 1);
	end

elseif isa(a,'uint16')
	if (~is_indexed)
		a = single(a) / 65535;
	else
		a = single(a) + 1;
	end
  
elseif isa(a,'int16')
	if (~is_indexed)
		a = (single(a) + 32768) / 65535;
	else
		error('Canny:invalidIndexedImage', 'An indexed image can be uint8, uint16, double, single, or logical.');
	end

end

	[m,n] = size(a);

	% The output edge map:
	eout = false(m,n);

	% Magic numbers
	GaussianDieOff = .0001;  
	PercentOfPixelsNotEdges = .7; % Used for selecting thresholds
	ThresholdRatio = .4;          % Low thresh is this fraction of the high.
	
	% Design the filters - a gaussian and its derivative
	
	pw = 1:30; % possible widths
	ssq = sigma^2;
	width = find(exp(-(pw.*pw)/(2*ssq)) > GaussianDieOff);
	if isempty(width)
		width = 1;  % the user entered a really small sigma
	else
		width = width(end);
	end
	
	t = (-width:width);
	gau = exp(-(t.*t)/(2*ssq))/(2*pi*ssq);     % the gaussian 1D filter
	
	% Find the directional derivative of 2D Gaussian (along X-axis)
	% Since the result is symmetric along X, we can get the derivative along
	% Y-axis simply by transposing the result for X direction.
	[x,y]  = meshgrid(-width:width,-width:width);
	dgau2D = -x.*exp(-(x.*x+y.*y)/(2*ssq))/(pi*ssq);
	
	% Convolve the filters with the image in each direction
	% The canny edge detector first requires convolution with
	% 2D gaussian, and then with the derivitave of a gaussian.
	% Since gaussian filter is separable, for smoothing, we can use 
	% two 1D convolutions in order to achieve the effect of convolving
	% with 2D Gaussian.  We convolve along rows and then columns.
	
	%smooth the image out
	aSmooth = img_fun('imfilter',a,gau,'conv','replicate');   % run the filter accross rows
	aSmooth = img_fun('imfilter',aSmooth,gau','conv','replicate'); % and then accross columns
	
	%apply directional derivatives
	ax = img_fun('imfilter',aSmooth, dgau2D, 'conv','replicate');
	ay = img_fun('imfilter',aSmooth, dgau2D', 'conv','replicate');
	clear aSmooth;
	
	mag = cvlib_mex('hypot', ax, ay);
	magmax = max(mag(:));
	if magmax > 0
		mag = cvlib_mex('CvtScale', mag, 1. / double(magmax));
	end
  
	% Select the thresholds
	if isempty(thresh) 
		counts = double(imhistc(mag, 64, 1, 1));   % Call MEX file to do work.
		highThresh = find(cumsum(counts) > PercentOfPixelsNotEdges*m*n);
		highThresh = highThresh(1) / 64;
		lowThresh = ThresholdRatio*highThresh;
		thresh = single([lowThresh highThresh]);
	elseif length(thresh)==1
		highThresh = thresh;
		if (thresh >= 1)
			error('Canny:thresholdMustBeLessThanOne','The threshold must be less than 1.');
		end
		lowThresh = ThresholdRatio*thresh;
		thresh = [lowThresh highThresh];
	elseif length(thresh)==2
		lowThresh = thresh(1);
		highThresh = thresh(2);
		if (lowThresh >= highThresh) || (highThresh >= 1)
			error('Canny:thresholdOutOfRange','Thresh must be [low high], where low < high < 1.');
		end
	end
  
	% The next step is to do the non-maximum supression.  
	% We will accrue indices which specify ON pixels in strong edgemap
	% The array e will become the weak edge map.
	idxStrong = [];  
	for (dir = 1:4)
		idxLocalMax = cannyFindLocalMaxima(dir,ax,ay,mag);
		idxWeak = idxLocalMax(mag(idxLocalMax) > lowThresh);
		eout(idxWeak)=1;
		idxStrong = [idxStrong; idxWeak(mag(idxWeak) > highThresh)];
	end
	
	if ~isempty(idxStrong) % result is all zeros if idxStrong is empty
		rstrong = rem(idxStrong-1, m)+1;
		cstrong = floor((idxStrong-1)/m)+1;
		eout = img_fun('bwselect',eout, cstrong, rstrong, 8);
		eout = img_fun('bwmorph',eout, 'thin', 1);  % Thin double (or triple) pixel wide contours
	end

% ---------------------------------------------------------------------
function idxLocalMax = cannyFindLocalMaxima(direction,ix,iy,mag)
%
% This sub-function helps with the non-maximum supression in the Canny
% edge detector.  The input parameters are:
% 
%   direction - the index of which direction the gradient is pointing, 
%               read from the diagram below. direction is 1, 2, 3, or 4.
%   ix        - input image filtered by derivative of gaussian along x 
%   iy        - input image filtered by derivative of gaussian along y
%   mag       - the gradient magnitude image
%
%    there are 4 cases:
%
%                         The X marks the pixel in question, and each
%         3     2         of the quadrants for the gradient vector
%       O----0----0       fall into two cases, divided by the 45 
%     4 |         | 1     degree line.  In one case the gradient
%       |         |       vector is more horizontal, and in the other
%       O    X    O       it is more vertical.  There are eight 
%       |         |       divisions, but for the non-maximum supression  
%    (1)|         |(4)    we are only worried about 4 of them since we 
%       O----O----O       use symmetric points about the center pixel.
%        (2)   (3)        


[m,n] = size(mag);

% Find the indices of all points whose gradient (specified by the 
% vector (ix,iy)) is going in the direction we're looking at.  

switch direction
	case 1
		cvlib_mex('CvtScale',iy,-1);      % Permanently make iy = -iy;
		idx = find((iy > 0 & ix > iy) | (iy < 0 & ix < iy));
	case 2
		idx = find((ix > 0 & iy >= ix)  | (ix < 0 & iy <= ix));
		% Undo the above iy = -iy. I have to do this due to the inplace nature of
		% operation done by cvlib_mex.
		cvlib_mex('CvtScale',iy,-1);
	case 3
		idx = find((ix<=0 & ix>iy) | (ix>=0 & ix<iy));
	case 4
		idx = find((iy<0 & ix<=iy) | (iy>0 & ix>=iy));
end

% Exclude the exterior pixels
if ~isempty(idx)
	v = mod(idx,m);
	extIdx = (v == 1 | v == 0 | idx <= m | (idx > (n-1)*m));
	idx(extIdx) = [];
end

ixv = ix(idx);  
iyv = iy(idx);   
gradmag = mag(idx);

% Do the linear interpolations for the interior pixels
switch direction
	case 1
		cvlib_mex('div',iyv, ixv)
		d = cvlib_mex('abs',iyv);  d_1 = cvlib_mex('subRS',d,1);     % 1 - d
		tmp = cvlib_mex('mul',mag(idx+m),d_1);  gradmag1 = cvlib_mex('mul',mag(idx+m-1),d);
		cvlib_mex('add',gradmag1,tmp);
		tmp = cvlib_mex('mul',mag(idx-m),d_1);  gradmag2 = cvlib_mex('mul',mag(idx-m+1),d);
		cvlib_mex('add',gradmag2,tmp);
	case 2
		cvlib_mex('div',ixv, iyv)
		d = cvlib_mex('abs',ixv);  d_1 = cvlib_mex('subRS',d,1);     % 1 - d
		tmp = cvlib_mex('mul',mag(idx-1),d_1);  gradmag1 = cvlib_mex('mul',mag(idx+m-1),d);
		cvlib_mex('add',gradmag1,tmp);
		tmp = cvlib_mex('mul',mag(idx+1),d_1);  gradmag2 = cvlib_mex('mul',mag(idx-m+1),d);
		cvlib_mex('add',gradmag2,tmp);
	case 3
		cvlib_mex('div',ixv, iyv)
		d = cvlib_mex('abs',ixv);  d_1 = cvlib_mex('subRS',d,1);     % 1 - d
		tmp = cvlib_mex('mul',mag(idx-1),d_1);  gradmag1 = cvlib_mex('mul',mag(idx-m-1),d);
		cvlib_mex('add',gradmag1,tmp);
		tmp = cvlib_mex('mul',mag(idx+1),d_1);  gradmag2 = cvlib_mex('mul',mag(idx+m+1),d);
		cvlib_mex('add',gradmag2,tmp);
	case 4
		cvlib_mex('div',iyv, ixv)
		d = cvlib_mex('abs',iyv);  d_1 = cvlib_mex('subRS',d,1);     % 1 - d
		tmp = cvlib_mex('mul',mag(idx-m),d_1);  gradmag1 = cvlib_mex('mul',mag(idx-m-1),d);
		cvlib_mex('add',gradmag1,tmp);
		tmp = cvlib_mex('mul',mag(idx+m),d_1);  gradmag2 = cvlib_mex('mul',mag(idx+m+1),d);
		cvlib_mex('add',gradmag2,tmp);
end
idxLocalMax = idx(gradmag >= gradmag1 & gradmag >= gradmag2); 

% ---------------------------------------------------------------------
function [I,Thresh,Sigma] = parse_inputs(varargin)

	error(nargchk(1,3,nargin));

	I = varargin{1};
	if (islogical(I)),	I = uint8(I);	end
	if (ndims(I) == 3),	I = cvlib_mex('color', I, 'rgb2gray');		end
	iptcheckinput(I,{'numeric'},{'nonsparse','2d'},mfilename,'I',1);

	Thresh=[];

	% Now parse the nargin-1 remaining input arguments

	Sigma = 1.0;          % Default Std dev of gaussian for canny
	threshSpecified = 0;  % Threshold is not yet specified
  
	for i = 2:nargin
		if numel(varargin{i})==2 && ~threshSpecified
			Thresh = varargin{i};
			threshSpecified = 1;
		elseif numel(varargin{i})==1 
			if ~threshSpecified
				Thresh = varargin{i};
				threshSpecified = 1;
			else
				Sigma = varargin{i};
			end
		elseif isempty(varargin{i}) && ~threshSpecified
			threshSpecified = 1;
		else
			error('Canny:invalidInputArguments', 'Invalid input arguments');
		end
	end
  
	if Sigma <= 0
		error('Canny:sigmaMustBePositive','Sigma must be positive');
	end
