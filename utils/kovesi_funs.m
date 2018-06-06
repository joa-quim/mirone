function  varargout = kovesi_funs(opt,varargin)

% Copyright (c) 2012-2014 Peter Kovesi
% Centre for Exploration Targeting
% The University of Western Australia
% peter.kovesi at uwa edu au
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in 
% all copies or substantial portions of the Software.
%
% The Software is provided "as is", without warranty of any kind.

	if (nargout)
		[varargout{1:nargout}] = feval(opt, varargin{:});
	else
		feval(opt, varargin{:});
	end

% ------------------------------------------------------------------------------------------
function [dim, mask] = ppdrc(im, wavelength, clip, n)
% PPDRC Phase Preserving Dynamic Range Compression
%
% Generates a series of dynamic range compressed images at different scales.
% This function is designed to reveal subtle features within high dynamic range
% images such as aeromagnetic and other potential field grids. Often this kind
% of data is presented using histogram equalisation in conjunction with a
% rainbow colourmap. A problem with histogram equalisation is that the contrast
% amplification of a feature depends on how commonly its data value occurs,
% rather than on the amplitude of the feature itself. In addition, the use of a
% rainbow colourmap can introduce undesirable perceptual distortions.
%
% Phase Preserving Dynamic Range Compression allows subtle features to be
% revealed without these distortions. Perceptually important phase information
% is preserved and the contrast amplification of anomalies in the signal is
% purely a function of their amplitude. It operates as follows: first a highpass
% filter is applied to the data, this controls the desired scale of analysis.
% The 2D analytic signal of the data is then computed to obtain local phase and
% amplitude at each point in the image. The amplitude is attenuated by adding 1
% and then taking its logarithm, the signal is then reconstructed using the
% original phase values.
%
% Usage: dim = ppdrc(im, wavelength, clip, savename, n)
%
% Arguments:      im - Image to be processed (can contain NaNs)
%         wavelength - Array of wavelengths, in pixels, of the  cut-in
%                      frequencies to be used when forming the highpass
%                      versions of the image.  Try a range of values starting
%                      with, say, a wavelength corresponding to half the size
%                      of the image and working down to something like 50
%                      grid units. 
%               clip - Percentage of output image histogram to clip.  Only a
%                      very small value should be used, say 0.01 or 0.02, but 
%                      this can be beneficial.  Defaults to 0.01%
%           savename - (optional) Basename of filname to be used when saving
%                      the output images.  Images are saved as
%                      'basename-n.png' where n is the highpass wavelength
%                      for that image .  You will be prompted to select a
%                      folder to save the images in. 
%                  n - Order of the Butterworth high pass filter.  Defaults
%                      to 2
%
% Returns:       dim - dynamic range reduced image.
%
% Note when specifying the array 'wavelength' it is suggested that you use
% wavelengths that increase in a geometric series.  You can use the function
% GEOSERIES to conveniently do this
% 
% Reference:
% Peter Kovesi, "Phase Preserving Tone Mapping of Non-Photographic High Dynamic
% Range Images".  Proceedings: Digital Image Computing: Techniques and
% Applications 2012 (DICTA 2012). Available via IEEE Xplore

	if (nargin == 1)
		wavelength = floor(max(size(im)) / 2);
		clip = 0.01;		n = 2;
	elseif (nargin == 2)
		clip = 0.01;		n = 2;
	elseif (nargin == 3)
		n = 2;
	end

	% Identify no-data regions in the image (assummed to be marked by NaN
	% values). These values are filled by a call to fillnan() when passing image
	% to highpassmonogenic.  While fillnan() is a very poor 'inpainting' scheme
	% it does keep artifacts at the boundaries of no-data regions fairly small.
	mask = isnan(im);
	if (any(mask(:)))
		[ph, lix, E] = highpassmonogenic(fillnan(im), wavelength, n);
		dim = histtruncate(sin(ph) .* log(1+E), clip, clip);	% Use log(1+x) instead of log1p(x) because it doesn't exist in R13
		dim(mask) = NaN;
	else
		[ph, lix, E] = highpassmonogenic(im, wavelength, n);
		dim = histtruncate(sin(ph) .* log(1+E), clip, clip);
	end
	
	if (isa(im, 'single')),			dim = single(dim);
	elseif (isa(im, 'int32')),		dim = int32(dim);
	elseif (isa(im, 'uint32')),		dim = uint32(dim);
	elseif (isa(im, 'int16')),		dim = int16(dim);
	elseif (isa(im, 'uint16')),		dim = uint16(dim);
	elseif (isa(im, 'uint8')),		dim = scaleto8(dim);
	end

% ------------------------------------------------------------------------------------------
function [newim, mask] = fillnan(im)
% FILLNAN  - fills NaN values in an image with closest non Nan value
%
% NaN values in an image are replaced with the value in the closest pixel that
% is not a NaN.  This can be used as a crude (but quick) 'inpainting' function
% to allow a FFT to be computed on an image containing NaN values.  While the
% 'inpainting' is very crude it is typically good enough to remove most of the
% edge effects one might get at the boundaries of the NaN regions.  The NaN
% regions should then be remasked out of the final processed image.
%
% Usage:  [newim, mask] = fillnan(im);
%
%   Argument:  im    - Image to be 'filled'
%   Returns:   newim - Filled image
%              mask  - Binary image indicating NaN regions in the original image

	% Generate distance transform from non NaN regions of the image. 
	% L will contain indices of closest non NaN points in the image
	mask = ~isnan(im);   

	[D,L] = img_fun('bwdist', mask);	clear D
	ind = find(isnan(im));  % Indices of points that are NaN

	% Fill NaN locations with value of closest non NaN pixel
	newim = im;
	newim(ind) = im(L(ind));

% ------------------------------------------------------------------------------------------
function [phase, orient, E, f, h1f, h2f] = highpassmonogenic(im, maxwavelength, n)
% HIGHPASSMONOGENIC Compute phase and amplitude on highpass images via monogenic filters
%
% Usage: [phase, orient, E, f, h1f, h2f] = highpassmonogenic(im, maxwavelength, n)
%
% Arguments:            im - Image to be processed.
%            maxwavelength - Wavelength in pixels of the  cut-in frequency
%                            of the Butterworth highpass filter. 
%                        n - The order of the Butterworth filter. This is an
%                            integer >= 1.  The higher the value the sharper
%                            the cutoff.  I generaly do not use a value
%                            higher than 2 to avoid ringing artifacts
%
% Returns:           phase - The local phase. Values are between -pi/2 and pi/2
%                   orient - The local orientation. Values between -pi and pi.
%                            Note that where the local phase is close to
%                            +-pi/2 the orientation will be poorly defined.
%                        E - Local energy, or amplitude, of the signal.
%
%   Note that maxwavelength can be an array in which case the outputs will all
%   be cell arrays with an element for each corresponding maxwavelength value.

	if (ndims(im) == 3),	im = img_fun('rgb2gray', im);	end

	[rows,cols] = size(im);    
	IM = fft2(double(im));
%    IM = perfft2(double(im));  % Periodic fft can be useful but I think we
                                % may want to add the smooth component of the
                                % decomposition back in after the
                                % filtering. Without this one can get edge artifcats
                                % The problem is calculating the periodic fft
                                % + computing the smooth component in the
                                % spatial domain involves 3 ffts

	% Generate horizontal and vertical frequency grids
	[radius, u1, u2] = filtergrid(rows,cols);

	% Get rid of the 0 radius value in the middle (at top left corner after
	% fftshifting) so that dividing by the radius, will not cause trouble.
	radius(1,1) = 1;

	H1 = 1i*u1 ./ radius;   % The two monogenic filters in the frequency domain
	H2 = 1i*u2 ./ radius;
	H1(1,1) = 0;	H2(1,1) = 0;
	radius(1,1) = 0;		% undo fudge
	clear('u1', 'u2');

	% High pass Butterworth filter
	H =  1 - 1 ./ (1 + (radius * maxwavelength).^(2*n));        
	clear('radius');
	IM = IM.*H;    clear('H');

	f = real(ifft2(IM));
	h1f = real(ifft2(H1.*IM));  clear('H1');
	h2f = real(ifft2(H2.*IM));  clear('H2', 'IM');

	phase = atan(f./sqrt(h1f.^2+h2f.^2 + eps));
	orient = atan2(h2f, h1f);
	E = sqrt(f.^2 + h1f.^2 + h2f.^2);

% ------------------------------------------------------------------------------------------
function [radius, u1, u2] = filtergrid(rows, cols)
% FILTERGRID Generates grid for constructing frequency domain filters
%
% Usage:  [radius, u1, u2] = filtergrid(rows, cols)
%         [radius, u1, u2] = filtergrid([rows, cols])
%
% Arguments:  rows, cols - Size of image/filter
%
% Returns:        radius - Grid of size [rows cols] containing normalised
%                          radius values from 0 to 0.5.  Grid is quadrant
%                          shifted so that 0 frequency is at radius(1,1)
%                 u1, u2 - Grids containing normalised frequency values
%                          ranging from -0.5 to 0.5 in x and y directions
%                          respectively. u1 and u2 are quadrant shifted.

	% Handle case where rows, cols has been supplied as a 2-vector
	if nargin == 1 && length(rows) == 2  
		tmp = rows;
		rows = tmp(1);
		cols = tmp(2);
	end

	% Set up X and Y spatial frequency matrices, u1 and u2 The following code
	% adjusts things appropriately for odd and even values of rows and columns
	% so that the 0 frequency point is placed appropriately.  See
	% https://blogs.uoregon.edu/seis/wiki/unpacking-the-matlab-fft/
	if mod(cols,2)
		u1range = (-(cols-1)/2:(cols-1)/2)/cols;
	else
		u1range = (-cols/2:(cols/2-1))/cols; 
	end

	if mod(rows,2)
		u2range = (-(rows-1)/2:(rows-1)/2)/rows;
	else
		u2range = (-rows/2:(rows/2-1))/rows; 
	end

	[u1,u2] = meshgrid(u1range, u2range);

	% Quadrant shift so that filters are constructed with 0 frequency at the corners
	u1 = ifftshift(u1);
	u2 = ifftshift(u2);

	% Construct spatial frequency values in terms of normalised radius from centre. 
	radius = sqrt(u1.^2 + u2.^2);

% ------------------------------------------------------------------------------------------
function [newim, sortv] = histtruncate(im, lHistCut, uHistCut, sortv)
% HISTTRUNCATE - Truncates ends of an image histogram.
%
% Function truncates a specified percentage of the lower and
% upper ends of an image histogram.
%
% This operation allows grey levels to be distributed across
% the primary part of the histogram.  This solves the problem
% when one has, say, a few very bright values in the image which
% have the overall effect of darkening the rest of the image after
% rescaling.
%
% Usage: 
%    [newim, sortv] = histtruncate(im, lHistCut, uHistCut)
%    [newim, sortv] = histtruncate(im, lHistCut, uHistCut, sortv)
%
% Arguments:
%    im          -  Image to be processed
%    lHistCut    -  Percentage of the lower end of the histogram to saturate.
%    uHistCut    -  Percentage of the upper end of the histogram
%                   to saturate. If omitted or empty defaults to the value for lHistCut.
%    sortv       -  Optional array of sorted image pixel values obtained
%                   from a previous call to histtruncate.  Supplying this
%                   data speeds the operation of histtruncate when one is
%                   repeatedly varying lHistCut and uHistCut.
%
% Returns:
%    newim       -  Image with values clipped at the specified histogram
%                   fraction values.  If the input image was colour the
%                   lightness values are clipped and stretched to the range
%                   0-1.  If the input image is greyscale no stretching is
%                   applied. You may want to use NORMAALISE to achieve this
%    sortv       -  Sorted image values for reuse in subsequent calls to histruncate.

	if ~exist('uHistCut', 'var') || isempty(uHistCut), uHistCut = lHistCut; end

	if lHistCut < 0 || lHistCut > 100 || uHistCut < 0 || uHistCut > 100
		error('Histogram truncation values must be between 0 and 100');
	end

	if ~exist('sortv', 'var'), sortv = []; end

	if ndims(im) == 3  % Assume colour image in RGB
		hsv = rgb2hsv(im);     % Convert to HSV 
		% Apply histogram truncation just to intensity component
		[hsv(:,:,3), sortv] = Ihisttruncate(hsv(:,:,3), lHistCut, uHistCut, sortv);

		% Stretch intensity component to 0-1
		hsv(:,:,3) = normalise(hsv(:,:,3));
		newim = hsv2rgb(hsv);  % Convert back to RGB
	else
		[newim, sortv] = Ihisttruncate(im, lHistCut, uHistCut, sortv);
	end

%-----------------------------------------------------------------------
% Internal function that does the work
%-----------------------------------------------------------------------
    
function [im, sortv] = Ihisttruncate(im, lHistCut, uHistCut, sortv)
    
	if ndims(im) > 2
		error('HISTTRUNCATE only defined for grey value images');
	end

	% Generate a sorted array of pixel values or use supplied values
	if isempty(sortv)
		sortv = sort(im(:));
	end

	% Any NaN values end up at the end of the sorted list. We need to
	% eliminate these
	sortv = sortv(~isnan(sortv));
	N = length(sortv(:));

	% Compute indicies corresponding to specified upper and lower fractions
	% of the histogram.
	lind = floor(1 + N*lHistCut/100);
	hind = ceil(N - N*uHistCut/100);

	low_in  = sortv(lind);
	high_in = sortv(hind);

	% Adjust image
	im(im < low_in) = low_in;
	im(im > high_in) = high_in;

