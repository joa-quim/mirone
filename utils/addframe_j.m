function [aviobj] = addframe_j(aviobj,varargin)
%ADDFRAME  Add video frame to AVI file.  
%   Copyright 1984-2002 The MathWorks, Inc.
%   $Revision: 1.9 $  $Date: 2002/06/05 20:04:44 $

numframes = nargin - 1;

for i = 1:numframes
    MovieLength = 1;
    mlMovie = 0;
    % Parse input arguments
    inputType = getInputType(aviobj,varargin{i});
    switch inputType
        case 'axes'
        case 'figure'
        case 'movie'
            mlMovie = 1;
            MovieLength = length(varargin{i});
            if (~isempty(varargin{i}(1).colormap) && isempty(aviobj.Bitmapheader.Colormap) && ...
                    aviobj.MainHeader.TotalFrames == 0 )
                aviobj.Colormap = varargin{i}(1).colormap;
            end
        case 'data'
            frame = varargin{i};
    end

    for j = 1:MovieLength
        if (mlMovie ),  frame = varargin{i}(j).cdata;   end
        frameClass = class(frame);
        if isempty(strmatch(frameClass,strvcat('double','uint8')))
            error('FRAME must be of either double or uint8 precision');
        end
        
        % Determine image dimensions
        height = size(frame,1); 
        width = size(frame,2);
        dims = size(frame,3);

        % Check requirements for the Intel Indeo codec
        % Intel Indeo requires images dimensions to be a multiple of four,
        % greater than 32, and no more than 4,194,304 pixels.
        if strmatch('iv',lower(aviobj.StreamHeader.fccHandler))
            if (aviobj.MainHeader.TotalFrames == 0) && (aviobj.Bitmapheader.biBitCount == 8) && ...
	            (aviobj.Bitmapheader.biClrUsed >236)
	            error('The colormap can not exceed 236 colors, as specified by the Intel Indeo compressor.');
            end
            hpad = rem(height,4);
            wpad = rem(width,4);
      
            if (width < 32) || (height < 32)
	            error('The minimum frame size for the Indeo compressor is 32x32.');
            end
            if width*height > 4194304
	            error('The Intel Indeo compressor can not compress frame sizes that exceed a maximum frame size of 4,194,304 pixels.');
            end
            if hpad
	            if aviobj.MainHeader.TotalFrames == 0
	                warning('The frame height has been padded to be a multiple of four as required by Intel Indeo.');
	            end
	            frame = [frame;zeros(4-hpad,size(frame,2),dims)];
            end
            if wpad
	            if aviobj.MainHeader.TotalFrames == 0
	                warning('The frame width has been padded to be a multiple of four as required by Intel Indeo.');
	            end
	            frame = [frame, zeros(size(frame,1),4-wpad,dims)];
            end

            % Determine adjusted image dimensions
            height = size(frame,1);
            width = size(frame,2);
            dims = size(frame,3);
        end
    
        % Truecolor images can not be compressed with RLE or MSVC compression 
        if dims == 3
            msg = 'Use a compression method other than RLE or MSVC for truecolor images.';
            if strmatch(lower(aviobj.StreamHeader.fccHandler),'mrle') 
	            error(msg);
            elseif strmatch(lower(aviobj.StreamHeader.fccHandler),'msvc')
	            error(msg);
            end
        end
    
        % If this is not the first frame, make sure it is consistent
        if aviobj.MainHeader.TotalFrames ~= 0
            ValidateFrame(aviobj,width, height,dims);
        end

        % Reshape image data
        frame = ReshapeImage(frame);

        % Compute memory requirements for frame storage
        numFrameElements = prod(size(frame));

        % If this is the first frame, set necessary fields
        if aviobj.MainHeader.TotalFrames==0
            aviobj.MainHeader.SuggestedBufferSize = numFrameElements;
            aviobj.StreamHeader.SuggestedBufferSize = numFrameElements;
            aviobj.MainHeader.Width = width;
            aviobj.MainHeader.Height = height;
            aviobj.Bitmapheader.biWidth = width;
            aviobj.Bitmapheader.biHeight = height;
            aviobj.Bitmapheader.biSizeImage = numFrameElements;
            if dims == 3 
	            aviobj.Bitmapheader.biBitCount = 24;
            else
	            aviobj.Bitmapheader.biBitCount = 8;
            end
        end

        % On Windows use Video for Windows to write the video stream
        if ispc
            % fps is calculated in avi.c by dividing the rate by the scale (100).
            % The scale of 100 is hard coded into avi.c
            rate = aviobj.StreamHeader.Rate; 
            avi('addframe',rot90(frame,-1), aviobj.Bitmapheader, ...
	        aviobj.MainHeader.TotalFrames,rate, ...
	        aviobj.StreamHeader.Quality,aviobj.FileHandle, ...
	        aviobj.StreamName,aviobj.KeyFrameEveryNth);
        end
        % Update the total frames
        aviobj.MainHeader.TotalFrames = aviobj.MainHeader.TotalFrames + 1;
    end
end

% ------------------------------------------------------------------------
function ValidateFrame(aviobj, width, height, dims)
% VALIDATEFRAME
%   Verify the frame is consistent with header information in AVIOBJ.  The
%   frame must have the same WIDTH, HEIGHT, and DIMS as the previous frames.

if width ~= aviobj.MainHeader.Width
    error(sprintf('Frame must be %d by %d.', aviobj.MainHeader.Width,aviobj.MainHeader.Height))
elseif height ~= aviobj.MainHeader.Height
    error(sprintf('Frame must be %d by %d.', aviobj.MainHeader.Width,aviobj.MainHeader.Height))
end

if (aviobj.Bitmapheader.biBitCount == 24) && (dims ~= 3)
    error('Frame must be a truecolor image.');
elseif (aviobj.Bitmapheader.biBitCount == 8) && (dims ~= 1)
    error('Frame must be an indexed image.')
end

% ------------------------------------------------------------------------
function X = ReshapeImage(X)
numcomps = size(X,3);

if (isa(X,'double'))
    if (numcomps == 3),  X = uint8(round(255*X));
    else                        X = uint8(X-1);    end
end

% Squeeze 3rd dimension into second
if (numcomps == 3)
    X = X(:,:,[3 2 1]);
    X = permute(X, [1 3 2]);
    X = reshape(X, [size(X,1) size(X,2)*size(X,3)]);
end

width = size(X,2);
tmp = rem(width,4);
if (tmp > 0)
    padding = 4 - tmp;
    X = cat(2, X, repmat(uint8(0), [size(X,1) padding]));
end

% ------------------------------------------------------------------------
function inputType = getInputType(aviobj,frame)
if ishandle(frame)
    inputType = get(frame,'type');
elseif isstruct(frame) && isfield(frame,'cdata')
    inputType = 'movie';
elseif isa(frame,'numeric')
    inputType = 'data';
else
    error('Invalid input argument.  Each frame must be a numeric matrix, a MATLAB movie structure, or a handle to a figure or axis.');
end
