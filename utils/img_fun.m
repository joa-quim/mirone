%#mex
function  varargout = img_fun(opt,varargin)

switch opt
    case 'adapthisteq'
        varargout{1} = adapthisteq(varargin{:});
    case 'bwareaopen'
        varargout{1} = bwareaopen(varargin{:});
    case 'bwboundaries'
        varargout{1} = bwboundaries(varargin{:});
	case 'bwdist'
        [varargout{1:nargout}] = bwdist(varargin{:});
    case 'bwmorph'
        varargout{1} = bwmorph(varargin{:});
    case 'bwselect'
        varargout{1} = bwselect(varargin{:});        
    case 'decorrstretch'
        varargout{1} = decorrstretch(varargin{:});
    case 'histeq_j'
        [out,T] = histeq_j(varargin{:});
        varargout{1} = out;  varargout{2} = T;
    case 'hough'
        [h,theta,rho] = hough(varargin{:});
        varargout{1} = h;  varargout{2} = theta;    varargout{3} = rho;
    case 'houghlines'
        varargout{1} = houghlines(varargin{:});
    case 'houghpeaks'
        varargout{1} = houghpeaks(varargin{:});
    case 'im2bw'
        [bw,level] = im2bw(varargin{:});
        varargout{1} = bw;  varargout{2} = level;
    case 'imfilter'
        varargout{1} = imfilter(varargin{:});
    case 'imadjust_j'
        varargout{1} = imadjust_j(varargin{:});
    case 'imresize'
        varargout{1} = imresize(varargin{:});
    case 'fspecial'
        varargout{1} = fspecial(varargin{:});
    case 'fwind1'
        varargout{1} = fwind1(varargin{:});
    case 'edge'
        varargout{1} = edge(varargin{:});
    case 'roifilt2'
        varargout{1} = roifilt2(varargin{:});
    case 'medfilt2'
        varargout{1} = medfilt2(varargin{:});
    case 'padarray'
        varargout{1} = padarray(varargin{:});
    case 'poly2mask'
        varargout{1} = poly2mask(varargin{:});
    case 'roipoly_j'
        varargout{1} = roipoly_j(varargin{:});
    case 'rgb2gray'
        varargout{1} = rgb2gray(varargin{:});
    case 'rgb2ind'
        [X,map] = rgb2ind(varargin{:});
        varargout{1} = X;  varargout{2} = map;
    case 'rangefilt'
        varargout{1} = rangefilt(varargin{:});
    case 'stdfilt'
        varargout{1} = stdfilt(varargin{:});
    case 'find_holes'
        varargout{1} = find_holes(varargin{:});
	otherwise
		if (nargout)
			[varargout{1:nargout}] = feval(opt, varargin{:});
		else
			feval(opt, varargin{:});
		end
end

%----------------------------------------------------------------------------------
function out = adapthisteq(varargin)
%ADAPTHISTEQ Performs Contrast-Limited Adaptive Histogram Equalization (CLAHE).
%   Copyright 1993-2003 The MathWorks, Inc.
%   $Revision: 1.8 $  $Date: 2003/03/25 20:01:48 $

[I, selectedRange, fullRange, numTiles, dimTile, clipLimit, numBins, ...
 noPadRect, distribution, alpha] = parseInputs_adapthisteq(varargin{:});

tileMappings = makeTileMappings(I, numTiles, dimTile, numBins, clipLimit, ...
                                selectedRange, fullRange, distribution, alpha);

%Synthesize the output image based on the individual tile mappings. 
out = makeClaheImage(I, tileMappings, numTiles, selectedRange, numBins, dimTile);

if ~isempty(noPadRect) %do we need to remove padding?
  out = out(noPadRect.ulRow:noPadRect.lrRow, noPadRect.ulCol:noPadRect.lrCol);
end

%----------------------------------------------------------------------------------
function f = adpmedian(g, Smax)
%ADPMEDIAN Perform adaptive median filtering.
%   F = ADPMEDIAN(G, SMAX) performs adaptive median filtering of
%   image G.  The median filter starts at size 3-by-3 and iterates up
%   to size SMAX-by-SMAX. SMAX must be an odd integer greater than 1.

%   Copyright 2002-2004 R. C. Gonzalez, R. E. Woods, & S. L. Eddins
%   Digital Image Processing Using MATLAB, Prentice-Hall, 2004
%   $Revision: 1.5 $  $Date: 2003/11/21 14:19:05 $

% SMAX must be an odd, positive integer greater than 1.
if (Smax <= 1) || (Smax/2 == round(Smax/2)) || (Smax ~= round(Smax))
   error('SMAX must be an odd integer > 1.')
end
[M, N] = size(g);

% Initial setup.
f = g;		f(:) = 0;
alreadyProcessed = false(size(g));

% Begin filtering.
for k = 3:2:Smax
   zmin = ordfilt2(g, 1, ones(k, k), 'symmetric');
   zmax = ordfilt2(g, k * k, ones(k, k), 'symmetric');
   zmed = medfilt2(g, [k k], 'symmetric');
   
   processUsingLevelB = (zmed > zmin) & (zmax > zmed) & ~alreadyProcessed; 
   zB = (g > zmin) & (zmax > g);
   outputZxy  = processUsingLevelB & zB;
   outputZmed = processUsingLevelB & ~zB;
   f(outputZxy) = g(outputZxy);
   f(outputZmed) = zmed(outputZmed);
   
   alreadyProcessed = alreadyProcessed | processUsingLevelB;
   if all(alreadyProcessed(:)),		break,		end
end

% Output zmed for any remaining unprocessed pixels. Note that this zmed was computed
% using a window of size Smax-by-Smax, which is the final value of k in the loop.
f(~alreadyProcessed) = zmed(~alreadyProcessed);

%-----------------------------------------------------------------------------
function tileMappings = makeTileMappings(I, numTiles, dimTile, numBins, clipLimit,...
                     selectedRange, fullRange, distribution, alpha)

numPixInTile = prod(dimTile);
% extract and process each tile
imgCol = 1;
for col=1:numTiles(2)
  imgRow = 1;
  for row=1:numTiles(1)
    tile = I(imgRow:imgRow+dimTile(1)-1,imgCol:imgCol+dimTile(2)-1);
    % for speed, call MEX file directly thus avoiding costly input parsing of imhist
    tileHist = imhistc(tile, numBins, 1, fullRange(2)); 
    tileHist = clipHistogram(tileHist, clipLimit, numBins);
    tileMapping = makeMapping(tileHist, selectedRange, fullRange, ...
                              numPixInTile, numBins, distribution, alpha);
    
    % assemble individual tile mappings by storing them in a cell array;
    tileMappings{row,col} = tileMapping;
    imgRow = imgRow + dimTile(1);    
  end
  imgCol = imgCol + dimTile(2); % move to the next column of tiles
end

%-----------------------------------------------------------------------------
% Calculate the equalized lookup table (mapping) based on cumulating the input 
% histogram.  Note: lookup table is rescaled in the selectedRange [Min..Max].
function mapping = makeMapping(imgHist, selectedRange, fullRange, ...
                               numPixInTile, numBins, distribution, alpha)

mapping = zeros(size(imgHist));
histSum = cumsum(imgHist);
valSpread  = selectedRange(2) - selectedRange(1);

switch distribution
 case 'uniform'
  scale =  valSpread/numPixInTile;
  mapping = min(selectedRange(1) + histSum*scale,selectedRange(2)); %limit to max
  
 case 'rayleigh'	% suitable for underwater imagery
  % pdf = (t./alpha^2).*exp(-t.^2/(2*alpha^2))*U(t)
  % cdf = 1-exp(-t.^2./(2*alpha^2))
  hconst = 2*alpha^2;
  vmax = 1 - exp(-1/hconst);
  val = vmax*(histSum/numPixInTile);
  val(val>=1) = 1-eps; % avoid log(0)
  temp = sqrt(-hconst*log(1-val));
  mapping = min(selectedRange(1)+temp*valSpread,selectedRange(2)); %limit to max
  
 case 'exponential'
  % pdf = alpha*exp(-alpha*t)*U(t)
  % cdf = 1-exp(-alpha*t)
  vmax = 1 - exp(-alpha);
  val = (vmax*histSum/numPixInTile);
  val(val>=1) = 1-eps;
  temp = -1/alpha*log(1-val);
  mapping = min(selectedRange(1)+temp*valSpread,selectedRange(2));
  
otherwise
    eid = sprintf(':%s:internalError', mfilename);
    error(eid, 'Unknown distribution type.');     % should never get here  
end

%rescale the result to be between 0 and 1 for later use by the GRAYXFORM 
%private mex function
mapping = mapping/fullRange(2);

%-----------------------------------------------------------------------------
function imgHist = clipHistogram(imgHist, clipLimit, numBins)
% This function clips the histogram according to the clipLimit and
% redistributes clipped pixels accross bins below the clipLimit

% total number of pixels overflowing clip limit in each bin
totalExcess = sum(max(imgHist - clipLimit,0));  

% clip the histogram and redistribute the excess pixels in each bin
avgBinIncr = floor(totalExcess/numBins);
upperLimit = clipLimit - avgBinIncr; % bins larger than this will be
                                     % set to clipLimit

% this loop should speed up the operation by putting multiple pixels
% into the "obvious" places first
for k=1:numBins
  if imgHist(k) > clipLimit
    imgHist(k) = clipLimit;
  else
    if imgHist(k) > upperLimit % high bin count
      totalExcess = totalExcess - (clipLimit - imgHist(k));
      imgHist(k) = clipLimit;
    else
      totalExcess = totalExcess - avgBinIncr;
      imgHist(k) = imgHist(k) + avgBinIncr;      
    end
  end
end

% this loops redistributes the remaining pixels, one pixel at a time
k = 1;
while (totalExcess ~= 0)
  %keep increasing the step as fewer and fewer pixels remain for
  %the redistribution (spread them evenly)
  stepSize = max(floor(numBins/totalExcess),1);
  for m=k:stepSize:numBins
    if imgHist(m) < clipLimit
      imgHist(m) = imgHist(m)+1;
      totalExcess = totalExcess - 1; %reduce excess
      if totalExcess == 0
        break;
      end
    end
  end
  
  k = k+1; %prevent from always placing the pixels in bin #1
  if k > numBins % start over if numBins was reached
    k = 1;
  end
end

%-----------------------------------------------------------------------------
function [c,map] = cmunique(varargin)
%CMUNIQUE Find unique colormap colors and corresponding image.

%   Copyright 1993-2003 The MathWorks, Inc.  
%   $Revision: 5.19 $  $Date: 2003/01/17 16:27:24 $

%   I/O Spec
%   ========
%   IN
%      X      - image of class uint8, uint16, or double
%      MAP    - M-by-3 array of doubles (colormap)
%   OUT
%      Y      - uint8 if NEWMAP has <= 256 entries, double 
%               if NEWMAP has > 256 entries.
%      NEWMAP - M-by-3 array of doubles (colormap)

checknargin(1,3,nargin,mfilename);
checkinput(varargin{1},{'double' 'uint8' 'uint16'},{'real' 'nonsparse'}, mfilename,'X',1);

% Convert all possible input arguments to an indexed image.
if nargin==1	% cmunique(I) or cmunique(RGB)
    arg1 = varargin{1};
    if ndims(arg1)==3	% cmunique(RGB)
        [c,map] = rgbToInd(arg1);
    else % cmunique(I)
        [c,map] = grayToInd(arg1);
    end
elseif nargin==2	% cmunique(a,cm)
    c = varargin{1}; map = varargin{2};
end

if ~isa(c, 'double')    % The promotion is necessary for the indexing into
    c = im2double(c, 'indexed');  % pos below --  ...loc(pos(c))...
end

tol = 1/1024;

map = round(map/tol)*tol;       % Quantize colormap entries to help matching below.

% Remove matching entries from colormap
 
% Sort colormap entries
[dum,ndx1] = sort(map(:,1));
[dum,ndx2] = sort(map(ndx1,2));
[dum,ndx3] = sort(map(ndx1(ndx2),3));
                % ndx maps from sorted cm to original cm
ndx = ndx1(ndx2(ndx3));
                % pos maps from original cm to sorted cm
pos = zeros(size(ndx)); pos(ndx) = 1:length(ndx);

% Find matching entries.   d indicates the location of matching entries
d = all(abs(diff(map(ndx,:)))'<tol)';

% Mapping from full cm to compressed cm. loc maps from sorted cm to compressed cm
loc = (1:length(ndx))' - [0;cumsum(d)]; 
c(:) = loc(pos(c));

% Remove matching entries (compress cm)
ndx(d) = [];
map = map(ndx,:);

% Remove colormap entries that are not used in c
[n,x] = imhist_j(c,map);
d = (n==0);             % Find unused colormap entries
loc = [1:length(map)]' - cumsum(d);

c(:) = loc(c);          % Update image values

% Remove unused entries (compress cm)
map(d,:) = [];

if max(c(:))<=256    % Output a uint8 array if we can
    c = uint8(c-1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,map] = rgbToInd(rgb)
% Convert rgb image to indexed by stuffing all pixel colors into a big colormap.
m = size(rgb,1);
n = size(rgb,2);
map = im2double(reshape(rgb,m*n,3));
x = reshape(1:m*n, m, n);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,map] = grayToInd(I)
% Convert intensity image to indexed by stuffing all pixel colors into a big colormap.
[m,n] = size(I);
map = im2double(repmat(I(:),1,3));
x = reshape(1:m*n, m, n);

%-----------------------------------------------------------------------------
function claheI = makeClaheImage(I, tileMappings, numTiles, selectedRange, numBins, dimTile)
% This function interpolates between neighboring tile mappings to produce a 
% new mapping in order to remove artificially induced tile borders.  
% Otherwise, these borders would become quite visible.  The resulting
% mapping is applied to the input image thus producing a CLAHE processed image.

%initialize the output image to zeros (preserve the class of the input image)
claheI = I;     claheI(:) = 0;

%compute the LUT for looking up original image values in the tile mappings,
%which we created earlier
if ~isa(I,'double')
  k = selectedRange(1)+1 : selectedRange(2)+1;
  aLut(k) = (k-1)-selectedRange(1);
  aLut = aLut/(selectedRange(2)-selectedRange(1));
else
  % remap from 0..1 to 0..numBins-1
  if numBins ~= 1
    binStep = 1/(numBins-1);
    start = ceil(selectedRange(1)/binStep);
    stop  = floor(selectedRange(2)/binStep);
    k = start+1:stop+1;
    aLut(k) = 0:1/(length(k)-1):1;
  else
    aLut(1) = 0; %in case someone specifies numBins = 1, which is just silly
  end
end

imgTileRow=1;
for k=1:numTiles(1)+1
  if k == 1  %special case: top row
    imgTileNumRows = dimTile(1)/2; %always divisible by 2 because of padding
    mapTileRows = [1 1];
  else 
    if k == numTiles(1)+1 %special case: bottom row      
      imgTileNumRows = dimTile(1)/2;
      mapTileRows = [numTiles(1) numTiles(1)];
    else %default values
      imgTileNumRows = dimTile(1); 
      mapTileRows = [k-1, k]; %[upperRow lowerRow]
    end
  end
  
  % loop over columns of the tileMappings cell array
  imgTileCol=1;
  for l=1:numTiles(2)+1
    if l == 1 %special case: left column
      imgTileNumCols = dimTile(2)/2;
      mapTileCols = [1, 1];
    else
      if l == numTiles(2)+1 % special case: right column
        imgTileNumCols = dimTile(2)/2;
        mapTileCols = [numTiles(2), numTiles(2)];
      else %default values
        imgTileNumCols = dimTile(2);
        mapTileCols = [l-1, l]; % right left
      end
    end
    
    % Extract four tile mappings
    ulMapTile = tileMappings{mapTileRows(1), mapTileCols(1)};
    urMapTile = tileMappings{mapTileRows(1), mapTileCols(2)};
    blMapTile = tileMappings{mapTileRows(2), mapTileCols(1)};
    brMapTile = tileMappings{mapTileRows(2), mapTileCols(2)};

    % Calculate the new greylevel assignments of pixels 
    % within a submatrix of the image specified by imgTileIdx. This 
    % is done by a bilinear interpolation between four different mappings 
    % in order to eliminate boundary artifacts.
    
    normFactor = imgTileNumRows*imgTileNumCols; %normalization factor  
    imgTileIdx = {imgTileRow:imgTileRow+imgTileNumRows-1, ...
                 imgTileCol:imgTileCol+imgTileNumCols-1};

    imgPixVals = grayxform(I(imgTileIdx{1},imgTileIdx{2}), aLut);
    
    % calculate the weights used for linear interpolation between the
    % four mappings
    rowW = repmat((0:imgTileNumRows-1)',1,imgTileNumCols);
    colW = repmat(0:imgTileNumCols-1,imgTileNumRows,1);
    rowRevW = repmat((imgTileNumRows:-1:1)',1,imgTileNumCols);
    colRevW = repmat(imgTileNumCols:-1:1,imgTileNumRows,1);
    
    claheI(imgTileIdx{1}, imgTileIdx{2}) = ...
        (rowRevW .* (colRevW .* double(grayxform(imgPixVals,ulMapTile)) + ...
                     colW    .* double(grayxform(imgPixVals,urMapTile)))+ ...
         rowW    .* (colRevW .* double(grayxform(imgPixVals,blMapTile)) + ...
                     colW    .* double(grayxform(imgPixVals,brMapTile))))...
        /normFactor;
    
    imgTileCol = imgTileCol + imgTileNumCols;    
  end %over tile cols
  imgTileRow = imgTileRow + imgTileNumRows;
end %over tile rows

%-----------------------------------------------------------------------------
function [I, selectedRange, fullRange, numTiles, dimTile, clipLimit,...
          numBins, noPadRect, distribution, alpha] = parseInputs_adapthisteq(varargin)

checknargin(1,13,nargin,mfilename);
I = varargin{1};
checkinput(I, {'uint8' 'uint16' 'double'}, {'real', '2d', 'nonsparse', 'nonempty'}, ...
           mfilename, 'I', 1);

if any(size(I) < 2)
  eid = sprintf('Images:%s:inputImageTooSmall', mfilename);
  msg = 'The input image width and height must be at least equal to 2.';
  error(eid, msg);
end

%Other options

%Set the defaults
distribution = 'uniform';
alpha   = 0.4;

if isa(I, 'double')
  fullRange = [0 1];
else
  fullRange(1) = I(1);         %copy class of the input image
  fullRange(1:2) = [-Inf Inf]; %will be clipped to min and max
  fullRange = double(fullRange);
end

selectedRange   = fullRange;

%Set the default to 256 bins regardless of the data type;
%the user can override this value at any time
numBins = 256;
normClipLimit = 0.01;
numTiles = [8 8];

checkAlpha = false;

validStrings = {'NumTiles','ClipLimit','NBins','Distribution','Alpha','Range'};

if nargin > 1
  done = false;
  idx = 2;
  while ~done
    input = varargin{idx};
    inputStr = checkstrs(input, validStrings,mfilename,'PARAM',idx);
    idx = idx+1; %advance index to point to the VAL portion of the input 
    if idx > nargin
      eid = sprintf('Images:%s:valFor%sMissing', mfilename, inputStr);
      msg = sprintf('Parameter ''%s'' must be followed by a value.', inputStr);
      error(eid, msg);        
    end
    
    switch inputStr
     case 'NumTiles'
       numTiles = varargin{idx};
       checkinput(numTiles, {'double'}, {'real', 'vector', ...
                           'integer', 'finite','positive'},...
                  mfilename, inputStr, idx);

       if (any(size(numTiles) ~= [1,2]))
         eid = sprintf('Images:%s:invalidNumTiles', mfilename);
         msg = sprintf(['Value of parameter,''%s'', must be a two '...
                        'element row vector.'], inputStr);
         error(eid, msg);
       end
       
       if any(numTiles < 2)
         eid = sprintf('Images:%s:invalidNumTiles', mfilename);
         msg = sprintf(['Both elements of ''%s'' value ', 'must be greater or equal to 2.'],inputStr);
         error(eid, msg);
       end
      
     case 'ClipLimit'
      normClipLimit = varargin{idx};
      checkinput(normClipLimit, {'double'},{'scalar','real','nonnegative'}, mfilename, inputStr, idx);
      
      if normClipLimit > 1
        eid = sprintf('Images:%s:invalidClipLimit', mfilename);
        msg = sprintf(['Value of parameter ''%s'' must be in the range ', 'from 0 to 1.'], inputStr);
        error(eid, msg);
      end
     
     case 'NBins'
      numBins = varargin{idx};      
      checkinput(numBins, {'double'}, {'scalar','real','integer', 'positive'}, mfilename, 'NBins', idx);
     
     case 'Distribution'
      validDist = {'rayleigh','exponential','uniform'};
      distribution = checkstrs(varargin{idx}, validDist, mfilename, 'Distribution', idx);
     case 'Alpha'
      alpha = varargin{idx};
      checkinput(alpha, {'double'},{'scalar','real',...
                          'nonnan','positive','finite'},mfilename, 'Alpha',idx);
      checkAlpha = true;
     case 'Range'
      validRangeStrings = {'original','full'};
      rangeStr = checkstrs(varargin{idx}, validRangeStrings,mfilename,'Range',idx);
      if strmatch(rangeStr,'original')
        selectedRange = double([min(I(:)), max(I(:))]);
      end
     otherwise
      eid = sprintf('Images:%s:internalError', mfilename);
      msg   = 'Unknown input string.'; %should never get here
      error(eid, msg);
    end
    if (idx >= nargin)
        done = true;
        break;
    end
    idx=idx+1;
  end
end

%% Pre-process the inputs
dimI = size(I);     dimTile = dimI ./ numTiles;

%check if tile size is reasonable
if any(dimTile < 1)
  eid = sprintf('Images:%s:inputImageTooSmall', mfilename);
  msg = sprintf(['The image I is too small to be split into [%s] ',...
                 'number of tiles.'], num2str(numTiles));
  error(eid, msg);
end

if checkAlpha
  if strcmp(distribution,'uniform')
    eid = sprintf('Images:%s:alphaShouldNotBeSpecified', mfilename);
    msg = sprintf(['Parameter ''Alpha'' cannot be specified for',...
                   ' ''%s'' distribution.'], distribution);
    error(eid, msg);
  end
end

%check if the image needs to be padded; pad if necessary;
%padding occurs if any dimension of a single tile is an odd number
%and/or when image dimensions are not divisible by the selected 
%number of tiles
rowDiv  = mod(dimI(1),numTiles(1)) == 0;
colDiv  = mod(dimI(2),numTiles(2)) == 0;

if rowDiv && colDiv
  rowEven = mod(dimTile(1),2) == 0;
  colEven = mod(dimTile(2),2) == 0;  
end

noPadRect = [];
if  ~(rowDiv && colDiv && rowEven && colEven)
  padRow = 0;
  padCol = 0;
  
  if ~rowDiv
    rowTileDim = floor(dimI(1)/numTiles(1)) + 1;
    padRow = rowTileDim*numTiles(1) - dimI(1);
  else
    rowTileDim = dimI(1)/numTiles(1);
  end
  
  if ~colDiv
    colTileDim = floor(dimI(2)/numTiles(2)) + 1;
    padCol = colTileDim*numTiles(2) - dimI(2);
  else
    colTileDim = dimI(2)/numTiles(2);
  end
  
  %check if tile dimensions are even numbers
  rowEven = mod(rowTileDim,2) == 0;
  colEven = mod(colTileDim,2) == 0;
  
  if (~rowEven),    padRow = padRow+numTiles(1);  end
  if (~colEven),    padCol = padCol+numTiles(2);  end
  
  padRowPre  = floor(padRow/2);  padRowPost = ceil(padRow/2);
  padColPre  = floor(padCol/2);  padColPost = ceil(padCol/2);
  
  I = padarray(I,[padRowPre  padColPre ],'symmetric','pre');
  I = padarray(I,[padRowPost padColPost],'symmetric','post');

  %UL corner (Row, Col), LR corner (Row, Col)
  noPadRect.ulRow = padRowPre+1;
  noPadRect.ulCol = padColPre+1;
  noPadRect.lrRow = padRowPre+dimI(1);
  noPadRect.lrCol = padColPre+dimI(2);
end

%redefine this variable to include the padding
dimI = size(I);

%size of the single tile
dimTile = dimI ./ numTiles;

%compute actual clip limit from the normalized value entered by the user
%maximum value of normClipLimit=1 results in standard AHE, i.e. no clipping;
%the minimum value minClipLimit would uniformly distribute the image pixels
%across the entire histogram, which would result in the lowest possible
%contrast value
numPixInTile = prod(dimTile);
minClipLimit = ceil(numPixInTile/numBins);
clipLimit = minClipLimit + round(normClipLimit*(numPixInTile-minClipLimit));

%------------------------------------------------------------------------------
function bw2 = bwareaopen(varargin)
%BWAREAOPEN Binary area open; remove small objects.
%   BW2 = BWAREAOPEN(BW,P) removes from a binary image all connected
%   components (objects) that have fewer than P pixels, producing another
%   binary image BW2.  The default connectivity is 8 for two dimensions,
%   26 for three dimensions, and CONNDEF(NDIMS(BW),'maximal') for higer
%   dimensions. 
%
%   BW2 = BWAREAOPEN(BW,P,CONN) specifies the desired connectivity.  CONN
%   may have the following scalar values:  
%
%       4     two-dimensional four-connected neighborhood
%       8     two-dimensional eight-connected neighborhood
%       6     three-dimensional six-connected neighborhood
%       18    three-dimensional 18-connected neighborhood
%       26    three-dimensional 26-connected neighborhood
%
%   Connectivity may be defined in a more general way for any dimension by
%   using for CONN a 3-by-3-by- ... -by-3 matrix of 0s and 1s.  The 1-valued
%   elements define neighborhood locations relative to the center element of
%   CONN.  CONN must be symmetric about its center element.
%
%   Class Support
%   -------------
%   BW can be a logical or numeric array of any dimension, 
%   and it must be nonsparse.
%
%   BW2 is logical.
%
%   Example
%   -------
%   Remove all objects in the image text.png containing fewer than 50
%   pixels.
%
%       bwOriginal = imread('text.png');
%       imview(bwOriginal)
%       bwAreaOpen_50pixels = bwareaopen(bwOriginal,50);
%       imview(bwAreaOpen_50pixels)
%
%   See also BWLABEL, BWLABELN, CONNDEF, REGIONPROPS.

%   Copyright 1993-2003 The MathWorks, Inc.
%   $Revision: 1.12 $  $Date: 2003/01/27 20:15:56 $

% Input/output specs
% ------------------
% BW:    N-D real full matrix
%        any numeric class
%        sparse not allowed
%        anything that's not logical is converted first using
%          bw = BW ~= 0
%        Empty ok
%        Inf's ok, treated as 1
%        NaN's ok, treated as 1
%
% P:     double scalar
%        nonnegative integer
%
% CONN:  connectivity
%
% BW2:   logical, same size as BW
%        contains only 0s and 1s.

[bw,p,conn] = parse_inputs_bwareaopen(varargin{:});

L = bwlabeln(bw,conn);
s = regionprops(L,'area');
area = [s.Area];
idx = find(area >= p);
bw2 = ismember(L,idx);

% ----------------------- parse_inputs ---
function [bw,p,conn] = parse_inputs_bwareaopen(varargin)

checknargin(2,3,nargin,mfilename)

bw = varargin{1};
checkinput(bw,{'numeric' 'logical'},{'nonsparse'},mfilename,'BW',1);
if (~islogical(bw)),	bw = bw ~= 0;	end

p = varargin{2};
checkinput(p,{'double'},{'scalar' 'integer' 'nonnegative'},mfilename,'P',2);

if (nargin >= 3),	conn = varargin{3};
else				conn = conndef(ndims(bw),'maximal');
end
checkconn(conn,mfilename,'CONN',3)

%------------------------------------------------------------------------------
function [B,L,N,A] = bwboundaries(varargin)
%BWBOUNDARIES Trace region boundaries in a binary image.
%   Copyright 1993-2003 The MathWorks, Inc.
%   $Revision: 1.2 $  $Date: 2003/03/05 22:29:20 $

[BW, conn, findholes] = parseInputs_bwboundaries(varargin{:});
[objs , L] = FindObjectBoundaries(BW, conn);
if (nargout == 1),	clear L;	end			% Bloody array of doubles to store short ints. Fck

if (findholes)
	if (~isa(BW,'logical'))
		warndlg('bwboundaries with holes works only with binary masks.','Warning');    return
	end
	[holes, LabeledHoles] = FindHoleBoundaries(BW, conn);
	% Generate combined holes+objects label matrix
	if (nargout > 1)
		L = L + (LabeledHoles~=0)*length(objs) + LabeledHoles;
	end
else
	holes = {};
end

B = [objs; holes];      % Create the output matrix
N = length(objs);       % Return number of object boundaries
if(nargout > 3)         % Produce an adjacency matrix showing parent-hole-child relationships
	A = CreateAdjMatrix(B, N);
end

%-----------------------------------------------------------------------------
function [BW, conn, findholes] = parseInputs_bwboundaries(varargin)
checknargin(1,4,nargin,mfilename);

BW = varargin{1};
checkinput(BW, {'numeric','logical'}, {'real','2d','nonsparse'}, mfilename, 'BW', 1);
if (~islogical(BW))     BW = BW ~= 0;   end
firstStringToProcess = 0;

if nargin < 2
  conn = 8;
else
  if ischar(varargin{2})
    firstStringToProcess = 2;
    conn = 8;
  else
    if (nargin > 2)     firstStringToProcess = 3;   end
    conn = varargin{2};
    checkinput(conn, {'double'}, {}, mfilename, 'CONN', 2);
    if (conn~=4 && conn~=8)
      eid = sprintf('Images:%s:badScalarConn', mfilename);
      msg = 'A scalar connectivity specifier CONN must either be 4 or 8';
      error(eid, msg);
    end
  end
end

findholes = true;
if firstStringToProcess
  validStrings = {'noholes', 'holes'};
  for k = firstStringToProcess:nargin
    % check for options
    string = checkstrs(varargin{k}, validStrings, mfilename, 'OPTION', k);
    switch string
     case 'noholes',    findholes = false;
     case 'holes',      findholes = true;
     otherwise
      error('Images:bwboundaries:unexpectedError', '%s', 'Unexpected logic error.')
    end      
  end
end

%-----------------------------------------------------------------------------
function varargout = bwdist(varargin)
%BWDIST Distance transform.
%   D = BWDIST(BW) computes the Euclidean distance transform of the
%   binary image BW. For each pixel in BW, the distance transform assigns
%   a number that is the distance between that pixel and the nearest
%   nonzero pixel of BW. BWDIST uses the Euclidean distance metric by
%   default.  BW can have any dimension.  D is the same size as BW.
%
%   [D,L] = BWDIST(BW) also computes the nearest-neighbor transform and
%   returns it as a label matrix, L.  L has the same size as BW and D.
%   Each element of L contains the linear index of the nearest nonzero
%   pixel of BW.
%
%   [D,L] = BWDIST(BW,METHOD) lets you compute an alternate distance
%   transform, depending on the value of METHOD.  METHOD can be
%   'cityblock', 'chessboard', 'quasi-euclidean', or 'euclidean'.  METHOD
%   defaults to 'euclidean' if not specified.  METHOD may be
%   abbreviated.
%
%   The different methods correspond to different distance metrics.  In
%   2-D, the cityblock distance between (x1,y1) and (x2,y2) is abs(x1-x2)
%   + abs(y1-y2).  The chessboard distance is max(abs(x1-x2),
%   abs(y1-y2)).  The quasi-Euclidean distance is:
%
%       abs(x1-x2) + (sqrt(2)-1)*abs(y1-y2),  if abs(x1-x2) > abs(y1-y2)
%       (sqrt(2)-1)*abs(x1-x2) + abs(y1-y2),  otherwise
%
%   The Euclidean distance is sqrt((x1-x2)^2 + (y1-y2)^2).
%
%   Note
%   ----
%   BWDIST uses fast algorithms to compute the true Euclidean distance
%   transform, especially in the 2-D case. The other methods are
%   provided primarily for pedagogical reasons. However, the alternative
%   distance transforms are sometimes significantly faster for
%   multidimensional input images, particularly those that have many
%   nonzero elements.
%
%   Class support
%   -------------
%   BW can be numeric or logical, and it must be nonsparse. D and L are
%   single and uint32 matrices with the same size as BW.
%
%   Examples
%   --------
%   Here is a simple example of the Euclidean distance transform:
%
%       bw = zeros(5,5); bw(2,2) = 1; bw(4,4) = 1;
%       [D,L] = bwdist(bw)

%   Copyright 1993-2004 The MathWorks, Inc.
%   $Revision: 1.7.4.4 $  $Date: 2004/08/10 01:38:38 $

	[BW,method] = parse_inputs_bwdist(varargin{:});

	% Computing the nearest-neighbor transform is expensive in memory, so we
	% only want to call the lower-level functions eucdist2, eucdistn, and ddist
	% with two output arguments if we have been called with two output arguments.
	if nargout <= 1
		varargout = cell(1,1);
	else
		varargout = cell(1,2);
	end

	iptcheckinput(BW,{'logical'},{'nonsparse'},mfilename,'BW',1);

	if strcmp(method,'euclidean')
		% Use a really fast method for 2-D Euclidean distance transforms, or
		% a reasonably fast kd-tree based method for multidimensional
		% Euclidean distance transforms.
		if ndims(BW) == 2
			[varargout{:}] = eucdist2(BW);
		else
			%[varargout{:}] = eucdistn(BW);
			error('Over 2D arrays is not implemented here.')
		end
	else
		% All methods other that Euclidean use the same algorithm,
		% implemented in private/ddist.  The only difference is in the
		% connectivity and weights used.
		switch method
		  case 'cityblock'
			conn = conndef(ndims(BW),'minimal');
			weights = ones(size(conn));

		  case 'chessboard'
			conn = conndef(ndims(BW),'maximal');
			weights = ones(size(conn));

		  case 'quasi-euclidean'
			conn = conndef(ndims(BW),'maximal');

			% For quasi-Euclidean, form a weights array whose values are the
			% distances between the corresponding elements and the center element.
			kk = cell(1,ndims(BW));
			[kk{:}] = ndgrid(-1:1);
			weights = zeros(size(conn));
			for p = 1:ndims(BW)
				% Although the Euclidean distance formula certainly involves
				% squaring each term, all terms here are either 0, 1, or -1,
				% so that's why the abs() term isn't squared below.
				weights = weights + abs(kk{p});
			end
			weights = sqrt(weights);

		  otherwise
			eid = 'Images:bwdist:unrecognizedMethodString';
			error(eid, 'Internal problem - unrecognized method string: %s', method);
		end

		% Postprocess the weights to make sure the center weight is
		% zero, and to keep only values corresponding to nonzero
		% connectivity values.
		weights = weights(:);
		weights((end+1)/2) = 0.0;
		weights(~conn) = [];

		% Call the dual-scan neighborhood based algorithm.
		[varargout{:}] = ddist(BW,conn,weights);
	end

%--------------------------------------------------
function [BW,method] = parse_inputs_bwdist(varargin)
	iptchecknargin(1,2,nargin,mfilename);
	iptcheckinput(varargin{1}, {'logical','numeric'}, {'nonsparse', 'real'}, mfilename, 'BW', 1);
	BW = varargin{1} ~= 0;

	if nargin < 2
		method = 'euclidean';
	else
		valid_methods = {'euclidean','cityblock','chessboard','quasi-euclidean'};
		method = checkstrs(varargin{2}, valid_methods, mfilename, 'METHOD', 2);
	end

%-----------------------------------------------------------------------------
function [B, L]= FindHoleBoundaries(BW, conn)
% Avoid topological errors.  If objects are 8 connected, then holes
% must be 4 connected and vice versa.
	if (conn == 4),		backgroundConn = 8;
	else				backgroundConn = 4;
	end

	% Turn holes into objects
	%BWcomplement = imcomplement(BW);
	BWcomplement = ~BW;			% Logical only

	% clear unwanted "hole" objects from the border
	% BWholes = imclearborder(BWcomplement, backgroundConn);	% NOT PORTED

	% get the holes!
	L = bwlabel(BWcomplement, backgroundConn);
	B = bwboundariesmex(L, backgroundConn);

% ----------------------------------------------------------------------------
function e = bweuler(a,n)
%BWEULER Compute the Euler number of binary image.
%   EUL = BWEULER(BW,N) returns the Euler number for the binary
%   image BW. EUL is a scalar whose value is the number of
%   objects in the image minus the total number of holes in those
%   objects.  N can have a value of either 4 or 8, where 4
%   specifies 4-connected objects and 8 specifies 8-connected
%   objects; if the argument is omitted, it defaults to 8. 
%
%   Class Support
%   -------------
%   BW can be numeric or logical and it must be real, nonsparse 
%   and two-dimensional.
%   EUL is of class double.
%
%   Example
%   -------
%       BW = imread('circles.png');
%       imview(BW)
%       bweuler(BW)
%
%   See also BWPERIM, BWMORPH.

%   Copyright 1993-2003 The MathWorks, Inc.  
%   $Revision: 5.25 $  $Date: 2003/01/27 20:15:56 $

% Reference: William Pratt, Digital Image Processing, John Wiley
% and Sons, 1991, pp. 630-634.
  
checknargin(1,2,nargin,'bweuler');

checkinput(a,{'numeric' 'logical'},{'nonsparse' 'real' '2d'}, 'bweuler', 'BW', 1);

if nargin < 2 
    n = 8; 
else
    checkinput(n,{'double'},{'scalar' 'real' 'integer'}, 'bweuler', 'N', 2);
end

if n~=8 && n~=4
    eid = 'Images:bweuler:invalidN';
    error(eid,'%s','N must be either 4 or 8.');
end

if n==4,
    lut = 4*[0 0.25 0.25 0 0.25 0  .5 -0.25 0.25  0.5  0 -0.25 0 -0.25 -0.25 0] + 2;
else
    lut = 4*[0 0.25 0.25 0 0.25 0 -.5 -0.25 0.25 -0.5  0 -0.25 0 -0.25 -0.25 0] + 2;
end

% Need to zero-pad the input
b = padarray(a,[1 1],'both');

weights = applylut(b,lut);
e = (sum(weights(:)) - 2*numel(b)) / 4;


%-----------------------------------------------------------------------------
function [A] = CreateAdjMatrix(B, numObjs)
A = sparse(false(length(B)));
levelCellArray = GroupBoundariesByTreeLevel(B, numObjs);

% scan through all the level pairs
for k = 1:length(levelCellArray)-1
  parentsIdx = levelCellArray{k};     % outside boundaries
  childrenIdx = levelCellArray{k+1};  % inside boundaries
  parents  = B(parentsIdx);
  children = B(childrenIdx);
  sampChildren = GetSamplePointsFromBoundaries(children);
  for m=1:length(parents)
    parent = parents{m};
    inside = inpolygon(sampChildren(:,2), sampChildren(:,1), parent(:,2), parent(:,1));
    % casting to logical is necessary because of the bug, see GECK #137394
    inside = logical(inside);
    A(childrenIdx(inside), parentsIdx(m)) = true;
  end
end

%-----------------------------------------------------------------------------
function points = GetSamplePointsFromBoundaries(B)
points = zeros(length(B),2);
for m = 1:length(B)
  boundary = B{m};
  points(m,:) = boundary(1,:);
end

%-----------------------------------------------------------------------------
function idxGroupedByLevel = GroupBoundariesByTreeLevel(B, numObjs)
% Produces a cell array of indices into the boundaries cell array B.  The
% first element of the output cell array holds a double array of indices
% of boundaries which are the outermost (first layer of an onion), the
% second holds the second layer, and so on.

processHoles = ~(length(B) == numObjs);
% parse the input
objIdx  = 1:numObjs;
objs  = B(objIdx);

if processHoles
  holeIdx = numObjs+1:length(B);
  holes = B(holeIdx);
else
  holes = {};
end

% initialize output and loop control variables
idxGroupedByLevel = {};
done     = false;
findHole = false; % start with an object boundary

while ~done
  if (findHole)
    I = FindOutermostBoundaries(holes);
    holes = holes(~I); % remove processed boundaries
    idxGroupedByLevel = [ idxGroupedByLevel, {holeIdx(I)} ];
    holeIdx = holeIdx(~I);   % remove indices of processed boundaries
  else
    I = FindOutermostBoundaries(objs);
    objs = objs(~I);
    idxGroupedByLevel = [ idxGroupedByLevel, {objIdx(I)} ];
    objIdx = objIdx(~I);
  end
  if(processHoles),  findHole = ~findHole;   end
  if (isempty(holes) && isempty(objs)),  done = true;    end
end  
  
%-----------------------------------------------------------------------------
function I = FindOutermostBoundaries(B)
% Returns a logical vector showing the locations of outermost boundaries 
% in the input vector (ie 1 for the boundaries that are outermost and
% 0 for all other boundaries)

% Look for parent boundaries
I = false(1,length(B));
for m = 1:length(B)
  boundary = B{m};
  x = boundary(1,2); % grab a sample point for testing
  y = boundary(1,1);
  surrounded = false;
  for n = [1:(m-1), (m+1):length(B)]	% exclude boundary under test
    boundary = B{n};
    if( inpolygon(x, y, boundary(:,2), boundary(:,1)) )
        surrounded=true;
        break;
    end
  end
  I(m) = ~surrounded;
end

%-----------------------------------------------------------------------------
function [B, L] = FindObjectBoundaries(BW, conn)
L = bwlabel(BW, conn);
B = bwboundariesmex(L, conn);

% ----------------------------------------------------------------------------
function [L,numComponents] = bwlabel(BW,mode)
%BWLABEL Label connected components in binary image.
%   Copyright 1993-2003 The MathWorks, Inc.  
%   $Revision: 1.31 $  $Date: 2003/03/05 22:29:18 $

checknargin(1,2,nargin,mfilename);
checkinput(BW, {'logical' 'numeric'}, {'real', '2d', 'nonsparse'}, mfilename, 'BW', 1);

if (nargin < 2)
    mode = 8;
else
    checkinput(mode, 'double', 'scalar', mfilename, 'N', 2);
end
if (~islogical(BW)),    BW = BW ~= 0;    end

[M,N] = size(BW);
% Compute run-length encoding and assign initial labels.
[sr,er,sc,labels,i,j] = bwlabel1(BW,mode);
if (isempty(labels))
    numLabels = 0;
else
    numLabels = max(labels);
end

% Create a sparse matrix representing the equivalence graph.
tmp = (1:numLabels)';
A = sparse([i;j;tmp], [j;i;tmp], 1, numLabels, numLabels);

% Determine the connected components of the equivalence graph
% and compute a new label vector.

% Find the strongly connected components of the adjacency graph
% of A.  dmperm finds row and column permutations that transform
% A into upper block triangular form.  Each block corresponds to
% a connected component; the original source rows in each block
% correspond to the members of the corresponding connected
% component.  The first two output% arguments (row and column
% permutations, respectively) are the same in this case because A
% is symmetric.  The vector r contains the locations of the
% blocks; the k-th block as indices r(k):r(k+1)-1.
[p,p,r] = dmperm(A);

% Compute vector containing the number of elements in each component.
sizes = diff(r);
numComponents = length(sizes);  % Number of components.

blocks = zeros(1,numLabels);
blocks(r(1:numComponents)) = 1;
blocks = cumsum(blocks);
blocks(p) = blocks;
labels = blocks(labels);

% Given label information, create output matrix.
L = bwlabel2(sr, er, sc, labels, M, N);

% ----------------------------------------------------------------------------
function [L,num] = bwlabeln(varargin)
%BWLABELN Label connected components in N-D binary image.
%   L = BWLABELN(BW) returns a label matrix, L, containing labels for the
%   connected components in BW.  BW can have any dimension; L is the same
%   size as BW.  The elements of L are integer values greater than or equal
%   to 0.  The pixels labeled 0 are the background.  The pixels labeled 1
%   make up one object, the pixels labeled 2 make up a second object, and so
%   on.  The default connectivity is 8 for two dimensions, 26 for three
%   dimensions, and CONNDEF(NDIMS(BW),'maximal') for higher dimensions.
%
%   [L,NUM] = BWLABELN(BW) returns the number of connected objects found
%   in BW.
%
%   [L,NUM] = BWLABELN(BW,CONN) specifies the desired connectivity.  CONN
%   may have the following scalar values:
%
%       4     two-dimensional four-connected neighborhood
%       8     two-dimensional eight-connected neighborhood
%       6     three-dimensional six-connected neighborhood
%       18    three-dimensional 18-connected neighborhood
%       26    three-dimensional 26-connected neighborhood
%
%   Connectivity may be defined in a more general way for any dimension by
%   using for CONN a 3-by-3-by- ... -by-3 matrix of 0s and 1s.  The 1-valued
%   elements define neighborhood locations relative to the center element of
%   CONN.  CONN must be symmetric about its center element.
%
%   Note: Comparing BWLABEL and BWLABELN
%   ------------------------------------
%   BWLABEL supports 2-D inputs only, whereas BWLABELN support any 
%   input dimension.  In some cases you might prefer to use BWLABELN even
%   for 2-D problems because it can be faster.  If you have a 2-D input
%   whose objects are relatively "thick" in the vertical direction,
%   BWLABEL will probably be faster; otherwise BWLABELN will probably be faster.
%
%   Class Support
%   -------------
%   BW can be numeric or logical, and it must be real and nonsparse. L is double.
%
%   Example
%   -------
%       BW = cat(3,[1 1 0; 0 0 0; 1 0 0],...
%                  [0 1 0; 0 0 0; 0 1 0],...
%                  [0 1 1; 0 0 0; 0 0 1])
%       bwlabeln(BW)

%   Copyright 1993-2003 The MathWorks, Inc.  
%   $Revision: 1.9 $  $Date: 2003/01/17 16:27:22 $

[A,conn] = parse_inputs_bwlabeln(varargin{:});

[L,num] = bwlabelnmex(A,conn);

% ----------------- parse_inputs ---
function [A,conn] = parse_inputs_bwlabeln(varargin)

checknargin(1,2,nargin,mfilename);

checkinput(varargin{1}, {'numeric', 'logical'}, {'real' 'nonsparse'}, mfilename, 'BW', 1);

A = varargin{1};
if (~islogical(A)),		A = A ~= 0;		end

if nargin < 2
    conn = conndef(ndims(A), 'maximal');
else
    conn = varargin{2};
    checkconn(conn,mfilename,'CONN',2);
end

%-----------------------------------------------------------------------------
function pixelx = axes2pix(dim, x, axesx)
%   Copyright 1993-2003 The MathWorks, Inc.  
%   $Revision: 5.13 $  $Date: 2003/01/17 16:28:14 $

error_str = nargchk(3, 3, nargin);
if (~ isempty(error_str));  error('There must be 3 input arguments.');  end
if (max(size(dim)) ~= 1)    error('First argument must be a scalar.');  end
if (min(size(x)) > 1)   error('X must be a vector.');   end
xfirst = x(1);      xlast = x(max(size(x)));

if (dim == 1)   pixelx = axesx - xfirst + 1;    return;     end
xslope = (dim - 1) / (xlast - xfirst);
if ((xslope == 1) && (xfirst == 1))
  pixelx = axesx;
else
  pixelx = xslope * (axesx - xfirst) + 1;
end

%-----------------------------------------------------------------------------
function h = fsamp2(f1,f2,hd,siz)
%FSAMP2 Design 2-D FIR filter using frequency sampling.
%   Copyright 1993-2003 The MathWorks, Inc.  
%   $Revision: 5.21 $  $Date: 2003/01/17 16:27:29 $

checknargin(1,4,nargin,mfilename);
if nargin==1, % Uniform spacing case (fast)
  hd = f1;
  hd = rot90(fftshift(rot90(hd,2)),2); % Inverse fftshift
  h = fftshift(ifft2(hd));
elseif nargin==2 || nargin==3,
  msg = 'Wrong number of input arguments.';
  eid = sprintf('Images:%s:expectedOneOrFourInputs',mfilename);
  error(eid, msg);
else, % Create filter of size SIZ to solve problem at the points (f1,f2,hd)
  if ~isa(hd, 'double')
    hd = double(hd);
  end  
  % Expand f1 and f2 if they are vectors.
  if min(size(f1))==1 && min(size(f2))==1 && any(size(hd)~=size(f1)),
    [f1,f2] = meshgrid(f1,f2);
  end
  if numel(hd)<prod(siz),
      msg = 'Not enough desired frequency points. Results may be inaccurate.';
      wid = sprintf('Images:%s:notEnoughFreqPoints',mfilename);
      warning(wid,msg);
  end
  % Convert frequency to radians.
  f1 = f1*pi; f2 = f2*pi;
  h = zeros(siz);
  [n1,n2] = meshgrid((0:siz(2)-1)-floor(siz(2)/2),(0:siz(1)-1)-floor(siz(1)/2));
  DFT = exp(-sqrt(-1)*f1(:)*n1(:)').*exp(-sqrt(-1)*f2(:)*n2(:)');
  h(:) = DFT\hd(:);
end  

% Convert to real if possible.
if all(max(abs(imag(h)))<sqrt(eps)), h = real(h); end
h = rot90(h,2); % Rotate for use with filter2

%----------------------------------------------------------------------------------
function h = fspecial(varargin)
%FSPECIAL Create 2-D special filters.
%   Copyright 1993-2003 The MathWorks, Inc.  
%   $Revision: 5.29 $  $Date: 2003/01/17 16:27:30 $

[type, p2, p3] = ParseInputs_fspecial(varargin{:});
switch type
  case 'average' % Smoothing filter
     siz = p2;
     h   = ones(siz)/prod(siz);
  case 'disk' % Disk filter
     rad   = p2;
     crad  = ceil(rad-0.5);
     [x,y] = meshgrid(-crad:crad,-crad:crad);
     maxxy = max(abs(x),abs(y));
     minxy = min(abs(x),abs(y));
     m1 = (rad^2 <  (maxxy+0.5).^2 + (minxy-0.5).^2).*(minxy-0.5) + ...
          (rad^2 >= (maxxy+0.5).^2 + (minxy-0.5).^2).* sqrt(rad^2 - (maxxy + 0.5).^2);
     m2 = (rad^2 >  (maxxy-0.5).^2 + (minxy+0.5).^2).*(minxy+0.5) + ...
          (rad^2 <= (maxxy-0.5).^2 + (minxy+0.5).^2).* sqrt(rad^2 - (maxxy - 0.5).^2);
     sgrid = (rad^2*(0.5*(asin(m2/rad) - asin(m1/rad)) + ...
             0.25*(sin(2*asin(m2/rad)) - sin(2*asin(m1/rad)))) - ...
             (maxxy-0.5).*(m2-m1) + (m1-minxy+0.5)) ... 
	          .*((((rad^2 < (maxxy+0.5).^2 + (minxy+0.5).^2) & ...
             (rad^2 > (maxxy-0.5).^2 + (minxy-0.5).^2)) | ...
	          ((minxy==0)&(maxxy-0.5 < rad)&(maxxy+0.5>=rad))));
     sgrid = sgrid + ((maxxy+0.5).^2 + (minxy+0.5).^2 < rad^2);
     sgrid(crad+1,crad+1) = min(pi*rad^2,pi/2);
     if ((crad>0) & (rad > crad-0.5) & (rad^2 < (crad-0.5)^2+0.25)) 
        m1  = sqrt(rad^2 - (crad - 0.5).^2);
	    m1n = m1/rad;
        sg0 = 2*(rad^2*(0.5*asin(m1n) + 0.25*sin(2*asin(m1n)))-m1*(crad-0.5));
        sgrid(2*crad+1,crad+1) = sg0;
        sgrid(crad+1,2*crad+1) = sg0;
        sgrid(crad+1,1)        = sg0;
        sgrid(1,crad+1)        = sg0;
        sgrid(2*crad,crad+1)   = sgrid(2*crad,crad+1) - sg0;
        sgrid(crad+1,2*crad)   = sgrid(crad+1,2*crad) - sg0;
        sgrid(crad+1,2)        = sgrid(crad+1,2)      - sg0;
        sgrid(2,crad+1)        = sgrid(2,crad+1)      - sg0;
     end
     sgrid(crad+1,crad+1) = min(sgrid(crad+1,crad+1),1);
     h = sgrid/sum(sgrid(:));
  case 'gaussian' % Gaussian filter
     siz   = (p2-1)/2;     std   = p3;
     [x,y] = meshgrid(-siz(2):siz(2),-siz(1):siz(1));
     arg   = -(x.*x + y.*y)/(2*std*std);
     h     = exp(arg);
     h(h<eps*max(h(:))) = 0;
     sumh = sum(h(:));
     if (sumh ~= 0)       h  = h/sumh;     end
  case 'laplacian' % Laplacian filter
     alpha = p2;
     alpha = max(0,min(alpha,1));
     h1    = alpha/(alpha+1); h2 = (1-alpha)/(alpha+1);
     h     = [h1 h2 h1;h2 -4/(alpha+1) h2;h1 h2 h1];
  case 'log' % Laplacian of Gaussian
     % first calculate Gaussian
     siz   = (p2-1)/2;
     std2   = p3^2;
     [x,y] = meshgrid(-siz(2):siz(2),-siz(1):siz(1));
     arg   = -(x.*x + y.*y)/(2*std2);
     h     = exp(arg);
     h(h<eps*max(h(:))) = 0;
     sumh = sum(h(:));
     if (sumh ~= 0)     h  = h/sumh;     end
     % now calculate Laplacian     
     h1 = h.*(x.*x + y.*y - 2*std2)/(std2^2);
     h     = h1 - sum(h1(:))/prod(p2); % make the filter sum to zero
  case 'motion' % Motion filter uses bilinear interpolation
     len = max(1,p2);
     half = (len-1)/2;% rotate half length around center
     phi = mod(p3,180)/180*pi;

     cosphi = cos(phi);
     sinphi = sin(phi);
     xsign = sign(cosphi);
     linewdt = 1;

     % define mesh for the half matrix, eps takes care of the right size
     % for 0 & 90 rotation
     sx = fix(half*cosphi + linewdt*xsign - len*eps);
     sy = fix(half*sinphi + linewdt - len*eps);
     [x, y] = meshgrid([0:xsign:sx],[0:sy]);

     % define shortest distance from a pixel to the rotated line 
     dist2line = (y*cosphi-x*sinphi);% distance perpendicular to the line

     rad = sqrt(x.^2 + y.^2);
     % find points beyond the line's end-point but within the line width
     lastpix = find((rad >= half)&(abs(dist2line)<=linewdt));
     %distance to the line's end-point parallel to the line 
     x2lastpix = half - abs((x(lastpix) + dist2line(lastpix)*sinphi)/cosphi);

     dist2line(lastpix) = sqrt(dist2line(lastpix).^2 + x2lastpix.^2);
     dist2line = linewdt + eps - abs(dist2line);
     dist2line(dist2line<0) = 0;% zero out anything beyond line width

     % unfold half-matrix to the full size
     h = rot90(dist2line,2);
     h(end+[1:end]-1,end+[1:end]-1) = dist2line;
     h = h./(sum(h(:)) + eps*len*len);
     if cosphi>0
       h = flipud(h);
     end
  case 'prewitt' % Prewitt filter
     h = [1 1 1;0 0 0;-1 -1 -1];
  case 'sobel' % Sobel filter
     h = [1 2 1;0 0 0;-1 -2 -1];
  case 'unsharp' % Unsharp filter
     alpha = p2;
     h     = [0 0 0;0 1 0;0 0 0] - fspecial('laplacian',alpha);
end

%---------------------------------------------------------------
function [type, p2, p3] = ParseInputs_fspecial(varargin)
p2 = [];    p3 = [];

% Check the number of input arguments.
checknargin(1,3,nargin,mfilename);
% Determine filter type from the user supplied string.
type = varargin{1};
type = checkstrs(type,{'gaussian','sobel','prewitt','laplacian','log',...
                       'average','unsharp','disk','motion'},mfilename,'TYPE',1);
  
% default values
switch type
	case 'average'
      p2 = [3 3];  % siz
   case 'disk'
      p2 = 5;      % rad
   case 'gaussian'
      p2 = [3 3];  % siz
      p3 = 0.5;    % std
   case {'laplacian', 'unsharp'}
      p2 = 1/5;    % alpha
   case 'log'
      p2 = [5 5];  % siz
      p3 = 0.5;    % std
   case 'motion'
      p2 = 9;     % len
      p3 = 0;      % theta
end

switch nargin
    case 1
        % Nothing to do here; the default values have already been assigned.        
    case 2
       p2 = varargin{2}; 
       switch type
          case {'sobel','prewitt'}
              msg = sprintf('%s: Too many arguments for this type of filter.', upper(mfilename));
              eid = sprintf('Images:%s:tooManyArgsForThisFilter', mfilename);
              error(eid,msg);
          case {'laplacian','unsharp'}
              checkinput(p2,{'double'},{'nonnegative','real','nonempty','finite','scalar'},...
                         mfilename,'ALPHA',2);
              if  p2 > 1
                  msg = sprintf('%s: ALPHA should be less than or equal 1 and greater than 0.', upper(mfilename));
                  eid = sprintf('Images:%s:outOfRangeAlpha', mfilename);
                  error(eid,msg);
              end
          case {'disk','motion'}
              checkinput(p2,{'double'},{'positive','finite','real','nonempty','scalar'},mfilename,'RADIUS or LEN',2);
          case {'gaussian','log','average'}
              checkinput(p2,{'double'},{'positive','finite','real','nonempty','integer'},mfilename,'HSIZE',2);
              if numel(p2) > 2
                  msg = 'HSIZE should have 1 or 2 elements.';
                  eid = sprintf('Images:%s:wrongSizeN', mfilename);
                  error(eid,msg);
              elseif (numel(p2)==1)
                  p2 = [p2 p2]; 
              end
       end       
    case 3
       p2 = varargin{2};       p3 = varargin{3};       
       switch type
          case 'motion'
              checkinput(p2,{'double'},{'positive','finite','real','nonempty','scalar'},mfilename,'LEN',2);
              checkinput(p3,{'double'},{'real','nonempty','finite','scalar'},mfilename,'THETA',3);
          case {'gaussian','log'}
              checkinput(p2,{'double'},{'positive','finite','real','nonempty','integer'},mfilename,'N',2);
              checkinput(p3,{'double'},{'positive','finite','real','nonempty','scalar'},mfilename,'SIGMA',3);
              if numel(p2) > 2
                  msg = sprintf('%s: size(N) should be less than or equal 2.', upper(mfilename));
                  eid = sprintf('Images:%s:wrongSizeN', mfilename);
                  error(eid,msg);
              elseif (numel(p2)==1)
                  p2 = [p2 p2]; 
              end
          otherwise   
              msg = sprintf('%s: Too many arguments for this type of filter.', upper(mfilename));
              eid = sprintf('Images:%s:tooManyArgsForThisFilter', mfilename);
              error(eid,msg);
      end
end

%----------------------------------------------------------------------------------
function h = fwind1(f1,f2,hd,w1,w2)
%FWIND1 Design 2-D FIR filter using 1-D window method.
%   Copyright 1993-2003 The MathWorks, Inc.  
%   $Revision: 5.20 $  $Date: 2003/01/17 16:27:30 $

checknargin(2,5,nargin,mfilename);
method1 = 'linear';     method2 = 'bilinear';

if nargin==2		% Uniform spacing case with Huang's method
  hd = f1; w1 = f2;
  n = length(w1); m = n;
elseif nargin==3	% Uniform spacing with separable window
  w2 = hd;  hd = f1; w1 = f2;
  n = length(w1); m = length(w2);
elseif nargin==4	% Non-uniform spacing with Huang's method
  n = length(w1); m = n;
else
  n = length(w1); m = length(w2);
end

% Create 2-D window using either Huang's method or separable method.
if nargin==2 || nargin==4 % Huang's method: Create 2-D circular window
  if any(abs(w1-rot90(w1,2))>sqrt(eps))
    msg = '1-D window must be symmetric to use Huang''s method.';
    eid = sprintf('Images:%s:oneDWindowMustBeSymmetric',mfilename);
    error(eid,msg);
  end
  if length(w1)<2
    msg = 'Length of window must be greater than 1.';
    eid = sprintf('Images:%s:windowLengthMustBeGreaterThanOne',mfilename);
    error(eid,msg);
  end

  t = (-(n-1)/2:(n-1)/2)*(2/(n-1));
  [t1,t2] = meshgrid(t,t);
  t12 = sqrt(t1.*t1 + t2.*t2);
  d = find(t12<t(1) | t12>t(length(w1)));
  if ~isempty(d), t12(d) = zeros(size(d)); end
  w = zeros(size(t12)); w(:) = interp1(t,w1,t12(:),method1);
  if ~isempty(d), w(d) = zeros(size(d)); end

else % Create separable window
  w = w2(:)*w1(:).';
end

% Design filter using fsamp2 and apply window
if nargin<4		% Uniformly spaced data
  % Interpolate Hd to be the same size as W, if necessary
  if any([m n]~=size(hd))
    if any(size(hd)<[2 2])
        eid = sprintf('Images:%s:hdMustHaveAtLeast2rowsAnd2cols',mfilename);
        error(eid,'Hd must have at least 2 rows and 2 columns.');
    end

    [f1,f2] = freqspace(size(hd));
    % Extrapolate hd so that interpolation is never out of range.
    [mh,nh] = size(hd);
    if floor(nh/2)==nh/2		% if even
      hd = [hd,hd(:,1)]; f1 = [f1 1];
    else
      hd = [zeros(mh,1) hd zeros(mh,1)]; 
      df = f1(2)-f1(1); f1 = [f1(1)-df f1 f1(nh)+df];
    end
    [mh,nh] = size(hd);
    if floor(mh/2)==mh/2		% if even
      hd = [hd;hd(1,:)]; f2 = [f2 1];
    else
      hd = [zeros(1,nh);hd;zeros(1,nh)]; 
      df = f2(2)-f2(1); f2 = [f2(1)-df f2 f2(mh)+df];
    end
    [t1,t2] = freqspace([m n],'meshgrid');
    
    % Promote to double for call to interp2
    if ~isa(hd,'double')
       hd = double(hd);
    end
    hd = interp2(f1,f2,hd,t1,t2,method2);
    d = find(isnan(hd)); if ~isempty(d), hd(d) = zeros(size(d)); end
  end
  h = fsamp2(hd) .* w;
else % Non-uniformly spaced data
  h = fsamp2(f1,f2,hd,size(w)) .* w;
end

%----------------------------------------------------------------------------------
function [out,T] = histeq_j(a,cm,hgram)
%HISTEQ Enhance contrast using histogram equalization.
%   Copyright 1993-2003 The MathWorks, Inc.  
%   $Revision: 5.20 $  $Date: 2003/01/27 20:15:58 $

checknargin(1,3,nargin,mfilename);
checkinput(a,'uint8 uint16 double','nonsparse 2d',mfilename,'I',1);
NPTS = 256;
if nargin==1		% Histogram equalization of intensity image
    n = 64; % Default n
    hgram = ones(1,n)*(numel(a)/n);
    n = NPTS;
    kind = 1;
elseif nargin==2
    if numel(cm)==1
        % histeq(I,N); Histogram equalization of intensity image
        m = cm;
        hgram = ones(1,m)*(numel(a)/m);
        n = NPTS;
        kind = 1;
    elseif size(cm,2)==3 && size(cm,1)>1
        % histeq(X,map); Histogram equalization of indexed image
        n = size(cm,1);
        hgram = ones(1,n)*(numel(a)/n);
        kind = 2;
        if isa(a, 'uint16')
            msg = 'Histogram equalization of UINT16 indexed images is not supported.';
            eid = sprintf('Images:%s:unsupportedUint16IndexedImages',mfilename);
            error(eid, msg);
        end
	else        % histeq(I, HGRAM); Histogram modification of intensity image
        hgram = cm;
        n = NPTS;
        kind = 1;
    end
else		% Histogram modification of indexed image
    n = size(cm,1);
    if length(hgram)~=n
        msg = 'HGRAM must be the same size as MAP.';
        eid = sprintf('Images:%s:HGRAMmustBeSameSizeAsMAP',mfilename);
        error(eid, msg);
    end
    kind = 2;
    if isa(a, 'uint16')
        msg = 'Histogram equalization of UINT16 indexed images is not supported.';
        eid = sprintf('Images:%s:unsupportedUint16IndexedImages',mfilename);
        error(eid, msg);
    end
end

if min(size(hgram))>1
    msg = 'HGRAM must be a vector.';
    eid = sprintf('Images:%s:hgramMustBeAVector',mfilename);
    error(eid, msg);
end

% Normalize hgram 
hgram = hgram*(numel(a)/sum(hgram));       % Set sum = prod(size(a))
m = length(hgram);

% Compute cumulative histograms
if kind==1
    nn = imhist_j(a,n)';
    cum = cumsum(nn);
else    % Convert image to equivalent gray image
    I = ind2gray(a,cm);
    nn = imhist_j(I,n)';
    cum = cumsum(nn);
end
cumd = cumsum(hgram*numel(a)/sum(hgram));

% Create transformation to an intensity image by minimizing the error
% between desired and actual cumulative histogram.
tol = ones(m,1)*min([nn(1:n-1),0;0,nn(2:n)])/2;
err = (cumd(:)*ones(1,n)-ones(m,1)*cum(:)')+tol;
d = find(err < -numel(a)*sqrt(eps));
if ~isempty(d), err(d) = numel(a)*ones(size(d)); end
[dum,T] = min(err);
T = (T-1)/(m-1); 

if kind == 1	% Modify intensity image
    b = grayxform(a, T);
else % Modify colormap by extending the (r,g,b) vectors.
  % Compute equivalent colormap luminance
  ntsc = rgb2ntsc(cm);
  % Map to new luminance using T, store in 2nd column of ntsc.
  ntsc(:,2) = T(floor(ntsc(:,1)*(n-1))+1)';
  % Scale (r,g,b) vectors by relative luminance change
  map = cm.*((ntsc(:,2)./max(ntsc(:,1),eps))*ones(1,3));
  % Clip the (r,g,b) vectors to the unit color cube
  map = map ./ (max(max(map')',1)*ones(1,3));
end

if (kind==1),	out = b;
else,			out = map;
end

% ----------------------------------------------------------------------------
function [h, theta, rho] = hough(varargin)
%HOUGH Hough transform.
%   HOUGH implements the Standard Hough Transform (SHT). HOUGH is designed
%   to detect lines. It uses the parametric representation of a line:
%
%                       rho = x*cos(theta) + y*sin(theta).
%   
%   The variable rho is the distance from the origin to the line along a 
%   vector perpendicular to the line.  Theta is the angle between
%   the x-axis and this vector. HOUGH generates parameter space matrix 
%   whose rows and columns correspond to rho and theta values respectively.
%   Peak values in this space represent potential lines in the input image.
% 
%   [H, THETA, RHO] = HOUGH(BW) computes the SHT of the binary image BW.
%   THETA (in degrees) and RHO are the arrays of rho and theta values over 
%   which the Hough transform matrix, H, was generated.
%   
%   [H, THETA, RHO] = HOUGH(BW,PARAM1,VAL1,PARAM2,VAL2) sets various
%   parameters.  Parameter names can be abbreviated, and case does not
%   matter. Each string parameter is followed by a value as indicated 
%   below:
%
%   'ThetaResolution' Real scalar between 0 and 90, exclusive.
%                     'ThetaResolution' specifies the spacing (in degrees)
%                     of the Hough transform bins along the theta axis.
%
%                     Default: 1.
%
%   'RhoResolution'   Real scalar between 0 and norm(size(BW)), exclusive.
%                     'RhoResolution' specifies the spacing of the Hough
%                     transform bins along the rho axis.
%
%                     Default: 1.
%
%   Notes
%   -----
%   The Hough transform matrix, H, is NRHO-by-NTHETA where 
%   NRHO = 2*ceil(norm(size(BW))/RhoResolution)-1, and 
%   NTHETA = 2*ceil(90/ThetaResolution). Theta angle values are in 
%   the range [-90, 90) degrees and rho values range from -DIAGONAL to 
%   DIAGONAL where DIAGONAL = RhoResolution*ceil(norm(size(BW))/
%   RhoResolution). Note that if 90/DTHETA is not an integer, 
%   the actual angle spacing will be 90/ceil(90/DTHETA).
%
%   Class Support
%   -------------
%   BW can be logical or numeric and it must be real, 2-D, and nonsparse.
%
%   Example
%   -------
%   Compute and display the Hough transform of the gantrycrane.png image
%
%      RGB = imread('gantrycrane.png');
%      I  = rgb2gray(RGB); % convert to intensity
%      BW = edge(I,'canny'); % extract edges
%      [H,T,R] = hough(BW,'RhoResolution',0.5,'ThetaResolution',0.5);
%
%      % display the original image
%      subplot(2,1,1);
%      imshow(RGB);
%      title('gantrycrane.png');
%
%      % display the hough matrix
%      subplot(2,1,2);
%      imshow(imadjust(mat2gray(H)),'XData',T,'YData',R,...
%             'InitialMagnification','fit');
%      title('Hough transform of gantrycrane.png');
%      xlabel('\theta'), ylabel('\rho');
%      axis on, axis normal, hold on;
%      colormap(hot);
%
%   See also HOUGHPEAKS and HOUGHLINES.

%   Copyright 1993-2004 The MathWorks, Inc.
%   $Revision: 1.1.8.1 $  $Date: 2004/08/10 01:39:43 $

%   References:
%   Rafael C. Gonzalez, Richard E. Woods, Steven L. Eddins, "Digital
%   Image Processing Using MATLAB", Prentice Hall, 2004

[bw, rho, theta] = parseInputs_hough(varargin{:});

h = houghmex(bw,rho,theta*pi/180);

%-----------------------------------------------------------------------------
function [bw, rho, theta] = parseInputs_hough(varargin)

iptchecknargin(1,5,nargin,mfilename);

bw     = varargin{1};
checkinput(bw, {'numeric','logical'},{'real', '2d', 'nonsparse', 'nonempty'}, mfilename, 'BW', 1);

if ~islogical(bw),      bw = bw~=0;     end

% Set the defaults
theta_res = 1;      rho_res = 1;

% Process parameter-value pairs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% validStrings = {'ThetaResolution','RhoResolution'};

if (nargin > 1)
    done = false;
    idx = 2;
    while ~done
        input = varargin{idx};
        %inputStr = iptcheckstrs(input, validStrings,mfilename,'PARAM',idx);
        inputStr = lower(input);
        idx = idx+1; %advance index to point to the VAL portion of the input 
    
        if idx > nargin
            msg = sprintf('Parameter ''%s'' must be followed by a value.', inputStr);
            error(sprintf('Images:%s:valFor%sMissing', mfilename, inputStr),'%s', msg);        
        end

        switch inputStr
            case 'thetaresolution'
                theta_res = varargin{idx};
                checkinput(theta_res, {'double'}, {'real', 'scalar', 'finite','positive'}, mfilename, inputStr, idx);

                if (theta_res >= 90)
                    msg = sprintf('Value of parameter,''%s'', must be between 0 and 90.', inputStr);
                    error(sprintf('Images:%s:invalidThetaRes', mfilename),'%s', msg);
                end
      
            case 'rhoresolution'
                rho_res = varargin{idx};
                checkinput(rho_res, {'double'}, {'real', 'scalar', 'finite','positive'}, mfilename, inputStr, idx);

                if (rho_res >= norm(size(bw)))
                    msg = sprintf('Value of parameter,''%s'', must be between 0 and norm(size(BW)).', inputStr);
                    error(sprintf('Images:%s:invalidThetaRes', mfilename),'%s', msg);
                end
     
            otherwise      %should never get here
                error(sprintf('Images:%s:internalError', mfilename),'%s', 'Unknown input string.');
        end
    
        if (idx >= nargin),     done = true;    end
        idx=idx+1;
    end
end

% Compute theta and rho
%%%%%%%%%%%%%%%%%%%%%%%
[M,N] = size(bw);
theta = linspace(-90, 0, ceil(90/theta_res) + 1);
theta = [theta -fliplr(theta(2:end - 1))];

D = sqrt((M - 1)^2 + (N - 1)^2);
q = ceil(D/rho_res);
nrho = 2*q - 1;
rho = linspace(-q*rho_res, q*rho_res, nrho);

%-----------------------------------------------------------------------------
function lines = houghlines(varargin)
%HOUGHLINES Extract line segments based on Hough transform.
%   LINES = HOUGHLINES(BW, THETA, RHO, PEAKS) extracts line segments
%   in the image BW associated with particular bins in a Hough 
%   transform.  THETA and RHO are vectors returned by function HOUGH.
%   Matrix PEAKS, which is returned by function HOUGHPEAKS,
%   contains the row and column coordinates of the Hough transform 
%   bins to use in searching for line segments. HOUGHLINES returns
%   LINES structure array whose length equals the number of merged
%   line segments found. Each element of the structure array has
%   these fields: 
%
%      point1  End-point of the line segment; two-element vector
%      point2  End-point of the line segment; two-element vector
%      theta   Angle (in degrees) of the Hough transform bin
%      rho     Rho-axis position of the Hough transform bin
%
%   The end-point vectors contain [X, Y] coordinates.
%
%   LINES = HOUGHLINES(...,PARAM1,VAL1,PARAM2,VAL2) sets various
%   parameters. Parameter names cannot be abbreviated, and case 
%   does not matter. Each string parameter is followed by a value 
%   as indicated below:
%
%   'FillGap'   Positive real scalar.
%               When HOUGHLINES finds two line segments associated
%               with the same Hough transform bin that are separated
%               by less than 'FillGap' distance, HOUGHLINES merges
%               them into a single line segment.
%
%               Default: 20
%
%   'MinLength' Positive real scalar.
%               Merged line segments shorter than 'MinLength'
%               are discarded.
%
%               Default: 40
%
%   Class Support
%   -------------
%   BW can be logical or numeric and it must be real, 2-D, and nonsparse.
%
%   Example
%   -------
%   Search for line segments corresponding to five peaks in the Hough 
%   transform of the rotated circuit.tif image. Additionally, highlight
%   the longest segment.
%
%      I  = imread('circuit.tif');
%      rotI = imrotate(I,33,'crop');
%      BW = edge(rotI,'canny');
%      [H,T,R] = hough(BW);
%      imshow(H,[],'XData',T,'YData',R,'InitialMagnification','fit');
%      xlabel('\theta'), ylabel('\rho');
%      axis on, axis normal, hold on;
%      P  = houghpeaks(H,5,'threshold',ceil(0.3*max(H(:))));
%      x = T(P(:,2)); y = R(P(:,1));
%      plot(x,y,'s','color','white');
%
%      % Find lines and plot them
%      lines = houghlines(BW,T,R,P,'FillGap',5,'MinLength',7);
%      figure, imshow(rotI), hold on
%      max_len = 0;
%      for k = 1:length(lines)
%        xy = [lines(k).point1; lines(k).point2];
%        plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
%
%        % plot beginnings and ends of lines
%        plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
%        plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');
%
%        % determine the endpoints of the longest line segment 
%        len = norm(lines(k).point1 - lines(k).point2);
%        if ( len > max_len)
%          max_len = len;
%          xy_long = xy;
%        end
%      end
%
%      % highlight the longest line segment
%      plot(xy_long(:,1),xy_long(:,2),'LineWidth',2,'Color','cyan');
%
%   See also HOUGH and HOUGHPEAKS.

%   Copyright 1993-2005 The MathWorks, Inc.
%   $Revision: 1.1.8.3 $  $Date: 2006/06/15 20:08:51 $

%   References:
%   Rafael C. Gonzalez, Richard E. Woods, Steven L. Eddins, "Digital
%   Image Processing Using MATLAB", Prentice Hall, 2003

[nonzeropix,theta,rho,peaks,fillgap,minlength] = parseInputs_houghlines(varargin{:});

minlength_sq = minlength^2;
fillgap_sq = fillgap^2;
numlines = 0; lines = struct;
for k = 1:size(peaks,1)
    % Get all pixels associated with Hough transform cell.
    [r, c] = houghpixels(nonzeropix, theta, rho, peaks(k,:));
    if (isempty(r)),      continue;       end
  
    % Compute distance^2 between the point pairs
    xy = [c r]; % x,y pairs in coordinate system with the origin at (1,1)
    diff_xy_sq = diff(xy,1,1).^2;
    dist_sq = sum(diff_xy_sq,2);
  
    % Find the gaps larger than the threshold.
    fillgap_idx = find(dist_sq > fillgap_sq);
    idx = [0; fillgap_idx; size(xy,1)];
    for p = 1:length(idx) - 1
        p1 = xy(idx(p) + 1,:); % offset by 1 to convert to 1 based index
        p2 = xy(idx(p + 1),:); % set the end (don't offset by one this time)

        linelength_sq = sum((p2-p1).^2);
        if linelength_sq >= minlength_sq
            numlines = numlines + 1;
            lines(numlines).point1 = p1;
            lines(numlines).point2 = p2;
            lines(numlines).theta = theta(peaks(k,2));
            lines(numlines).rho = rho(peaks(k,1));
        end
    end
end

%-----------------------------------------------------------------------------
function [r, c] = houghpixels(nonzeropix, theta, rho, peak)
%HOUGHPIXELS Compute image pixels belonging to Hough transform bin.
%   [R, C] = HOUGHPIXELS(NONZEROPIX, THETA, RHO, PEAK) computes the
%   row-column indices (R, C) for nonzero pixels NONZEROPIX that map
%   to a particular Hough transform bin, PEAK which is a two element
%   vector [RBIN CBIN].  RBIN and CBIN are scalars indicating the 
%   row-column bin location in the Hough transform matrix returned by
%   function HOUGH.  THETA and RHO are the second and third output 
%   arguments from the HOUGH function.

x = nonzeropix(:,1);        y = nonzeropix(:,2);

theta_c = theta(peak(2)) * pi / 180;
rho_xy = x*cos(theta_c) + y*sin(theta_c);
nrho = length(rho);
slope = (nrho - 1)/(rho(end) - rho(1));
rho_bin_index = round(slope*(rho_xy - rho(1)) + 1);

idx = find(rho_bin_index == peak(1));
r = y(idx) + 1;     c = x(idx) + 1;

%-----------------------------------------------------------------------------
function [nonzeropix,theta,rho,peaks,fillgap,minlength] = parseInputs_houghlines(varargin)

iptchecknargin(1,8,nargin,mfilename);

idx = 1;
bw = varargin{idx};
checkinput(bw, {'numeric','logical'}, {'real', '2d', 'nonsparse', 'nonempty'}, mfilename, 'BW', idx);

idx = idx+1;
theta = varargin{idx};
checkinput(theta, {'double'}, {'real','vector','finite','nonsparse','nonempty'}, mfilename, 'THETA', idx);

idx = idx+1;
rho = varargin{idx};
checkinput(rho, {'double'}, {'real','vector','finite','nonsparse','nonempty'}, mfilename, 'RHO', idx);

idx = idx+1;
peaks = varargin{idx};
checkinput(peaks, {'double'}, {'real','2d','nonsparse','integer'}, mfilename, 'PEAKS', idx);

if size(peaks,2) ~= 2
  msg = sprintf('PEAKS must be a Q-by-2 matrix');
  error(sprintf('Images:%s:invalidPEAKS', mfilename),'%s',msg);
end

% Set the defaults
fillgap = 20;       minlength = 40; 

% Process parameter-value pairs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% validStrings = {'FillGap','MinLength'};
idx = idx+1;

if (nargin > idx-1) % we have parameter/value pairs
    done = false;

    while ~done
        input = varargin{idx};
        %inputStr = iptcheckstrs(input, validStrings,mfilename,'PARAM',idx);
        inputStr = lower(input);
        idx = idx+1;    % advance index to point to the VAL portion of the input 
    
        if (idx > nargin)
            msg = sprintf('Parameter ''%s'' must be followed by a value.', inputStr);
            error(sprintf('Images:%s:valFor%sMissing', mfilename, inputStr),'%s',msg);        
        end
    
        switch inputStr
            case 'fillgap'
                fillgap = varargin{idx};
                checkinput(fillgap, {'double'}, {'finite','real', 'scalar', 'positive'}, mfilename, inputStr, idx);
            case 'minlength'
                minlength = varargin{idx};
                checkinput(minlength, {'double'}, {'finite','real', 'scalar', 'positive'}, mfilename, inputStr, idx);
            otherwise      %should never get here
                error(sprintf('Images:%s:internalError', mfilename),'%s','Unknown input string.');
        end
    
        if (idx >= nargin),     done = true;    end
        idx=idx+1;
    end
end

% Compute the required parameters
[y, x] = find(bw);
nonzeropix = [x, y] - 1;

%-----------------------------------------------------------------------------
function peaks = houghpeaks(varargin)
%HOUGHPEAKS Identify peaks in Hough transform.
%   PEAKS = HOUGHPEAKS(H,NUMPEAKS) locates peaks in the Hough 
%   transform matrix, H, generated by the HOUGH function. NUMPEAKS 
%   specifies the maximum number of peaks to identify. PEAKS is 
%   a Q-by-2 matrix, where Q can range from 0 to NUMPEAKS. Q holds
%   the row and column coordinates of the peaks. If NUMPEAKS is 
%   omitted, it defaults to 1.
%
%   PEAKS = HOUGHPEAKS(...,PARAM1,VAL1,PARAM2,VAL2) sets various 
%   parameters. Parameter names can be abbreviated, and case 
%   does not matter. Each string parameter is followed by a value 
%   as indicated below:
%
%   'Threshold' Nonnegative scalar.
%               Values of H below 'Threshold' will not be considered
%               to be peaks. Threshold can vary from 0 to Inf.
%   
%               Default: 0.5*max(H(:))
%
%   'NHoodSize' Two-element vector of positive odd integers: [M N].
%               'NHoodSize' specifies the size of the suppression
%               neighborhood. This is the neighborhood around each 
%               peak that is set to zero after the peak is identified.
%
%               Default: smallest odd values greater than or equal to
%                        size(H)/50.
%
%   Class Support
%   -------------
%   H is the output of the HOUGH function. NUMPEAKS is a positive
%   integer scalar.
%
%   Example
%   -------
%   Locate and display two peaks in the Hough transform of the 
%   rotated circuit.tif image.
%
%      I  = imread('circuit.tif');
%      BW = edge(imrotate(I,50,'crop'),'canny');
%      [H,T,R] = hough(BW);
%      P  = houghpeaks(H,2);
%      imshow(H,[],'XData',T,'YData',R,'InitialMagnification','fit');
%      xlabel('\theta'), ylabel('\rho');
%      axis on, axis normal, hold on;
%      plot(T(P(:,2)),R(P(:,1)),'s','color','white');
%
%   See also HOUGH and HOUGHLINES.

%   Copyright 1993-2005 The MathWorks, Inc.
%   $Revision: 1.1.8.5 $  $Date: 2006/06/15 20:08:52 $

%   References:
%   Rafael C. Gonzalez, Richard E. Woods, Steven L. Eddins, "Digital
%   Image Processing Using MATLAB", Prentice Hall, 2004

[h, numpeaks, threshold, nhood] = parseInputs_houghpeaks(varargin{:});

% initialize the loop variables
done = false;
hnew = h;
nhood_center = (nhood-1)/2;
peaks = [];

while ~done
  [dummy, max_idx] = max(hnew(:));
  [p, q] = ind2sub(size(hnew), max_idx);
  
  p = p(1); q = q(1);
  if hnew(p, q) >= threshold
    peaks = [peaks; [p q]]; % add the peak to the list
    
    % Suppress this maximum and its close neighbors.
    p1 = p - nhood_center(1); p2 = p + nhood_center(1);
    q1 = q - nhood_center(2); q2 = q + nhood_center(2);
    % Throw away neighbor coordinates that are out of bounds in
    % the rho direction.
    [qq, pp] = meshgrid(q1:q2, max(p1,1):min(p2,size(h,1)));
    pp = pp(:); qq = qq(:);
    
    % For coordinates that are out of bounds in the theta
    % direction, we want to consider that H is antisymmetric
    % along the rho axis for theta = +/- 90 degrees.
    theta_too_low = find(qq < 1);
    qq(theta_too_low) = size(h, 2) + qq(theta_too_low);
    pp(theta_too_low) = size(h, 1) - pp(theta_too_low) + 1;
    theta_too_high = find(qq > size(h, 2));
    qq(theta_too_high) = qq(theta_too_high) - size(h, 2);
    pp(theta_too_high) = size(h, 1) - pp(theta_too_high) + 1;
    
    % Convert to linear indices to zero out all the values.
    hnew(sub2ind(size(hnew), pp, qq)) = 0;
    done = size(peaks,1) == numpeaks;
  else
    done = true;
  end
end

%-----------------------------------------------------------------------------
function [h, numpeaks, threshold, nhood] = parseInputs_houghpeaks(varargin)

iptchecknargin(1,6,nargin,mfilename);

h = varargin{1};
checkinput(h, {'double', 'single'}, {'real', '2d', 'nonsparse', 'nonempty',...
                   'finite', 'integer'}, mfilename, 'H', 1);
numpeaks = 1;   % set default value for numpeaks

idx = 2;
if nargin > 1
    if ~ischar(varargin{2})
        numpeaks = varargin{2};
        checkinput(numpeaks, {'double', 'single'}, {'real', 'scalar', 'integer', ...
                            'positive', 'nonempty'}, mfilename, 'NUMPEAKS', 2);
        idx = 3;
    end
end

% Initialize to empty
nhood = [];     threshold = [];

% Process parameter-value pairs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% validStrings = {'Threshold','NHoodSize'};

if nargin > idx-1 % we have parameter/value pairs
    done = false;
    while ~done
        input = varargin{idx};
        %inputStr = iptcheckstrs(input, validStrings,mfilename,'PARAM',idx);
        inputStr = lower(input);
        idx = idx+1; %advance index to point to the VAL portion of the input 
    
        if idx > nargin
            msg = sprintf('Parameter ''%s'' must be followed by a value.', inputStr);
            error(sprintf('Images:%s:valFor%sMissing', mfilename, inputStr),'%s', msg);        
        end
    
        switch inputStr
            case 'threshold'
                threshold = varargin{idx};
                checkinput(threshold, {'double'}, {'real', 'scalar','nonnan' ...
                                  'nonnegative'}, mfilename, inputStr, idx);
            case 'nhoodsize'
                nhood = varargin{idx};
                checkinput(nhood, {'double'}, {'real', 'vector', 'finite', ...
                                  'integer','positive','odd'}, mfilename, inputStr, idx);
      
                if (any(size(nhood) ~= [1,2]))
                    msg = sprintf(['Value of parameter,''%s'', must be a two element row vector.'], inputStr);
                    error(sprintf('Images:%s:invalidNHoodSize', mfilename),'%s', msg);
                end
      
                if (any(nhood > size(h)))
                    msg = sprintf(['Dimensions specified by ''%s'', '...
                               'must be smaller than '...
                               'size of Hough matrix H.'], inputStr);
                    error(sprintf('Images:%s:tooBigNHoodSize', mfilename),'%s', msg);
                end     
            otherwise      %should never get here
                    error(sprintf('Images:%s:internalError', mfilename),'%s', 'Unknown input string.');
        end
    
        if (idx >= nargin),     done = true;    end
        idx=idx+1;
    end
end

% Set the defaults if necessary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(nhood)
    nhood = size(h)/50; 
    nhood = max(2*ceil(nhood/2) + 1, 1); % Make sure the nhood size is odd.
end

if (isempty(threshold)),    threshold = 0.5 * max(h(:));    end

%----------------------------------------------------------------------------------
function [BW,level] = im2bw(varargin)
%IM2BW Convert image to binary image by thresholding.
%   IM2BW produces binary images from indexed, intensity, or RGB images. To do this, it converts the 
%   input image to grayscale format (if it is not already an intensity image), and then converts 
%   this grayscale image to binary by thresholding. The output binary image BW has values of 0 (black)
%   for all pixels in the input image with luminance less than LEVEL and 1 (white) for all other pixels.
%   (Note that you specify LEVEL in the range [0,1], regardless of the class of the input image.)
%  
%   BW = IM2BW(I,LEVEL) converts the intensity image I to black
%   and white.
%
%   BW = IM2BW(X,MAP,LEVEL) converts the indexed image X with
%   colormap MAP to black and white.
%
%   BW = IM2BW(RGB,LEVEL) converts the RGB image RGB to black and
%   white.
%
%   [BW,LEVEL] = IM2BW(...) also returns the LEVEL corresponding to the input class
%
%   If LEVEL is not provided it will be computed from the function GRAYTHRESH 
%
%   Class Support
%   -------------
%   The input image can be of class uint8, uint16, or double and
%   it must be nonsparse. The output image BW is of class logical.
%

%   Copyright 1993-2003 The MathWorks, Inc.  
%   $Revision: 5.26 $  $Date: 2003/01/27 20:15:58 $

[A,map,level] = parse_inputs_im2bw(varargin{:});

if ndims(A)==3              % RGB is given
    A = rgb2gray(A);
elseif ~isempty(map)        % indexed image is given
    A = ind2gray(A,map);
end                         % nothing to do for intensity image

if (isempty(level))
    level = graythresh(A);
end

if isa(A, 'uint8')
    A = (A >= 255*level);       level = 255*level;          % In case we want it to ouput
elseif isa(A, 'uint16')
    A = (A >= 65535*level);     level = 65535*level;
elseif isa(A, 'logical')
    % A is already a binary image and does not require thresholding 
else    % must be double
    A = (A >= level);
end
  
BW = A;

% --------------------------------------------------------------------------
function [A,map,level] = parse_inputs_im2bw(varargin)
% Defaults:
map = [];       level = [];

checknargin(1,4,nargin,mfilename);

checkinput(varargin{1}, {'uint8', 'uint16', 'logical', 'double'},...
    {'real', 'nonsparse'},mfilename,'I, X or RGB',1);

switch nargin
	case 1                          % im2bw(RGB) | im2bw(I)
        A = varargin{1};
	case 2
        A = varargin{1};            % im2bw(RGB,level) | im2bw(I,level)
        level = varargin{2};        % im2bw(X,MAP)
	case 3                          % im2bw(X,MAP,level)
        A = varargin{1};
        map = varargin{2};
        level = varargin{3};
end

% Check validity of the input parameters 
if (ndims(A) == 3) && (size(A,3)~=3)        % check RGB image array
    error('Truecolor RGB image has to be an M-by-N-by-3 array.');
end

if (nargin == 2) && (ndims(A) == 2) && (size(level,2) == 3)       % it is a colormap
    map = level;                    % and we assume that image given is an indexed image X
    level = [];
end

if ~isempty(map)                    % check colormap if given
    if (size(map,2) ~= 3) || ndims(map) > 2
        error('Input colormap has to be a 2D array with exactly 3 columns.');
    elseif (min(map(:))<0) || (max(map(:))>1)
        error('All colormap intensities must be between 0 and 1.');
    end
end

if ( ~isempty(level) && ((numel(level) ~= 1) || (max(level(:)) > 1) || (min(level(:)) < 0)) )
    error('Threshold luminance LEVEL has to be a non-negative number between 0 and 1');
end

%----------------------------------------------------------------------------------
function level = graythresh(I)
%GRAYTHRESH Compute global image threshold using Otsu's method.
%   LEVEL = GRAYTHRESH(I) computes a global threshold (LEVEL) that can be
%   used to convert an intensity image to a binary image with IM2BW. LEVEL
%   is a normalized intensity value that lies in the range [0, 1].
%   GRAYTHRESH uses Otsu's method, which chooses the threshold to minimize
%   the intraclass variance of the thresholded black and white pixels.
%
%   Class Support
%   -------------
%   The input image I can be of class uint8, uint16, or double and it
%   must be nonsparse.  LEVEL is a double scalar.

%   Copyright 1993-2003 The MathWorks, Inc.  
%   $Revision: 1.11 $  $Date: 2003/01/27 20:15:58 $

% One input argument required.
checknargin(1,1,nargin,mfilename);
checkinput(I,'uint8 uint16 double','nonsparse',mfilename,'I',1);

% Convert all N-D arrays into a single column.  Convert to uint8 for
% fastest histogram computation.
if (~isa(I,'uint8')),   I = im2uint8(I(:));     end

num_bins = 256;
counts = imhist_j(I,num_bins);

% Variables names are chosen to be similar to the formulas in the Otsu paper.
p = counts / sum(counts);
omega = cumsum(p);
mu = cumsum(p .* (1:num_bins)');
mu_t = mu(end);

% Save the warning state and disable warnings to prevent divide-by-zero warnings.
state = warning;
warning off;
sigma_b_squared = (mu_t * omega - mu).^2 ./ (omega .* (1 - omega));

warning(state);         % Restore the warning state.

% Find the location of the maximum value of sigma_b_squared.
% The maximum may extend over several bins, so average together the locations.
% If maxval is NaN, meaning that sigma_b_squared is all NaN, then return 0.
maxval = max(sigma_b_squared);
if isfinite(maxval)
    idx = mean(find(sigma_b_squared == maxval));
    % Normalize the threshold to the range [0, 1].
    level = (idx - 1) / (num_bins - 1);
else
    level = 0.0;
end

%----------------------------------------------------------------------------------
function d = im2double(img, typestr)
%IM2DOUBLE Convert image to double precision.
%   Copyright 1993-2003 The MathWorks, Inc.  
%   $Revision: 1.17 $  $Date: 2003/01/17 16:27:33 $

checknargin(1,2,nargin,mfilename);
checkinput(img,{'double','logical','uint8','uint16'},{},mfilename,'Image',1);
if isa(img, 'double')
   d = img;
elseif isa(img, 'logical')
   d = double(img);
elseif isa(img, 'uint8') || isa(img, 'uint16')
   if nargin==1
      if isa(img, 'uint8')
          d = double(img)/255;
      else
          d = double(img)/65535;
      end
   elseif nargin==2
      checkstrs(typestr,{'indexed'},mfilename,'type',2);
      d = double(img)+1;
   end
end

%----------------------------------------------------------------------------------
function u = im2uint8(varargin)
%IM2UINT8 Convert image to eight-bit unsigned integers.
%   Copyright 1993-2003 The MathWorks, Inc.  
%   $Revision: 1.21 $  $Date: 2003/01/17 16:27:34 $

checknargin(1,2,nargin,mfilename);
img = varargin{1};
checkinput(img,{'double','logical','uint8','uint16'},{},mfilename,'Image',1);
if isa(img, 'uint8')
    u = img; 
elseif isa(img, 'logical')
    u=uint8(img);
    u(img)=255;
elseif isa(img, 'double') || isa(img, 'uint16')
    if nargin==1         % intensity image; call MEX-file
         u = grayto8(img);
    elseif nargin==2
       typestr = varargin{2};
       checkstrs(typestr,{'indexed'},mfilename,'type',2);
       if (isa(img, 'uint16'))
            if (max(img(:)) > 255)
                eid = sprintf('Images:%s:tooManyColorsFor8bitStorage');
                error(eid,'Too many colors for 8-bit integer storage.');
            else
                u = uint8(img);
            end
        else          % img is double
            if max(img(:))>=257 
                eid = sprintf('Images:%s:tooManyColorsFor8bitStorage');
                error(eid,'Too many colors for 8-bit integer storage.');
            elseif min(img(:))<1
                eid = sprintf('Images:%s:invalidIndexedImage');
                error(eid,'Invalid indexed image: an index was less than 1.');
            else
                u = uint8(img-1);
            end
        end
    end
end

%----------------------------------------------------------------------------------
function u = im2uint16(img, typestr)
%IM2UINT16 Convert image to sixteen-bit unsigned integers.
%   Copyright 1993-2003 The MathWorks, Inc.  
%   $Revision: 1.13 $  $Date: 2003/01/17 16:27:34 $

checknargin(1,2,nargin,mfilename);
if isa(img, 'uint16')
    u = img; 
elseif isa(img, 'double') || isa(img, 'uint8')
    if (nargin==1)        % intensity image; call MEX-file
        u = grayto16(img);
    elseif nargin==2
        if ~ischar(typestr) || (typestr(1) ~= 'i')
            eid = 'Images:im2uint16:invalidInput';
            msg = 'Invalid input arguments';
            error(eid,'%s',msg);
        else 
            if (isa(img, 'uint8'))
                u = uint16(img);
            else          % img is double
                if max(img(:))>=65537 
                    eid = 'Images:im2uint16:tooManyColors';
                    msg = 'Too many colors for 16-bit integer storage.';
                    error(eid,'%s',msg);
                elseif min(img(:))<1
                    eid = 'Images:im2uint16:invalidIndex';
                    msg = 'Invalid indexed image: an index was less than 1.';
                    error(eid,'%s',msg);
                else
                    u = uint16(img-1);
                end
            end
        end
    end
    
elseif islogical(img)
    u = uint16(img);    u(img) = 65535;
else
    error('im2uint16:unsupportedInputClass','Unsupported input class.');
end

%----------------------------------------------------------------------------------
function Z = imadd(X,Y,output_class)
%IMADD Add two images, or add constant to image.
%   Copyright 1993-2003 The MathWorks, Inc.
%   $Revision: 1.12 $  $Date: 2003/03/05 19:42:51 $

checknargin(2, 3, nargin, mfilename);
if (nargin < 3)
    if (islogical(X)),	output_class = 'double';
	else,				output_class = class(X);
    end
else
    valid_strings = {'uint8' 'uint16' 'uint32' 'int8' 'int16' 'int32' 'single' 'double'};
    output_class = checkstrs(output_class, valid_strings, mfilename,'OUTPUT_CLASS', 3);
end
if (numel(Y) == 1) && isa(Y, 'double')
    Z = imlincomb(1.0, X, Y, output_class);
else
    Z = imlincomb(1.0, X, 1.0, Y, output_class);
end

%----------------------------------------------------------------------------------
function B = imdilate(A,se,varargin)
	checknargin(2,4,nargin,'imdilate');
	B = morphop(A,se,'dilate','imdilate',varargin{:});

%----------------------------------------------------------------------------------
function B = imerode(A,se,varargin)
	checknargin(2,5,nargin,'imerode');
	B = morphop(A,se,'erode','imerode',varargin{:});

%----------------------------------------------------------------------------------
function out = imadjust_j(varargin)
%IMADJUST Adjust image intensity values or colormap.

%Parse inputs and initialize variables
[img,grayFlag,rgbFlag,low_in,high_in,low_out,high_out,gamma] = ParseInputs_imadjust_j(varargin{:});
if ( isa(img,'uint8') || ( isa(img,'uint16') && numel(img) > 65536 ) )
	out = AdjustWithLUT(img,grayFlag,rgbFlag,low_in,high_in,low_out,high_out,gamma);
else
	out = AdjustGeneric(img,grayFlag,rgbFlag,low_in,high_in,low_out,high_out, gamma);
end

%------------------------------------
function out=AdjustWithLUT(img,grayFlag,rgbFlag,low_in,high_in,low_out,high_out,gamma)
out = img;      grayFlag = true;    rgbFlag = false;
for p=1:size(img,3)
  if isa(img,'uint8')
    lut = linspace(0,1,256);
    lut = AdjustGeneric(lut,grayFlag,rgbFlag,low_in(p),high_in(p),low_out(p), high_out(p),gamma(p));
    lut = im2uint8(lut);
  else    %uint16
    lut = linspace(0,1,65536);
    lut = AdjustGeneric(lut,grayFlag,rgbFlag,low_in(p),high_in(p),low_out(p), high_out(p),gamma(p));
    lut = im2uint16(lut);
  end
  out(:,:,p) = uintlut(img(:,:,p),lut);
end

%---------------------------------------------
function out = AdjustGeneric(img,grayFlag,rgbFlag,low_in,high_in,low_out,high_out,gamma)
%Promote image data to double if it is uint8 or uint16
classin = class(img);
classChanged = false;
if ~isa(img,'double')
  classChanged = true;
  img = im2double(img);
end

%Initialize output
out = zeros(size(img));

if ~rgbFlag
  % Transform colormap or Grayscale intensity image (v1 functionality)
  if (grayFlag),	 n = 1; 
  else,				n = size(img,1); 
  end
    
  % Make sure img is in the range [low_in,high_in]
  img(:) = max(low_in(ones(n,1),:),min(high_in(ones(n,1),:),img));
  
  % The transformation
  d = ones(n,1);
  out = ( (img-low_in(d,:))./(high_in(d,:)-low_in(d,:)) ).^(gamma(d,:));
  out(:) = out.*(high_out(d,:)-low_out(d,:)) + low_out(d,:);
else 
  % Loop over image planes and perform transformation
  for p=1:3    % Make sure img is in the range [low_in,high_in]
    img(:,:,p) =  max(low_in(p),min(high_in(p),img(:,:,p)));
    % The transformation
    out(:,:,p) = ( (img(:,:,p)-low_in(p))./(high_in(p)-low_in(p)) ).^(gamma(p));
    out(:,:,p) = out(:,:,p).*(high_out(p)-low_out(p)) + low_out(p);
  end
end

if classChanged
  out = changeclass(classin,out);
end
    
%--------------------------------------------------------------------
% Subfunction ParseInputs
function [img,grayFlag,rgbFlag,low_in,high_in,low_out,high_out,gamma] = ParseInputs_imadjust_j(varargin)
checknargin(1,4,nargin,mfilename);
% Default values
lowhigh_in  = [0; 1];
lowhigh_out = [0; 1];
gamma = 1;

img = varargin{1};
[rgbFlag, grayFlag] = getImType(img);

switch nargin
  case 1
    if grayFlag
      lowhigh_in = stretchlim_j(img);
    else
      error('imadjust_j:oneArgOnlyGrayscale','IMADJUST(I) is only supported for grayscale images.');
    end    
  case 2
    if ~isempty(varargin{2})
        lowhigh_in = varargin{2};
    end
  case 3
    if ~isempty(varargin{2})
        lowhigh_in = varargin{2};
    end
    if ~isempty(varargin{3})
        lowhigh_out = varargin{3};
    end
case 4
    if ~isempty(varargin{2})
        lowhigh_in = varargin{2};
    end
    if ~isempty(varargin{3})
        lowhigh_out = varargin{3};
    end
    if ~isempty(varargin{4})
        gamma = varargin{4};
    end
end

[low_in, high_in]   = range_split(lowhigh_in, grayFlag,2,'[LOW_IN; HIGH_IN]');
[low_out, high_out] = range_split(lowhigh_out,grayFlag,3,'[LOW_OUT; HIGH_OUT]');

check_lowhigh(low_in,high_in,low_out,high_out)
gamma = check_gamma(gamma,grayFlag);

%--------------------------------------------------------------------
function [rgbFlag, grayFlag] = getImType(img)
if (ndims(img)==3 && size(img,3)==3)  % RGB Intensity Image
  checkinput(img, {'double' 'uint8' 'uint16'}, '', mfilename, 'RGB1', 1);
  rgbFlag = true;
  grayFlag = false;  
elseif size(img,2)~=3  % Grayscale Intensity Image
  checkinput(img, {'double' 'uint8' 'uint16'}, '2d', mfilename, 'I', 1);
  grayFlag = true;
  rgbFlag = false;  
else  % Colormap
  checkmap(img,mfilename,'MAP',1);
  grayFlag = false;
  rgbFlag = false;  
end

%--------------------------------------------------------------------
function [range_min, range_max] = range_split(range,grayFlag,argument_position,variable_name)

msg1 = sprintf('Function %s expected its %s input argument, %s', 'range_split',...
               num2ordinal(argument_position),variable_name);

if grayFlag
  if numel(range) ~= 2
    eid = sprintf('Images:%s:InputMustBe2ElVec','range_split');
    error(eid, '%s\n%s', msg1, 'to be a two-element vector.');
  end
    range_min = range(1); range_max = range(2);
else
  if (numel(range) ~= 2) && ~isequal(size(range),[2 3])
    eid = sprintf('Images:%s:InputMustBe2ElVecOr2by3Matrix','range_split');
    error(eid, '%s\n%s', msg1, 'to be a two-element vector or a 2-by-3 matrix.');
  end
  if min(size(range))==1
    % Create triples for RGB image or Colormap
    range_min = range(1)*ones(1,3);    range_max = range(2)*ones(1,3);
  else
    range_min = range(1,:);    range_max = range(2,:);
  end
end

%--------------------------------------------------------------------
function check_lowhigh(low_in,high_in,low_out,high_out)
if any(low_in>=high_in)
    msg = sprintf('%s: LOW_IN must be less than HIGH_IN.',upper(mfilename));
    eid = sprintf('Images:%s:lowMustBeSmallerThanHigh',mfilename);
    error(eid,msg);
end

if min(low_in)<0 || max(low_in)>1 || min(high_in)<0 || max(high_in)>1 || ...
        min(low_out)<0 || max(low_out)>1 || min(high_out)<0  || max(high_out)>1
    msg = sprintf('%s: %s %s', upper(mfilename),...
                  'LOW_IN, HIGH_IN, LOW_OUT and HIGH_OUT',...
                  'must be in the range [0.0, 1.0].');
    eid = sprintf('Images:%s:parametersAreOutOfRange',mfilename);
    error(eid,msg);
end

%--------------------------------------------------------------------
function gamma = check_gamma(gamma,grayFlag)
if grayFlag
    checkinput(gamma,'double',{'scalar', 'nonnegative'}, mfilename, 'GAMMA', 4)
else
    checkinput(gamma,'double',{'nonnegative','2d'}, mfilename, 'GAMMA', 4)
    if (numel(gamma) == 1),    gamma = gamma*ones(1,3);   end
end

%----------------------------------------------------------------------------------
function b = imfilter(varargin)
%IMFILTER Multidimensional image filtering.

[a,h,boundary,flags] = parse_inputs_imfilter(varargin{:});  
rank_a = ndims(a);  rank_h = ndims(h);

% Pad dimensions with ones if filter and image rank are different
size_h = [size(h) ones(1,rank_a-rank_h)];
size_a = [size(a) ones(1,rank_h-rank_a)];

if bitand(flags,8)      %Full output
  im_size = size_a+size_h-1;
  pad = [size(h)-1 ones(1,rank_a-rank_h)];
else                    %Same output
  im_size = size(a);
  %Calculate the number of pad pixels
  filter_center = floor((size(h)+1)/2);
  pad = [size(h)-filter_center ones(1,rank_a-rank_h)];
end

%Empty Inputs
% 'Same' output then size(b) = size(a)
% 'Full' output then size(b) = size(h)+size(a)-1 
if isempty(a)
  if bitand(flags,4) %Same
    b = a;
  else %Full
    if all(im_size>0)
      b = a;      b = b(:);
      b(prod(im_size)) = 0;
      b = reshape(b,im_size);
    elseif all(im_size>=0)
      b = feval(class(a),zeros(im_size));
    else
      msg = ['Error in size of B.  At least one dimension is negative. ',...
             '\n''Full'' output size calculation is: size(B) = size(A) ',...
             '+ size(H) - 1.'];
      error('imfilter:negativeDimensionBadSizeB',msg);
    end
  end
  return
end

if  isempty(h)
  if bitand(flags,4) %Same
    b = a;    b(:) = 0;
  else %Full
    if all(im_size>0)
      b = a;
      if im_size < size_a  %Output is smaller than input
        b(:) = [];
      else %Grow the array, is this a no-op?
        b(:) = 0;        b = b(:);
      end
      b(prod(im_size)) = 0;
      b = reshape(b,im_size);
    elseif all(im_size>=0)
      b = feval(class(a),zeros(im_size));
    else
      eid = sprintf('imsize:negativeDimensionBadSizeB');
      msg = ['Error in size of B.  At least one dimension is negative. ',...
             '\n''Full'' output size calculation is: size(B) = size(A) +',...
             ' size(H) - 1.'];
      error('imfilter:negativeDimensionBadSizeB',msg);
    end
  end
  return
end

im_size = int32(im_size);

%Starting point in padded image, zero based.
start = int32(pad);

%Pad image
a = padarray(a,pad,boundary,'both');

%Create connectivity matrix.  Only use nonzero values of the filter.
conn_logical = h~=0;
conn = double( conn_logical );  %input to the mex file must be double
nonzero_h = h(conn_logical);

% Seperate real and imaginary parts of the filter (h) in the M-code and
% filter imaginary and real parts of the image (a) in the mex code. 
checkMexFileInputs(a,im_size,real(h),real(nonzero_h),conn,start,flags);
b1 = imfilter_mex(a,im_size,real(h),real(nonzero_h),conn,start,flags);

if ~isreal(h)
  checkMexFileInputs(a,im_size,imag(h),imag(nonzero_h),conn,start,flags);
  b2 = imfilter_mex(a,im_size,imag(h),imag(nonzero_h),conn,start,flags);
end

%If input is not complex, the output should not be complex. COMPLEX always
%creates an imaginary part even if the imaginary part is zeros.
if isreal(h)        % b will always be real
  b = b1;
elseif isreal(a)    % b1 and b2 will always be real. b will always be complex
  b = complex(b1,b2);
else                % b1 and/or b2 may be complex.  b will always be complex
  b = complex(imsubtract(real(b1),imag(b2)),imadd(imag(b1),real(b2)));
end

%======================================================================
function [a,h,boundary,flags] = parse_inputs_imfilter(a,h,varargin)
checknargin(2,5,nargin,mfilename);
checkinput(a,{'numeric' 'logical'},{'nonsparse'},mfilename,'A',1);
checkinput(h,{'double'},{'nonsparse'},mfilename,'H',2);

%Assign defaults
flags = 0;
boundary = 0;  %Scalar value of zero
output = 'same';
do_fcn = 'corr';

allStrings = {'replicate', 'symmetric', 'circular', 'conv', 'corr', 'full','same'};

for k = 1:length(varargin)
  if ischar(varargin{k})
    string = checkstrs(varargin{k}, allStrings,mfilename, 'OPTION',k+2);
    switch string
     case {'replicate', 'symmetric', 'circular'}
      boundary = string;
     case {'full','same'}
      output = string;
     case {'conv','corr'}
      do_fcn = string;
    end
  else
    checkinput(varargin{k},{'numeric'},{'nonsparse'},mfilename,'OPTION',k+2);
    boundary = varargin{k};
  end %else
end

if strcmp(output,'full'),       flags = bitor(flags,8);
elseif strcmp(output,'same'),   flags = bitor(flags,4);
end

if strcmp(do_fcn,'conv'),       flags = bitor(flags,2);
elseif strcmp(do_fcn,'corr'),   flags = bitor(flags,0);
end

%--------------------------------------------------------------
function checkMexFileInputs(varargin)
% a
a = varargin{1};
checkinput(a,{'numeric' 'logical'},{'nonsparse'},mfilename,'A',1);
% im_size
im_size = varargin{2};
if ~strcmp(class(im_size),'int32') || issparse(im_size)
  displayInternalError('im_size');
end
% h
h = varargin{3};
if ~isa(h,'double') || ~isreal(h) || issparse(h)
	displayInternalError('h');
end
% nonzero_h
nonzero_h = varargin{4};
if ~isa(nonzero_h,'double') || ~isreal(nonzero_h) || issparse(nonzero_h)
	displayInternalError('nonzero_h');
end
% start
start = varargin{6};
if ~strcmp(class(start),'int32') || issparse(start)
	displayInternalError('start');
end
% flags
flags = varargin{7};
if ~isa(flags,'double') ||  any(size(flags) ~= 1)
	displayInternalError('flags');
end

%--------------------------------------------------------------
function displayInternalError(string)
eid = sprintf(':%s:internalError',mfilename);
error(eid,'%s',sprintf('Internal error: %s is not valid.',upper(string)));

%----------------------------------------------------------------------------------
function im = dither(varargin)
%DITHER Convert image using dithering.

%   Copyright 1993-2003 The MathWorks, Inc.  
%   $Revision: 5.25 $  $Date: 2003/01/17 16:27:28 $

%   References: 
%      R. W. Floyd and L. Steinberg, "An Adaptive Algorithm for
%         Spatial Gray Scale," International Symposium Digest of Technical
%         Papers, Society for Information Displays, 36.  1975.
%      Spencer W. Thomas, "Efficient Inverse Color Map Computation",
%         Graphics Gems II, (ed. James Arvo), Academic Press: Boston.
%         1991. (includes source code)

[X,m,qm,qe] = dither_parse_inputs(varargin{:});

if (ndims(X) == 2)% Convert intensity image to binary by dithering
    im = logical(ditherc(X,m,qm,qe));
else % Create an indexed image from RGB 
    im = ditherc(X,m,qm,qe);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,m,qm,qe] = dither_parse_inputs(varargin)
% Outputs:  X  the input RGB (3D) or intensity image (2D)
%           m  colormap (:,3)
%           qm number of quantization bits for colormap
%           qe number of quantization bits for errors, qe>qm

checknargin(1,6,nargin,'dither');

% Default values:
qm = 5;     qe = 8;

switch nargin
	case 1                        % dither(I)
        X = varargin{1};
        m = gray(2); 
	case 2                        % dither(RGB,m)
        X = varargin{1};
        m = varargin{2};
	case 4
        X = varargin{1};
        m = varargin{2};
        qm = varargin{3};
        qe = varargin{4};
	otherwise
        eid = sprintf(':%s:invalidInput','dither');
        error(eid,'Invalid input arguments in function %s.','dither');
end

% Check validity of the input parameters 
if (ndims(X)==3) && (nargin==1)
    eid = sprintf(':%s:imageMustBe2D','dither');  
    error(eid,'DITHER(I): the intensity image I has to be a two-dimensional array.');
elseif (ndims(X)==2) && (nargin==2)
    eid = sprintf(':%s:imageMustBe3D','dither');  
    error(eid,'DITHER(RGB,map): the RGB image has to be a three-dimensional array.');
end

X = im2uint8(X);
 
if ((size(m,2) ~= 3) || (size(m,1) == 1) || ndims(m) > 2)
  eid = sprintf(':%s:colormapMustBe2D','dither');  
  error(eid,['In function %s, input colormap has to be a ',...
             '2D array with at least 2 rows and exactly 3 columns.'], 'dither');
end

%----------------------------------------------------------------------------------
function [yout,x] = imhist_j(varargin)
%IMHIST_J returns histogram of image data.
%   Copyright 1993-2003 The MathWorks, Inc.
%   $Revision: 5.25 $  $Date: 2003/01/17 16:27:37 $

[a, n, isScaled, top, map] = parse_inputs_imhist_j(varargin{:});
if islogical(a)
    if (n ~= 2)
        error(':imhist:invalidParameterForLogical', '%s', 'N must be set to two for a logical image.');
    end
    y(2) = sum(a(:));
    y(1) = numel(a) - y(2);
    y = y';
else
    y = imhistc(a, n, isScaled, top); % Call MEX file to do work.
end

switch class(a)
case 'uint8',    x = linspace(0,255,n)';
case 'uint16',   x = linspace(0,65535,n)';
case 'double',   x = linspace(0,1,n)';
case 'logical',  x = [0,1]';
otherwise
    message2 = ['The input image must be uint8, uint16, double, or' ' logical.'];
    error(':imhist:invalidImageClass','%s', message2);             
end
yout = y;

%-------------------------------------------------------------------------
function [a, n, isScaled, top, map] = parse_inputs_imhist_j(varargin)
checknargin(1,2,nargin,mfilename);
a = varargin{1};
checkinput(a, 'double uint8 logical uint16', '2d', mfilename, 'I or X', 1);
n = 256;

if (isa(a,'double'))
    isScaled = 1;    top = 1;    map = []; 
elseif (isa(a,'uint8'))
    isScaled = 1;    top = 255;  map = [];
elseif (islogical(a))
    n = 2;    isScaled = 1;      top = 1;    map = [];
else % a must be uint16
    isScaled = 1;    top = 65535;    map = [];
end
    
if (nargin ==2)
    if (numel(varargin{2}) == 1)        % IMHIST(I, N)
        n = varargin{2};
        checkinput(n, 'numeric', 'real positive integer', mfilename, 'N', 2);
    elseif (size(varargin{2},2) == 3)   % IMHIST(X,MAP) or invalid second argument
        n = size(varargin{2},1);        isScaled = 0;
        top = n;        map = varargin{2};
    else
        messageId = sprintf('Images:%s:invalidSecondArgument', mfilename);
        error(messageId, '%s','Second argument must be a colormap or a positive integer.'); 
    end    
end

%----------------------------------------------------------------------------------
function Z = imlincomb(varargin)
%IMLINCOMB Compute linear combination of images.
[images, scalars, output_class] = ParseInputs_imlincomb(varargin{:});
Z = imlincombc(images, scalars, output_class);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [images, scalars, output_class] = ParseInputs_imlincomb(varargin)
checknargin(2, Inf, nargin, 'imlincomb');

if ischar(varargin{end})
  valid_strings = {'uint8' 'uint16' 'uint32' 'int8' 'int16' 'int32' 'single' 'double'};
  output_class = checkstrs(varargin{end}, valid_strings, 'imlincomb', 'OUTPUT_CLASS', 3);
  varargin(end) = [];
else
  if islogical(varargin{2})
    output_class = 'double';
  else
    output_class = class(varargin{2});
  end
end

%check images
images = varargin(2:2:end);
if ~iscell(images) || isempty(images)
  displayInternalError('images');
end

% assign and check scalars
for p = 1:2:length(varargin)
  checkinput(varargin{p}, 'double', 'real nonsparse scalar', 'imlincomb', sprintf('K%d', (p+1)/2), p);
end
scalars = [varargin{1:2:end}];

%make sure it is a vector
if ( ndims(scalars)~=2 || (all(size(scalars)~=1) && any(size(scalars)~=0)) )
  displayInternalError('scalars');
end

%----------------------------------------------------------------------------------
function Z = imsubtract(X,Y)
%IMSUBTRACT Subtract two images, or subtract constant from image.
%   Copyright 1993-2003 The MathWorks, Inc. 
%   $Revision: 1.12 $  $Date: 2003/03/05 19:42:51 $

error(nargchk(2,2,nargin))
if (numel(Y) == 1) && strcmp(class(Y),'double')
    Z = imlincomb(1.0, X, -Y);
else
    Z = imlincomb(1.0, X, -1.0, Y);
end

%----------------------------------------------------------------------------------
function I = ind2gray(varargin)
%IND2GRAY Convert indexed image to intensity image.
%   Copyright 1993-2003 The MathWorks, Inc.  
%   $Revision: 5.17 $  $Date: 2003/01/27 20:16:03 $

[a,cm] = parse_inputs_ind2gray(varargin{:});
%initialize output matrix
I = a;

% calculate gray colormap
graycm = rgb2gray(cm);  graycm = graycm(:,1); 

% do transformation
if isa(a,'double')
  % Make sure A is in the range from 1 to size(cm,1)
  a = max(1,min(a,length(graycm)));  
  I(:) = graycm(a);
else
  graycm = changeclass(class(a),graycm);
  if isa(a,'uint8')
    lut = uint8(zeros(1,256));
  else
    lut = uint16(zeros(1,65536));
  end
  lut(1:length(graycm)) = graycm;
  lut(length(graycm):end) = graycm(end);
  I(:) = uintlut(a,lut);
end

%----------------------------------------------------------------------------------
function [ind,map] = parse_inputs_ind2gray(varargin)
  
checknargin(2,2,nargin,'ind2gray');
% For backward compatability, this function handles an indexed image that is
% logical. This usage will be removed in a future release.
checkinput(varargin{1},{'uint8', 'logical','double', 'uint16'},{'nonempty'}, ...
           'ind2gray','X',1);
ind = varargin{1};
if islogical(ind)
  wid = sprintf('Images:%s:invalidType','ind2gray');
  msg = ['X should be a double, uint8, or uint16 array.  Convert your image to ' ...
         'double using IM2DOUBLE(X,''INDEXED'').'];
  warning(wid,'%s',msg);
  ind = im2double(ind,'indexed');
end

% For backward compatability, this function handles colormaps that are not
% double. This usage will be removed in a future release.
checkinput(varargin{2},{'double','uint8','uint16'},{'nonempty','2d'},'ind2gray','MAP',2);
if ( size(varargin{2},2) ~=3 || size(varargin{2},1) < 1 )
  eid = sprintf('Images:%s:invalidSizeForColormap','ind2gray');
  error(eid,'%s','MAP must be a m x 3 array.');
else
  map = varargin{2};
end
if ~isa(map,'double')
  msg = ['MAP should be a double m x 3 array with values in the range [0,1].'...
         'Convert your map to double using IM2DOUBLE.'];
  warning('ind2gray:notAValidColormap','%s',msg);
  map = im2double(map);
end

%----------------------------------------------------------------------------------
function B = ordfilt2(varargin)
%ORDFILT2 Perform 2-D order-statistic filtering.
%   Copyright 1993-2003 The MathWorks, Inc.  
%   $Revision: 5.22 $  $Date: 2003/03/05 22:29:18 $

[A,order,domain,s,padopt] = ParseInputs_ordfilt2(varargin{:});
domainSize = size(domain);
center = floor((domainSize + 1) / 2);
[r,c] = find(domain);
r = r - center(1);
c = c - center(2);
padSize = max(max(abs(r)), max(abs(c)));
originalSize = size(A);
if (strcmp(padopt, 'zeros'))
    A = padarray(A, padSize * [1 1], 0, 'both');
elseif (strcmp(padopt, 'ones'))
    % padopt of 'ones' is for support of medfilt2; it is undocumented
    A = padarray(A, padSize * [1 1], 1, 'both');
else
    A = padarray(A, padSize * [1 1], 'symmetric', 'both');
end
Ma = size(A,1);
offsets = c*Ma + r;

% make sure that offsets are valid
if ~isreal(offsets) || any(floor(offsets) ~= offsets) || any(~isfinite(offsets))
    %should never get here
    eid = sprintf('Images:%s:internalError', mfilename);
    error(eid,'%s','Internal error: bad OFFSETS.');
end

if isempty(s)  %ORDFILT2(A,ORDER,DOMAIN)
  B = ordf(A, order, offsets, [padSize padSize] + 1, originalSize, domainSize);
else          %ORDFILT2(A,ORDER,DOMAIN,S,PADOPT)
  B = ordf(A, order, offsets, [padSize padSize] + 1, originalSize, domainSize, s);
end

%%% ParseInputs
function [A,order,domain,s,padopt,msg] = ParseInputs_ordfilt2(varargin)
s = [];
padopt = 'zeros';   msg = '';
checknargin(3,5,nargin,'ordfilt2');
A = varargin{1};
order = varargin{2};
domain = varargin{3};
options = {'zeros', 'ones', 'symmetric'};
% padopt of 'ones' is for supporting medfilt2; it is undocumented.

if (nargin == 4)
  if (ischar(varargin{4}))
    padopt = checkstrs(varargin{4},options,'ordfilt2','PADOPT',4);
  else
    s = varargin{4};
  end
elseif (nargin == 5)
  s = varargin{4};
  padopt = checkstrs(varargin{5},options,mfilename,'PADOPT',5);  
end

% make sure that arguments are valid
checkinput(order,'double',{'real','scalar','integer'},'ordfilt2', 'ORDER',2);

if ~isempty(s)
	if (~isa(A, 'double')),	A = double(A);	end
	checkinput(A, 'double', {'2d','real'}, 'ordfilt2', 'A', 1);
	s = s(domain);
	checkinput(s, 'double', 'real', 'ordfilt2', 'S', 4);
else
	checkinput(A, {'numeric','logical'}, {'2d','real'}, 'ordfilt2', 'A', 1);
end

%----------------------------------------------------------------------------------
function b = padarray(varargin)
%PADARRAY Pad an array.
%   Copyright 1993-2003 The MathWorks, Inc.  
%   $Revision: 1.12 $  $Date: 2003/01/17 16:27:48 $

[a, method, padSize, padVal, direction] = ParseInputs_padarray(varargin{:});
if isempty(a)	% treat empty matrix similar for any method
   if strcmp(direction,'both')
      sizeB = size(a) + 2*padSize;
   else
      sizeB = size(a) + padSize;
   end
   b = mkconstarray(class(a), padVal, sizeB);   
else
  switch method
    case 'constant'
        b = ConstantPad(a, padSize, padVal, direction);
    case 'circular'
        b = CircularPad(a, padSize, direction);
    case 'symmetric'
        b = SymmetricPad(a, padSize, direction);
    case 'replicate'
        b = ReplicatePad(a, padSize, direction);
  end      
end

if (islogical(a)),    b = logical(b);    end

%----------------------------------------------------------------------------------
function b = ConstantPad(a, padSize, padVal, direction)
numDims = numel(padSize);

% Form index vectors to subsasgn input array into output array.
% Also compute the size of the output array.
idx   = cell(1,numDims);
sizeB = zeros(1,numDims);
for k = 1:numDims
    M = size(a,k);
    switch direction
        case 'pre'
            idx{k}   = (1:M) + padSize(k);
            sizeB(k) = M + padSize(k);
        case 'post'
            idx{k}   = 1:M;
            sizeB(k) = M + padSize(k);
        case 'both'
            idx{k}   = (1:M) + padSize(k);
            sizeB(k) = M + 2*padSize(k);
    end
end

% Initialize output array with the padding value.  Make sure the
% output array is the same type as the input.
b         = mkconstarray(class(a), padVal, sizeB);
b(idx{:}) = a;

%----------------------------------------------------------------------------------
function b = CircularPad(a, padSize, direction)

numDims = numel(padSize);
% Form index vectors to subsasgn input array into output array.
% Also compute the size of the output array.
idx   = cell(1,numDims);
for k = 1:numDims
  M = size(a,k);
  dimNums = [1:M];
  p = padSize(k);
    
  switch direction
    case 'pre'
       idx{k}   = dimNums(mod([-p:M-1], M) + 1);
    case 'post'
      idx{k}   = dimNums(mod([0:M+p-1], M) + 1);
    case 'both'
      idx{k}   = dimNums(mod([-p:M+p-1], M) + 1);
  end
end
b = a(idx{:});

%----------------------------------------------------------------------------------
function b = SymmetricPad(a, padSize, direction)

numDims = numel(padSize);
% Form index vectors to subsasgn input array into output array.
% Also compute the size of the output array.
idx   = cell(1,numDims);
for k = 1:numDims
  M = size(a,k);
  dimNums = [1:M M:-1:1];
  p = padSize(k);
    
  switch direction
    case 'pre',      idx{k} = dimNums(mod([-p:M-1], 2*M) + 1);
    case 'post',     idx{k} = dimNums(mod([0:M+p-1], 2*M) + 1);
    case 'both',     idx{k} = dimNums(mod([-p:M+p-1], 2*M) + 1);
  end
end
b = a(idx{:});

%----------------------------------------------------------------------------------
function b = ReplicatePad(a, padSize, direction)

numDims = numel(padSize);
% Form index vectors to subsasgn input array into output array.
% Also compute the size of the output array.
idx   = cell(1,numDims);
for k = 1:numDims
  M = size(a,k);
  p = padSize(k);
  onesVector = ones(1,p);
    
  switch direction
    case 'pre'
      idx{k}   = [onesVector 1:M];
    case 'post'
      idx{k}   = [1:M M*onesVector];
    case 'both'
      idx{k}   = [onesVector 1:M M*onesVector];
  end
end
 b = a(idx{:});

%----------------------------------------------------------------
function [a, method, padSize, padVal, direction] = ParseInputs_padarray(varargin)
% default values
method  = 'constant';     padVal = 0;
direction = 'both';

checknargin(2,4,nargin,mfilename);
a = varargin{1};

padSize = varargin{2};
checkinput(padSize, {'double'}, {'real' 'vector' 'nonnan' 'nonnegative' ...
                    'integer'}, mfilename, 'PADSIZE', 2);

% Preprocess the padding size
if (numel(padSize) < ndims(a))
    padSize           = padSize(:);
    padSize(ndims(a)) = 0;
end

if nargin > 2
    firstStringToProcess = 3;
    if ~ischar(varargin{3})        % Third input must be pad value.
        padVal = varargin{3};
        checkinput(padVal, {'numeric' 'logical'}, {'scalar'}, mfilename, 'PADVAL', 3);
        firstStringToProcess = 4;
    end
    
    for k = firstStringToProcess:nargin
        validStrings = {'circular' 'replicate' 'symmetric' 'pre' 'post' 'both'};
        string = checkstrs(varargin{k}, validStrings, mfilename,'METHOD or DIRECTION', k);
        switch string
         case {'circular' 'replicate' 'symmetric'}
          method = string;
         case {'pre' 'post' 'both'}
          direction = string;
         otherwise
          error('Images:padarray:unexpectedError', '%s','Unexpected logic error.')
        end
    end
end
    
% Check the input array type
if strcmp(method,'constant') && ~(isnumeric(a) || islogical(a))
    id = sprintf('Images:%s:badTypeForConstantPadding', mfilename);
    msg1 = sprintf('Function %s expected A (argument 1)',mfilename);
    msg2 = 'to be numeric or logical for constant padding.';
    error(id,'%s\n%s',msg1,msg2);
end

%----------------------------------------------------------------------------------
function BW = poly2mask(x,y,M,N)
%POLY2MASK Convert region polygon to region mask.
%   Copyright 1993-2003 The MathWorks, Inc.
%   $Revision: 1.1 $  $Date: 2003/01/17 16:27:49 $

checknargin(4,4,nargin,mfilename);
checkinput(x,{'double'},{'real','vector','finite'},mfilename,'X',1);
checkinput(y,{'double'},{'real','vector','finite'},mfilename,'Y',2);
checkinput(M,{'double'},{'real','integer','nonnegative'},mfilename,'M',3);
checkinput(N,{'double'},{'real','integer','nonnegative'},mfilename,'N',4);
if length(x) ~= length(y)
    error('Images:poly2mask:vectorSizeMismatch','%s',...
          'Function POLY2MASK expected its first two inputs, X and Y, to be vectors with the same length.');
end

if isempty(x)
    BW = false(M,N);
    return;
end

if x(end) ~= x(1),    x(end+1) = x(1);  end
if y(end) ~= y(1),    y(end+1) = y(1);  end

[xe,ye] = poly2edgelist(x,y);
BW = edgelist2mask(M,N,xe,ye);

%----------------------------------------------------------------------------------
function I = rgb2gray(varargin)
%RGB2GRAY Convert RGB image or colormap to grayscale.
%   Copyright 1993-2003 The MathWorks, Inc.  
%   $Revision: 5.21 $  $Date: 2003/01/17 16:27:50 $

X = parse_inputs_rgb2gray(varargin{:});
origSize = size(X);

% Determine if input includes a 3-D array 
threeD = (ndims(X)==3);

% Calculate transformation matrix
T = inv([1.0 0.956 0.621; 1.0 -0.272 -0.647; 1.0 -1.106 1.703]);
coef = T(1,:)';

if threeD  %RGB
	% Shape input matrix so that it is a n x 3 array and initialize output matrix  
	X = reshape(X(:),origSize(1)*origSize(2),3);
	sizeOutput = [origSize(1), origSize(2)];
	
	% Do transformation
	if isa(X, 'double')
	    I = X*coef;
	    I = min(max(I,0),1);
	else
	    %uint8 or uint16
	    I = imlincomb(coef(1),X(:,1),coef(2),X(:,2),coef(3),X(:,3),class(X));
	end
	%Make sure that the output matrix has the right size
	I = reshape(I,sizeOutput);   
else
  % For backward compatability, this function handles uint8 and uint16
  % colormaps. This usage will be removed in a future release.
  if ~isa(X,'double')
    I = imlincomb(coef(1),X(:,1),coef(2),X(:,2),coef(3),X(:,3), class(X));
  else
    I = X*coef;
    I = min(max(I,0),1);
  end
  I = im2double(I);
  I = [I,I,I];
end

%---------------------------------------------------------------
function X = parse_inputs_rgb2gray(varargin)

checknargin(1,3,nargin,mfilename);
switch nargin
 case 1
  if ndims(varargin{1})==2
    if (size(varargin{1},2) ~=3 || size(varargin{1},1) < 1)
      eid = sprintf('Images:%s:invalidSizeForColormap',mfilename);
      error(eid,'%s','MAP must be a m x 3 array.');
    end
    if ~isa(varargin{1},'double')
      wid = sprintf('Images:%s:notAValidColormap',mfilename);
      msg = ['MAP should be a double m x 3 array with values in the range [0,1].'...
             'Convert your map to double using IM2DOUBLE.'];
      warning(wid,'%s',msg);
    end
  elseif (ndims(varargin{1})==3)
    if ((size(varargin{1},3) ~= 3))
      eid = sprintf('Images:%s:invalidInputSize',mfilename);
      error(eid,'%s','RGB must be a m x n x 3 array.');
    end
  else
    msg = ['RGB2GRAY only accepts a 2-D input for MAP or a 3-D input for RGB.'];
    error(sprintf('Images:%s:invalidInputSize',mfilename),'%s',msg);
  end
  X = varargin{1};  
 case 2
  eid = sprintf('Images:%s:invalidNumberOfArguments',mfilename);
  error(eid,'%s','Two input arguments are not allowed.');
 case 3
  msg = ['RGB2GRAY(R,G,B) is an obsolete syntax.',...
         'Use a three-dimensional array to represent RGB image.'];
  warning(sprintf('Images:%s:obsoleteSyntax',mfilename),'%s',msg);
  if ( any(size(varargin{1})~=size(varargin{2})) || any(size(varargin{1})~=size(varargin{3})) )
    eid = sprintf('Images:%s:inputsDontHaveSameSize',mfilename);
    error(eid,'%s','R, G, and B must all be the same size.');
  end
  X = cat(3,varargin{1},varargin{2},varargin{3});
end

%no logical arrays
if islogical(X)
  eid = sprintf('Images:%s:invalidType',mfilename);
  error(eid,'%s','RGB2GRAY does not accept logical arrays as inputs.');
end

%----------------------------------------------------------------------------------
function varargout = rgb2ntsc(varargin)
%RGB2NTSC Convert RGB values to NTSC colorspace.
%   Copyright 1993-2003 The MathWorks, Inc.  
%   $Revision: 5.18 $  $Date: 2003/01/17 16:27:50 $

A = parse_inputs_rgb2ntsc(varargin{:});
T = [1.0 0.956 0.621; 1.0 -0.272 -0.647; 1.0 -1.106 1.703].';
[so(1), so(2), thirdD] = size(A);
if thirdD == 1	% A is RGBMAP, M-by-3 colormap
  A = A/T;
else % A is truecolor image RBG
  A = reshape(reshape(A,so(1)*so(2),thirdD)/T,so(1),so(2),thirdD);
end

% Output
if nargout < 2				% YIQMAP = RGB2NTSC(RGBMAP)
  varargout{1} = A;
elseif  nargout == 3		% [Y,I,Q] = RGB2NTSC(RGB) OBSOLETE
  if thirdD == 1
    eid = sprintf('Images:%s:wrongNumberOfOutputArguments', mfilename);
    error(eid,'%s','YIQMAP = RGB2NTSC(RGBMAP) must have one output argument.');
  else
    msg = ['[Y,I,Q] = RGB2NTSC(RGB) is an obsolete output syntax. ',...
           'Use a three-dimensional array to represent YIQ image.'];
    warning(sprintf('Images:%s:obsoleteOutputSyntax', mfilename),'%s',msg);
  end
  for k = 1:3
    varargout{k} = A(:,:,k);
  end
else 
  eid = sprintf('Images:%s:wrongNumberOfOutputArguments', mfilename);
  error(eid,'%s',sprintf('RGB2NTSC cannot return %d output arguments.', nargout));
end

%--------------------------------------------------------
function A = parse_inputs_rgb2ntsc(varargin)
checknargin(1,3,nargin,mfilename);
switch nargin
 case 1  % rgb2ntsc(RGB) or rgb2ntsc(RGBMAP)
  A = varargin{1};
 case 2
  eid = sprintf('Images:%s:wrongNumberOfInputs',mfilename);
  error(eid,'%s','Two input arguments are not allowed');
 case 3  % rgb2ntsc(R,G,B) OBSOLETE
  msg = ['RGB2NTSC(R,G,B) is an obsolete syntax. ',...
         'Use a three-dimensional array to represent RGB image.'];
  warning(sprintf('Images:%s:obsoleteInputSyntax',mfilename),'%s',msg);
  %check class
  if ( strcmp(class(varargin{1}),class(varargin{2})) == 0 || ...
        strcmp(class(varargin{1}),class(varargin{3})) == 0 )
    eid = sprintf('Images:%s:inputsHaveDifferentClass',mfilename);
    error(eid,'%s','R, G, and B arrays must have the same class.');
  end
 
  if isequal(size(varargin{1}),size(varargin{2}),size(varargin{3}))
    A = cat(3,varargin{1},varargin{2},varargin{3});
  else
    eid = sprintf('Images:%s:unequalSizes',mfilename);
    error(eid,'%s','R, G, and B arrays must have the same size.');
  end
end

%no logical
if islogical(A)
  eid = sprintf('Images:%s:invalidType',mfilename);
  error(eid,'%s','A truecolor image cannot be logical.');
end

% Check validity of the input parameters 
if ndims(A)==2 
  % Check colormap 
  id = sprintf('Images:%s:invalidColormap',mfilename);
  if ( size(A,2)~=3 || size(A,1) < 1 ) 
    error(id,'%s','RGBMAP must be an M-by-3 array.');
  end
  if ~isa(A,'double')
    msg = ['MAP should be a double m x 3 array with values in the range [0,1].'...
           'Convert your map to double using IM2DOUBLE.'];
    warning(id,'%s',msg);
    A = im2double(A);
  end
elseif ndims(A)==3
  % Check RGB
  if size(A,3)~=3
    eid = sprintf('Images:%s:invalidTruecolorImage',mfilename);
    error(eid,'%s','RGB image must be an M-by-N-by-3 array.');
  end
  A = im2double(A);
else
  eid = sprintf('Images:%s:invalidSize',mfilename);
  msg = 'RGB2NTSC only accepts a 2-D input for RGBMAP and a 3-D input for RGB.';
  error(eid,'%s',msg);
end

% A has to be double because YIQ colorspace can contain negative values.
if ~isa(A, 'double'),   A = im2double(A);   end

%----------------------------------------------------------------------------------
function J = roifilt2(varargin)
%ROIFILT2 Filters a region of interest in an image.
%   Copyright 1993-2003 The MathWorks, Inc.
%   $Revision: 1.22 $  $Date: 2003/01/27 20:16:05 $

[H, J, BW, params, fcnflag] = parse_inputs_roifilt2(varargin{:});
% Assigning default value to minrow.  In the case when J = ROIFILT2(H, I,
% BW), minrow will be set to at least one.
minrow = 0;

if (fcnflag == 1 )    %case when J = ROIFILT2(I, BW, 'fun', P1,...)
    filtI = feval(H, J, params{:});
    if ~isa(J, class(filtI))
        J = feval(class(filtI), J);  
    end
else 
    % case when J = ROIFILT2(H, I, BW).  Determine rectangle that encloses
    % the non-zero elements in BW.  The rectangle vector is chosen so that
    % no non-zero element in J is considered to be a boundary pixel by
    % imfilter.  In other words, the row and column padding should be equal
    % to the row size and column size of H, respectively.  Also, rectangle
    % cannot be bigger than size of original image.
    [row, col] = find(BW==1);
    colpad = ceil(size(H, 2) / 2);
    rowpad = ceil(size(H, 1) / 2);
    mincol = max(1, min(col(:)) - colpad);
    minrow = max(1, min(row(:)) - rowpad); 
    maxcol = min(size(J, 2), max(col(:)) + colpad);
    maxrow = min(size(J, 1), max(row(:)) + rowpad);
    
    % perform filtering on y that is cropped to the rectangle. 
    I = J;
    J = J(minrow:maxrow, mincol:maxcol);
    BW = BW(minrow:maxrow, mincol:maxcol);
    filtI = imfilter(J, H);
end

if ~isequal(size(filtI), size(J)) 
    error('The filtered image is not the same size as the original image.');
end

J(BW) = filtI(BW);
if (minrow ~= 0)    % case when J = ROIFILT2(H ,I, BW).
    I(minrow: maxrow, mincol: maxcol) = J;
    J = I;
end

%------------------------------------------------------------------------
function [filter, I, mask, param, flag] = parse_inputs_roifilt2(varargin)
error(nargchk(3, inf, nargin));
% initialize elements
flag = 0;
param = varargin(4:end);

% [fun,fcnchk_msg] = fcnchk(varargin{3}, length(param));
% if isempty(fcnchk_msg)    % J = ROIFILT2(I,BW,'fun',P1,...)
%     I = varargin{1};    mask = varargin{2};
%     filter = fun;    flag = 1;
% else                      % J = ROIFILT2(H, I, BW)
%     filter = varargin{1};    I = varargin{2};    mask = varargin{3};
% end
% The above  does not compile and seems never being used in Mirone
filter = varargin{1};    I = varargin{2};    mask = varargin{3};

if (ndims(I) ~=  2),    error('I must be a two dimensional array.'); end
if ~islogical(mask),    mask = mask ~= 0;   end
if ~isequal(size(mask),size(I))
    error('Image and binary mask must be the same size.');
end

%----------------------------------------------------------------------------------
function varargout = roipoly_j(varargin)
%ROIPOLY Select polygonal region of interest.
%   Copyright 1993-2003 The MathWorks, Inc.  
%   $Revision: 5.29 $  $Date: 2003/01/27 20:16:05 $
% Hacked version for not using options that use imshow

[xdata,ydata,num_rows,num_cols,xi,yi] = parse_inputs_roipoly_j(varargin{:});
if length(xi)~=length(yi)
  error('XI and YI must be the same length.'); 
end

% Make sure polygon is closed.
if (~isempty(xi))
    if ((xi(1) ~= xi(end)) || (yi(1) ~= yi(end)))
        xi = [xi;xi(1)]; yi = [yi;yi(1)];
    end
end
% Transform xi,yi into pixel coordinates.
roix = axes2pix(num_cols, xdata, xi);
roiy = axes2pix(num_rows, ydata, yi);

d = poly2mask(roix, roiy, num_rows, num_cols);

switch nargout
case 1
    varargout{1} = d;
case 2
    varargout{1} = d;    varargout{2} = xi;
case 3
    varargout{1} = d;    varargout{2} = xi;    varargout{3} = yi;
case 4
    varargout{1} = xdata;    varargout{2} = ydata;
    varargout{3} = d;        varargout{4} = xi;
case 5
    varargout{1} = xdata;    varargout{2} = ydata;
    varargout{3} = d;        varargout{4} = xi;
    varargout{5} = yi;
otherwise
    error('Too many output arguments');
end

%-----------------------------------------------------------------
function [x,y,nrows,ncols,xi,yi] = parse_inputs_roipoly_j(varargin)

switch nargin
% Hacked for not using options (nargin == 0-2) that use imshow
case 3    % SYNTAX: roipoly_j(A,xi,yi)
    a = varargin{1};
    nrows = size(a,1);    ncols = size(a,2);
    xi = varargin{2}(:);  yi = varargin{3}(:);
    x = [1 ncols]; y = [1 nrows];
case 4    % SYNTAX: roipoly_j(m,n,xi,yi)
    nrows = varargin{1};    ncols = varargin{2};
    xi = varargin{3}(:);    yi = varargin{4}(:);
    x = [1 ncols]; y = [1 nrows];
case 5    % SYNTAX: roipoly_j(x,y,A,xi,yi)
    x = varargin{1};        y = varargin{2}; 
    a = varargin{3};
    xi = varargin{4}(:);    yi = varargin{5}(:);
    nrows = size(a,1);      ncols = size(a,2);
    x = [x(1) x(numel(x))];
    y = [y(1) y(numel(y))];
case 6    % SYNTAX: roipoly_j(x,y,m,n,xi,yi)
    x = varargin{1};        y = varargin{2}; 
    nrows = varargin{3};    ncols = varargin{4};
    xi = varargin{5}(:);    yi = varargin{6}(:);
    x = [x(1) x(numel(x))];
    y = [y(1) y(numel(y))];
	otherwise
    error('Invalid input arguments.');
end

if isa(xi, 'uint8'), xi=double(xi); end
if isa(yi, 'uint8'), yi=double(yi); end
if isa(x, 'uint8'), x=double(x); end
if isa(y, 'uint8'), y=double(y); end
if isa(nrows, 'uint8'), nrows=double(nrows); end
if isa(ncols, 'uint8'), ncols=double(ncols); end

%----------------------------------------------------------------------------------
function [X,map] = rgb2ind(varargin)
%RGB2IND Convert RGB image to indexed image.

%   Copyright 1993-2003 The MathWorks, Inc.  
%   $Revision: 5.24 $  $Date: 2003/01/29 19:57:29 $

[RGB,m,dith] = rgb2ind_parse_inputs(varargin{:});

[so(1), so(2), thirdD] = size(RGB);

% Converts depending on what is m:
if (isempty(m))      % Convert RGB image to an indexed image.
    X = reshape([1:so(1)*so(2)],so(1),so(2));
    if so(1)*so(2) <= 256
        X = uint8(X-1);
    elseif so(1)*so(2) <= 65536
        X = uint16(X-1);
    end
    map = im2double(reshape(RGB,so(1)*so(2),3));

elseif (length(m) == 1)    % TOL or N is given
    RGB = im2uint8(RGB);

    if (m < 1)     % tol is given. Use uniform quantization
        max_colors = 65536;
        max_N = floor(max_colors^(1/3)) - 1;
        N = round(1/m);
        if (N > max_N)
            N = max_N;
            warning(sprintf('Too many colors; increasing tolerance to %g', 1/N));
        end
        
        [x,y,z] = meshgrid([0:N]/N);
        map = [x(:),y(:),z(:)];

        if dith(1) == 'n'
            RGB = round(im2double(RGB)*N);
            X = RGB(:,:,3)*((N+1)^2)+RGB(:,:,1)*(N+1)+RGB(:,:,2)+1;
        else
            X = dither(RGB,map);
        end
        [X,map] = cmunique(X,map);

    else    % N is given. Use variance minimization quantization
        [map,X] = cq(RGB,m);
        map = double(map) / 255;
        if (dith(1) == 'd') && (size(map,1) > 1)
            % Use standalone dither if map is an approximation.
            X = dither(RGB,map);
        end
    end

else        % MAP is given
    RGB = im2uint8(RGB);
    map = m;
    if (dith(1) == 'n')          % 'nodither'
        X = dither(RGB,map,5,4); % Use dither to do inverse colormap mapping.
    else
        X = dither(RGB,map);
    end
end

if isa(map, 'uint8')    % Make sure that the colormap is doubles
    map = double(map)/255;
elseif isa(map, 'uint16')
    map = double(map)/65535;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [RGB,m,dith] = rgb2ind_parse_inputs(varargin)
% Outputs:  RGB     image
%           m       colormap
%           dith    dithering option
% Defaults:
dith = 'dither';    m = [];

error(nargchk(1,5,nargin));
switch nargin
	case 1               % rgb2ind(RGB)
        RGB = varargin{1};
        warning(sprintf('%s\n%s', 'RGB2IND(RGB) is an obsolete syntax.', ...
                      'Specify number of colors, tolerance, or colormap.'));
	case 2               % rgb2ind(RGB,x) where x = MAP | N | TOL
        RGB = varargin{1};
        m = varargin{2};
    case 3               % rgb2ind(R,G,B) OBSOLETE
        RGB = varargin{1};  %              where x = MAP | N | TOL
        m = varargin{2};
        dith = varargin{3};
	otherwise
        error('Invalid input arguments.');
end

% Check validity of the input parameters
if (ndims(RGB)==3) && (size(RGB,3) ~= 3) || (ndims(RGB) > 3)
    error('RGB image has to be M-by-N-by-3 array.');
end

if any(m(:) < 0)
    error('Colormap, Number of colors, or Tolerance have to be non-negative.');
elseif ((length(m)==1) && (m~=round(m)) && (m>1))
    error('Number of colors in the colormap has to be a non-negative integer.');
elseif ((length(m) > 1) && (ndims(m) == 2) && (size(m,1) == 1) && (size(m,2) ~= 3))    % MAP
    error('Input colormap has to be a 2D array with at least 2 rows and exactly 3 columns.');
elseif (length(m)>1) && (max(m(:))>1)
    error('All colormap intensities must be between 0 and 1.');
end

if ischar(dith)	% dither option
	strings = {'dither','nodither'};
	idx = strmatch(lower(dith),strings);
	if isempty(idx)
	    error(sprintf('Unknown dither option: %s',dith));
	elseif (length(idx) > 1)
	    error(sprintf('Ambiguous dither option: %s',dith));
	else
	    dith = strings{idx};
	end  
else
    error(sprintf('Dither option has to be a string.'));  
end

%----------------------------------------------------------------------------------
function lowhigh = stretchlim_j(varargin)
%STRETCHLIM Find limits to contrast stretch an image.
%   Copyright 1993-2003 The MathWorks, Inc.  
%   $Revision: 1.9 $ $Date: 2003/01/27 20:16:05 $

[img,tol] = ParseInputs_stretchlim_j(varargin{:});
nbins = 256;
tol_low = tol(1);
tol_high = tol(2);
[m,n,p] = size(img);

for i = 1:p                          % Find limits, one plane at a time
    N = imhist_j(img(:,:,i),nbins);
    cdf = cumsum(N)/sum(N);
    ilow = min(find(cdf > tol_low));
    ihigh = min(find(cdf >= tol_high));   
    ilowhigh(:,i) = [ilow;ihigh];

    if ilow==ihigh               
      ilowhigh(:,i) = [1; nbins]; % limits are same, use default   
    end
end   
lowhigh = (ilowhigh - 1)/(nbins-1);  % convert to range [0 1]

%-----------------------------------------------------
function [img,tol] = ParseInputs_stretchlim_j(varargin)

tol = [.01 .99];  % default 
checknargin(1, 2, nargin, mfilename);

img = varargin{1};
if nargin > 1
    tol = varargin{2};
    switch numel(tol)
      case 1
        tol(2) = 1 - tol;
      case 2
        if (tol(1) >= tol(2))
            msgId = 'Images:stretchlim:invalidTolOrder';
            error(msgId,'%s','TOL(1) must be less than TOL(2).');
        end
      otherwise
        msgId = 'Images:stretchlim:invalidTolSize';
        error(msgId,'%s','TOL must have 1 or 2 elements.');
    end
end

if ( any(tol < 0) || any(tol>1) || any(isnan(tol)) )
    msgId = 'Images:stretchlim:tolOutOfRange';
    error(msgId,'%s','TOL must be in the range [0 1].');
end

checkinput(img, {'uint8', 'uint16', 'double'}, {'real', 'nonsparse'}, mfilename, 'I or RGB', 1);

if (ndims(img) > 3) 
    msgId = 'Images:stretchlim:dimTooHigh';
    error(msgId,'%s','STRETCHLIM only supports individual images.');
end

%----------------------------------------------------------------------------------
function B = uintlut(varargin)
%UINTLUT computes new values of A based on lookup table LUT.
% Used to be uintlutc
B = intlutc(varargin{:});

%----------------------------------------------------------------------------------
function image = changeclass(class, varargin)
%CHANGECLASS will change the storage class of an image.
%   Copyright 1993-2003 The MathWorks, Inc.  
%   $Revision: 1.9 $  $Date: 2003/01/17 16:28:15 $
switch class
	case 'uint8',	image = im2uint8(varargin{:});
	case 'uint16',	image = im2uint16(varargin{:});
	case 'double',	image = im2double(varargin{:});
	otherwise
		error('Unsupported IPT data class.');
end

%----------------------------------------------------------------------------------
function checkinput(A, classes, attributes, function_name, variable_name, argument_position)
%CHECKINPUT Check validity of array.
%   Copyright 1993-2003 The MathWorks, Inc.
%   $Revision: 1.8 $  $Date: 2003/03/03 14:21:21 $

% Input arguments are not checked for validity.
check_classes(A, classes, function_name, variable_name, argument_position);
check_attributes(A, attributes, function_name, variable_name, argument_position);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tf = is_numeric(A)
numeric_classes = {'double' 'uint8' 'uint16' 'uint32' 'int8' 'int16' 'int32' 'single'};

tf = false;
for p = 1:length(numeric_classes)
    if isa(A, numeric_classes{p})
        tf = true;		break
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = expand_numeric(in)
% Converts the string 'numeric' to the equivalent cell array containing the
% names of the numeric classes.

out = in(:);
idx = strmatch('numeric', out, 'exact');
if (length(idx) == 1)
    out(idx) = [];
    numeric_classes = {'uint8', 'int8', 'uint16', 'int16', 'uint32', 'int32', 'single', 'double'}';
    out = [out; numeric_classes];
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function check_classes(A, classes, function_name, variable_name, argument_position)
if isempty(classes),    return;     end

if ischar(classes)
    if isempty(classes)        % Work around bug in strread.
        classes = {};
    else
        classes = strread(classes, '%s');
    end
end

is_valid_type = false;
for k = 1:length(classes)
    if strcmp(classes{k}, 'numeric') && is_numeric(A)
        is_valid_type = true;
        break;
    else        
        if isa(A, classes{k})
            is_valid_type = true;
            break;
        end
    end
end

if ~is_valid_type
    messageId = sprintf('Images:%s:%s', function_name, 'invalidType');
    classes = expand_numeric(classes);
    validTypes = '';
    for k = 1:length(classes)
        validTypes = [validTypes, classes{k}, ', '];
    end
    validTypes(end-1:end) = [];
    message1 = sprintf('Function %s expected its %s input argument, %s,', ...
                       upper(function_name), ...
                       num2ordinal(argument_position), variable_name);
    message2 = 'to be one of these types:';
    message3 = sprintf('Instead its type was %s.', class(A));
    error(messageId, '%s\n%s\n\n  %s\n\n%s', message1, message2, validTypes, message3);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function check_attributes(A, attributes, function_name, variable_name, argument_position)
if ischar(attributes)
    if isempty(attributes)        % Work around bug in strread.
        attributes = {};
    else
        attributes = strread(attributes, '%s');
    end
end

table = init_table;
for k = 1:length(attributes)
    if strcmp(attributes{k}, '2d')
        tableEntry = table.twod;
    else
        tableEntry = table.(attributes{k});
    end
    if ~feval(tableEntry.checkFunction, A)
        messageId = sprintf('Images:%s:%s', function_name, tableEntry.mnemonic);
        message1 = sprintf('Function %s expected its %s input argument, %s,', ...
                           upper(function_name), num2ordinal(argument_position), variable_name);
        message2 = sprintf('to be %s.', tableEntry.endOfMessage);
        error(messageId, '%s\n%s', message1, message2);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tf = check_real(A)
try		tf = isreal(A);
catch,	tf = false;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tf = check_even(A)
try		tf = ~any(rem(double(A(:)),2));
catch,	tf = false;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tf = check_vector(A)
try		tf = (ndims(A) == 2) && (any(size(A) == 1) || all(size(A) == 0));
catch,	tf = false;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tf = check_row(A)
try		tf = (ndims(A) == 2) && ((size(A,1) == 1) || isequal(size(A), [0 0]));
catch,	tf = false;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tf = check_column(A)
try		tf = (ndims(A) == 2) && ((size(A,2) == 1) || isequal(size(A), [0 0]));
catch,	tf = false;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tf = check_scalar(A)
try		tf = all(size(A) == 1);
catch,	tf = false;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tf = check_2d(A)
try		tf = ndims(A) == 2;
catch,	tf = false;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tf = check_nonsparse(A)
try		tf = ~issparse(A);
catch,	tf = false;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tf = check_nonempty(A)
try		tf = ~isempty(A);
catch,	tf = false;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tf = check_integer(A)
try
    A = A(:);
    switch class(A)
      case {'double','single'}
        % In MATLAB 6.5 and earlier, floor(A) isn't support for single
        % A, so convert to double.
        A = double(A);
        tf = all(floor(A) == A) & all(isfinite(A));

      case {'uint8','int8','uint16','int16','uint32','int32','logical'}
        tf = true;
      otherwise
        tf = false;
    end
catch
    tf = false;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tf = check_nonnegative(A)
try     tf = all(A(:) >= 0);
catch,  tf = false;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tf = check_positive(A)
try     tf = all(A(:) > 0);
catch,  tf = false;     end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tf = check_nonnan(A)
try
    tf = ~any(isnan(A(:)));
catch
    % if isnan isn't defined for the class of A, then we'll end up here.  If isnan isn't
    % defined then we'll assume that A can't contain NaNs.
    tf = true;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tf = check_finite(A)
try
    tf = all(isfinite(A(:)));
catch
    % if isfinite isn't defined for the class of A,
    % then we'll end up here.  If isfinite isn't
    % defined then we'll assume that A is finite.
    tf = true;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tf = check_nonzero(A)
try     tf = ~all(A(:) == 0);
catch,  tf = false;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = init_table
persistent table
if isempty(table)
    table.real.checkFunction        = @check_real;
    table.real.mnemonic             = 'expectedReal';
    table.real.endOfMessage         = 'real';
    
    table.vector.checkFunction      = @check_vector;
    table.vector.mnemonic           = 'expectedVector';
    table.vector.endOfMessage       = 'a vector';
    
    table.row.checkFunction         = @check_row;
    table.row.mnemonic              = 'expectedRow';
    table.row.endOfMessage          = 'a row vector';
    
    table.column.checkFunction      = @check_column;
    table.column.mnemonic           = 'expectedColumn';
    table.column.endOfMessage       = 'a column vector';
    
    table.scalar.checkFunction      = @check_scalar;
    table.scalar.mnemonic           = 'expectedScalar';
    table.scalar.endOfMessage       = 'a scalar';
    
    table.twod.checkFunction        = @check_2d;
    table.twod.mnemonic             = 'expected2D';
    table.twod.endOfMessage         = 'two-dimensional';
    
    table.nonsparse.checkFunction   = @check_nonsparse;
    table.nonsparse.mnemonic        = 'expectedNonsparse';
    table.nonsparse.endOfMessage    = 'nonsparse';
    
    table.nonempty.checkFunction    = @check_nonempty;
    table.nonempty.mnemonic         = 'expectedNonempty';
    table.nonempty.endOfMessage     = 'nonempty';
    
    table.integer.checkFunction     = @check_integer;
    table.integer.mnemonic          = 'expectedInteger';
    table.integer.endOfMessage      = 'integer-valued';
        
    table.nonnegative.checkFunction = @check_nonnegative;
    table.nonnegative.mnemonic      = 'expectedNonnegative';
    table.nonnegative.endOfMessage  = 'nonnegative';
    
    table.positive.checkFunction    = @check_positive;
    table.positive.mnemonic         = 'expectedPositive';
    table.positive.endOfMessage     = 'positive';
    
    table.nonnan.checkFunction      = @check_nonnan;
    table.nonnan.mnemonic           = 'expectedNonNaN';
    table.nonnan.endOfMessage       = 'non-NaN';
    
    table.finite.checkFunction      = @check_finite;
    table.finite.mnemonic           = 'expectedFinite';
    table.finite.endOfMessage       = 'finite';
    
    table.nonzero.checkFunction     = @check_nonzero;
    table.nonzero.mnemonic          = 'expectedNonZero';
    table.nonzero.endOfMessage      = 'non-zero';
    
    table.even.checkFunction        = @check_even;
    table.even.mnemonic             = 'expectedEven';
    table.even.endOfMessage         = 'even'; 
end
out = table;

%----------------------------------------------------------------------------------
function out = checkstrs(in, valid_strings, function_name, variable_name, argument_position)
%CHECKSTRS Check validity of option string.
%   Copyright 1993-2003 The MathWorks, Inc.
%   $Revision: 1.4 $  $Date: 2003/01/17 16:28:16 $

% Except for IN, input arguments are not checked for validity.
checkinput(in, 'char', 'row', function_name, variable_name, argument_position);
start = regexpi(valid_strings, ['^' in]);
if ~iscell(start),	start = {start};   end
matches = ~cellfun('isempty',start);
idx = find(matches);
num_matches = length(idx);

if num_matches == 1
    out = valid_strings{idx};
else
    out = substringMatch(valid_strings(idx));
    if isempty(out)
        % Convert valid_strings to a single string containing a space-separated list
        % of valid strings.
        list = '';
        for k = 1:length(valid_strings)
            list = [list ', ' valid_strings{k}];
        end
        list(1:2) = [];
        msg1 = sprintf('Function %s expected its %s input argument, %s,', ...
                       upper(function_name), num2ordinal(argument_position),variable_name);
        msg2 = 'to match one of these strings:';
        
        if num_matches == 0
            msg3 = sprintf('The input, ''%s'', did not match any of the valid strings.', in);
            id = sprintf('Images:%s:unrecognizedStringChoice', function_name);
        else
            msg3 = sprintf('The input, ''%s'', matched more than one valid string.', in);
            id = sprintf('Images:%s:ambiguousStringChoice', function_name);
        end
        error(id,'%s\n%s\n\n  %s\n\n%s', msg1, msg2, list, msg3);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function str = substringMatch(strings)
%   STR = substringMatch(STRINGS) looks at STRINGS (a cell array of
%   strings) to see whether the shortest string is a proper substring of
%   all the other strings.  If it is, then substringMatch returns the
%   shortest string; otherwise, it returns the empty string.
if isempty(strings)
    str = '';
else
    len = cellfun('prodofsize',strings);
    [tmp,sortIdx] = sort(len);
    strings = strings(sortIdx);
    start = regexpi(strings(2:end), ['^' strings{1}]);
    if isempty(start) || (iscell(start) && any(cellfun('isempty',start)))
        str = '';
    else
        str = strings{1};
    end
end

%----------------------------------------------------------------------------------
function checkmap(map, function_name, variable_name, argument_position)
%CHECKMAP Check colormap for validity.
%   Copyright 1993-2003 The MathWorks, Inc.
%   $Revision: 1.1 $  $Date: 2003/01/17 16:28:15 $

if ~isa(map,'double') || isempty(map) || ...
      (ndims(map) ~= 2) || (size(map,2) ~= 3) || ...
      issparse(map)
    msgId = sprintf('Images:%s:badMapMatrix',function_name);
    error(msgId,'Function %s expected its %s input argument, %s, to be a valid colormap.\nValid colormaps must be nonempty, double, 2-D matrices with 3 columns.', ...
          upper(function_name), num2ordinal(argument_position), variable_name);
end

if any(map(:) < 0) || any(map(:) > 1)
    msgId = sprintf('Images:%s:badMapValues',function_name);
    error(msgId,'Function %s expected its %s input argument, %s, to be a valid colormap.\nValid colormaps cannot have values outside the range [0,1].', ...
          upper(function_name), num2ordinal(argument_position), variable_name);
end

%----------------------------------------------------------------------------------
function checknargin(low, high, numInputs, function_name)
%CHECKNARGIN Check number of input arguments.
%   Copyright 1993-2003 The MathWorks, Inc.
%   $Revision: 1.3 $  $Date: 2003/01/17 16:28:15 $

% Input arguments are not checked for validity.
if numInputs < low
  if low == 1
    msg1 = sprintf('Function %s expected at least 1 input argument', upper(function_name));
  else
    msg1 = sprintf('Function %s expected at least %d input arguments', upper(function_name), low);
  end
  
  if numInputs == 1
    msg2 = 'but was called instead with 1 input argument.';
  else
    msg2 = sprintf('but was called instead with %d input arguments.', numInputs);
  end
  
  error(sprintf('Images:%s:tooFewInputs', function_name), '%s\n%s', msg1, msg2);
elseif numInputs > high
  if high == 1
    msg1 = sprintf('Function %s expected at most 1 input argument', upper(function_name));
  else
    msg1 = sprintf('Function %s expected at most %d input arguments', upper(function_name), high);
  end
  
  if numInputs == 1
    msg2 = 'but was called instead with 1 input argument.';
  else
    msg2 = sprintf('but was called instead with %d input arguments.', numInputs);
  end
  error(sprintf('Images:%s:tooManyInputs', function_name), '%s\n%s', msg1, msg2);
end

%----------------------------------------------------------------------------------
function BW = edgelist2mask(M,N,xe,ye,scale)
%EDGELIST2MASK Convert horizontal edge list to a region mask.
%   Copyright 1993-2003 The MathWorks, Inc.
%   $Revision: 1.1 $  $Date: 2003/01/17 16:28:30 $

if (nargin < 5)    scale = 5;   end
shift = (scale - 1)/2;

% Scale x values, throwing away edgelist points that aren't on a pixel's
% center column. 
xe = (xe+shift)/5;  idx = xe == floor(xe);
xe = xe(idx);       ye = ye(idx);

% Scale y values.
ye = ceil((ye + shift)/scale);

% Throw away horizontal edges that are too far left, too far right, or below the image.
bad_indices = find((xe < 1) | (xe > N) | (ye > M));
xe(bad_indices) = [];       ye(bad_indices) = [];

% Treat horizontal edges above the top of the image as they are along the upper edge.
ye = max(1,ye);
S = sparse(ye,xe,1,M,N);
BW = parityscan(full(S));

%----------------------------------------------------------------------------------
function [x,y] = intline(x1, x2, y1, y2)
%INTLINE Integer-coordinate line drawing algorithm.
%   Copyright 1993-2003 The MathWorks, Inc.  
%   $Revision: 5.12 $  $Date: 2003/01/17 16:28:32 $

dx = abs(x2 - x1);      dy = abs(y2 - y1);

% Check for degenerate case.
if ((dx == 0) && (dy == 0)),  x = x1;  y = y1;  return;   end

flip = 0;
if (dx >= dy)
  if (x1 > x2)    % Always "draw" from left to right.
    t = x1; x1 = x2; x2 = t;
    t = y1; y1 = y2; y2 = t;
    flip = 1;
  end
  m = (y2 - y1)/(x2 - x1);
  x = (x1:x2).';
  y = round(y1 + m*(x - x1));
else
  if (y1 > y2)    % Always "draw" from bottom to top.
    t = x1; x1 = x2; x2 = t;
    t = y1; y1 = y2; y2 = t;
    flip = 1;
  end
  m = (x2 - x1)/(y2 - y1);
  y = (y1:y2).';
  x = round(x1 + m*(y - y1));
end
if (flip)   x = flipud(x);  y = flipud(y);  end

%----------------------------------------------------------------------------------
function out = mkconstarray(class, value, size)
%MKCONSTARRAY creates a constant array of a specified numeric class.
%   Copyright 1993-2003 The MathWorks, Inc.  
out = repmat(feval(class, value), size);

%----------------------------------------------------------------------------------
function b = medfilt2(varargin)
%MEDFILT2 Perform 2-D median filtering.
%   Copyright 1993-2003 The MathWorks, Inc.  
%   $Revision: 5.22 $  $Date: 2003/01/27 20:16:04 $

[a, mn, padopt] = parse_inputs_medfilt2(varargin{:});
domain = ones(mn);
if (rem(prod(mn), 2) == 1)
    order = (prod(mn)+1)/2;
    b = ordfilt2(a, order, domain, padopt);
else
    order1 = prod(mn)/2;
    order2 = order1+1;
    b = ordfilt2(a, order1, domain, padopt);
    b2 = ordfilt2(a, order2, domain, padopt);
    idx = find(b ~= b2);
    b(idx) = (double(b(idx)) + double(b2(idx)))/2;
end

%--------------------------------------------------------
function [a, mn, padopt] = parse_inputs_medfilt2(varargin)
checknargin(1,4,nargin,mfilename);
a = varargin{1};    charLocation = [];
for k = 2:nargin
    if (ischar(varargin{k})),    charLocation = [charLocation k];    end
end

if (length(charLocation) > 1)    % More than one string in input list
    eid = 'Images:medfilt2:tooManyStringInputs';
    error(eid,'%s','Too many input string arguments.');
elseif (isempty(charLocation))    % No string specified
    padopt = 'zeros';
else
    options = {'indexed', 'zeros', 'symmetric'};
    padopt = checkstrs(varargin{charLocation}, options, mfilename, 'PADOPT', charLocation);
    varargin(charLocation) = [];
end

if (strcmp(padopt, 'indexed'))
    if (isa(a,'double')),    padopt = 'ones';
	else,                    padopt = 'zeros';
    end
end

if numel(varargin) == 1
  mn = [3 3];% default
elseif numel(varargin) >= 2
  mn = varargin{2}(:).';
  if size(mn,2)~=2
    msg = 'MEDFILT2(A,[M N]): Second argument must consist of two integers.';
    error('medfilt2:secondArgMustConsistOfTwoInts', msg);
  elseif length(varargin) > 2
    msg = ['MEDFILT2(A,[M N],[Mb Nb],...) is an obsolete syntax. [Mb Nb] argument is ignored.'];
    warning('medfilt2:obsoleteSyntax', msg);
  end
end

% ---------------------------------------------------------------------------------
function B = morphop(varargin)
%MORPHOP Dilate or erode image.
%   B = MORPHOP(OP_TYPE,A,SE,...) computes the erosion or dilation of A,
%   depending on whether OP_TYPE is 'erode' or 'dilate'.  SE is a
%   STREL array or an NHOOD array.  MORPHOP is intended to be called only
%   by IMDILATE or IMERODE.  Any additional arguments passed into IMDILATE
%   or IMERODE should be passed into MORPHOP following SE.  See the help
%   entries for IMDILATE and IMERODE for more details about the allowable syntaxes.

%   Copyright 1993-2003 The MathWorks, Inc.
%   $Revision: 1.10 $  $Date: 2003/01/17 16:28:38 $

[A,se,pre_pad, pre_pack,post_crop,post_unpack,op_type,is_packed, unpacked_M,mex_method] = ParseInputs_morphop(varargin{:});

num_strels = length(se);

if pre_pad
    % Find the array offsets and heights for each structuring element
    % in the sequence.
    offsets = cell(1,num_strels);
    for k = 1:num_strels
        offsets{k} = getneighbors(se(k));
    end
    
    % Now compute how padding is needed based on the strel offsets.
    [pad_ul, pad_lr] = PadSize_morphop(offsets,op_type);
    P = length(pad_ul);
    Q = ndims(A);
    if P < Q
        pad_ul = [pad_ul zeros(1,Q-P)];
        pad_lr = [pad_lr zeros(1,Q-P)];
    end
    
    if (is_packed)			% Input is packed binary.  Adjust padding appropriately.
        pad_ul(1) = ceil(pad_ul(1) / 32);
        pad_lr(1) = ceil(pad_lr(1) / 32);
    end

    if strcmp(op_type, 'dilate')
        pad_val = -Inf;
    else
        pad_val = Inf;
    end
    if islogical(A)			% Use 0s and 1s instead of plus/minus Inf.
        pad_val = max(min(pad_val, 1), 0);
    end
    A = padarray(A,pad_ul,pad_val,'pre');
    A = padarray(A,pad_lr,pad_val,'post');
end

if pre_pack
    unpacked_M = size(A,1);
    A = bwpack(A);
end

% Apply the sequence of dilations/erosions.
B = A;
for k = 1:num_strels
    B = morphmex(mex_method, B, double(getnhood(se(k))), getheight(se(k)), unpacked_M);
end

% Image postprocessing steps.
if post_unpack
    B = bwunpack(B,unpacked_M);
end

if post_crop
    % Extract the "middle" of the result; it should be the same size as
    % the input image.
    idx = cell(1,ndims(B));
    for k = 1:ndims(B)
        P = size(B,k) - pad_ul(k) - pad_lr(k);
        first = pad_ul(k) + 1;
        last = first + P - 1;
        idx{k} = first:last;
    end
    B = B(idx{:});
end

%%%%%%%%%% ParseInputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,se,pre_pad,pre_pack, post_crop,post_unpack,op_type,input_is_packed, unpacked_M,mex_method] = ...
	ParseInputs_morphop(A,se,op_type,func_name,varargin)

checknargin(2,5,nargin-2,func_name);

% Get the required inputs and check them for validity.
se = strelcheck(se,func_name,'SE',2);
A = CheckInputImage_morphop(A, func_name);

% Process optional arguments.
[padopt,packopt,unpacked_M] = ProcessOptionalArgs_morphop(func_name, varargin{:});
if strcmp(packopt,'ispacked')
    CheckUnpackedM_morphop(unpacked_M, size(A,1));
end

% Figure out the appropriate image preprocessing steps, image 
% postprocessing steps, and MEX-file method to invoke.

% First, find out the values of all the necessary predicates.
se = getsequence(se);
num_strels = length(se);
strel_is_all_flat = all(isflat(se));
input_numdims = ndims(A);
strel_is_single = num_strels == 1;
class_A = class(A);
input_is_uint32 = strcmp(class_A,'uint32');
input_is_packed = strcmp(packopt,'ispacked');
input_is_logical = islogical(A);
input_is_2d = ndims(A) == 2;
output_is_full = strcmp(padopt,'full');

strel_is_all_2d = true;
for k = 1:length(se)
    if (ndims(getnhood(se(k))) > 2)
        strel_is_all_2d = false;
        break
    end
end

% Check for error conditions related to packing
if input_is_packed && strcmp(op_type, 'erode') && (unpacked_M < 1)
    error(sprintf('Images:%s:missingPackedM', func_name), '%s', 'M must be provided for packed erosion.');
end
if input_is_packed && ~strel_is_all_2d
    error(sprintf('Images:%s:packedStrelNot2D', func_name), ...
          '%s', 'Cannot perform packed erosion or dilation unless structuring element is 2-D.');
end
if input_is_packed && ~input_is_uint32
    error(sprintf('Images:%s:invalidPackedInputType', func_name), ...
          '%s', 'Input image must be uint32 for packed erosion or dilation.');
end
if input_is_packed && ~strel_is_all_flat
    error(sprintf('Images:%s:nonflatStrelPacked', func_name), ...
          '%s', 'Structuring element must be flat for packed erosion or dilation.');
end
if input_is_packed && (input_numdims > 2)
    error(sprintf('Images:%s:packedImageNot2D', func_name), ...
          '%s', 'Cannot perform packed erosion or dilation unless input image is 2-D.');
end
if input_is_packed && output_is_full
    error(sprintf('Images:%s:packedFull', func_name), ...
          '%s', 'Cannot perform packed erosion or dilation with the ''full'' option.');
end

% Next, use predicate values to determine the necessary
% preprocessing and postprocessing steps.

% If the user has asked for full-size output, or if there are multiple
% and/or decomposed strels, then pre-pad the input image.
pre_pad = ~strel_is_single || output_is_full;

% If the input image is logical, then the strel must be flat.
if input_is_logical && ~strel_is_all_flat
    msgId = sprintf('Images:%s:binaryWithNonflatStrel', func_name);
    error(msgId,'Function %s cannot perform dilate a binary image with a nonflat structuring element.', func_name);
end

% If the input image is logical and not packed, and if there are multiple
% all-flat strels, the prepack the input image.
pre_pack = ~strel_is_single && input_is_logical && input_is_2d && strel_is_all_flat && strel_is_all_2d;

% If we had to pre-pad the input but the user didn't specify the 'full'
% option, then crop the image before returning it.
post_crop = pre_pad && ~output_is_full;

% If this function pre-packed the image, unpack it before returning it.
post_unpack = pre_pack;

% Finally, determine the appropriate MEX-file method to invoke.
if pre_pack || strcmp(packopt,'ispacked')
    mex_method = sprintf('%s_binary_packed',op_type);
elseif input_is_logical
    mex_method = sprintf('%s_binary',op_type);
elseif strel_is_all_flat
    mex_method = sprintf('%s_gray_flat',op_type);
else
    mex_method = sprintf('%s_gray_nonflat',op_type);
end

%%%%%%%%%% ProcessOptionalArgs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [padopt,packopt,unpacked_M] = ProcessOptionalArgs_morphop(func_name, varargin)

% Default values
padopt = 'same';	packopt = 'notpacked';	unpacked_M = -1;	check_M = false;
allowed_strings = {'same','full','ispacked','notpacked'};

for k = 1:length(varargin)
	if ischar(varargin{k})
		string = checkstrs(varargin{k}, allowed_strings, func_name, 'OPTION', k+2);
		switch string
			case {'full','same'},			padopt = string;
			case {'ispacked','notpacked'}, packopt = string;
		end
	else
		unpacked_M = varargin{k};
		check_M = true;
		M_pos = k+2;
	end
end

if check_M
    checkinput(unpacked_M, {'double'}, {'real' 'nonsparse' 'scalar' 'integer' 'nonnegative'}, func_name, 'M', M_pos);
end

%%%%%%%%%% CheckInputImage %%%%%%%%%%%%%%%%%%%%%%%
function B = CheckInputImage_morphop(A,op_function)
	B = A;
	checkinput(A, {'numeric' 'logical'}, {'real' 'nonsparse'}, op_function, 'IM', 1);

%%%%%%%%%% CheckUnpackedM %%%%%%%%%%%%%%%%%%%%
function CheckUnpackedM_morphop(unpacked_M, M)
	if unpacked_M >= 0
		d = 32*M - unpacked_M;
		if (d < 0) || (d > 31)
			error('imerode:inconsistentUnpackedM','M is not consistent with the row dimension of the image.');
		end
	end

%%%%%%%%%% PadSize %%%%%%%%%%%%%%%%%%%%%%%%%%
function [pad_ul, pad_lr] = PadSize_morphop(offsets,op_type)

if isempty(offsets)
    pad_ul = zeros(1,2);
    pad_lr = zeros(1,2);
else
    num_dims = size(offsets{1},2);
    for k = 2:length(offsets)
        num_dims = max(num_dims, size(offsets{k},2));
    end
    for k = 1:length(offsets)
        offsets{k} = [offsets{k} zeros(size(offsets{k},1), num_dims - size(offsets{k},2))];
    end
    
    pad_ul = zeros(1,num_dims);
    pad_lr = zeros(1,num_dims);
    
    for k = 1:length(offsets)
        offsets_k = offsets{k};
        if ~isempty(offsets_k)
            pad_ul = pad_ul + max(0, -min(offsets_k,[],1));
            pad_lr = pad_lr + max(0, max(offsets_k,[],1));
        end
    end
    
    if strcmp(op_type,'erode')        % Swap
        tmp = pad_ul;
        pad_ul = pad_lr;
        pad_lr = tmp;
    end
end

% ---------------------------------------------------------------------------------
function tf = isflat(se)
%   ISFLAT(SE) returns true (1) if the structuring element SE is flat;
%   otherwise it returns false (0).  If SE is a STREL array, then TF is
%   the same size as SE.
% tf:           double logical array, same size as se, containing 0s and 1s. 

	tf = false(size(se));
	for k = 1:numel(se)
		tf(k) = ~any(se(k).height(:));
	end

% ---------------------------------------------------------------------------------
function se = strelcheck(in,func_name,arg_name,arg_position)
%STRELCHECK Check validity of STREL object, or convert neighborhood to STREL.
%   SE = STREL(IN) returns IN if it is already a STREL; otherwise it
%   assumes IN is a neighborhood-style array and tries to convert it to a STREL.

%   Copyright 1993-2003 The MathWorks, Inc.  $Revision: 1.7 $
  
msg = sprintf('%s\n%s', 'Expected a structuring element, which can be either',...
              'a STREL object or an array containing 0s and 1s.');

if isa(in, 'strel')
    se = in;
else
    if ~( isnumeric(in) || islogical(in) )
        msgId = sprintf('Images:%s:invalidStrelType', func_name);
        error(msgId, 'Function %s expected its %s input argument, %s, to be either numeric or logical.', ...
                      func_name, num2ordinal(arg_position), arg_name);
    else
        if (issparse(in)),	in = full(in);	end
        in = double(in);
        if ~isempty(in)
            bad_elements = (in ~= 0) & (in ~= 1);
            if any(bad_elements(:))
                msgId = sprintf('Images:%s:invalidStrelValues', func_name);
                error(msgId, '%s, the %s input argument to function %s, contained values other than 0 or 1.', ...
                      arg_name, num2ordinal(arg_position), func_name);
            end
        end
        se = Localstrel(in);
    end
end

% ---------------------------------------------------------------------------------
function [offsets,heights] = getneighbors(se)
%   [OFFSETS,HEIGHTS] = GETNEIGHBORS(SE) returns the relative locations
%   and corresponding heights for each of the neighbors in the
%   structuring element object SE.  OFFSETS is a P-by-N array where P is
%   the number of neighbors in the structuring element and N is the
%   dimensionality of the structuring element.  Each row of OFFSETS
%   contains the location of the corresponding neighbor, relative to the
%   center of the structuring element.  HEIGHTS is a P-element column
%   vector containing the height of each structuring element neighbor.

%   Copyright 1993-2003 The MathWorks, Inc.
%   $Revision: 1.5 $  $Date: 2003/01/17 16:28:42 $

if length(se) ~= 1
	error('SE must be a 1-by-1 STREL array.');
end

num_dims = ndims(se.nhood);
idx = find(se.nhood);
heights = se.height(idx);
size_nhood = size(se.nhood);
center = floor((size_nhood+1)/2);
subs = cell(1,num_dims);
[subs{:}] = ind2sub(size_nhood,idx);
offsets = [subs{:}];
offsets = reshape(offsets,length(idx),num_dims);
offsets = offsets - repmat(center, size(offsets,1), 1);

%----------------------------------------------------------------------------------
function seq = getsequence(se)
%   SEQ = GETSEQUENCE(SE), where SE is a structuring element array,
%   returns another structuring element array SEQ containing the
%   individual structuring elements that form the decomposition of SE.
%   SEQ is equivalent to SE, but the elements of SEQ have no decomposition. 
	if length(se) > 1
		se = se(:);
		seq = getsequence(se(1));
		for k = 2:length(se)
			seq = [seq; getsequence(se(k))];
		end
	elseif isempty(se)
		% A bit of a hack here to return a 1-by-0 strel array.
		seq.nhood = [];
		seq.height = [];
		seq.decomposition = [];
		seq.version = 1;
		seq(1) = [];
	else
		if isempty(se.decomposition)
			seq = se;
		else
			seq = getsequence(se.decomposition(1));
			for k = 2:length(se.decomposition)
				seq = [seq; getsequence(se.decomposition(k))];
			end
		end
	end

% ---------------------------------------------------------------------------------
function nhood = getnhood(se)
%   NHOOD = GETNHOOD(SE) returns the neighborhood associated with the structuring element SE.
	if length(se) ~= 1
		error('SE must be a 1-by-1 STREL array.');
	end
	nhood = se.nhood;

% ---------------------------------------------------------------------------------
function height = getheight(se)
%   H = GETHEIGHT(SE) returns an array the same size as GETNHOOD(SE)
%   containing the height associated with each of the structuring element
%   neighbors.  H is all zeros for a flat structuring element.
% se:      a 1-by-1 strel array; it may have no neighbors.
% height:  a double array with the same size as NHOOD = GETNHOOD(SE).

	if length(se) ~= 1
		error('SE must be a 1-by-1 STREL array.');
	end
	height = se.height;

% ---------------------------------------------------------------------------------	
function se = Localstrel(varargin)
% A STREL limited to MakeArbitraryStrel strel and no object shit

if (nargin == 0)		% No input arguments --- return empty strel
	se.nhood = [];
	se.height = [];
	se.decomposition = [];
	se.version = 1;
elseif ((nargin == 1) && isa(varargin{1}, 'strel'))
    % One strel input --- return it unchanged
    se = varargin{1};
    return
else
	[type,params] = ParseInputs_strel(varargin{:});
	if (~strcmp(type,'arbitrary'))
        error('This limited strel only supports "arbitrary" strel type.');
	end
	se = MakeArbitraryStrel(params{:});
end

%%% MakeArbitraryStrel
function se = MakeArbitraryStrel(nhood,height)

se.decomposition = [];
se.version = 1;
se.nhood = nhood ~= 0;
se.height = height;

if (~isempty(nhood) && all(nhood(:)) && ~any(height(:)))
    % Strel is flat with an all-ones neighborhood.  Decide whether to decompose it.
    size_nhood = size(nhood);
    % Heuristic --- if theoretical computation advantage is
    % at least a factor of two, then assume that the advantage
    % is worth the overhead cost of performing dilation or erosion twice.
    advantage = prod(size_nhood) / sum(size_nhood);
    if (advantage >= 2)
        num_dims = ndims(nhood);
        se.decomposition = Localstrel;
        for k = 1:ndims(nhood)
            size_k = ones(1,num_dims);
            size_k(k) = size(nhood,k);
            if k == 1
                se.decomposition = Localstrel(ones(size_k));
            else
                se.decomposition(k) = Localstrel(ones(size_k));
            end
        end
    end
end

% ---------------------------------------------------------------------------------	
function [type,params] = ParseInputs_strel(varargin)

checknargin(1, 4, nargin, 'strel');
type = 'arbitrary';
params = varargin;
num_params = numel(params);

switch type
  case 'arbitrary'
    if num_params < 1
        error('Images:strel:tooFewInputsForArbitrary','%s','Too few inputs.');
    end
    
    % Check validity of the NHOOD argument.
    nhood = params{1};
    checkinput(nhood, {'numeric', 'logical'}, {'real'}, 'strel', 'NHOOD', 2);
    
    % Check validity of the HEIGHT argument.
    if num_params >= 2
        height = params{2};
        checkinput(height, {'double'}, {'real', 'nonnan'}, 'strel', 'HEIGHT', 3);
        if ~isequal(size(height), size(nhood))
            msg = 'For arbitrary strels, the HEIGHT input must be a real double matrix with the same size as the NHOOD input.';
            error('Images:strel:sizeMismatch','%s',msg);
        end
    else
        params{2} = zeros(size(nhood));
    end
end

% ---------------------------------------------------------------------------------
function string = num2ordinal(number)
%NUM2ORDINAL Convert positive integer to ordinal string.
if number <= 20
  table1 = {'first' 'second' 'third' 'fourth' 'fifth' 'sixth' 'seventh' ...
            'eighth' 'ninth' 'tenth' 'eleventh' 'twelfth' 'thirteenth' ...
            'fourteenth' 'fifteenth' 'sixteenth' 'seventeenth' ...
            'eighteenth' 'nineteenth' 'twentieth'};
	string = table1{number};
else
	table2 = {'th' 'st' 'nd' 'rd' 'th' 'th' 'th' 'th' 'th' 'th'};
	ones_digit = rem(number, 10);
	string = sprintf('%d%s',number,table2{ones_digit + 1});
end

%----------------------------------------------------------------------------------
function [xe, ye] = poly2edgelist(x,y,scale)
%POLY2EDGELIST Computes list of horizontal edges from polygon.
%   Copyright 1993-2003 The MathWorks, Inc.
%   $Revision: 1.1 $  $Date: 2003/01/17 16:28:39 $

if (nargin < 3)    scale = 5;   end
% Scale and quantize (x,y) locations to the higher resolution grid.
x = round(scale*(x - 0.5) + 1);
y = round(scale*(y - 0.5) + 1);

num_segments = length(x) - 1;
x_segments = cell(num_segments,1);
y_segments = cell(num_segments,1);
for k = 1:num_segments
    [x_segments{k},y_segments{k}] = intline(x(k),x(k+1),y(k),y(k+1));
end

% Concatenate segment vertices.
x = cat(1,x_segments{:});   y = cat(1,y_segments{:});

% Horizontal edges are located where the x-value changes.
d = diff(x);
edge_indices = find(d);
xe = x(edge_indices);

% Wherever the diff is negative, the x-coordinate should be x-1 instead of x.
shift = find(d(edge_indices) < 0);
xe(shift) = xe(shift) - 1;

% In order for the result to be the same no matter which direction we are
% tracing the polynomial, the y-value for a diagonal transition has to be
% biased the same way no matter what.  We'll always chooser the smaller
% y-value associated with diagonal transitions.
ye = min(y(edge_indices), y(edge_indices+1));

%----------------------------------------------------------------------------------
function [eout,thresh] = edge(varargin)
%EDGE Find edges in intensity image.
%   Copyright 1993-2003 The MathWorks, Inc.  
%   $Revision: 5.28 $  $Date: 2003/01/27 20:15:58 $

[a,method,thresh,sigma,H,kx,ky] = parse_inputs_edge(varargin{:});
% Transform to a double precision intensity image if necessary
if ~isa(a, 'double')   a = im2double(a);    end
m = size(a,1);
n = size(a,2);
rr = 2:m-1; cc=2:n-1;

% The output edge map:
e = false(m, n);

if strcmp(method,'canny')
   % Magic numbers
   GaussianDieOff = .0001;  
   PercentOfPixelsNotEdges = .7; % Used for selecting thresholds
   ThresholdRatio = .4;          % Low thresh is this fraction of the high.
   
   % Design the filters - a gaussian and its derivative
   pw = 1:30; % possible widths
   ssq = sigma*sigma;
   width = max(find(exp(-(pw.*pw)/(2*sigma*sigma))>GaussianDieOff));
   if (isempty(width))      width = 1;  end % the user entered a really small sigma

   t = (-width:width);
   gau = exp(-(t.*t)/(2*ssq))/(2*pi*ssq);     % the gaussian 1D filter

   % Find the directional derivative of 2D Gaussian (along X-axis)
   % Since the result is symmetric along X, we can get the derivative along
   % Y-axis simply by transposing the result for X direction.
   [x,y]=meshgrid(-width:width,-width:width);
   dgau2D=-x.*exp(-(x.*x+y.*y)/(2*ssq))/(pi*ssq);
      
   % Convolve the filters with the image in each direction
   % The canny edge detector first requires convolution with
   % 2D gaussian, and then with the derivitave of a gaussian.
   % Since gaussian filter is separable, for smoothing, we can use 
   % two 1D convolutions in order to achieve the effect of convolving
   % with 2D Gaussian.  We convolve along rows and then columns.

   %smooth the image out
   aSmooth=imfilter(a,gau,'conv','replicate');         % run the filter accross rows
   aSmooth=imfilter(aSmooth,gau','conv','replicate');  % and then accross columns
     
   %apply directional derivatives
   ax = imfilter(aSmooth, dgau2D, 'conv','replicate');
   ay = imfilter(aSmooth, dgau2D', 'conv','replicate');

   mag = sqrt((ax.*ax) + (ay.*ay));
   magmax = max(mag(:));
   if (magmax>0)      mag = mag / magmax;   end % normalize
   
   % Select the thresholds                                                                      
   if isempty(thresh) 
      [counts,x]=imhist_j(mag, 64);
      highThresh = min(find(cumsum(counts) > PercentOfPixelsNotEdges*m*n)) / 64;
      lowThresh = ThresholdRatio*highThresh;
      thresh = [lowThresh highThresh];
   elseif length(thresh)==1
      highThresh = thresh;
      if (thresh>=1),    error('The threshold must be less than 1.');      end
      lowThresh = ThresholdRatio*thresh;
      thresh = [lowThresh highThresh];
   elseif length(thresh)==2
      lowThresh = thresh(1);
      highThresh = thresh(2);
      if (lowThresh >= highThresh) || (highThresh >= 1)
         error('Thresh must be [low high], where low < high < 1.');
      end
   end
   
   % The next step is to do the non-maximum supression.  
   % We will accrue indices which specify ON pixels in strong edgemap
   % The array e will become the weak edge map.
   idxStrong = [];  
   for dir = 1:4
      idxLocalMax = cannyFindLocalMaxima(dir,ax,ay,mag);
      idxWeak = idxLocalMax(mag(idxLocalMax) > lowThresh);
      e(idxWeak)=1;
      idxStrong = [idxStrong; idxWeak(mag(idxWeak) > highThresh)];
   end
   
   rstrong = rem(idxStrong-1, m)+1;
   cstrong = floor((idxStrong-1)/m)+1;
   e = bwselect(e, cstrong, rstrong, 8);
   e = bwmorph(e, 'thin', 1);  % Thin double (or triple) pixel wide contours
   
elseif any(strcmp(method, {'log','marr-hildreth','zerocross'}))
   % We don't use image blocks here
   if isempty(H)
      fsize = ceil(sigma*3) * 2 + 1;  % choose an odd fsize > 6*sigma;
      op = fspecial('log',fsize,sigma); 
   else 
      op = H; 
   end
   
   op = op - sum(op(:))/numel(op); % make the op to sum to zero
   b = filter2(op,a);
   
   if (isempty(thresh)),      thresh = .75*mean2(abs(b(rr,cc)));   end
  
   % Look for the zero crossings:  +-, -+ and their transposes 
   % We arbitrarily choose the edge to be the negative point
   [rx,cx] = find( b(rr,cc) < 0 & b(rr,cc+1) > 0 & abs( b(rr,cc)-b(rr,cc+1) ) > thresh );   % [- +]
   e((rx+1) + cx*m) = 1;
   [rx,cx] = find( b(rr,cc-1) > 0 & b(rr,cc) < 0 & abs( b(rr,cc-1)-b(rr,cc) ) > thresh );   % [+ -]
   e((rx+1) + cx*m) = 1;
   [rx,cx] = find( b(rr,cc) < 0 & b(rr+1,cc) > 0 & abs( b(rr,cc)-b(rr+1,cc) ) > thresh);   % [- +]'
   e((rx+1) + cx*m) = 1;
   [rx,cx] = find( b(rr-1,cc) > 0 & b(rr,cc) < 0 & abs( b(rr-1,cc)-b(rr,cc) ) > thresh);   % [+ -]'
   e((rx+1) + cx*m) = 1;
   
   % Most likely this covers all of the cases.   Just check to see if there
   % are any points where the LoG was precisely zero:
   [rz,cz] = find( b(rr,cc)==0 );
   if ~isempty(rz)
      % Look for the zero crossings: +0-, -0+ and their transposes
      % The edge lies on the Zero point
      zero = (rz+1) + cz*m;   % Linear index for zero points
      zz = b(zero-1) < 0 & b(zero+1) > 0 & abs( b(zero-1)-b(zero+1) ) > 2*thresh;     % [- 0 +]'
      e(zero(zz)) = 1;
      zz = b(zero-1) > 0 & b(zero+1) < 0 & abs( b(zero-1)-b(zero+1) ) > 2*thresh;     % [+ 0 -]'
      e(zero(zz)) = 1;
      zz = b(zero-m) < 0 & b(zero+m) > 0 & abs( b(zero-m)-b(zero+m) ) > 2*thresh;     % [- 0 +]
      e(zero(zz)) = 1;
      zz = b(zero-m) > 0 & b(zero+m) < 0 & abs( b(zero-m)-b(zero+m) ) > 2*thresh;     % [+ 0 -]
      e(zero(zz)) = 1;
   end

else  % one of the easy methods (roberts,sobel,prewitt)
   % Determine edges in blocks for easy methods 
   nr = length(rr); nc = length(cc);
   blk = bestblk([nr nc]);
   nblks = floor([nr nc]./blk); nrem = [nr nc] - nblks.*blk;
   mblocks = nblks(1); nblocks = nblks(2);
   mb = blk(1); nb = blk(2);
   
   if strcmp(method,'sobel')
      op = [-1 -2 -1;0 0 0;1 2 1]/8; % Sobel approximation to derivative
      bx = abs(filter2(op',a)); by = abs(filter2(op,a));
      b = kx*bx.*bx + ky*by.*by;
      if isempty(thresh)	% Determine cutoff based on RMS estimate of noise
         cutoff = 4*sum(sum(b(rr,cc)))/numel(b(rr,cc)); thresh = sqrt(cutoff);
      else                   % Use relative tolerance specified by the user
         cutoff = (thresh).^2;
      end
      rows = 1:blk(1);
      for i=0:mblocks
         if i==mblocks, rows = (1:nrem(1)); end
         for j=0:nblocks
            if j==0, cols = 1:blk(2); elseif j==nblocks, cols=(1:nrem(2)); end
            if ~isempty(rows) && ~isempty(cols)
               r = rr(i*mb+rows); c = cc(j*nb+cols);
               e(r,c) = (b(r,c)>cutoff) & ...
               ( ( (bx(r,c) >= (kx*by(r,c)-eps*100)) & ...
               (b(r,c-1) <= b(r,c)) & (b(r,c) > b(r,c+1)) ) | ...
               ( (by(r,c) >= (ky*bx(r,c)-eps*100 )) & ...
               (b(r-1,c) <= b(r,c)) & (b(r,c) > b(r+1,c))));
            end
         end
      end
      
   elseif strcmp(method,'prewitt')
      op = [-1 -1 -1;0 0 0;1 1 1]/6; % Prewitt approximation to derivative
      bx = abs(filter2(op',a)); by = abs(filter2(op,a));
      b = kx*bx.*bx + ky*by.*by;
      if isempty(thresh)	% Determine cutoff based on RMS estimate of noise
         cutoff = 4*sum(sum(b(rr,cc)))/numel(b(rr,cc)); thresh = sqrt(cutoff);
      else                   % Use relative tolerance specified by the user
         cutoff = (thresh).^2;
      end
      rows = 1:blk(1);
      for i=0:mblocks
         if i==mblocks, rows = (1:nrem(1)); end
         for j=0:nblocks
            if j==0, cols = 1:blk(2); elseif j==nblocks, cols=(1:nrem(2)); end
            if ~isempty(rows) && ~isempty(cols)
               r = rr(i*mb+rows); c = cc(j*nb+cols);
               e(r,c) = (b(r,c)>cutoff) & ...
               ( ( (bx(r,c) >= (kx*by(r,c)-eps*100) ) & ...
               (b(r,c-1) <= b(r,c)) & (b(r,c) > b(r,c+1)) ) | ...
               ((by(r,c) >= (ky*bx(r,c)-eps*100) )  & ...
               (b(r-1,c) <= b(r,c)) & (b(r,c) > b(r+1,c)) ) );
            end
         end
      end
      
   elseif strcmp(method, 'roberts')
      op = [1 0;0 -1]/sqrt(2); % Roberts approximation to diagonal derivative
      bx = abs(filter2(op,a)); by = abs(filter2(rot90(op),a));
      b = kx*bx.*bx + ky*by.*by;
      if isempty(thresh)	% Determine cutoff based on RMS estimate of noise
         cutoff = 6*sum(sum(b(rr,cc)))/numel(b(rr,cc)); thresh = sqrt(cutoff);
      else                   % Use relative tolerance specified by the user
         cutoff = (thresh).^2;
      end
      rows = 1:blk(1);
      for i=0:mblocks
         if i==mblocks, rows = (1:nrem(1)); end
         for j=0:nblocks
            if j==0, cols = 1:blk(2); elseif j==nblocks, cols=(1:nrem(2)); end
            if ~isempty(rows) && ~isempty(cols)
               r = rr(i*mb+rows); c = cc(j*nb+cols);
               e(r,c) = (b(r,c)>cutoff) & ...
               ( ( (bx(r,c) >= (kx*by(r,c)-eps*100)) & ...
               (b(r-1,c-1) <= b(r,c)) & (b(r,c) > b(r+1,c+1)) ) | ...
               ( (by(r,c) >= (ky*bx(r,c)-eps*100)) & ...
               (b(r-1,c+1) <= b(r,c)) & (b(r,c) > b(r+1,c-1)) ) );
            end
         end
      end
   else
      error([method,' is not a valid method.']);
   end
end
eout = e;

%----------------------------------------------------------------------------------
%   Local Function : cannyFindLocalMaxima
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

[m,n,o] = size(mag);

% Find the indices of all points whose gradient (specified by the 
% vector (ix,iy)) is going in the direction we're looking at.  

switch direction
case 1,   idx = find((iy<=0 & ix>-iy)  | (iy>=0 & ix<-iy));
case 2,   idx = find((ix>0 & -iy>=ix)  | (ix<0 & -iy<=ix));
case 3,   idx = find((ix<=0 & ix>iy) | (ix>=0 & ix<iy));
case 4,   idx = find((iy<0 & ix<=iy) | (iy>0 & ix>=iy));
end

% Exclude the exterior pixels
if ~isempty(idx)
   v = mod(idx,m);
   extIdx = v==1 | v==0 | idx<=m | (idx>(n-1)*m);
   idx(extIdx) = [];
end

ixv = ix(idx);  
iyv = iy(idx);   
gradmag = mag(idx);

% Do the linear interpolations for the interior pixels
switch direction
case 1
   d = abs(iyv./ixv);
   gradmag1 = mag(idx+m).*(1-d) + mag(idx+m-1).*d; 
   gradmag2 = mag(idx-m).*(1-d) + mag(idx-m+1).*d; 
case 2
   d = abs(ixv./iyv);
   gradmag1 = mag(idx-1).*(1-d) + mag(idx+m-1).*d; 
   gradmag2 = mag(idx+1).*(1-d) + mag(idx-m+1).*d; 
case 3
   d = abs(ixv./iyv);
   gradmag1 = mag(idx-1).*(1-d) + mag(idx-m-1).*d; 
   gradmag2 = mag(idx+1).*(1-d) + mag(idx+m+1).*d; 
case 4
   d = abs(iyv./ixv);
   gradmag1 = mag(idx-m).*(1-d) + mag(idx-m-1).*d; 
   gradmag2 = mag(idx+m).*(1-d) + mag(idx+m+1).*d; 
end
idxLocalMax = idx(gradmag>=gradmag1 & gradmag>=gradmag2); 

%----------------------------------------------------------------------------------
function [I,Method,Thresh,Sigma,H,kx,ky] = parse_inputs_edge(varargin)
% OUTPUTS:
%   I      Image Data
%   Method Edge detection method
%   Thresh Threshold value
%   Sigma  standard deviation of Gaussian
%   H      Filter for Zero-crossing detection
%   kx,ky  From Directionality vector

error(nargchk(1,5,nargin));
I = varargin{1};
checkinput(I,{'double','logical','uint8','uint16'},{'nonsparse','2d'},mfilename,'I',1);

% Defaults
Method='sobel';     Thresh=[];  Direction='both';   Sigma=2;    H=[];   K=[1 1];

methods = {'canny','prewitt','sobel','marr-hildreth','log','roberts','zerocross'};
directions = {'both','horizontal','vertical'};

% Now parse the nargin-1 remaining input arguments

% First get the strings - we do this because the intepretation of the 
% rest of the arguments will depend on the method.
nonstr = [];   % ordered indices of non-string arguments
for i = 2:nargin
   if ischar(varargin{i})
      str = lower(varargin{i});
      j = strmatch(str,methods);
      k = strmatch(str,directions);
      if ~isempty(j)
         Method = methods{j(1)};
         if strcmp(Method,'marr-hildreth')  
            warning('''Marr-Hildreth'' is an obsolete syntax, use ''LoG'' instead.');
         end
      elseif ~isempty(k)
         Direction = directions{k(1)};
      else
         error(['Invalid input string: ''' varargin{i} '''.']);
      end
   else
      nonstr = [nonstr i];
   end
end

% Now get the rest of the arguments 
switch Method
case {'prewitt','sobel','roberts'}
   threshSpecified = 0;  % Threshold is not yet specified
   for i = nonstr
      if numel(varargin{i})<=1 && ~threshSpecified % Scalar or empty
         Thresh = varargin{i};
         threshSpecified = 1;
      elseif numel(varargin{i})==2  % The dreaded K vector
         warning(['BW = EDGE(... , K) is an obsolete syntax. '...
           'Use BW = EDGE(... , DIRECTION), where DIRECTION is a string.']);
         K=varargin{i};     
      else
         error('Invalid input arguments');
      end
   end
case 'canny'
   Sigma = 1.0;          % Default Std dev of gaussian for canny
   threshSpecified = 0;  % Threshold is not yet specified
   for i = nonstr
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
         % Thresh = [];
         threshSpecified = 1;
      else
         error('Invalid input arguments');
      end
   end
case 'log'
   threshSpecified = 0;  % Threshold is not yet specified
   for i = nonstr
      if numel(varargin{i})<=1  % Scalar or empty
         if ~threshSpecified
            Thresh = varargin{i};
            threshSpecified = 1;
         else
            Sigma = varargin{i};
         end
      else
         error('Invalid input arguments');
      end
   end
case 'zerocross'
   threshSpecified = 0;  % Threshold is not yet specified
   for i = nonstr
      if numel(varargin{i})<=1 && ~threshSpecified % Scalar or empty
         Thresh = varargin{i};
         threshSpecified = 1;
      elseif numel(varargin{i}) > 1 % The filter for zerocross
         H = varargin{i};
      else
         error('Invalid input arguments');
      end
   end
case 'marr-hildreth'
   for i = nonstr
      if numel(varargin{i})<=1  % Scalar or empty
         Thresh = varargin{i};
      elseif numel(varargin{i})==2  % The dreaded K vector 
         warning('The [kx ky] direction factor has no effect for ''Marr-Hildreth''.');
      elseif numel(varargin{i}) > 2 % The filter for zerocross
         H = varargin{i};
      else
         error('Invalid input arguments');
      end
   end
otherwise
   error('Invalid input arguments');
end   

if (Sigma<=0)   error('Sigma must be positive');    end

switch Direction
case 'both',        kx = K(1); ky = K(2); 
case 'horizontal',  kx = 0; ky = 1; % Directionality factor
case 'vertical',    kx = 1; ky = 0; % Directionality factor
otherwise
   error('Unrecognized direction string');
end

%--------------------------------------------------------------
function [mb,nb] = bestblk(siz,k)
%BESTBLK Choose block size for block processing.
%   Copyright 1993-2003 The MathWorks, Inc.  
%   $Revision: 5.16 $  $Date: 2003/01/17 16:27:19 $

if nargin==1, k = 100; end  % Default block size

% Find possible factors of siz that make good blocks
% Define acceptable block sizes
m = floor(k):-1:floor(min(ceil(siz(1)/10),ceil(k/2)));
n = floor(k):-1:floor(min(ceil(siz(2)/10),ceil(k/2)));

% Choose that largest acceptable block that has the minimum padding.
[dum,ndx] = min(ceil(siz(1)./m).*m-siz(1)); blk(1) = m(ndx);
[dum,ndx] = min(ceil(siz(2)./n).*n-siz(2)); blk(2) = n(ndx);

if nargout==2,
    mb = blk(1); nb = blk(2);
else
    mb = blk;
end

%--------------------------------------------------------------
function y = mean2(x)
%MEAN2 Compute mean of matrix elements.
	y = sum(x(:))/numel(x);

%--------------------------------------------------------------
function varargout = bwselect(varargin)
%BWSELECT Select objects in binary image.
%   Copyright 1993-2003 The MathWorks, Inc.  
%   $Revision: 1.25 $  $Date: 2003/01/27 20:15:57 $

[xdata,ydata,BW,xi,yi,r,c,n] = ParseInputs_bwselect(varargin{:});

seed_indices = sub2ind(size(BW), r(:), c(:));
BW2 = imfill(~BW, seed_indices, n);
BW2 = BW2 & BW;

switch nargout
case 1    % BW2 = BWSELECT(...)
    varargout{1} = BW2;
case 2    % [BW2,IDX] = BWSELECT(...)
    varargout{1} = BW2;
    varargout{2} = find(BW2);
otherwise    % [X,Y,BW2,...] = BWSELECT(...)
    varargout{1} = xdata;
    varargout{2} = ydata;
    varargout{3} = BW2;
    if (nargout >= 4)        % [X,Y,BW2,IDX,...] = BWSELECT(...)
        varargout{4} = find(BW2);
    end
    if (nargout >= 5)        % [X,Y,BW2,IDX,Xi,...] = BWSELECT(...)
        varargout{5} = xi;
    end
    if (nargout >= 6)        % [X,Y,BW2,IDX,Xi,Yi] = BWSELECT(...)
        varargout{6} = yi;
    end
end

% ----- Subfunction ParseInputs -----------------------------------------------
function [xdata,ydata,BW,xi,yi,r,c,style,newFig] = ParseInputs_bwselect(varargin)
% Hacked for cases 0-2. Those depended on imshow and other non-wanted stuff
style = 8;
check_style = false;
xdata = []; ydata = [];
BW = [];
check_BW = false;
xi = [];    yi = [];
r = [];     c = [];
newFig = 0;

checknargin(0,6,nargin,mfilename);

switch nargin
case 3    % BWSELECT(BW,Xi,Yi)
    BW = varargin{1};       BW_position = 1;
    check_BW = true;
    xdata = [1 size(BW,2)]; ydata = [1 size(BW,1)];
    xi = varargin{2};       yi = varargin{3};
    r = round(yi);          c = round(xi);
case 4    % BWSELECT(BW,Xi,Yi,N)
    BW = varargin{1};       BW_position = 1;
    check_BW = true;
    xdata = [1 size(BW,2)]; ydata = [1 size(BW,1)];
    xi = varargin{2};       yi = varargin{3};
    r = round(yi);          c = round(xi);
    style = varargin{4};    style_position = 4;
    check_style = true;
case 5    % BWSELECT(X,Y,BW,Xi,Yi)
    xdata = varargin{1};    ydata = varargin{2};
    BW = varargin{3};       BW_position = 3;
    check_BW = true;
    xi = varargin{4};       yi = varargin{5};
    r = round(axes2pix(size(BW,1), ydata, yi));
    c = round(axes2pix(size(BW,2), xdata, xi));
case 6    % BWSELECT(X,Y,BW,Xi,Yi,N)
    xdata = varargin{1};    ydata = varargin{2};
    BW = varargin{3};       BW_position = 3;
    check_BW = true;
    xi = varargin{4};       yi = varargin{5};
    style = varargin{6};    style_position = 6;
    check_style = true;
    r = round(axes2pix(size(BW,1), ydata, yi));
    c = round(axes2pix(size(BW,2), xdata, xi));
end

if check_BW
    checkinput(BW,{'logical' 'numeric'},{'2d' 'nonsparse'}, mfilename, 'BW', BW_position);
end

if (~islogical(BW))    BW = BW ~= 0;    end

if (check_style)
    checkinput(style, {'numeric'}, {'scalar'}, mfilename, 'N', 2);
end

badPix = find((r < 1) | (r > size(BW,1)) | (c < 1) | (c > size(BW,2)));
if (~isempty(badPix))
    msgId = sprintf('Images:%s:outOfRange', mfilename);
    warning(msgId, '%s', 'Ignoring out-of-range input coordinates');
    r(badPix) = [];
    c(badPix) = [];
end 

% ----------------------------------------------------------------------------------
function [I2,locations] = imfill(varargin)
%IMFILL Fill image regions.
%   Copyright 1993-2003 The MathWorks, Inc.  
%   $Revision: 1.13 $  $Date: 2003/01/27 23:59:37 $

[I,locations,conn,do_fillholes] = parse_inputs_imfill(varargin{:});
if do_fillholes
    if (islogical(I))   mask = uint8(I);
    else                mask = I;    end
    mask = padarray(mask, ones(1,ndims(mask)), -Inf, 'both');

    marker = mask;
    idx = cell(1,ndims(I));
    for k = 1:ndims(I)
        idx{k} = 2:(size(marker,k) - 1);
    end
    marker(idx{:}) = Inf;

    mask = imcomplement(mask);
    marker = imcomplement(marker);
    I2 = imreconstruct(marker, mask, conn);
    I2 = imcomplement(I2);
    I2 = I2(idx{:});
    if (islogical(I))   I2 = I2 ~= 0;    end
else    
    mask = imcomplement(I);
    marker = mask;
    marker(:) = 0;
    marker(locations) = mask(locations);
    marker = imreconstruct(marker, mask, conn);
    I2 = I | marker;
end

%----- Subfunction ParseInputs ---------------------------------------------------
function [IM,locations,conn,do_fillholes] = parse_inputs_imfill(varargin)
    
checknargin(1,3,nargin,mfilename);
IM = varargin{1};
checkinput(IM, {'numeric' 'logical'}, {'nonsparse' 'real'}, mfilename, 'I1 or BW1', 1);
do_interactive = false;
do_fillholes = false;
conn = conndef(ndims(IM),'minimal');
do_conn_check = false;

locations = [];
do_location_check = false;

switch nargin
  case 1
    if islogical(IM)        % IMFILL(BW1)
        do_interactive = true;
    else        % IMFILL(I1)
        do_fillholes = true;
    end
  case 2
    if islogical(IM)
        if ischar(varargin{2})            % IMFILL(BW1, 'holes')
            checkstrs(varargin{2}, {'holes'}, mfilename, 'OPTION', 2);
            do_fillholes = true;
        else            % IMFILL(BW1, LOCATIONS)
            locations = varargin{2};
            do_location_check = true;
        end
    else
        if ischar(varargin{2})            % IMFILL(I1, 'holes')
            checkstrs(varargin{2}, {'holes'}, mfilename, 'OPTION', 2);
            do_fillholes = true;
        else            % IMFILL(I1, CONN)
            conn = varargin{2};
            do_conn_check = true;
            conn_position = 2;
        end
    end
  case 3
    if islogical(IM)
        if ischar(varargin{3})  % IMFILL(BW1,CONN,'holes')
            checkstrs(varargin{3}, {'holes'}, mfilename, 'OPTION', 3);
            do_fillholes = true;
            conn = varargin{2};
            do_conn_check = true;
            conn_position = 2;
        else
            if isequal(varargin{2}, 0)  % IMFILL(BW1,0,CONN)
                do_interactive = true;
                conn = varargin{3};
                do_conn_check = true;
                conn_position = 2;
            else                % IMFILL(BW1,LOCATIONS,CONN)
                locations = varargin{2};
                do_location_check = true;
                conn = varargin{3};
                do_conn_check = true;
                conn_position = 3;
            end
        end
    else        % IMFILL(I1,CONN,'holes')
        checkstrs(varargin{3}, {'holes'}, mfilename, 'OPTION', 3);
        do_fillholes = true;
        conn = varargin{2};
        do_conn_check = true;
        conn_position = 2;
    end
end

if do_conn_check
    checkconn(conn, mfilename, 'CONN', conn_position);
end

if do_location_check
    checkinput(locations, {'double'}, {'real' 'positive' 'integer' '2d'}, mfilename, 'LOCATIONS', 2);
    locations = check_locations(locations, size(IM));
elseif do_interactive
    warndlg('Interactive get locations not available','Warning')
end

% Convert to linear indices if necessary.
if ~do_fillholes && (size(locations,2) ~= 1)
    idx = cell(1,ndims(IM));
    for k = 1:ndims(IM)
        idx{k} = locations(:,k);
    end
    locations = sub2ind(size(IM), idx{:});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function locations = check_locations(locations, image_size)
%   Checks validity of LOCATIONS.  Converts LOCATIONS to linear index
%   form.  Warns if any locations are out of range.

checkinput(locations, {'double'}, {'real' 'positive' 'integer' '2d'}, mfilename, 'LOCATIONS', 2);
num_dims = length(image_size);
if (size(locations,2) ~= 1) && (size(locations,2) ~= num_dims)
    msgId = sprintf('Images:%s:badLocationSize', mfilename);
    error(msgId,'Function %s expected its %s input argument, LOCATIONS, to have either 1 or NDIMS(IM) columns.', ...
          mfilename, num2ordinal(2));
end

if size(locations,2) == 1
    bad_pix = (locations < 1) | (locations > prod(image_size));
else
    bad_pix = zeros(size(locations,1),1);
    for k = 1:num_dims
        bad_pix = bad_pix | ((locations(:,k) < 1) | (locations(:,k) > image_size(k)));
    end
end
    
if any(bad_pix)
    msgId = sprintf('Images:%s:outOfRange', mfilename);
    warning(msgId, '%s', 'Ignoring out-of-range locations.');
    locations(bad_pix,:) = [];
end

% --------------------------------------------------------------------------------
function conn = conndef(num_dims,type)
%CONNDEF Default connectivity array.
%   Copyright 1993-2003 The MathWorks, Inc.
%   $Revision: 1.6 $  $Date: 2003/01/17 16:27:24 $

error(nargchk(2,2,nargin))

if (~ischar(type))    error('TYPE must be a string.');  end

if ~isnumeric(num_dims) || (numel(num_dims) ~= 1)
    error('NUM_DIMS must be a scalar integer >= 2.');
end
num_dims = double(num_dims);
if (num_dims ~= round(num_dims)) || (num_dims < 2)
    error('NUM_DIMS must be a scalar integer >= 2.');
end

allowed_strings = {'minimal','maximal'};
idx = strmatch(lower(type),allowed_strings);
if isempty(idx)
    error(sprintf('Unrecognized TYPE string: %s', type));
elseif length(idx) > 1
    error(sprintf('Ambiguous TYPE string: %s', type));
else
    type = allowed_strings{idx};
end

switch type
  case 'minimal'
    conn = zeros(repmat(3,1,num_dims));
    conn((end+1)/2) = 1;
    idx = repmat({2},1,num_dims);
    for k = 1:num_dims
        idx{k} = 1;
        conn(idx{:}) = 1;
        idx{k} = 3;
        conn(idx{:}) = 1;
        idx{k} = 2;
    end
  case 'maximal'
    conn = ones(repmat(3,1,num_dims));
  otherwise
    error('Internal error: unexpected TYPE string.');
end

% --------------------------------------------------------------------------------
function checkconn(conn,function_name,variable_name,arg_position)
%CHECKCONN Check validity of connectivity argument.
%    CHECKCONN(CONN) issues an error message if CONN is not a valid connectivity array.
%    Copyright 1993-2003 The MathWorks, Inc.
%    $Revision: 1.7 $  $Date: 2003/01/17 16:28:15 $

checkinput(conn,{'double' 'logical'},{'real' 'nonsparse'},function_name,variable_name,arg_position);

if all(size(conn) == 1)
    if (conn ~= 1) && (conn ~= 4) && (conn ~= 8) && (conn ~= 6) && (conn ~= 18) && (conn ~= 26)
        msg1 = first_line(variable_name, function_name, arg_position);
        msg2 = 'A scalar connectivity specifier must be 1, 4, 6, 8, 18, or 26.';
        error(sprintf('Images:%s:badScalarConn', function_name),'%s\n%s',msg1,msg2);
    end
else
    if any(size(conn) ~= 3)
        msg1 = first_line(variable_name, function_name, arg_position);
        msg2 = 'A nonscalar connectivity specifier must be 3-by-3-by- ... -by-3.';
        error(sprintf('Images:%s:badConnSize', function_name),'%s\n%s',msg1,msg2);
    end
    
    if any((conn(:) ~= 1) & (conn(:) ~= 0))
        msg1 = first_line(variable_name, function_name, arg_position);
        msg2 = 'A nonscalar connectivity specifier must contain only 0s and 1s.';
        error(sprintf('Images:%s:badConnValue', function_name),'%s\n%s',msg1,msg2);
    end
    
    if conn((end+1)/2) == 0
        msg1 = first_line(variable_name, function_name, arg_position);
        msg2 = 'The central element of a connectivity specifier must be nonzero.';
        error(sprintf('Images:%s:badConnCenter', function_name),'%s\n%s',msg1,msg2);
    end
    
    if ~isequal(conn(1:end), conn(end:-1:1))
        msg1 = first_line(variable_name, function_name, arg_position);
        msg2 = 'A connectivity specifier must be symmetric about its center.';
        error(sprintf('Images:%s:nonsymmetricConn', function_name),'%s\n%s',msg1,msg2);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function str = first_line(variable_name, function_name, arg_position)
str = sprintf('Function %s expected its %s input argument, %s,\nto be a valid connectivity specifier.', ...
              upper(function_name), num2ordinal(arg_position), variable_name);

% -----------------------------------------------------------------------------
function im = imreconstruct(varargin)
%IMRECONSTRUCT Perform morphological reconstruction.
%   Copyright 1993-2003 The MathWorks, Inc.
%   $Revision: 1.11 $  $Date: 2003/03/05 22:29:19 $

[marker,mask,conn] = ParseInputs_imreconstruct(varargin{:});
if nargin == 3
  im = imreconstructmex(marker,mask,conn);
else
  im = imreconstructmex(marker,mask);
end

%---------------------------------------------------
function [Marker,Mask,Conn] = ParseInputs_imreconstruct(varargin)
Conn = [];
checknargin(2,3,nargin,mfilename);
checkinput(varargin{1},{'numeric','logical'},{'real','nonsparse'},mfilename, 'MARKER', 1);
checkinput(varargin{2},{'numeric','logical'},{'real','nonsparse'},mfilename, 'MASK', 2);

Marker = varargin{1};   Mask = varargin{2};

if nargin==3
  checkconn(varargin{3},mfilename,'CONN',3);
  Conn = varargin{3};
end

% ----------------------------------------------------------------
function im2 = imcomplement(im)
%IMCOMPLEMENT Complement image.
%   Copyright 1993-2003 The MathWorks, Inc.
%   $Revision: 1.14 $  $Date: 2003/01/27 20:16:00 $

switch class(im)
 case 'logical'
  im2 = ~im;
 case 'double'
  im2 = 1 - im;
 case 'single'
  im2 = imlincomb(-1,im,1);
 case 'uint8'
  lut = uint8(255:-1:0);
  im2 = uintlut(im,lut);
 case 'uint16'
  lut = uint16(65535:-1:0);
  im2 = uintlut(im,lut);
 case 'uint32'
  im2 = imlincomb(-1,im,double(uint32(inf)));
 case 'int8'
  im2 = imlincomb(-1,im,double(int8(inf)));
 case 'int16'
  im2 = imlincomb(-1,im,double(int16(inf)));
 case 'int32'
  im2 = imlincomb(-1,im,double(int32(inf)));
 otherwise
  error('Invalid input image class.')
end

% ----------------------------------------------------------------
function [cout,lut] = bwmorph(a,op,n)
%BWMORPH Perform morphological operations on binary image.
%   Copyright 1993-2003 The MathWorks, Inc.  
%   $Revision: 5.30 $  $Date: 2003/01/27 20:15:56 $
% Hacked to use only the THIN, BRIDGE, CLEAN, PERIM4, PERIM8, SHRINK, DILATE, ERODE & SKEL methods.
% When more are needed, add them.
error(nargchk(2,3,nargin));
if (nargin < 3)    n = 1;   end

checkinput(a,{'numeric' 'logical'},{'real' 'nonsparse' '2d'}, mfilename, 'BW', 1);
if (~islogical(a))    a = a ~= 0;   end

if isstr(op),    % BWMORPH(A, 'op', n)
    % Find out what operation has been requested
    opString = lower(op);
    matchStrings = [ ...
            'bothat  ' 
            'bridge  '
            'clean   '
            'close   '
            'diag    '
            'dilate  '
            'erode   '
            'fatten  '
            'fill    '
            'hbreak  '
            'majority'
            'perim4  '
            'perim8  '
            'open    '
            'remove  '
            'shrink  '
            'skeleton'
            'spur    '
            'thicken '
            'thin    '
            'tophat  '];
    functionStrings = [
            'bothat  '
            'bridge  '
            'clean   '
            'close   '
            'diag_   '
            'dilate  '
            'erode   '
            'fatten  '
            'fill_   '
            'hbreak  '
            'majority'
            'perim4  '
            'perim8  '
            'open    '
            'remove  '
            'shrink  '
            'skel    '
            'spur    '
            'thicken '
            'thin    '
            'tophat  '];
    idx = strmatch(opString, matchStrings);
    if (isempty(idx))
        error(sprintf('Unknown operation "%s"', opString));
    elseif (length(idx) > 1)
        error(sprintf('Input string "%s" matches multiple operations', opString));
    end
   
    % Call the appropriate subfunction
    fcn = deblank(functionStrings(idx,:));
    c = a;    iter = 1;    done = n == 0;
    while (~done)
        lastc = c;
        [c,lut] = feval(fcn, c);
        done = ((iter >= n) | isequal(lastc, c));
        iter = iter + 1;
    end
else    % BWMORPH(A, lut, n)
    % Pass on the call to applylut
    lut = op;
    if (isempty(lut))   error('LUT can''t be empty');   end
    c = a;    done = n == 0;    iter = 1;
    while (~done)
        lastc = c;
        c = applylut(c, lut);
        done = ((iter >= n) | isequal(lastc, c));
        iter = iter + 1;
    end
    
end

cout = c;
if ((nargout == 2) && isempty(lut))
    message = ['LUT output argument is no longer supported', ...
                ' for the "', deblank(matchStrings(idx,:)), '" operation'];
    warning(message);
end

% ------------------- Function THIN
function [c,lut] = thin(a)
% Louisa Lam, Seong-Whan Lee, and Ching Y. Wuen, "Thinning Methodologies-A
% Comprehensive Survey," IEEE TrPAMI, vol. 14, no. 9, pp. 869-885, 1992.  The
% algorithm is described at the bottom of the first column and the top of the
% second column on page 879.

lut = [];
if (isempty(a))
    c = zeros(size(a));    return;
end

G1 = uint8(luts('lutthin1'));   G2 = uint8(luts('lutthin2'));
G3 = uint8(luts('lutthin3'));   G4 = uint8(luts('lutthin4'));

% Make a lookup table that will produce a lookup table indices.  This is avoid
% doing more work in calling applylut multiple times with the same input than
% we really need to.

lutlut = 1:512;
lookup = applylut(a, lutlut);
% Preallocate a uint8 zeros matrix
d = uint8(0);
[m,n] = size(a);
d(m,n) = 0;

% First subiteration
d(:) = G1(lookup) & G2(lookup) & G3(lookup);
c = a & ~d;

% Second subiteration
lookup = applylut(c, lutlut);
d(:) = G1(lookup) & G2(lookup) & G4(lookup);
c = c & ~d;

% ------------------------------------------------------------------------
% Function BRIDGE
function [c,lut] = bridge(a)
    lut = lutbridge;
    c = applylut(a, lut);

% ------------------------------------------------------------------------
% Function CLEAN
function [c,lut] = clean(a)
    lut = lutclean;
    c = applylut(a, lut);

% ------------------------------------------------------------------------
% Function PERIM4
function [c,lut] = perim4(a)
    lut = lutper4;
    c = applylut(a, lut);
    
% ------------------------------------------------------------------------
% Function DILATE
function [c,lut] = dilate(a)
    lut = ones(512,1);      lut(1) = 0;     % This is what 'lutdilate' does
    c = applylut(a, lut);
    
% ------------------------------------------------------------------------
% Function ERODE
function [c,lut] = erode(a)
    lut = zeros(512,1);      lut(512) = 1;  % This is what 'luterode' does
    c = applylut(a, lut);

% ------------------------------------------------------------------------
% Function SHRINK
function [c,lut] = shrink(a)
	lut = [];
	table = uint8(lutshrink);
	c = a;
	
	% First subiteration
	m = applylut(c, table);
	sub = c & ~m;
	c(1:2:end,1:2:end) = sub(1:2:end,1:2:end);
	
	% Second subiteration
	m = applylut(c, table);
	sub = c & ~m;
	c(2:2:end,2:2:end) = sub(2:2:end,2:2:end);
	
	% Third subiteration
	m = applylut(c, table);
	sub = c & ~m;
	c(1:2:end,2:2:end) = sub(1:2:end,2:2:end);
	
	% Fourth subiteration
	m = applylut(c, table);
	sub = c & ~m;
	c(2:2:end,1:2:end) = sub(2:2:end,1:2:end);

% ------------------------------------------------------------------------
% Function SKEL
function [c,lut] = skel(a)

	lut = [];
    skel1 = zeros(512,1);   skel1([90 92 218 220]) = 1;
    skel2 = zeros(512,1);   skel2([153 154 217 218 409 410 473 474]) = 1;
    skel3 = zeros(512,1);   skel3([24 32 56 64]) = 1;
    skel4 = zeros(512,1);   skel4([27 28 31 32 91 92 95 96]) = 1;
    skel5 = zeros(512,1);   skel5([465 473 497 505]) = 1;
    skel6 = zeros(512,1);   skel6([51 52 55 56 307 308 311 312]) = 1;
    skel7 = zeros(512,1);   skel7([309 311 437 439]) = 1;
    skel8 = zeros(512,1);   skel8([177 181 241 245 433 437 497 501]) = 1;
    
	c = a;
	c = c & ~applylut(c, skel1);
	c = c & ~applylut(c, skel2);
	c = c & ~applylut(c, skel3);
	c = c & ~applylut(c, skel4);
	c = c & ~applylut(c, skel5);
	c = c & ~applylut(c, skel6);
	c = c & ~applylut(c, skel7);
	c = c & ~applylut(c, skel8);

% ------------------------------------------------------------------------
function B = applylut(varargin)
%APPLYLUT Perform neighborhood operations using lookup tables.
%   Copyright 1993-2003 The MathWorks, Inc.
%   $Revision: 1.22 $  $Date: 2003/03/05 22:29:20 $
    [A,lut] = ParseInputs_applylut(varargin{:});
    B = applylutc(A,lut);

%---------------------------------
function [A,LUT] = ParseInputs_applylut(varargin)
checknargin(2,2,nargin,mfilename);
checkinput(varargin{1}, {'numeric','logical'},{'real','nonsparse','2d'}, mfilename, 'A', 1);
checkinput(varargin{2}, {'numeric','logical'},{'real','vector'}, mfilename, 'LUT', 2); 

% force A to be logical
A = varargin{1};
if (~islogical(A))    A = A ~= 0;   end

% force LUT to be logical
LUT = varargin{2};
if (~isa(LUT,'double'))  LUT = double(LUT); end

% ------------------------------------------------------------------------
function J = stdfilt(varargin)
%STDFILT Local standard deviation of image.
%   J = STDFILT(I) returns the array J, where each output pixel contains the
%   standard deviation value of the 3-by-3 neighborhood around the corresponding
%   pixel in the input image I. I can have any dimension.  The output image J is
%   the same size as the input image I.
%
%   For pixels on the borders of I, STDFILT uses symmetric padding.  In
%   symmetric padding, the values of padding pixels are a mirror reflection
%   of the border pixels in I.
%  
%   J = STDFILT(I,NHOOD) performs standard deviation filtering of the input
%   image I where you specify the neighborhood in NHOOD.  NHOOD is a
%   multidimensional array of zeros and ones where the nonzero elements specify
%   the neighbors.  NHOOD's size must be odd in each dimension. 
%
%   By default, STDFILT uses the neighborhood ones(3). STDFILT determines the
%   center element of the neighborhood by FLOOR((SIZE(NHOOD) + 1)/2). For
%   information about specifying neighborhoods, see Notes.
%
%   Class Support
%   -------------    
%   I can be logical or numeric and must be real and nonsparse.  NHOOD can be
%   logical or numeric and must contain zeros and/or ones.  I and NHOOD can have
%   any dimension. J is uint8, single or double.
%
%   Notes
%   -----    
%   To specify the neighborhoods of various shapes, such as a disk, use the
%   STREL function to create a structuring element object and then use the
%   GETNHOOD function to extract the neighborhood from the structuring element object.

%   Copyright 1993-2005 The MathWorks, Inc. $Revision.2 $  

[I, h] = ParseInputs_stdfilt(varargin{:});
c = class(I);
if (~isa(I,'double')),	I = double(I);	end
n = sum(h(:));

% If n = 1 then return default J (all zeros) to avoid the divideByZero warning.
% Otherwise, calculate standard deviation. The formula for standard deviation
% can be rewritten in terms of the theoretical definition of
% convolution. However, in practise, use correlation in IMFILTER to avoid a
% flipped answer when NHOOD is asymmetric.
% conv1 = imfilter(I.^2,h,'symmetric') / (n-1); 
% conv2 = imfilter(I,h,'symmetric').^2 / (n*(n-1));
% std = sqrt(conv1-conv2).  
% These equations can be further optimized for speed.

n1 = n - 1;
if n ~= 1
	conv1 = imfilter(I.^2, h/n1 , 'symmetric');
	conv2 = imfilter(I, h, 'symmetric').^2 / (n*n1);
	clear I;
	J = sqrt(max((conv1 - conv2),0));
else
	J = zeros(size(I));
end

% I only antecipate to send in uint8 or singles but we should have a generic function here
if (strcmp(c,'uint8') || strcmp(c,'logical'))
	fac = 255 / max(J(:));
	J = uint8(round(J * fac));
elseif (strcmp(c,'single'))				% Grid?
	J = single(J);
end

%%%%%%%%%%%%%%%ParseInputs%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [I,H] = ParseInputs_stdfilt(varargin)

iptchecknargin(1,2,nargin,'stdfilt');
iptcheckinput(varargin{1},{'numeric','logical'},{'real','nonsparse'}, 'stdfilt', 'I',1);
I = varargin{1};

if nargin == 2
	iptcheckinput(varargin{2},{'logical','numeric'},{'nonsparse'}, 'stdfilt','NHOOD',2);
	H = varargin{2};
	
	eid = 'stdfilt:invalidNeighborhood';
	% H must contain zeros and/or ones.
	bad_elements = (H ~= 0) & (H ~= 1);
	if any(bad_elements(:))
		error(eid,'%s','NHOOD must be a matrix that contains zeros and/or ones.');
	end
	
	% H's size must be a factor of 2n-1 (odd).
	sizeH = size(H);
	if any(floor(sizeH/2) == (sizeH/2) )
		error(eid,'%s','NHOOD must have a size that is odd in each dimension.');
	end
	if (~isa(H,'double')),	H = double(H);		end
else
	H = ones(3);
end

% ------------------------------------------------------------------------
function J = rangefilt(varargin)
%RANGEFILT Local range of image.  
%   J = RANGEFILT(I) returns the array J, where each output pixel contains the
%   range value (maximum value - minimum value) of the 3-by-3 neighborhood
%   around the corresponding pixel in the input image I. I can have any
%   dimension.  The output image J is the same size as the input image I.
%  
%   J = RANGEFILT(I,NHOOD) performs range filtering of the input image I where
%   you specify the neighborhood in NHOOD.  NHOOD is a multidimensional array
%   of zeros and ones where the nonzero elements specify the neighborhood for
%   the range filtering operation.  NHOOD's size must be odd in each dimension.
%  
%   By default, RANGEFILT uses the neighborhood true(3). RANGEFILT determines
%   the center element of the neighborhood by FLOOR((SIZE(NHOOD) + 1)/2). For
%   information about specifying neighborhoods, see Notes.
%
%   Class Support
%   -------------      
%   I can be logical or numeric and must be real and nonsparse.  NHOOD can be
%   logical or numeric and must contain zeros and/or ones.
%  
%   The output image J is the same class as I, except for signed integer data
%   types. The output class for signed integer data types is the corresponding
%   unsigned integer data type.  For example, if the class of I is int8, then
%   the class of J is uint8.
%
%   Notes
%   -----    
%   RANGEFILT uses the morphological functions IMDILATE and IMERODE to
%   determine the maximum and minimum values in the specified neighborhood.
%   Consequently, RANGEFILT uses the padding behavior of these morphological
%   functions.  
%  
%   To specify the neighborhoods of various shapes, such as a disk, use the
%   STREL function to create a structuring element object and then use the
%   GETNHOOD function to extract the neighborhood from the structuring element object.

%   Copyright 1993-2005 The MathWorks, Inc. $Revision.1 $

[I, h] = ParseInputs_rngfilt(varargin{:});

% NHOOD is reflected across its origin in order for IMDILATE
% to return the local maxima of I in NHOOD if it is asymmetric. A symmetric NHOOD
% is naturally unaffected by this reflection.
reflectH = h(:);
reflectH = flipud(reflectH);
reflectH = reshape(reflectH, size(h));

dilateI = imdilate(I,reflectH);

% IMERODE returns the local minima of I in NHOOD.
erodeI = imerode(I,h);  

% Set the output classes for signed integer data types.
class_in = class(I);
class_out = class_in;
switch class_in
	case 'int8'
		class_out = 'uint8';
	case 'int16'
		class_out = 'uint16';
	case 'int32'
		class_out = 'uint32';
end

% Calculate the range with imlincomb instead of imsubtract so that you can
% specify the output class.  Use the relational operator to calculate the
% range for a logical image to be efficient.
if strcmp(class_in, 'logical')
	J = dilateI > erodeI;
else
	J = imlincomb(1, dilateI, -1, erodeI, class_out);
end
         
%%%%%%%%%%%%%%%ParseInputs%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [I,H] = ParseInputs_rngfilt(varargin)

iptchecknargin(1,2,nargin,'rangefilt');
iptcheckinput(varargin{1},{'numeric' 'logical'}, {'real','nonsparse','nonnan'}, mfilename,'I',1);
I = varargin{1};

if nargin == 2
	iptcheckinput(varargin{2},{'logical','numeric'},{'nonempty','nonsparse'}, 'rangefilt','NHOOD',2);
	H = varargin{2};
	
	eid = 'rangefilt:invalidNeighborhood';
	% H must contain zeros and or ones.
	bad_elements = (H ~= 0) & (H ~= 1);
	if any(bad_elements(:))
		error(eid,'%s','NHOOD must be a matrix that contains zeros and/or ones.');
	end
	
	% H's size must be odd.
	sizeH = size(H);
	if any(floor(sizeH/2) == (sizeH/2) )
		error(eid,'%s','NHOOD''s size must be odd in each dimension.');
	end
	H = H ~= 0;		% Convert to logical
else
	H = true(3);
end

% ------------------------------------------------------------------------
function B = find_holes(a)
	% Track ==>   a = bwperim(a,8);   B = bwboundaries(a,8,'noholes');
	lut = luts('lutper8');
	BW = applylutc(a, lut);
	L = bwlabel(BW, 8);
	O = bwboundariesmex(L, 8);
	B = [O; {}];

%---------------------------------------------------------------------
function  varargout = luts(opt,varargin)

switch opt
    case 'lutthin1'
        varargout{1} = lutthin1(varargin{:});
    case 'lutthin2'
        varargout{1} = lutthin2(varargin{:});
    case 'lutthin3'
        varargout{1} = lutthin3(varargin{:});
    case 'lutthin4'
        varargout{1} = lutthin4(varargin{:});
    case 'lutper8'
        varargout{1} = lutper8(varargin{:});
end

%---------------------------------------------------------------------
function lut = lutthin1()
%LUTTHIN1 Compute "thin1" look-up table.
%   Copyright 1993-2003 The MathWorks, Inc.  
%   $Revision: 1.11 $  $Date: 2003/01/17 16:28:37 $
lut = [ ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     0     1     1     1     1     0     1     1 ...
     1     1     1     1     0     0     1     1     0     0     0     0 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     1     0     1     1     1     0     1     1     0     0     1     1 ...
     0     0     1     1     0     0     0     0     0     0     0     0 ...
     0     0     0     0     0     0     0     0     1     0     0     0 ...
     0     0     0     0     1     1     1     1     0     0     1     1 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     0     0     1     1     0     0     1     1     0     0     0     0 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     1     0     0     0     0     0     0     0     1     1     1     1 ...
     0     0     1     1     0     0     0     0     0     0     0     0 ...
     0     0     0     0     0     0     0     0     1     0     1     1 ...
     1     0     1     1     1     1     0     0     1     1     0     0 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     1     0     0     0     0     0     0     0 ...
     1     1     1     1     0     0     1     1     0     0     0     0 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     1     0     1     1     1     0     1     1     1     1     0     0 ...
     1     1     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     0     0     0     0     1     0     0     0 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     1     0     1     1     1     0     1     1 ...
     0     0     1     1     0     0     1     1     0     0     0     0 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     0     0     1     1     0     0     1     1 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     1     0     0     0     0     0     0     0 ...
     1     1     1     1     0     0     1     1     0     0     0     0 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     1     0     1     1     1     0     1     1     1     1     0     0 ...
     1     1     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     0     0     0     0     1     0     0     0 ...
     0     0     0     0     1     1     1     1     0     0     1     1 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     1     0     1     1     1     0     1     1 ...
     1     1     0     0     1     1     0     0];
lut = lut(:);

%---------------------------------------------------------------------
function lut = lutthin2()
%LUTTHIN2 Compute "thin2" look-up table.
%   Copyright 1993-2003 The MathWorks, Inc.  
%   $Revision: 1.6 $  $Date: 2003/01/17 16:28:37 $
lut = [ ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     0     0     0     0     0     1     0     1 ...
     0     0     1     1     1     1     1     1     0     0     0     0 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     0     1     1     1     0     1     1     1     1     1     1     1 ...
     1     1     1     1     0     0     0     0     0     0     0     0 ...
     0     0     0     0     0     0     0     0     0     1     1     1 ...
     1     1     1     1     0     1     1     1     1     1     1     1 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     1     1     1     1     1     1     1     1 ...
     1     1     1     1     1     1     1     1     0     0     0     0 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     0     1     1     1     1     1     1     1     1     1     1     1 ...
     1     1     1     1     0     0     0     0     0     0     0     0 ...
     0     0     0     0     0     0     0     0     1     1     1     1 ...
     1     1     1     1     1     1     0     0     1     0     0     0 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     0     1     1     1     1     1     1     1 ...
     1     1     1     1     1     1     1     1     0     0     0     0 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     1     1     1     0     1     0     1     0     1     1     0     0 ...
     1     0     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     0     0     0     0     0     1     1     1 ...
     1     1     1     1     1     1     1     1     1     1     1     1 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     0     1     1     1     1     1     1     1 ...
     1     1     1     1     1     1     1     1     0     0     0     0 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     1     1     1     1     1     0     1     0     1     1     1     1 ...
     1     0     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     0     0     0     0     1     1     1     0 ...
     1     0     1     0     1     1     0     0     1     0     0     0 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     0     1     1     1     1     1     1     1 ...
     1     1     1     1     1     0     0     0     0     0     0     0 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     1     1     1     1     1     1     1     1     1     1     0     0 ...
     1     0     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     0     0     0     0     1     1     1     1 ...
     1     0     1     0     1     1     1     1     1     0     0     0 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     1     1     1     0     1     0     1     0 ...
     1     1     0     0     1     0     0     0];
lut = lut(:);

%---------------------------------------------------------------------
function lut = lutthin3()
%LUTTHIN3 Compute "thin3" look-up table.
%   Copyright 1993-2003 The MathWorks, Inc.  
%   $Revision: 1.11 $  $Date: 2003/01/17 16:28:38 $
lut = [ ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     1     1     1     1     1     1     1     1 ...
     1     1     1     1     1     1     1     1     0     0     0     0 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     1     1     1     1     1     1     1     1     1     1     1     1 ...
     1     1     1     1     0     0     0     0     0     0     0     0 ...
     0     0     0     0     0     0     0     0     1     1     1     1 ...
     1     1     1     1     1     1     1     1     1     1     1     1 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     1     1     1     1     1     1     1     1 ...
     1     1     1     1     1     1     1     1     0     0     0     0 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     0     0     0     0     1     1     1     1 ...
     1     1     1     1     1     1     1     1     1     1     1     1 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     1     1     1     1     1     1     1     1 ...
     1     1     1     1     1     1     1     1     0     0     0     0 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     1     1     1     1     1     1     1     1     1     1     1     1 ...
     1     1     1     1     0     0     0     0     0     0     0     0 ...
     0     0     0     0     0     0     0     0     1     1     1     1 ...
     1     1     1     1     1     1     1     1     1     1     1     1 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     1     1     1     1     1     1     1     1 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     1     1     1     1     1     1     1     1     0     0     0     0 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     0     0     0     0];
lut = lut(:);

%---------------------------------------------------------------------
function lut = lutthin4()
%LUTTHIN4 Compute "thin4" look-up table.
%   Copyright 1993-2003 The MathWorks, Inc.  
%   $Revision: 1.6 $  $Date: 2003/01/17 16:28:38 $
lut = [ ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     1     1     0     1     1     1     0     0 ...
     1     1     0     1     1     1     0     0     0     0     0     0 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     1     1     0     0     1     1     0     0     1     1     0     0 ...
     1     1     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     0     0     0     0     1     1     0     1 ...
     1     1     0     0     1     1     0     1     1     1     0     0 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     1     1     0     0     1     1     0     0 ...
     1     1     0     0     1     1     0     0     0     0     0     0 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     1     1     0     1     1     1     0     0     1     1     0     1 ...
     1     1     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     0     0     0     0     1     1     0     0 ...
     1     1     0     0     1     1     0     0     1     1     0     0 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     1     1     0     1     1     1     0     0 ...
     1     1     0     1     1     1     0     0     0     0     0     0 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     1     1     0     0     1     1     0     0     1     1     0     0 ...
     1     1     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     0     0     0     0     1     1     0     1 ...
     1     1     0     0     1     1     0     1     1     1     0     0 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     1     1     0     0     1     1     0     0 ...
     1     1     0     0     1     1     0     0     0     0     0     0 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     1     1     0     1     1     1     0     0     1     1     0     1 ...
     1     1     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     0     0     0     0     1     1     0     0 ...
     1     1     0     0     1     1     0     0     1     1     0     0 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     1     1     0     1     1     1     0     0 ...
     1     1     0     1     1     1     0     0     0     0     0     0 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     1     1     0     0     1     1     0     0     1     1     0     0 ...
     1     1     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     0     0     0     0     1     1     0     1 ...
     1     1     0     0     1     1     0     1     1     1     0     0 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     1     1     0     0     1     1     0     0 ...
     1     1     0     0     1     1     0     0];
lut = lut(:);

%---------------------------------------------------------------------
function lut = lutper4()
%LUTPER4 Compute "perim4" look-up table.
%   Copyright 1993-2003 The MathWorks, Inc.  
%   $Revision: 1.11 $  $Date: 2003/01/17 16:28:35 $
lut = [ ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     1     1     1     1     1     1     1     1 ...
     1     1     1     1     1     1     1     1     0     0     0     0 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     1     1     1     1     1     1     1     1     1     1     1     1 ...
     1     1     1     1     0     0     0     0     0     0     0     0 ...
     0     0     0     0     0     0     0     0     1     1     1     1 ...
     1     1     1     1     1     1     1     1     1     1     1     1 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     1     1     1     1     1     1     1     1 ...
     1     1     1     1     1     1     1     1     0     0     0     0 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     1     1     1     1     1     1     1     1     1     1     1     1 ...
     1     1     1     1     0     0     0     0     0     0     0     0 ...
     0     0     0     0     0     0     0     0     1     1     1     1 ...
     1     1     1     1     1     1     0     0     1     1     0     0 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     1     1     1     1     1     1     1     1 ...
     1     1     1     1     1     1     1     1     0     0     0     0 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     1     1     1     1     1     1     1     1     1     1     0     0 ...
     1     1     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     0     0     0     0     1     1     1     1 ...
     1     1     1     1     1     1     1     1     1     1     1     1 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     1     1     1     1     1     1     1     1 ...
     1     1     1     1     1     1     1     1     0     0     0     0 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     1     1     1     1     1     1     1     1     1     1     1     1 ...
     1     1     1     1     0     0     0     0     0     0     0     0 ...
     0     0     0     0     0     0     0     0     1     1     1     1 ...
     1     1     1     1     1     1     1     1     1     1     1     1 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     1     1     1     1     1     1     1     1 ...
     1     1     1     1     1     1     1     1     0     0     0     0 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     1     1     1     1     1     1     1     1     1     1     0     0 ...
     1     1     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     0     0     0     0     1     1     1     1 ...
     1     1     1     1     1     1     1     1     1     1     1     1 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     1     1     1     1     1     1     1     1 ...
     1     1     0     0     1     1     0     0];
lut = lut(:);

%---------------------------------------------------------------------
function lut = lutper8()
%LUTPER8 Compute "perim8" look-up table.
%   Copyright 1993-2003 The MathWorks, Inc.  
lut = [ ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     1     1     1     1     1     1     1     1 ...
     1     1     1     1     1     1     1     1     0     0     0     0 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     1     1     1     1     1     1     1     1     1     1     1     1 ...
     1     1     1     1     0     0     0     0     0     0     0     0 ...
     0     0     0     0     0     0     0     0     1     1     1     1 ...
     1     1     1     1     1     1     1     1     1     1     1     1 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     1     1     1     1     1     1     1     1 ...
     1     1     1     1     1     1     1     1     0     0     0     0 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     1     1     1     1     1     1     1     1     1     1     1     1 ...
     1     1     1     1     0     0     0     0     0     0     0     0 ...
     0     0     0     0     0     0     0     0     1     1     1     1 ...
     1     1     1     1     1     1     1     1     1     1     1     1 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     1     1     1     1     1     1     1     1 ...
     1     1     1     1     1     1     1     1     0     0     0     0 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     1     1     1     1     1     1     1     1     1     1     1     1 ...
     1     1     1     1     0     0     0     0     0     0     0     0 ...
     0     0     0     0     0     0     0     0     1     1     1     1 ...
     1     1     1     1     1     1     1     1     1     1     1     1 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     1     1     1     1     1     1     1     1 ...
     1     1     1     1     1     1     1     1     0     0     0     0 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     1     1     1     1     1     1     1     1     1     1     1     1 ...
     1     1     1     1     0     0     0     0     0     0     0     0 ...
     0     0     0     0     0     0     0     0     1     1     1     1 ...
     1     1     1     1     1     1     1     1     1     1     1     1 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     1     1     1     1     1     1     1     1 ...
     1     1     1     1     1     1     1     1     0     0     0     0 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     1     1     1     1     1     1     1     1     1     1     1     1 ...
     1     1     1     1     0     0     0     0     0     0     0     0 ...
     0     0     0     0     0     0     0     0     1     1     1     1 ...
     1     1     1     1     1     1     1     1     1     1     1     1 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     1     1     1     1     1     1     1     1 ...
     1     1     1     1     1     1     1     0];
lut = lut(:);

%---------------------------------------------------------------------
function lut = lutclean()
%LUTCLEAN Compute "clean" look-up table.
%   Copyright 1993-2003 The MathWorks, Inc.  
%   $Revision: 1.11 $  $Date: 2003/01/17 16:28:34 $
lut = [ ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     0     1     1     1     1     1     1     1 ...
     1     1     1     1     1     1     1     1     0     0     0     0 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     1     1     1     1     1     1     1     1     1     1     1     1 ...
     1     1     1     1     0     0     0     0     0     0     0     0 ...
     0     0     0     0     0     0     0     0     1     1     1     1 ...
     1     1     1     1     1     1     1     1     1     1     1     1 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     1     1     1     1     1     1     1     1 ...
     1     1     1     1     1     1     1     1     0     0     0     0 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     1     1     1     1     1     1     1     1     1     1     1     1 ...
     1     1     1     1     0     0     0     0     0     0     0     0 ...
     0     0     0     0     0     0     0     0     1     1     1     1 ...
     1     1     1     1     1     1     1     1     1     1     1     1 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     1     1     1     1     1     1     1     1 ...
     1     1     1     1     1     1     1     1     0     0     0     0 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     1     1     1     1     1     1     1     1     1     1     1     1 ...
     1     1     1     1     0     0     0     0     0     0     0     0 ...
     0     0     0     0     0     0     0     0     1     1     1     1 ...
     1     1     1     1     1     1     1     1     1     1     1     1 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     1     1     1     1     1     1     1     1 ...
     1     1     1     1     1     1     1     1     0     0     0     0 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     1     1     1     1     1     1     1     1     1     1     1     1 ...
     1     1     1     1     0     0     0     0     0     0     0     0 ...
     0     0     0     0     0     0     0     0     1     1     1     1 ...
     1     1     1     1     1     1     1     1     1     1     1     1 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     1     1     1     1     1     1     1     1 ...
     1     1     1     1     1     1     1     1     0     0     0     0 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     1     1     1     1     1     1     1     1     1     1     1     1 ...
     1     1     1     1     0     0     0     0     0     0     0     0 ...
     0     0     0     0     0     0     0     0     1     1     1     1 ...
     1     1     1     1     1     1     1     1     1     1     1     1 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     1     1     1     1     1     1     1     1 ...
     1     1     1     1     1     1     1     1];
lut = lut(:);

%---------------------------------------------------------------------
function lut = lutbridge()
%LUTBRIDGE Compute "bridge" look-up table.
%   Copyright 1993-2003 The MathWorks, Inc.  
%   $Revision: 1.11 $  $Date: 2003/01/17 16:28:34 $
lut = [ ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     1     1     0     0     1     1     1     1     1     1     1     1 ...
     1     1     1     1     1     1     1     1     0     1     0     0 ...
     0     1     0     0     1     1     0     0     1     1     0     0 ...
     1     1     1     1     1     1     1     1     1     1     1     1 ...
     1     1     1     1     0     0     1     1     1     1     1     1 ...
     0     0     0     0     1     1     0     0     1     1     1     1 ...
     1     1     1     1     1     1     1     1     1     1     1     1 ...
     1     1     1     1     1     1     1     1     1     1     0     0 ...
     1     1     0     0     1     1     1     1     1     1     1     1 ...
     1     1     1     1     1     1     1     1     0     1     1     1 ...
     1     1     1     1     0     0     0     0     1     1     0     0 ...
     1     1     1     1     1     1     1     1     1     1     1     1 ...
     1     1     1     1     0     1     0     0     0     1     0     0 ...
     0     0     0     0     0     0     0     0     1     1     1     1 ...
     1     1     1     1     1     1     1     1     1     1     1     1 ...
     0     1     1     1     1     1     1     1     0     0     0     0 ...
     1     1     0     0     1     1     1     1     1     1     1     1 ...
     1     1     1     1     1     1     1     1     0     1     0     0 ...
     0     1     0     0     0     0     0     0     0     0     0     0 ...
     1     1     1     1     1     1     1     1     1     1     1     1 ...
     1     1     1     1     0     1     1     1     0     1     1     1 ...
     1     1     1     1     1     1     1     1     1     1     1     1 ...
     1     1     1     1     1     1     1     1     1     1     1     1 ...
     0     1     0     0     0     1     0     0     1     1     0     0 ...
     1     1     0     0     1     1     1     1     1     1     1     1 ...
     1     1     1     1     1     1     1     1     0     1     1     1 ...
     1     1     1     1     1     1     1     1     1     1     1     1 ...
     1     1     1     1     1     1     1     1     1     1     1     1 ...
     1     1     1     1     1     1     1     1     1     1     1     1 ...
     1     1     0     0     1     1     0     0     1     1     1     1 ...
     1     1     1     1     1     1     1     1     1     1     1     1 ...
     0     1     1     1     1     1     1     1     0     0     0     0 ...
     1     1     0     0     1     1     1     1     1     1     1     1 ...
     1     1     1     1     1     1     1     1     0     1     0     0 ...
     0     1     0     0     0     0     0     0     0     0     0     0 ...
     1     1     1     1     1     1     1     1     1     1     1     1 ...
     1     1     1     1     0     1     1     1     1     1     1     1 ...
     0     0     0     0     1     1     0     0     1     1     1     1 ...
     1     1     1     1     1     1     1     1     1     1     1     1 ...
     0     1     0     0     0     1     0     0     0     0     0     0 ...
     0     0     0     0     1     1     1     1     1     1     1     1 ...
     1     1     1     1     1     1     1     1];
lut = lut(:);

%---------------------------------------------------------------------
function lut = lutshrink()
%LUTSHRINK Compute "shrink" look-up table.
%   Copyright 1993-2003 The MathWorks, Inc.  
%   $Revision: 1.11 $  $Date: 2003/01/17 16:28:36 $
lut = [ ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     0     1     1     1     1     0     1     1 ...
     1     1     1     1     0     0     1     1     0     0     0     0 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     1     0     1     1     1     0     1     1     0     0     1     1 ...
     0     0     1     1     0     0     0     0     0     0     0     0 ...
     0     0     0     0     0     0     0     0     1     0     0     0 ...
     0     0     0     0     1     1     1     1     0     0     1     1 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     0     0     1     1     0     0     1     1     0     0     0     0 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     1     0     0     0     0     0     0     0     1     1     1     1 ...
     0     0     1     1     0     0     0     0     0     0     0     0 ...
     0     0     0     0     0     0     0     0     1     0     1     1 ...
     1     0     1     1     1     1     0     0     1     1     0     0 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     1     0     0     0     0     0     0     0 ...
     1     1     1     1     0     0     1     1     0     0     0     0 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     1     0     1     1     1     0     1     1     1     1     0     0 ...
     1     1     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     0     0     0     0     1     0     0     0 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     1     0     1     1     1     0     1     1 ...
     0     0     1     1     0     0     1     1     0     0     0     0 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     0     0     1     1     0     0     1     1 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     1     0     0     0     0     0     0     0 ...
     1     1     1     1     0     0     1     1     0     0     0     0 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     1     0     1     1     1     0     1     1     1     1     0     0 ...
     1     1     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     0     0     0     0     1     0     0     0 ...
     0     0     0     0     1     1     1     1     0     0     1     1 ...
     0     0     0     0     0     0     0     0     0     0     0     0 ...
     0     0     0     0     1     0     1     1     1     0     1     1 ...
     1     1     0     0     1     1     0     0]';

%---------------------------------------------------------------------
function S = decorrstretch(varargin)
%DECORRSTRETCH Apply a decorrelation stretch to a multichannel image.
%   Copyright 1993-2003 The MathWorks, Inc.
%   $Revision: 1.4 $  $Date: 2003/03/05 16:40:31 $

% Parse and validate input arguments.
[A, mode, targetMean, targetSigma, tol, rowsubs, colsubs] = parseInputs_decorrstretch(varargin{:});

% Convert to double, if necessary.
inputClass = class(A);
if ~strcmp(inputClass,'double')
    A = im2double(A);
    targetMean  = im2double(feval(inputClass,targetMean));
    targetSigma = im2double(feval(inputClass,targetSigma));
end

% Apply decorrelation stretch.
S = decorr(A, strcmp(mode,'correlation'), targetMean, targetSigma, rowsubs, colsubs);

% Apply optional contrast stretch.
if ~isempty(tol)
    low_high = stretchlim_j(S,tol);
    S = LocalImadjust(S,low_high);
    S(S < 0) = 0;
    S(S > 1) = 1;
end

% Restore input class.
S = changeclass(inputClass,S);

%--------------------------------------------------------------------------
function S = decorr(A, useCorr, targetMean, targetSigma, rowsubs, colsubs)
% Decorrelation stretch for a multiband image of class double.

[r c nbands] = size(A);        % Save the shape
npixels = r * c;                % Number of pixels
A = reshape(A,[npixels nbands]);     % Reshape to numPixels-by-numBands

if isempty(rowsubs)
    B = A;
else
    ind = sub2ind([r c], rowsubs, colsubs);
    B = A(ind,:);
end

meanB = mean(B,1);        % Mean pixel value in each spectral band
n = size(B,1);            % Equals npixels if rowsubs is empty
if n == 1
    cov = zeros(nbands);
else
    cov = (B' * B - n * meanB' * meanB)/(n - 1);  % Sample covariance matrix
end

[T, offset]  = fitdecorrtrans(meanB, cov, useCorr, targetMean, targetSigma);

S = A * T + repmat(offset,[npixels 1]);
S = reshape(S, [r c nbands]);

%--------------------------------------------------------------------------
function out = LocalImadjust(img,low_high)
% A short, specialized version of IMADJUST that works with
% an arbitrary number of image planes.

low  = low_high(1,:);
high = low_high(2,:);
out = zeros(size(img));

% Loop over image planes and perform transformation.
for p = 1:size(img,3)
    % Make sure img is in the range [low,high].
    img(:,:,p) =  max(low(p),min(high(p),img(:,:,p)));

    % Transform.
    out(:,:,p) = (img(:,:,p)-low(p))./(high(p)-low(p));
end

%--------------------------------------------------------------------------
function [A, mode, targetMean, targetSigma, tol, rowsubs, colsubs] = ...
    parseInputs_decorrstretch(varargin)

% Number of arguments passed to DECORRSTRETCH
checknargin(1, 12, nargin, 'decorrstretch');

% Defaults.
A = [];
mode = 'correlation';
targetMean = [];    targetSigma = [];
rowsubs = [];       colsubs = [];
tol = [];

% Validate the image array.
A = varargin{1};

checkinput(A, {'double','uint8','uint16'},{'nonempty','real','nonnan','finite'},...
              'decorrstretch', 'A', 1);
          
if ndims(A) > 3
    error('Images:decorrstretch:expected2Dor3D','A has more than three dimensions.');
end

% Validate the property name-value pairs.
nbands = size(A,3);
validPropertyNames = {'Mode','TargetMean','TargetSigma','Tolerance','SampleSubscripts'};
for k = 2:2:nargin
    propName = checkstrs(varargin{k}, validPropertyNames, mfilename, 'PARAM', k);
    switch propName
      case 'Mode'
        checkExistence(k, nargin, mfilename, 'mode string', propName);
        mode = checkstrs(varargin{k+1}, {'correlation','covariance'},...
                         mfilename, 'mode', k+1);
      case 'TargetMean'
        checkExistence(k, nargin, mfilename, 'target mean vector', propName);
        targetMean = checkTargetMean(varargin{k+1}, nbands, k+1);
      case 'TargetSigma'
        checkExistence(k, nargin, mfilename, 'target sigma vector', propName);
        targetSigma = checkTargetSigma(varargin{k+1}, nbands, k+1);
      case 'Tolerance'
        checkExistence(k, nargin, mfilename, 'linear stretch fraction(s)', propName);
        tol = checkTol(varargin{k+1}, k+1);
      case 'SampleSubscripts'
        checkExistence(k, nargin, mfilename, 'pixel subscript arrays', propName);
        [rowsubs, colsubs] = checkSubs(varargin{k+1}, size(A,1), size(A,2), k+1);
      otherwise
        error('Images:decorrstretch:internalProblem',...
           'Internal problem: unrecognized property name: %s',propName);
    end
end

%--------------------------------------------------------------------------
function checkExistence(position, nargs, fcnName, propertyDescription, propertyName)
% Error if missing the property value following a property name.

if (position + 1 > nargs)
    error(sprintf('Images:%s:missingParameterValue',fcnName),...
          'Expected %s to follow parameter name ''%s''.',...
          propertyDescription, propertyName);
end

%--------------------------------------------------------------------------
function tol = checkTol(tol, position)

% Validate the linear-stretch tolerance.
checkinput(tol, {'double'}, {'nonempty','real','nonnan','nonnegative','finite'},...
                mfilename, 'TOL', position);

if any(tol < 0) || any(tol > 1)
    error('Images:decorrstretch:tolOutOfRange','Elements of TOL must be in the range [0 1].');
end

n = numel(tol);
if n > 2
    error('Images:decorrstretch:tolHasTooManyElements','TOL must have 1 or 2 elements.');
end

if (n == 2) && ~(tol(1) < tol(2))
    error('Images:decorrstretch:tolNotIncreasing','TOL(1) must be less than TOL(2).');
end

if (n == 1) && (tol(1) >= 0.5)
    error('Images:decorrstretch:tolOutOfRange','Scalar TOL must be in the range [0 0.5).');
end

if (n == 1),    tol = [tol 1-tol];  end

%--------------------------------------------------------------------------
function targetMean = checkTargetMean(targetMean, nbands, position)

checkinput(targetMean, {'double'}, {'nonempty','real','nonnan','finite','vector'},...
                mfilename, 'TARGETMEAN', position);

targetMean = targetMean(:)';  % Make sure it's a row vector.

if (numel(targetMean) > 1) && (size(targetMean,2) ~= nbands)
    error('Images:decorrstretch:targetMeanWrongSize',...
          'Expected TARGETMEAN to be a scalar or vector of length %d.', nbands);
end
  
%--------------------------------------------------------------------------
function targetSigma = checkTargetSigma(targetSigma, nbands, position)

checkinput(targetSigma, {'double'}, {'nonnegative','nonempty','real','nonnan','finite','vector'},...
                mfilename, 'TARGETSIGMA', position);

targetSigma = targetSigma(:)';  % Make sure it's a row vector.

if (numel(targetSigma) > 1) && (size(targetSigma,2) ~= nbands)
    error('Images:decorrstretch:targetSigmaWrongSize',...
          'Expected TARGETSIGMA to be a scalar or vector of length %d.', nbands, nbands);
end

% Convert to a diagonal matrix for convenient computation.
targetSigma = diag(targetSigma);

%--------------------------------------------------------------------------
function [rowsubs, colsubs] = checkSubs(subscell, nrows, ncols, position)

if ~iscell(subscell) || numel(subscell) ~= 2
    error('Images:decorrstretch:sampleSubsNotTwoElementCell',...
          'SAMPLESUBS must be a two-element cell array.');
end
    
rowsubs = subscell{1}(:);   colsubs = subscell{2}(:);

checkinput(rowsubs, {'double'},{'nonempty','integer','positive'},...
                mfilename, 'ROWSUBS', position);

checkinput(colsubs, {'double'},{'nonempty','integer','positive'},...
                mfilename, 'COLSUBS', position);

if any(rowsubs > nrows)
    error('Images:decorrstretch:subscriptsOutOfRange',...
          'Sample row subscripts must be in the range [1 SIZE(A,1)].');
end

if any(colsubs > ncols)
    error('Images:decorrstretch:subscriptsOutOfRange',...
          'Sample column subscripts must be in the range [1 SIZE(A,2)].');
end

if numel(rowsubs) ~= numel(colsubs)
    error('Images:decorrstretch:subscriptArraySizeMismatch',...
          'Sample row and column subscript arrays must be the same size.');
end

% ----------------------------------------------------------------------------------
function [T, offset] = fitdecorrtrans(means, Cov, useCorr, targetMean, targetSigma)
% FITDECORRTRANS   Fit decorrelating transformation to image statistics.
% Copyright 2002 The MathWorks, Inc.

% Square-root variances in a diagonal matrix.
S = diag(sqrt(diag(Cov)));  

if isempty(targetSigma)    % Restore original sample variances.
    targetSigma = S;
end

if useCorr
    Corr = pinv(S) * Cov * pinv(S);
    Corr(logical(eye(size(Corr,1)))) = 1;
    [V D] = eig(Corr);
    T = pinv(S) * V * decorrWeight(D) * V' * targetSigma;
else
    [V D] = eig(Cov);
    T = V * decorrWeight(D) * V' * targetSigma;
end

% Get the output variances right even for correlated bands, except
% for zero-variance bands---which can't be stretched at all.
T = T * pinv(diag(sqrt(diag(T' * Cov * T)))) * targetSigma;

if isempty(targetMean)    % Restore original sample means.
    targetMean = means;
end

offset = targetMean - means * T;

%--------------------------------------------------------------------------
function W = decorrWeight(D)
% Given the diagonal eigenvalue matrix D, compute the decorrelating
% weights W.  In the full rank, well-conditioned case, decorrWeight(D)
% returns the same result as sqrt(inv(D)).  In addition, it provides
% a graceful way to handle rank-deficient or near-rank-deficient
% (ill-conditioned) cases resulting from situations of perfect or
% near-perfect band-to-band correlation and/or bands with zero variance.

D(D < 0) = 0;   W = sqrt(pinv(D));

%--------------------------------------------------------------------------
function S = pinv(D)
% Pseudoinverse of a diagonal matrix, with a larger-than-standard
% tolerance to help in handling edge cases.  We've provided our
% own in order to: (1) Avoid replacing all calls to PINV with calls to
% PINV(...,TOL) and (2) Take advantage of the fact that our input is
% always diagonal so we don't need to call SVD.

d = diag(D);
tol =length(d) * max(d) * sqrt(eps);
keep = d > tol;
s = ones(size(d));
s(keep) = s(keep) ./ d(keep);
s(~keep) = 0;
S = diag(s);

% --------------------------------------------------------------------------
function varargout = imresize(varargin)
%IMRESIZE Resize image.
%   IMRESIZE resizes an image of any type using the specified
%   interpolation method. Supported interpolation methods
%   include:
%
%        'nearest'  (default) nearest neighbor interpolation
%
%        'bilinear' bilinear interpolation
%
%        'bicubic'  bicubic interpolation
%
%   B = IMRESIZE(A,M,METHOD) returns an image that is M times the
%   size of A. If M is between 0 and 1.0, B is smaller than A. If
%   M is greater than 1.0, B is larger than A. If METHOD is
%   omitted, IMRESIZE uses nearest neighbor interpolation.
%
%   B = IMRESIZE(A,[MROWS MCOLS],METHOD) returns an image of size
%   MROWS-by-MCOLS. If the specified size does not produce the
%   same aspect ratio as the input image has, the output image is
%   distorted.
%
%   When the specified output size is smaller than the size of
%   the input image, and METHOD is 'bilinear' or 'bicubic',
%   IMRESIZE applies a lowpass filter before interpolation to
%   reduce aliasing. The default filter size is 11-by-11.
%
%   You can specify a different length for the default filter
%   using:
%
%        [...] = IMRESIZE(...,METHOD,N)
%
%   N is an integer scalar specifying the size of the filter,
%   which is N-by-N. If N is 0, IMRESIZE omits the filtering
%   step.
%
%   You can also specify your own filter H using:
%
%        [...] = IMRESIZE(...,METHOD,H)
%
%   H is any two-dimensional FIR filter (such as those returned
%   by FTRANS2, FWIND1, FWIND2, or FSAMP2).
%
%   Class Support
%   -------------
%   The input image A can be numeric or logical and it must be
%   nonsparse. The output image is of the same class as the
%   input image.
%
%   Example
%   -------
%        I = imread('rice.png');
%        J = imresize(I,.5);
%        figure, imshow(I), figure, imshow(J)
%
%   See also IMROTATE, IMTRANSFORM, TFORMARRAY.

%   Obsolete Syntaxes:
%
%   [R1,G1,B1] = IMRESIZE(R,G,B,M,'method') or
%   [R1,G1,B1] = IMRESIZE(R,G,B,[MROWS NCOLS],'method') resizes
%   the RGB image in the matrices R,G,B.  'bilinear' is the
%   default interpolation method.

%   Copyright 1992-2005 The MathWorks, Inc.
%   $Revision: 5.30.4.7 $  $Date: 2005/06/20 03:07:00 $

[A,m,method,h] = parse_inputs_imresize(varargin{:});

% Preserve classes
inputClass = class(A);      classChanged = 0;   logicalIn = islogical(A);

% Define old and new image sizes, and actual scaling
% sc is a two-element vector: [vert_scale_factor, horiz_scale_factor].
[so(1),so(2),thirdD] = size(A); % old image size
if (length(m) == 1)     % m is the scale factor.
    sn = max(floor(m*so(1:2)),1); % new image size=(integer>0)
    sc = [m m];
else                    % m is new image size
    sn = m;
    sc = sn ./ so;
end

if switch_to_nearest_method(sn, so, method)
    warning('Images:imresize:inputTooSmall', ...
            'Input is too small for bilinear or bicubic method; using nearest-neighbor method instead.');
    method = 'nearest';
end

% Filtering is under the following conditions
bi_interp = (method(1)=='b'); % non-default interpolation only
defflt_reducedim=(length(h)<2)&any(sn<so);%default filter & reduced image
if length(h)==1,
    nonzero_odr = (h~=0);       % non-zero filter order
else
    nonzero_odr = 1;
end;
custm_flt = (length(h)>1);      %custom supplied filter H

if bi_interp && nonzero_odr && any([defflt_reducedim,custm_flt]),
    if (~isa(A,'double')),      %change format to double to perform imfilter
        A = im2double(A);
        classChanged = 1;
    end

    if defflt_reducedim,        % Design anti-aliasing filter for reduced image
        drec = find(sn<so);     % find direction of filtering
        for k = drec,           % create filter for drec-direction
            if isempty(h),      % make filter order corresponding to scale
                h = 11;
            end
            hh(k,:) = DesignFilter(h,sn(k)/so(k));
        end
        if length(drec)==1,%filters in one direction only
            % first direction is column, second is row
            h = reshape(hh(k,:),(h-1)*(k==1)+1,(h-1)*(k==2)+1);
        else % filters in both directions
            for k=1:thirdD,%loop if A matrix is 3D
                A(:,:,k) = imfilter(imfilter(A(:,:,k), hh(2,:),'replicate'), hh(1,:).','replicate');
            end
        end
    end
    if custm_flt || (defflt_reducedim && (length(drec)==1)), % filters in one direction
        for k=1:thirdD,%loop if A matrix is 3D
            A(:,:,k) = imfilter(A(:,:,k),h,'replicate');
        end
    end
end

% Construct an affine tform that:
%   *  maps (u,v) = (0.5,0.5) in input space to (x,y) = (0.5,0.5) in output space.
%
%   *  maps (u,v) = (1.5,1.5) in input space to (x,y) =
%      (0.5+sc(2),0.5+sc(1)) in output space.

a = [sc(2),         0,                0
     0,             sc(1),            0
     0.5*(1-sc(2)), 0.5*(1-sc(1)),    1];
T = transform_fun('maketform','affine', a);

% Interpolation
if (method(1)=='n')     % nearest neighbor (default)
    subscripts = repmat({':'}, [1 ndims(A)]);
    X = [(1:sn(2)).', ones(sn(2), 1)];
    U = transform_fun('tforminv',T, X);
    c = min(round(U(:,1)), so(2));
    
    X = [ones(sn(1), 1), (1:sn(1)).'];
    U = transform_fun('tforminv',T, X);
    r = min(round(U(:,2)), so(1));
    
    subscripts{1} = r;
    subscripts{2} = c;
    A = A(subscripts{:});
else            % bilinear or bicubic
    if strcmp(method,'bicubic')
        R = transform_fun('makeresampler','cubic','replicate');
    else
        R = transform_fun('makeresampler','linear','replicate');
    end

    % In the construction of the affine transform matrix above, the first coordinate
    % lies along the horizontal dimensions, and the second coordinate lies along
    % the vertical dimension.  However, the default coordinate convention for
    % tformarray follows MATLAB's array indexing ordering. That is, the first
    % spatial transform dimension corresponds to the rows, the second spatial
    % transform dimension corresponds to the columns, etc. That's the reason for
    % specifying [2 1], [2 1], and [sn(2) sn(1)] in the call to tformarray below.
    A = transform_fun('tformarray',A, T, R, [2 1], [2 1], [sn(2) sn(1)], [], []);
end

% Change format from double back to the original
if logicalIn,  % output should be logical (i.e. binary image)
    if ~islogical(A) % A became double because of imfilter, turn it back to logical
        A = A > .5;
    end
elseif classChanged,
    A = changeclass(inputClass, A);
end

% Output
varargout{1} = A;

% ----------------------------------------------------------------------
function [A,m,method,h] = parse_inputs_imresize(varargin)
% Outputs:  A       the input image
%           m       the resize scaling factor or the new size
%           method  interpolation method (nearest,bilinear,bicubic)
%           h       if 0, skip filtering; if non-zero scalar, use filter
%                   of size h; if empty, use filter of size 11;
%                   otherwise h is the anti-aliasing filter provided by user
% Defaults:
method = 'nearest';     h = [];

error(nargchk(2,6,nargin));
switch nargin
    case 2,                     % imresize(A,m)
        A = varargin{1};        m = varargin{2};
    case 3,                     % imresize(A,m,method)
        A = varargin{1};        m = varargin{2};
        method = varargin{3};
    case 4,                     % imresize(A,m,method,h)
        A = varargin{1};        m = varargin{2};
        method = varargin{3};   h = varargin{4};
    otherwise,
        error('Images:imresize:invalidInputs', '%s', 'Invalid input arguments.');
end

checkinput(A,{'numeric', 'logical'},{'nonsparse'},mfilename,'A',1);

% Check validity of the input parameters
if isempty(m) || (ndims(m)>2) || any(m<=0) || length(m(:))>2,
    eid = 'Images:imresize:invalidScaleFactor';
    error(eid, '%s', 'M must be either a scalar multiplier or a 1-by-2 size vector.');
elseif length(m)==2,% make sure that m is a row of non-negative integers
    m = ceil(m(:).');
end

if ischar(method)
    strings = {'nearest','bilinear','bicubic'};
    idx = find(strncmpi(method, strings, numel(method)));
    if isempty(idx)
        eid = 'Images:imresize:unrecognizedInterpolationMethod';
        error(eid, 'Unknown interpolation method: %s', method);
    elseif (length(idx) > 1)
        eid = 'Images:imresize:ambiguousInterpolationMethod';
        error(eid, '%s', 'Ambiguous interpolation method: %s',method);
    else
        method = strings{idx};
    end
else
    eid = 'Images:imresize:expectedString';
    error(eid, '%s', 'Interpolation method has to be a string.');
end

if (length(h) == 1)        % represents filter order
    if (h<0) || (h~=round(h)),
        eid = 'Images:imresize:invalidFilterOrder';
        error(eid, 'Filter order has to be a non-negative integer, not %g',h);
    end
elseif (length(h)>1) && (ndims(h)>2)    % custom supplied filter
    eid = 'Images:imresize:expected2DFilter';
    error(eid, '%s', 'Filter has to be a 2-D array.');
end

% ----------------------------------------------------------------------
function b = DesignFilter(N,Wn)
% Modified from SPT v3 fir1.m and hanning.m
% first creates only first half of the filter
% and later mirrows it to the other half

odd = rem(N,2);
vec = 1:floor(N/2);
vec2 = pi*(vec-(1-odd)/2);

wind = .54-.46*cos(2*pi*(vec-1)/(N-1));
b = [fliplr(sin(Wn*vec2)./vec2).*wind Wn];      % first half is ready
b = b([vec floor(N/2)+(1:odd) fliplr(vec)]);    % entire filter
b = b/abs(polyval(b,1));                        % norm

% ----------------------------------------------------------------------
function tf = switch_to_nearest_method(new_size, old_size, method)
% Returns true if method is not 'nearest' but 'nearest' must be used
% because of the input and output image sizes.

tf = any(new_size < 4) && any(new_size < old_size) && (method(1) ~= 'n');

% -----------------------------------------------------------------------
function outstats = regionprops(varargin)
%REGIONPROPS Measure properties of image regions.
%   STATS = REGIONPROPS(L,PROPERTIES) measures a set of properties for each
%   labeled region in the label matrix L. Positive integer elements of L
%   correspond to different regions. For example, the set of elements of L
%   equal to 1 corresponds to region 1; the set of elements of L equal to 2
%   corresponds to region 2; and so on. STATS is a structure array of length
%   max(L(:)). The fields of the structure array denote different properties
%   for each region, as specified by PROPERTIES.
%
%   PROPERTIES can be a comma-separated list of strings, a cell array
%   containing strings, the string 'all', or the string 'basic'. The set of
%   valid measurement strings includes:
%
%     'Area'              'ConvexHull'    'EulerNumber'
%     'Centroid'          'ConvexImage'   'Extrema'       
%     'BoundingBox'       'ConvexArea'    'EquivDiameter' 
%     'SubarrayIdx'       'Image'         'Solidity'      
%     'MajorAxisLength'   'PixelList'     'Extent'        
%     'MinorAxisLength'   'PixelIdxList'  'FilledImage'  
%     'Orientation'                       'FilledArea'                   
%     'Eccentricity'                       
%                                                         
%   Property strings are case insensitive and can be abbreviated.
%
%   If PROPERTIES is the string 'all', then all of the above measurements
%   are computed. If PROPERTIES is not specified or if it is the string
%   'basic', then these measurements are computed: 'Area', 'Centroid', and
%   'BoundingBox'.
%
%   Note - REGIONPROPS and binary images
%   ------------------------------------
%   REGIONPROPS does not accept a binary image as its first input.  There
%   are two common ways to convert a binary image to a label matrix:
%
%       1.  L = bwlabel(BW);
%
%       2.  L = double(BW);
%
%   Suppose that BW were a logical matrix containing these values:
%
%       1 1 0 0 0 0
%       1 1 0 0 0 0
%       0 0 0 0 0 0
%       0 0 0 0 1 1
%       0 0 0 0 1 1
%
%   The first method of forming a label matrix, L = bwlabel(BW), results
%   in a label matrix containing two contiguous regions labeled by the
%   integer values 1 and 2.  The second method of forming a label matrix,
%   L = double(BW), results in a label matrix containing one
%   discontiguous region labeled by the integer value 1.  Since each
%   result is legitimately desirable in certain situations, REGIONPROPS
%   does not accept binary images and convert them using either method.
%   You should convert a binary image to a label matrix using one of
%   these methods (or another method if appropriate) before calling
%   REGIONPROPS.
%
%   Example
%   -------
%   Label the connected pixel components in the text.png image, compute
%   their centroids, and superimpose the centroid locations on the
%   image.
%
%       bw = imread('text.png');
%       L = bwlabel(bw);
%       s  = regionprops(L, 'centroid');
%       centroids = cat(1, s.Centroid);
%       imshow(bw)
%       hold on
%       plot(centroids(:,1), centroids(:,2), 'b*')
%       hold off
%
%   Class Support
%   -------------
%   The input label matrix L can have any numeric class.
%
%   See also BWLABEL, BWLABELN, ISMEMBER, WATERSHED.

%   Copyright 1993-2003 The MathWorks, Inc.  
%   $Revision: 1.10 $  $Date: 2003/03/03 14:26:34 $

officialStats = {'Area'
                 'Centroid'
                 'BoundingBox'
                 'SubarrayIdx'
                 'MajorAxisLength'
                 'MinorAxisLength'
                 'Eccentricity'
                 'Orientation'
                 'ConvexHull'
                 'ConvexImage'
                 'ConvexArea'
                 'Image'
                 'FilledImage'
                 'FilledArea'
                 'EulerNumber'
                 'Extrema'
                 'EquivDiameter'
                 'Solidity'
                 'Extent'
                 'PixelIdxList'
                 'PixelList'};
tempStats = {'PerimeterCornerPixelList'};
allStats = [officialStats ; tempStats];

[L, requestedStats, msg] = ParseInputs_regionprops(officialStats, varargin{:});
if (~isempty(msg))
    error(msg);
end

if ndims(L) > 2
    % Remove stats that aren't supported for N-D input and issue warning messages as appropriate.
    requestedStats = PreprocessRequestedStats(requestedStats);
end

if isempty(requestedStats)
    error('No input properties')
end

if (isempty(L))
    numObjs = 0;
else
    numObjs = round(double(max(L(:))));
end

% Initialize the stats structure array.
numStats = length(allStats);
empties = cell(numStats, numObjs);
stats = cell2struct(empties, allStats, 1);

% Initialize the computedStats structure array.
zz = cell(numStats, 1);
for k = 1:numStats
    zz{k} = 0;
end
computedStats = cell2struct(zz, allStats, 1);

for k = 1:length(requestedStats)
	switch requestedStats{k}
		case 'Area',		[stats, computedStats] = ComputeArea(L, stats, computedStats);
		case 'FilledImage',	[stats, computedStats] = ComputeFilledImage(L,stats,computedStats);
		case 'FilledArea',	[stats, computedStats] = ComputeFilledArea(L,stats,computedStats);
		case 'ConvexArea',	[stats, computedStats] = ComputeConvexArea(L, stats, computedStats);
		case 'Centroid',	[stats, computedStats] = ComputeCentroid(L, stats, computedStats);
		case 'EulerNumber',	[stats, computedStats] = ComputeEulerNumber(L,stats,computedStats);
		case 'EquivDiameter',[stats, computedStats] = ComputeEquivDiameter(L, stats, computedStats);
		case 'Extrema',		[stats, computedStats] = ComputeExtrema(L, stats, computedStats);
		case 'BoundingBox',	[stats, computedStats] = ComputeBoundingBox(L, stats, computedStats);
		case 'SubarrayIdx',	[stats, computedStats] = ComputeSubarrayIdx(L, stats, computedStats);
		case {'MajorAxisLength', 'MinorAxisLength', 'Orientation', 'Eccentricity'}
            [stats, computedStats] = ComputeEllipseParams(L, stats, computedStats);
		case 'Solidity',	[stats, computedStats] = ComputeSolidity(L, stats, computedStats);
		case 'Extent',		[stats, computedStats] = ComputeExtent(L, stats, computedStats);
		case 'ConvexImage',	[stats, computedStats] = ComputeConvexImage(L, stats, computedStats);
		case 'ConvexHull',	[stats, computedStats] = ComputeConvexHull(L, stats, computedStats);
		case 'Image',		[stats, computedStats] = ComputeImage(L, stats, computedStats);
		case 'PixelIdxList',[stats, computedStats] = ComputePixelIdxList(L, stats, computedStats);
		case 'PixelList',	[stats, computedStats] = ComputePixelList(L, stats, computedStats);
	end
end

% Initialize the output stats structure array.
numStats = length(requestedStats);
empties = cell(numStats, numObjs);
outstats = cell2struct(empties, requestedStats, 1);

% Initialize the subsref structure.
s(1).type = '()';
s(1).subs = {};
s(2).type = '.';
s(2).subs = '';

% Copy only the requested stats into the output.
v = version;	v = str2double(v(1));
if (v < 9)		% Later ML fck buggued with the subsref call bellow
	for k = 1:numObjs
		for p = 1:length(requestedStats)
			s(1).subs = {k};
			s(2).subs = {requestedStats{p}};

			% In normal MATLAB syntax, the line below is the same as:
			%    outstats(k).fieldname = stats(k).fieldname
			% where fieldname is the string contained in
			% requestedStats{p}.  If you don't give subsasgn an
			% output argument it changes its first input argument in-place.
			outstats = subsasgn(outstats, s, subsref(stats, s));
		end
	end
else
	for k = 1:numObjs
		for p = 1:length(requestedStats)
			outstats(k).(requestedStats{p}) = stats(k).(requestedStats{p});
		end
	end
end
    
%%%
%%% ComputeArea
%%%
function [stats, computedStats] = ComputeArea(L, stats, computedStats)
%   The area is defined to be the number of pixels belonging to
%   the region.

if ~computedStats.Area
    computedStats.Area = 1;

    [stats, computedStats] = ComputePixelIdxList(L, stats, computedStats);

    for k = 1:length(stats)
        stats(k).Area = size(stats(k).PixelIdxList, 1);
    end
end

%%%
%%% ComputeEquivDiameter
%%%
function [stats, computedStats] = ComputeEquivDiameter(L, stats, computedStats)
%   Computes the diameter of the circle that has the same area as
%   the region.
%   Ref: Russ, The Image Processing Handbook, 2nd ed, 1994, page
%   511.

if ~computedStats.EquivDiameter
    computedStats.EquivDiameter = 1;
    
    if ndims(L) > 2
		NoNDSupport('EquivDiameter');		return
    end
        
    [stats, computedStats] = ComputeArea(L, stats, computedStats);

    factor = 2/sqrt(pi);
    for k = 1:length(stats)
        stats(k).EquivDiameter = factor * sqrt(stats(k).Area);
    end
end

%%%
%%% ComputeFilledImage
%%%
function [stats, computedStats] = ComputeFilledImage(L,stats,computedStats)
%   Uses imfill to fill holes in the region.

if ~computedStats.FilledImage
    computedStats.FilledImage = 1;
    
    [stats, computedStats] = ComputeImage(L, stats, computedStats);
    
    conn = conndef(ndims(L),'minimal');
    
    for k = 1:length(stats)
        stats(k).FilledImage = imfill(stats(k).Image,conn,'holes');
    end
end

%%%
%%% ComputeConvexArea
%%%
function [stats, computedStats] = ComputeConvexArea(L, stats, computedStats)
%   Computes the number of "on" pixels in ConvexImage.

if ~computedStats.ConvexArea
    computedStats.ConvexArea = 1;
    
    if ndims(L) > 2
        NoNDSupport('ConvexArea');		return
    end
        
    [stats, computedStats] = ComputeConvexImage(L, stats, computedStats);
    
    for k = 1:length(stats)
        stats(k).ConvexArea = sum(stats(k).ConvexImage(:));
    end
end

%%%
%%% ComputeFilledArea
%%%
function [stats, computedStats] = ComputeFilledArea(L,stats,computedStats)
%   Computes the number of "on" pixels in FilledImage.

if ~computedStats.FilledArea
    computedStats.FilledArea = 1;
    
    [stats, computedStats] = ComputeFilledImage(L,stats,computedStats);

    for k = 1:length(stats)
        stats(k).FilledArea = sum(stats(k).FilledImage(:));
    end
end

%%%
%%% ComputeConvexImage
%%%
function [stats, computedStats] = ComputeConvexImage(L, stats, computedStats)
%   Uses ROIPOLY to fill in the convex hull.

if ~computedStats.ConvexImage
    computedStats.ConvexImage = 1;
    
    if ndims(L) > 2
        NoNDSupport('ConvexImage');
        return
    end
        
    [stats, computedStats] = ComputeConvexHull(L, stats, computedStats);
    [stats, computedStats] = ComputeBoundingBox(L, stats, computedStats);
    
    for k = 1:length(stats)
        M = stats(k).BoundingBox(4);
        N = stats(k).BoundingBox(3);
        hull = stats(k).ConvexHull;
        if (isempty(hull))
            stats(k).ConvexImage = false(M,N);
        else
            firstRow = stats(k).BoundingBox(2) + 0.5;
            firstCol = stats(k).BoundingBox(1) + 0.5;
            r = hull(:,2) - firstRow + 1;
            c = hull(:,1) - firstCol + 1;
            stats(k).ConvexImage = roipoly_j(M, N, c, r);
        end
    end
end

%%%
%%% ComputeCentroid
%%%
function [stats, computedStats] = ComputeCentroid(L, stats, computedStats)
%   [mean(r) mean(c)]

if ~computedStats.Centroid
    computedStats.Centroid = 1;
    
    [stats, computedStats] = ComputePixelList(L, stats, computedStats);

    
    % Save the warning state and disable warnings to prevent divide-by-zero
    % warnings.
    state = warning;
    warning off;
    
    for k = 1:length(stats)
        stats(k).Centroid = mean(stats(k).PixelList,1);
    end
    
    % Restore the warning state.
    warning(state);
end

%%%
%%% ComputeEulerNumber
%%%
function [stats, computedStats] = ComputeEulerNumber(L,stats,computedStats)
%   Calls BWEULER on 'Image' using 8-connectivity

if ~computedStats.EulerNumber
    computedStats.EulerNumber = 1;
    
    if ndims(L) > 2
        NoNDSupport('EulerNumber');
        return
    end
    
    [stats, computedStats] = ComputeImage(L, stats, computedStats);
    
    for k = 1:length(stats)
        stats(k).EulerNumber = bweuler(stats(k).Image,8);
    end
end

%%%
%%% ComputeExtrema
%%%
function [stats, computedStats] = ComputeExtrema(L, stats, computedStats)
%   A 8-by-2 array; each row contains the x and y spatial
%   coordinates for these extrema:  leftmost-top, rightmost-top,
%   topmost-right, bottommost-right, rightmost-bottom, leftmost-bottom,
%   bottommost-left, topmost-left. 
%   reference: Haralick and Shapiro, Computer and Robot Vision
%   vol I, Addison-Wesley 1992, pp. 62-64.

if ~computedStats.Extrema
    computedStats.Extrema = 1;
    
    if ndims(L) > 2
        NoNDSupport('Extrema');
        return
    end
        
    [stats, computedStats] = ComputePixelList(L, stats, computedStats);
    
    for k = 1:length(stats)
        pixelList = stats(k).PixelList;
        if (isempty(pixelList))
            stats(k).Extrema = zeros(8,2) + 0.5;
        else
            r = pixelList(:,2);
            c = pixelList(:,1);
            
            minR = min(r);
            maxR = max(r);
            minC = min(c);
            maxC = max(c);
            
            minRSet = r==minR;
            maxRSet = r==maxR;
            minCSet = c==minC;
            maxCSet = c==maxC;

            % Points 1 and 2 are on the top row.
            r1 = minR;
            r2 = minR;
            % Find the minimum and maximum column coordinates for top-row pixels.
            tmp = c(minRSet);
            c1 = min(tmp);
            c2 = max(tmp);
            
            % Points 3 and 4 are on the right column.
            % Find the minimum and maximum row coordinates for right-column pixels.
            tmp = r(maxCSet);
            r3 = min(tmp);
            r4 = max(tmp);
            c3 = maxC;
            c4 = maxC;

            % Points 5 and 6 are on the bottom row.
            r5 = maxR;
            r6 = maxR;
            % Find the minimum and maximum column coordinates for bottom-row pixels.
            tmp = c(maxRSet);
            c5 = max(tmp);
            c6 = min(tmp);
            
            % Points 7 and 8 are on the left column.
            % Find the minimum and maximum row coordinates for left-column pixels.
            tmp = r(minCSet);
            r7 = max(tmp);
            r8 = min(tmp);
            c7 = minC;
            c8 = minC;
            
            stats(k).Extrema = [c1-0.5 r1-0.5
                c2+0.5 r2-0.5
                c3+0.5 r3-0.5
                c4+0.5 r4+0.5
                c5+0.5 r5+0.5
                c6-0.5 r6+0.5
                c7-0.5 r7+0.5
                c8-0.5 r8-0.5];
        end
    end
    
end
        
%%%
%%% ComputeBoundingBox
%%%
function [stats, computedStats] = ComputeBoundingBox(L, stats, computedStats)
%   [minC minR width height]; minC and minR end in .5.

if ~computedStats.BoundingBox
    computedStats.BoundingBox = 1;
    
    [stats, computedStats] = ComputePixelList(L, stats, computedStats);
    
    num_dims = ndims(L);
    
    for k = 1:length(stats)
        list = stats(k).PixelList;
        if (isempty(list))
            stats(k).BoundingBox = [0.5*ones(1,num_dims) zeros(1,num_dims)];
        else
            min_corner = min(list,[],1) - 0.5;
            max_corner = max(list,[],1) + 0.5;
            stats(k).BoundingBox = [min_corner (max_corner - min_corner)];
        end
    end
end

%%%
%%% ComputeSubarrayIdx
%%%
function [stats, computedStats] = ComputeSubarrayIdx(L, stats, computedStats)
%   Find a cell-array containing indices so that L(idx{:}) extracts the
%   elements of L inside the bounding box.

if ~computedStats.SubarrayIdx
    computedStats.SubarrayIdx = 1;
    
    [stats, computedStats] = ComputeBoundingBox(L, stats, computedStats);
    num_dims = ndims(L);
    idx = cell(1,num_dims);
    for k = 1:length(stats)
        boundingBox = stats(k).BoundingBox;
        left = boundingBox(1:(end/2));
        right = boundingBox((1+end/2):end);
        left = left(1,[2 1 3:end]);
        right = right(1,[2 1 3:end]);
        for p = 1:num_dims
            first = left(p) + 0.5;
            last = first + right(p) - 1;
            idx{p} = first:last;
        end
        stats(k).SubarrayIdx = idx;
    end
end

%%%
%%% ComputeEllipseParams
%%%
function [stats, computedStats] = ComputeEllipseParams(L, stats, computedStats)
%   Find the ellipse that has the same 2nd-order moments as the region.
%   Compute the axes lengths, orientation, and eccentricity of the ellipse.
%   Ref: Haralick and Shapiro, Computer and Robot Vision
%   vol I, Addison-Wesley 1992, Appendix A.


if ~(computedStats.MajorAxisLength & computedStats.MinorAxisLength & ...
            computedStats.Orientation & computedStats.Eccentricity)
    computedStats.MajorAxisLength = 1;
    computedStats.MinorAxisLength = 1;
    computedStats.Eccentricity = 1;
    computedStats.Orientation = 1;
    
    if ndims(L) > 2
        NoNDSupport({'MajorAxisLength', 'MinorAxisLength', 'Eccentricity', 'Orientation'});
        return
    end
        
    [stats, computedStats] = ComputePixelList(L, stats, computedStats);
    [stats, computedStats] = ComputeCentroid(L, stats, computedStats);

    % Disable divide-by-zero warning
    warning_state = warning;
    warning('off');
    
    for k = 1:length(stats)
        list = stats(k).PixelList;
        if (isempty(list))
            stats(k).MajorAxisLength = 0;
            stats(k).MinorAxisLength = 0;
            stats(k).Eccentricity = 0;
            stats(k).Orientation = 0;
            
        else
            % Assign X and Y variables so that we're measuring orientation
            % counterclockwise from the horizontal axis.
            
            xbar = stats(k).Centroid(1);
            ybar = stats(k).Centroid(2);
            
            x = list(:,1) - xbar;
            y = -(list(:,2) - ybar);
            
            N = length(x);
            uxx = sum(x.^2)/N + 1/12;
            uyy = sum(y.^2)/N + 1/12;
            uxy = sum(x.*y)/N;
            
            common = sqrt((uxx - uyy)^2 + 4*uxy^2);
            stats(k).MajorAxisLength = 2*sqrt(2)*sqrt(uxx + uyy + common);
            stats(k).MinorAxisLength = 2*sqrt(2)*sqrt(uxx + uyy - common);
            stats(k).Eccentricity = 2*sqrt((stats(k).MajorAxisLength/2)^2 - ...
                    (stats(k).MinorAxisLength/2)^2) / ...
                    stats(k).MajorAxisLength;
            
            if (uyy > uxx)
                num = uyy - uxx + sqrt((uyy - uxx)^2 + 4*uxy^2);
                den = 2*uxy;
            else
                num = 2*uxy;
                den = uxx - uyy + sqrt((uxx - uyy)^2 + 4*uxy^2);
            end
            if (num == 0) && (den == 0)
                stats(k).Orientation = 0;
            else
                stats(k).Orientation = (180/pi) * atan(num/den);
            end
        end
    end
    
    % Restore warning state.
    warning(warning_state);
end
    
%%%
%%% ComputeSolidity
%%%
function [stats, computedStats] = ComputeSolidity(L, stats, computedStats)
%   Area / ConvexArea

if ~computedStats.Solidity
    computedStats.Solidity = 1;
    
    if ndims(L) > 2
        NoNDSupport('Solidity');
        return
    end
        
    [stats, computedStats] = ComputeArea(L, stats, computedStats);
    [stats, computedStats] = ComputeConvexArea(L, stats, computedStats);
    
    for k = 1:length(stats)
        if (stats(k).ConvexArea == 0)
            stats(k).Solidity = NaN;
        else
            stats(k).Solidity = stats(k).Area / stats(k).ConvexArea;
        end
    end
end

%%%
%%% ComputeExtent
%%%
function [stats, computedStats] = ComputeExtent(L, stats, computedStats)
%   Area / (BoundingBox(3) * BoundingBox(4))

if ~computedStats.Extent
    computedStats.Extent = 1;
    
    if ndims(L) > 2
        NoNDSupport('Extent');
        return
    end
        
    [stats, computedStats] = ComputeArea(L, stats, computedStats);
    [stats, computedStats] = ComputeBoundingBox(L, stats, computedStats);
    
    for k = 1:length(stats)
        if (stats(k).Area == 0)
            stats(k).Extent = NaN;
        else
            stats(k).Extent = stats(k).Area / prod(stats(k).BoundingBox(3:4));
        end
    end
end

%%%
%%% ComputeImage
%%%
function [stats, computedStats] = ComputeImage(L, stats, computedStats)
%   Binary image containing "on" pixels corresponding to pixels
%   belonging to the region.  The size of the image corresponds
%   to the size of the bounding box for each region.

if ~computedStats.Image
    computedStats.Image = 1;

    [stats, computedStats] = ComputeSubarrayIdx(L, stats, computedStats);

    for k = 1:length(stats)
        subarray = L(stats(k).SubarrayIdx{:});
        if ~isempty(subarray)
            stats(k).Image = (subarray == k);
        else
            stats(k).Image = logical(subarray);
        end
    end
end

%%%
%%% ComputePixelIdxList
%%%
function [stats, computedStats] = ComputePixelIdxList(L, stats, computedStats)
%   A P-by-1 matrix, where P is the number of pixels belonging to
%   the region.  Each element contains the linear index of the corresponding pixel.

if ~computedStats.PixelIdxList
    computedStats.PixelIdxList = 1;
    
    % Form a sparse matrix containing one column per region.  In
    % column P, the location of nonzero values correspond to the
    % linear indices of pixels in L that have value P.  For
    % example, S(100,5) is nonzero if and only L(100) equals 5.
    idx = find(L);
    elementValues = L(idx);
    S = sparse(idx, double(elementValues), 1);

    for k = 1:length(stats)
        stats(k).PixelIdxList = find(S(:,k));
    end
end

%%%
%%% ComputePixelList
%%%
function [stats, computedStats] = ComputePixelList(L, stats, computedStats)
%   A P-by-2 matrix, where P is the number of pixels belonging to
%   the region.  Each row contains the row and column coordinates of a pixel.

if ~computedStats.PixelList
    computedStats.PixelList = 1;
    
    [stats, computedStats] = ComputePixelIdxList(L, stats, computedStats);
    
    % Loop over each column of the sparse matrix.  Finding the
    % row indices of the nonzero entries in S(:,P) is equivalent
    % to finding the linear indices of pixels in L that equal P.
    % Convert the linear indices to subscripts and store
    % the results in the pixel list.  Reverse the order of the first
    % two subscripts to form x-y order.
    In = cell(1,ndims(L));
    for k = 1:length(stats)
        if ~isempty(stats(k).PixelIdxList)
            [In{:}] = ind2sub(size(L), stats(k).PixelIdxList);
            stats(k).PixelList = [In{:}];
            stats(k).PixelList = stats(k).PixelList(:,[2 1 3:end]);
        else
            stats(k).PixelList = zeros(0,ndims(L));
        end
    end
end

%%%
%%% ComputePerimeterCornerPixelList
%%%
function [stats, computedStats] = ComputePerimeterCornerPixelList(L, stats, computedStats)
%   Find the pixels on the perimeter of the region; make a list
%   of the coordinates of their corners; sort and remove duplicates.

if ~computedStats.PerimeterCornerPixelList
    computedStats.PerimeterCornerPixelList = 1;
    
    if ndims(L) > 2
        NoNDSupport('PerimeterCornerPixelList');
        return
    end
    
    [stats, computedStats] = ComputeImage(L, stats, computedStats);
    [stats, computedStats] = ComputeBoundingBox(L, stats, computedStats);

    for k = 1:length(stats)
        perimImage = bwmorph(stats(k).Image, 'perim8');
        firstRow = stats(k).BoundingBox(2) + 0.5;
        firstCol = stats(k).BoundingBox(1) + 0.5;
        [r,c] = find(perimImage);
        % Force rectangular empties.
        r = r(:) + firstRow - 1;
        c = c(:) + firstCol - 1;
        rr = [r-.5 ; r    ; r+.5 ; r   ];
        cc = [c    ; c+.5 ; c    ; c-.5];
        stats(k).PerimeterCornerPixelList = [cc rr];
    end
    
end

%%%
%%% ComputeConvexHull
%%%
function [stats, computedStats] = ComputeConvexHull(L, stats, computedStats)
%   A P-by-2 array representing the convex hull of the region.
%   The first column contains row coordinates; the second column
%   contains column coordinates.  The resulting polygon goes
%   through pixel corners, not pixel centers.

if ~computedStats.ConvexHull
    computedStats.ConvexHull = 1;
    
    if ndims(L) > 2
        NoNDSupport('ConvexHull');
        return
    end
    
    [stats, computedStats] = ComputePerimeterCornerPixelList(L, stats, computedStats);
    [stats, computedStats] = ComputeBoundingBox(L, stats, computedStats);

    for k = 1:length(stats)
        list = stats(k).PerimeterCornerPixelList;
        if (isempty(list))
            stats(k).ConvexHull = zeros(0,2);
        else
            rr = list(:,2);
            cc = list(:,1);
            hullIdx = convhull(rr, cc);
            stats(k).ConvexHull = list(hullIdx,:);
        end
    end
end

%%%
%%% ParseInputs
%%%
function [L,reqStats,msg] = ParseInputs_regionprops(officialStats, varargin)

L = [];
reqStats = [];
msg = '';

if (length(varargin) < 1)
    msg = 'Too few input arguments.';
    return
end

L = varargin{1};

if islogical(L)
    msg = 'Use bwlabel(BW) or double(BW) convert binary image to label matrix before calling regionprops.';
    error('Images:regionprops:binaryInput', '%s', msg);
end

checkinput(L, {'numeric'}, {'real', 'integer', 'nonnegative'}, mfilename, 'L', 1);

list = varargin(2:end);
if (~isempty(list) && ~iscell(list{1}) && strcmpi(list{1}, 'all'))
    reqStats = officialStats;
    reqStatsIdx = 1:length(officialStats);
    
elseif (isempty(list) || (~iscell(list{1}) && strcmpi(list{1},'basic')))
    % Default list
    reqStats = {'Area' 'Centroid' 'BoundingBox'};
else
    
    if (iscell(list{1})),	list = list{1};		end
    list = list(:);

    officialStatsL = lower(officialStats);

    reqStatsIdx = [];
    for k = 1:length(list)
        if (~isstr(list{k}))
            msg = 'Invalid input argument.';
            return
        end
        
        idx = strmatch(lower(list{k}), officialStatsL);
        if (isempty(idx))
            msg = sprintf('Unknown measurement: "%s".', list{k});
            return
        elseif (length(idx) > 1)
            msg = sprintf('Ambiguous measurement: "%s".', list{k});
            return
        else
            reqStatsIdx = [reqStatsIdx; idx];
        end
    end
    
    reqStats = officialStats(reqStatsIdx);
end

%%%
%%% NoNDSupport
%%%
function NoNDSupport(str)
%   Issue a warning message about lack of N-D support for a given  measurement or measurements.

if iscell(str)
    warn_str = sprintf('%s: %s ', 'These measurements are not supported if ndims(L) > 2', sprintf('%s ', str{:}));
else
    warn_str = sprintf('%s: %s', 'This measurement is not supported if ndims(L) > 2', str);
end
warning(warn_str)

%%% PreprocessRequestedStats
function requestedStats = PreprocessRequestedStats(requestedStats)
%   Remove any requested stats that are not supported for N-D input and issue an appropriate warning.

no_nd_measurements = {'MajorAxisLength'
                    'MinorAxisLength'
                    'Eccentricity'
                    'Orientation'
                    'ConvexHull'
                    'ConvexImage'
                    'ConvexArea'
                    'EulerNumber'
                    'Extrema'
                    'EquivDiameter'
                    'Solidity'
                    'Extent'
                    'PerimeterCornerPixelList'};

bad_stats = find(ismember(requestedStats, no_nd_measurements));
if ~isempty(bad_stats)
    NoNDSupport(requestedStats(bad_stats));
end

requestedStats(bad_stats) = [];

% ------------------------------------------------------------------------------------------------
function iptchecknargin(low, high, numInputs, function_name)
%IPTCHECKNARGIN Check number of input arguments.
%   IPTCHECKNARGIN(LOW,HIGH,NUM_INPUTS,FUNC_NAME) checks whether
%   the number of input arguments NUM_INPUTS is in the range specified
%   by LOW and HIGH. If NUM_INPUTS is not in this range, IPTCHECKNARGIN
%   issues a formatted error message.
%
%   LOW must be a scalar nonnegative integer.
%
%   HIGH must be a scalar nonnegative integer or Inf.
%
%   FUNC_NAME is a string that specifies the name used in the formatted
%   error message to identify the function checking its input
%   arguments.
%
%   Example
%   -------
%
%   Create a function and use IPTCHECKNARGIN to check that the 
%   number of arguments passed to the function is within the 
%   expected range.
%
%       function test_function(varargin)
%       iptchecknargin(1,3,nargin,mfilename);
%
%   Trigger the error message by executing the function at 
%   the MATLAB command line, specifying more the expected 
%   number of arguments.
%   
%       test_function(1,2,3,4)
%  
%   See also IPTCHECKHANDLE, IPTCHECKINPUT, IPTCHECKMAP, IPTCHECKSTRS,
%            IPTNUM2ORDINAL.

%   Copyright 1993-2004 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2004/11/19 17:36:50 $

% Input arguments are not checked for validity.

if numInputs < low
  msgId = sprintf('Images:%s:tooFewInputs', function_name);
  if low == 1
    msg1 = sprintf('Function %s expected at least 1 input argument', upper(function_name));
  else
    msg1 = sprintf('Function %s expected at least %d input arguments', upper(function_name), low);
  end
  
  if numInputs == 1
    msg2 = 'but was called instead with 1 input argument.';
  else
    msg2 = sprintf('but was called instead with %d input arguments.', numInputs);
  end
  
  error(msgId, '%s\n%s', msg1, msg2);
  
elseif numInputs > high
  msgId = sprintf('Images:%s:tooManyInputs', function_name);

  if high == 1
    msg1 = sprintf('Function %s expected at most 1 input argument', upper(function_name));
  else
    msg1 = sprintf('Function %s expected at most %d input arguments', upper(function_name), high);
  end
  
  if numInputs == 1
    msg2 = 'but was called instead with 1 input argument.';
  else
    msg2 = sprintf('but was called instead with %d input arguments.', numInputs);
  end
  
  error(msgId, '%s\n%s', msg1, msg2);
end
