function writegif(varargin)
% Extract of tests only bits of original imwrite to call the true writegif function

[data, map, filename, format, paramPairs, msg] = parse_inputs(varargin{:});

if (~isempty(msg)),    error(msg);      end
if (isempty(data)),    error('Image data can not be empty.');   end

if (isempty(format))
    format = get_format_from_filename(filename);
    if (isempty(format))
        error('Unable to determine the file format from the filename.');
    end
end

% Test if we can write to this directory.
fid = fopen(filename, 'a');
if (fid < 0)
    error('Could not open file for writing.  Check directory or file permissions.')
end
fclose(fid);

% Remove test file if not appending.
if (~any(strcmp('append',paramPairs))),    delete(filename);    end

% Call the writing function if it exists.
LocalWritegif(data, map, filename, paramPairs{:});

% ---------------------------------------------------------------------------
%%% Function parse_inputs
%%%
function [data, map, filename, format, paramPairs, msg] = parse_inputs(varargin)

data = [];          map = [];
filename = '';      format = '';
paramPairs = {};    msg = '';

if (nargin < 2)
    msg = 'Too few input arguments';
    return
end

firstString = [];
for k = 1:length(varargin)
    if (ischar(varargin{k}))
        firstString = k;
        break;
    end
end

if (isempty(firstString))
    msg = 'Invalid input arguments: missing filename';
    return
end

switch firstString
case 1
    msg = 'Invalid input arguments: first argument should not be a string';
    return
    
case 2      % imwrite(data, filename, ...)
    data = varargin{1};
    filename = varargin{2};
    
case 3      % imwrite(data, map, filename, ...)
    data = varargin{1};
    map = varargin{2};
    filename = varargin{3};
    if (size(map,2) ~= 3)
        msg = 'Invalid colormap';
        return
    end
    
    if ((min(map(:)) < 0) || (max(map(:)) > 1))
        warning('Colormap data should be in the range 0 to 1.')
    end

otherwise
    msg = 'Invalid input arguments: missing filename';
    return
end

if (length(varargin) > firstString)
    % There are additional arguments after the filename.
    if (~ischar(varargin{firstString + 1}))
        msg = 'Invalid input arguments.';
        return
    end
    
    % Is the argument after the filename a format specifier?
    fmt_s = imformats(varargin{firstString + 1});
    
    if (~isempty(fmt_s))        % imwrite(..., filename, fmt, ...)
        format = varargin{firstString + 1};
        paramPairs = varargin((firstString + 2):end);
        
    else                        % imwrite(..., filename, prop1, val1, prop2, val2, ...)
        paramPairs = varargin((firstString + 1):end);
    end
    
end

% Do some validity checking on param-value pairs
if (rem(length(paramPairs), 2) ~= 0)
    msg = 'Invalid input syntax';
    return
end

for k = 1:2:length(paramPairs)
    if (~ischar(paramPairs{k}))
        msg = 'Parameter names must be strings';
        return
    end
end

% ------------------------------------------------------------------------------
function format = get_format_from_filename(filename)

format = '';
idx = find(filename == '.');

if (~isempty(idx))
    ext = filename((idx(end) + 1):end);
    fmt_s = imformats(ext);
    if (~isempty(fmt_s)),   format = ext;    end
end

% ------------------------------------------------------------------------------
function LocalWritegif(mat,map,filename,varargin)
%WRITEGIF Write a GIF (Graphics Interchange Format) file to disk.
%	WRITEGIF(X,MAP,'filename') writes the indexed image X,MAP
%   to the file specified by the string FILENAME.
%
%   X can be a single image (M-by-N) or a series of
%   frames(M-by-N-by-1-by-P).  X must be of type uint8, logical, or double.
%   If X is uint8 or logical, then colormap indexing starts at zero.  If X
%   is double, then colormap indexing starts at one.
%
%   MAP can be a single colormap (M-by-3) that is applied to all frames, or a series of colormaps
%   (M-by-3-by-P) where P is the number of frames specified in X.  
%
%   WRITEGIF(X,[],FILENAME) writes the grayscale image GRAY
%   to the file.
%
%   WRITEGIF(...,'writemode','append') appends a single image to a file existing
%   on disk.
%
%   WRITEGIF(...,'comment',TEXT) writes an image to a file with the
%   comment specified in TEXT.  The comment is placed immediately before the image.
%   TEXT can be either a character arroy or a cell array of strings.  If
%   TEXT is a cell array of strings, then a carriage return is added after
%   each row.
%
%   WRITEGIF(...,'disposalmethod',DMETH) specifies the disposal
%   method for the image.  DMETH must be one of
%   'leaveInPlace','restoreBG','restorePrevious', or 'doNotSpecify'.
%
%   WRITEGIF(...,'delaytime',TIME) specifies the delay before displaying
%   the next image.  TIME must be a scalar value measured in seconds
%   between 0 and 655 inclusive.
%   
%   WRITEGIF(...,'transparentcolor',COLOR) specifies the transparent color
%   for the image.  COLOR must be a scalar index into the colormap.  If X
%   is uint8 or logical, then indexing starts at 0.  If X is double, then
%   indexing starts at 1.
%
%   WRITEGIF(...,'backgroundcolor',COLOR) specifies the background color
%   for the image.  COLOR must be a scalar index into the colormap.  If X
%   is uint8 or logical, then indexing starts at 0.  If X is double, then
%   indexing starts at 1.
%
%   WRITEGIF(...,'loopcount',COUNT) specifies the number of times to repeat
%   the animation.  If COUNT is Inf, the animation will be continuously
%   looped.  If COUNT is 0, the animation will be played once.  If
%   COUNT is 1, the animation will be played twice, etc.  The maximum value
%   of COUNT (other than Inf) is 65535.

%   WRITEGIF(...,'screensize',SIZE) specifies the screensize for the
%   image.  SIZE must be a two element vector where the first element is
%   the screen height and the second element is the screen width.
%
%   WRITEGIF(...,'location',LOC) specifies the offset of the top left
%   corner of the image relative to the top left corner of the screen.  LOC
%   must be a two element vector where the first element is the offset from
%   the top and the second element is the offset from the left.
%

%	The extension '.gif' will be added to 'filename' if it doesn't
%	already have an extension.
%

%	Rick Paxson 10-22-92
%	Copyright 1993-2004 The MathWorks, Inc.
%	$Revision: 1.1.6.5 $  $Date: 2005/11/15 01:13:26 $

%minimum of two arguments
error(nargchk(2, Inf, nargin));

%check filename string
if ~ischar(filename) || isempty(filename)
  error('Filename must be a non-empty string.');
end

%check for number of dimensions greater than four
nd = ndims(mat);

if nd > 4
  error(sprintf('%d-D data not supported for GIF files. Data must be 2-D or 4-D.', nd));
end

%RGB output not supported
if (nd == 3)
    if (size(mat,3) == 3)
        error('RGB output not supported for GIF files');        
    else
        error(sprintf('%d-D data not supported for GIF files. Data must be 2-D or 4-D.', nd));
    end
end

%check image data type
if (~isa(mat, 'logical') && ~isa(mat,'uint8') && ~isa(mat,'double'))
    error('X must be double, uint8, or logical');
end

ndcmap = ndims(map);
if ndcmap > 3
  error(sprintf('%d-D colormap not supported for GIF files. Colormap must be 2-D or 3-D.', ndcmap));
elseif (ndcmap == 3) && ~(nd == 4)
  error('Colormap may be 3-D only if 4-D image data is provided.');
end

if size(map,1)>256
    warning('GIF supports colormaps with a maximum of 256 colors. Colormap will be truncated.');
    map = map(1:256,:,:);
end
    
%this variable determines whether or not color index params are zero based
%or one based (this is determined by the data type of mat)
zerobased = 1;

%convert to zero based indexing if data is one based (type double)
if(isa(mat,'double'))
    mat = mat-1;
    zerobased = 0;
end

%%%%%%%%%%%%%
%Commented code dealing with userinput and aspectratio
%are commented because the user does not need control
%over these GIF features.
%%%%%%%%%%%%%

%parameter defaults
writemode = 0;  %default writemode is replace (0)
userinput = 0;
disposalmethod = 0;
delaytime = 50;
transparentcolor = -1;
comment = '';  
%bitdepth = ceil(log2(size(map,1)));
backgroundcolor = 0;
aspectratio = 1;
loopcount = -1;  %default value from user's perspective is actually 0
location = [0 0];
screensize = [size(mat,2) size(mat,1)];


% Process param/value pairs
propStrings = {'writemode','comment','disposalmethod',...
        'delaytime','transparentcolor',...
        'loopcount','location','screensize','backgroundcolor'};
%        'userinput','aspectratio'...

%Keep track of which file scope parameters
%are requested so that a warning can be generated if writemode is append.
commentrequested = false;
loopcountrequested = false;
screensizerequested = false;
backgroundcolorrequested = false;

for k = 1:2:length(varargin)
    prop = lower(varargin{k});
    value = varargin{k+1};
    if (~isstr(prop)),  error('Parameter name must be a string');    end
    idx = strmatch(prop, propStrings);
    if (isempty(idx))
        error(sprintf('Unrecognized parameter name "%s"', prop));
    elseif (length(idx) > 1)
        error(sprintf('Ambiguous parameter name "%s"', prop));
    end
    
    prop = deblank(propStrings{idx});
    
    switch prop
        case 'writemode'
            value = lower(value);
            switch value;
                case 'append',      writemode = 1;
                case 'overwrite',   writemode = 0;
                otherwise
                    error('''WriteMode'' must be either ''append'' or ''overwrite''');
            end
        case 'comment'
            comment = value;
            if (~ischar(comment) && ~iscellstr(comment))
                error(['''Comment'' value must be a column vector cell array of' ...
                     ' strings or a char matrix']);
            elseif(iscellstr(comment))
                tmp = '';
                for str = comment      %TODO - make sure newline characters are correct
                    tmp = [tmp str{1} char(13) char(10)];
                end
                comment = tmp;
            end
            commentrequested = true;
        case 'disposalmethod'
            value = lower(value);
            switch value;
                case 'leaveinplace',        disposalmethod=1;
                case 'restorebg',           disposalmethod=2;    
                case 'restoreprevious',     disposalmethod=3;   
                case 'donotspecify',        disposalmethod=0;    
                otherwise
                    error('''DisposalMethod'' must be ''leaveInPlace'',''restoreBG'',''restorePrevious'', or ''doNotSpecify''');
            end
        case 'delaytime'
            if(isnumeric(value) && (value<=655 && value>=0))
                %convert seconds to the nearest hundredth of a second
                delaytime = round(value*100);
            else
                error('DelayTime must be a number between 0 and 655 inclusive.');
            end
        case 'transparentcolor'
            if(isnumeric(value))
                transparentcolor = round(double(value));
                transparentcolor = transparentcolor - double(not(zerobased));  %if one-based image data, then this is also treated as one-based
            else
                error('TransparentColor must be an integer.');
            end
        case 'loopcount'
            if(isnumeric(value) && (value==inf || (value<=65535 && value>=0)))
                loopcount = round(value);
                switch(loopcount)
                    case 0,     loopcount = -1;
                    case Inf,   loopcount = 0;
                end
            else
                error('LoopCount must be an integer between 0 and 65535 inclusive or equal to Inf');
            end
            loopcountrequested = true;
        case 'location'
            if(isnumeric(value) && ndims(value)==2 && length(value)==2 && min(size(value))==1)
                location = round(value);
            else
                error('Location must be a two element vector');
            end
        case 'screensize'
            if(isnumeric(value) && ndims(value)==2 && length(value)==2 && min(size(value))==1)
                screensize = round(value);
            else
                error('ScreenSize must be a two element vector');
            end
            screensizerequested=true;
        case 'backgroundcolor'
            if(isnumeric(value))
                backgroundcolor = double(value);
                if(not(zerobased))
                    backgroundcolor = backgroundcolor - 1;  %if one-based image data, then this is also treated as one-based
                end
            end
            backgroundcolorrequested = true;
%         case 'aspectratio'
%             tmp = (value*64)/15;
%             if(isnumeric(value) && value>=0 && value<=255)
%                 aspectratio = round(tmp);
%             else
%                 error('AspectRatio must be an integer between 0 and the 255 inclusive');
%             end
%             if(writemode)
%                 warning('The AspectRatio parameter has no effect when appending to a file');
%             end
%         case 'userinput'
%             switch value;
%                 case 0
%                     userinput=0;
%                 case 1
%                     userinput=1;    
%                 otherwise
%                     error('UserInput must be either 0 or 1');
%             end            
    end
end

%add .gif extension to filename if necessary
if isempty(strfind(filename,'.')),    filename = [filename '.gif'];     end

%warn if file scope parameters requested in append mode
if(writemode == 1)
    if(commentrequested)
        warning('The ''Comment'' parameter has no effect in append mode')
        comment = '';
    end
    if(loopcountrequested)
        warning('The ''LoopCount'' parameter has no effect in append mode')
        loopcount = -1;
    end
    if(screensizerequested)
        warning('The ''ScreenSize'' parameter has no effect in append mode')
        screensize = [size(mat,2) size(mat,1)];
    end
    if(backgroundcolorrequested)
        warning('The ''BackgroundColor'' parameter has no effect in append mode')
        backgroundcolor = 0;
    end
end

%generate grayscale colormap if map is empty and we are not appending 
%(will use global colormap if appending)
if(isempty(map) && not(writemode==1))
    if(islogical(mat)),     map = gray(2);
    else                    map = gray(256);
    end
end

% assume normal orientation of the colormap and image data passed in
if nd==2
    mat = mat';
elseif nd ==4
    mat = permute(mat,[2,1,3,4]);
end

if ndcmap==2
    map = map';
elseif ndcmap==3
    map = permute(map,[2,1,3]);
end

if(isempty(map)),   cmaplength = 256;
else                cmaplength = size(map,2);
end

greaterthancmap = find(mat>=cmaplength);
lessthancmap = mat<0;

if(~isempty(greaterthancmap) || ~isempty(find(mat<0)))
    warning('Image data contains values that are out of range.  Out of range values will be given the nearest valid value.');
    mat(greaterthancmap) = cmaplength-1;
    mat(lessthancmap) = 0;
end

if(transparentcolor>=cmaplength || transparentcolor < -1)
    warning('''TransparentColor'' property out of range and will have no effect.');
    transparentcolor = -1;
end
if(backgroundcolor>=cmaplength || backgroundcolor < -1)
    warning('''BackgroundColor'' property out of range and will have no effect.');
    backgroundcolor = -1;
end

%make sure data is uint8 format
if (~isa(mat,'uint8')),     mat = uint8(mat);   end

wgifc(mat, map, filename,writemode,userinput,disposalmethod,delaytime,...
    transparentcolor,comment,backgroundcolor,aspectratio,loopcount,location,screensize);
