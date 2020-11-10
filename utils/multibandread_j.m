function im = multibandread_j(filename,dataSize,precision,offset,interleave,byteOrder,varargin)
%MULTIBANDREAD Read band interleaved data from a binary file
%   X = MULTIBANDREAD(FILENAME,SIZE,PRECISION,
%                      OFFSET,INTERLEAVE,BYTEORDER)
%   reads band-sequential (BSQ), band-interleaved-by-line (BIL), or
%   band-interleaved-by-pixel (BIP) data from a binary file, FILENAME.  X is
%   a 2-D array if only one band is read, otherwise it is 3-D. X is returned
%   as an array of data type double by default.  Use the PRECISION argument
%   to map the data to a different data type.
%
%   X = MULTIBANDREAD(FILENAME,SIZE,PRECISION,OFFSET,INTERLEAVE,
%                    BYTEORDER,SUBSET,SUBSET,SUBSET)
%   reads a subset of the data in the file. Up to 3 SUBSET parameters may be
%   used to subset independently along the Row, Column, and Band dimensions.
%
%   Parameters:
%
%     FILENAME: A string containing the name of the file to be read.
%
%     SIZE: A 3 element vector of integers consisting of
%     [HEIGHT, WIDTH, N]. HEIGHT is the total number of rows, WIDTH is
%     the total number of elements in each row, and N is the total number
%     of bands. This will be the dimensions of the data if it read in its
%     entirety.
%
%     PRECISION: A string to specify the format of the data to be read. For
%     example, 'uint8', 'double', 'integer*4'. By default X is returned as
%     an array of class double. Use the PRECISION parameter to format the
%     data to a different class.  For example, a precision of 'uint8=>uint8'
%     (or '*uint8') will return the data as a UINT8 array.  'uint8=>single'
%     will read each 8 bit pixel and store it in MATLAB in single
%     precision. See the help for FREAD for a more complete description of
%     PRECISION.
%
%     OFFSET: The zero-based location of the first data element in the file.
%     This value represents number of bytes from the beginning of the file
%     to where the data begins.
%
%     INTERLEAVE: The format in which the data is stored.  This can be
%     either 'bsq','bil', or 'bip' for Band-Seqential,
%     Band-Interleaved-by-Line or Band-Interleaved-by-Pixel respectively.
%
%     BYTEORDER: The byte ordering (machine format) in which the data is
%     stored. This can be 'ieee-le' for little-endian or 'ieee-be' for
%     big-endian.  All other machine formats described in the help for FOPEN
%     are also valid values for BYTEORDER.
%
%     SUBSET: (optional) A cell array containing either {DIM,INDEX} or
%     {DIM,METHOD,INDEX}. DIM is one of three strings: 'Column', 'Row', or
%     'Band' specifying which dimension to subset along.  METHOD is 'Direct'
%     or 'Range'. If METHOD is omitted, then the default is 'Direct'. If
%     using 'Direct' subsetting, INDEX is a vector specifying the indices to
%     read along the Band dimension.  If METHOD is 'Range', INDEX is a 2 or
%     3 element vector of [START, INCREMENT, STOP] specifying the range and
%     step size to read along the dimension. If INDEX is 2 elements, then
%     INCREMENT is assumed to be one.
%
%   Examples:
%
%   % Setup initial parameters for a dataset.
%   rows=3; cols=3; bands=5;
%   filename = tempname;
%
%   % Define the dataset.
%   fid = fopen(filename, 'w', 'ieee-le');
%   fwrite(fid, 1:rows*cols*bands, 'double');
%   fclose(fid);
%
%   % Read the every other band of the data using the Band-Seqential format.
%   im1 = multibandread_j(filename, [rows cols bands], ...
%             'double', 0, 'bsq', 'ieee-le', ...
%             {'Band', 'Range', [1 2 bands]} )
%
%   % Read the first two rows and columns of data using
%   % Band-Interleaved-by-Pixel format.
%   im2 = multibandread_j(filename, [rows cols bands], ...
%             'double', 0, 'bip', 'ieee-le', ...
%             {'Row', 'Range', [1 2]}, ...
%             {'Column', 'Range', [1 2]} )
%
%   % Read the data using Band-Interleaved-by-Line format.
%   im3 = multibandread_j(filename, [rows cols bands], ...
%             'double', 0, 'bil', 'ieee-le')
%
%   % Delete the file that we created.
%        delete(filename);
%
%   % The FITS file 'tst0012.fits' contains int16 BIL data starting at
%   % byte 74880.
%   im4 = multibandread_j( 'tst0012.fits', [31 73 5], ...
%             'int16', 74880, 'bil', 'ieee-be', ...
%             {'Band', 'Range', [1 3]} );
%   im5 = double(im4)/max(max(max(im4)));
%   imagesc(im5);
%
%   See also MULTIBANDWRITE, FWRITE, FREAD.

%   Copyright 2001-2005 The MathWorks, Inc.
%   $Revision: 1.1.6.6 $  $Date: 2005/11/15 01:10:46 $

error(nargchk(6,9,nargin));
error(nargchk(0,4,nargout));

% Get any subsetted dimensions
info = parseInputs(filename, dataSize, precision, offset, byteOrder, varargin{:});

% Make sure that the file is large enough for the requested operation
fileInfo = dir(info.filename);
fileSize = fileInfo.bytes;
if info.bitPrecision
    specifiedSize = info.offset + prod(info.dataSize) * (info.eltsize/8);
else
    specifiedSize = info.offset + prod(info.dataSize) * info.eltsize;
end
if fileSize < specifiedSize
    error('MATLAB:multibandread_j:badFileSize', ...
          ['The file is too small to contain the specified data.' ...
          'Check the size, offset, and precision arguments.']);
end

% Create a cell array of the dimension indices
ndx = {info.rowIndex info.colIndex info.bandIndex};

% Take care of the file ordering
switch lower(interleave)
    case 'bil'
        readOrder = [2 3 1];
        permOrder = [3 1 2];
    case 'bip'
        readOrder = [3 2 1];
        permOrder = readOrder;
    case 'bsq'
        readOrder = [2 1 3];
        permOrder = readOrder;
    otherwise
        error('MATLAB:multibandread_j:badInterleave', ...
            'INTERLEAVE must be ''bsq'', ''bil'', or ''bip''.');
end

% Decide which reading algorithm to use
% if isMappableType(info.inputClass)      % Read from a memory mapped file
%     im = readMemFile(filename, info, ndx, readOrder);
%     im = permute(im, permOrder);
% elseif ( strcmpi(interleave, 'bip') && (isempty(info.subset) || isequal(info.subset, 'b')) )
if ( strcmpi(interleave, 'bip') && (isempty(info.subset) || isequal(info.subset, 'b')) )
    % Special optimization for BIP cases
    im = readDiskFileBip(filename, info);
else                                    % Use the general-purpose routine.
    im = readDiskFile(filename, info, ndx, readOrder);
    im = permute(im, permOrder);
end

% %==========================================================================
% function im = readMemFile(filename, info, ndx, readOrder)
% % Memory map the file
% m = memmapfile(filename, 'offset', info.offset, 'repeat', 1, ...
%     'format', {info.inputClass info.dataSize(readOrder) 'x'});
% 
% % Permute the indices so that they are in read order
% ndx = ndx(readOrder);
% 
% % Do any necessary subsetting.
% im = m.data.x(ndx{1}, ndx{2}, ndx{3});
% 
% [str,maxsize,endian] = computer;
% 
% % Change the endianness, if necessary
% if strcmpi(endian, 'l')
%     if strcmpi(info.byteOrder, 'ieee-be') || strcmpi(info.byteOrder, 'b')
%         im = swapbytes(im);
%     end
% else
%     if strcmpi(info.byteOrder, 'ieee-le') || strcmpi(info.byteOrder, 'l')
%         im = swapbytes(im);
%     end
% end
% 
% % Change the type of the output, if necessary.
% if ~strcmp(info.inputClass, info.outputClass)
%     im = feval(info.outputClass, im);
% end

%==========================================================================
function im = readDiskFile(filename, info, srcNdx, readOrder)
% A general-purpose routine to read from the disk
% We use fread, which will handle non-integral (bit) types.
info.fid = fopen(filename, 'r');
lastReadPos = 0;
skip(info,info.offset,lastReadPos);

% Do permutation of sizes and indices
srcNdx = srcNdx(readOrder);
dim = info.dataSize(readOrder);

% Preallocate image output array
outputSize = [length(srcNdx{1}), length(srcNdx{2}), length(srcNdx{3})];
%im = zeros(outputSize(1), outputSize(2), outputSize(3), info.outputClass);
im = alloc_mex(outputSize(1), outputSize(2), outputSize(3), info.outputClass);

% Determine the start and ending read positions
kStart = srcNdx{1}(1);
kEnd = srcNdx{1}(end);

% srcNdx is a vector which contains the desired row, column, and band subsets
% of the input.  destNdx contains the destination in the output matrix.
destNdx(3) = 1;
for i=srcNdx{3}
    pos(1) = (i-1)*dim(1)*dim(2);
    destNdx(2) = 1;
    for j=srcNdx{2}
        pos(2) = (j-1)*dim(1);

        % Determine what to read
        posStart = pos(1) + pos(2) + kStart;
        posEnd = pos(1) + pos(2) + kEnd;
        readAmt = posEnd - posStart + 1;

        % Read the entire dimension
        skipNum = (posStart-1)-lastReadPos;
        if skipNum
            fread(info.fid, skipNum, info.precision);
        end
        [data, count] = fread(info.fid, readAmt, info.precision);
        lastReadPos = posEnd;
        if count ~= readAmt
            readError(info.fid);
        end

        % Assign the specified subset of what was read to the output matrix
        im(:,destNdx(2),destNdx(3)) = data(srcNdx{1}-kStart+1);
        destNdx(2) = destNdx(2) + 1;
    end
    destNdx(3) = destNdx(3) + 1;
end

fclose(info.fid);

%==========================================================================
function im = readDiskFileBip(filename, info)
% Read a file from disk, using optimizations applicable to the BIP case
% when we read bands at a time.
info.fid = fopen(filename, 'r');
lastReadPos = 0;
skip(info,info.offset,lastReadPos);

% extract the dataSize into meaningful terms
height = info.dataSize(1);
width  = info.dataSize(2);
bands  = info.dataSize(3);
%im = zeros(height, width, bands, info.outputClass);
im = alloc_mex(height, width, bands, info.outputClass);

% Read the file, one band at a time.
plane = 1;
for i=info.bandIndex
    skip(info,info.offset,i-1);
    [data,count] = fread(info.fid, height*width, ['1*' info.precision], (bands-1)*info.eltsize);
    if count ~= height*width
        readError(info.fid);
    end
    im(:,:,plane) = reshape(data,size(im(:,:,plane)));
    plane = plane + 1;
end

im = permute(im,[2 1 3]);
fclose(info.fid);

%==========================================================================
function readError(fid)
% Throw an error based on the fid status
if feof(fid)
    fclose(fid);
    error('MATLAB:multibandread_j:unexpectedEOF', ...
        'Error reading data. End of file reached unexpectedly.');
else
    fclose(fid);
    error('MATLAB:multibandread_j:readProblem', 'Error reading data. Unknown error.');
end

%==========================================================================
function skip(info,offset,skipSize)
% Skip to a specified position in the file
if info.bitPrecision
    fseek(info.fid,offset,'bof');
    fread(info.fid,skipSize,info.precision);
else
    fseek(info.fid,offset+skipSize*info.eltsize,'bof');
end

%==========================================================================
function info = parseInputs(filename, dataSize, precision, offset, byteOrder, varargin)
% Open the file. Determine pixel width and input/output classes.
fid = fopen(filename,'r',byteOrder);
if fid == -1
    error('MATLAB:multibandread_j:fileOpen', 'Unable to open %s for reading.',filename);
end
info = getPixelInfo(fid,precision);
info.filename = fopen(fid);
fclose(fid);

% Assign sizes
info.offset = offset;
if ~isnumeric(offset) || (offset < 0)
    error('MATLAB:multibandread_j:badOffset','OFFSET must be a number greater than zero.')
end
info.dataSize = dataSize;
if length(dataSize) ~= 3
    error('MATLAB:multibandread_j:badSize', ...
        'SIZE must be a 3 element vector of integers containing [HEIGHT, WIDTH, N]');
end
info.byteOrder  = byteOrder;

% Get the default indices.
% if isMappableType(info.inputClass)
%     info.rowIndex  = ':';
%     info.colIndex  = ':';
%     info.bandIndex = ':';
% else
    info.rowIndex  = 1:info.dataSize(1);
    info.colIndex  = 1:info.dataSize(2);
    info.bandIndex = 1:info.dataSize(3);
%end

% 'subset' is a string with 0-3 characters (r, c, b). The string
% represents which dimension is being subset
info.subset = '';

% Analyze the parameters that specify the subset of interest
if numel(varargin) > 0
    for i=1:length(varargin)
        % Determine the subsetting method
        methods = {'direct', 'range'};
        if (length(varargin{i}) == 2)       %Default to 'Direct'
            method = 'direct';
            n = 2;
        else
            match = strmatch(lower(varargin{i}{2}),methods);
            if ~isempty(match)
                method = methods{match};
            else
                error('MATLAB:multibandread_j:badSubset', ...
                    'Unrecogized subset string ''%s'' for METHOD.', ...
                    varargin{i}{2});
            end
            n = 3;
        end
        % Determine the orientation of the subset
        dimensions = {'row', 'column', 'band'};
        match = strmatch(lower(varargin{i}{1}),dimensions);
        if ~isempty(match)
            dim = dimensions{match};
        else
            error('MATLAB:multibandread_j:badSubset', ...
                'Unrecogized subset string ''%s'' for DIM.',varargin{i}{1});
        end
        % Get the indicies for the subset
        info.subset(i) = dimensions{match}(1); %build a string 'rcb'
        switch dim
            case 'row',     info.rowIndex = getIndices(method, i, n, varargin{:});
            case 'column',  info.colIndex = getIndices(method, i, n, varargin{:});
            case 'band',    info.bandIndex = getIndices(method, i, n, varargin{:});
        end
    end
end

%==========================================================================
function ndx = getIndices(method, i, n, varargin)
% Use a subsetting method
if strcmp(method,'direct')
    ndx = varargin{i}{n};
else
    switch length(varargin{i}{n})
        case 2
            ndx = feval('colon',varargin{i}{n}(1),varargin{i}{n}(2));
        case 3
            ndx = feval('colon',varargin{i}{n}(1), varargin{i}{n}(2), varargin{i}{n}(3));
        otherwise
            error('MATLAB:multibandread_j:badSubset', ...
                'Value for ''Range'' in the SUBSET parameter must be a 2 or 3 element vector.')
    end
end

%==========================================================================
function info = getPixelInfo(fid, precision)
% Returns size of each pixel.  Size is in bytes unless precision is
% ubitN or bitN, in which case width is in bits.

% Reformat the precision string
info.precision = precision(~isspace(precision));
if strncmp(info.precision(1),'*',1)
    info.precision(1) = [];
    info.precision = [info.precision '=>' info.precision];
end

% Determine the input and output types (classes)
lastInputChar = strfind(info.precision,'=>')-1;
if isempty(lastInputChar)
    lastInputChar=length(info.precision);
end
info.inputClass = precision(1:lastInputChar);
p = ftell(fid);
tmp = fread(fid, 1, info.precision);
info.eltsize = ftell(fid)-p;
info.outputClass = class(tmp);

% If it is a bit precision, parse the precision string to determine eltsize.
if isempty(strfind(info.precision,'bit'))
    info.bitPrecision = false;
else
    info.bitPrecision = true;
    info.eltsize = sscanf(info.inputClass(~isletter(info.inputClass)), '%d');
    if isempty(info.eltsize)
        error('MATLAB:multibandread_j:badPrecision', ...
            'Unrecognized PRECISION, ''%s''.',info.precision);
    end
end

% %==========================================================================
% function rslt = isMappableType(type)
% switch(type)
%     case {'int8' 'uint8' 'int16' 'uint16' 'int32' 'uint32' ...
%             'int64' 'uint64' 'single' 'double'}
%         rslt = true;
%     otherwise
%         rslt = false;
% end
