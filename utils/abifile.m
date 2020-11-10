function aviobj = abifile(filename,varargin)
%AVIFILE Create a new AVI file
%   Copyright 1984-2002 The MathWorks, Inc.
       
	if nargin == 0
		error('You must provide a filename as the input.');
	elseif ~isa(filename,'char')
		error('The first input argument must be a filename.');
	end

	FramesPerSec = 10;      Quality = 7500;

	%Global defaults
	aviobj.Filename = '';
	if (isunix),    aviobj.Compression = 0;
	else            aviobj.Compression = 1;
	end
	aviobj.FileHandle = [];
	aviobj.TempDataFile = [tempname '.tmp'];
	aviobj.PrevFrame = [];
	aviobj.StreamName = filename;
	% 2 Key frames per second is the default
	aviobj.KeyFrameEveryNth = floor(FramesPerSec/2); 
	aviobj.CurrentState = 'Open';

	% Other defaults 
	aviobj.MainHeader.Fps= 1/FramesPerSec*10^6; % Microseconds per frame, 15 fps
	aviobj.MainHeader.MaxBytesPerSec=1000000;
	aviobj.MainHeader.Reserved=0;
	aviobj.MainHeader.Flags= 16; % AVIF_HASINDEX from vfw.h
	aviobj.MainHeader.TotalFrames=0;
	aviobj.MainHeader.InitialFrames= 0;
	aviobj.MainHeader.Streams= 1;
	aviobj.MainHeader.SuggestedBufferSize=0;
	aviobj.MainHeader.Width=0;
	aviobj.MainHeader.Height=0;
	aviobj.MainHeader.Scale=100; % Rate/Scale = Samples per second
	aviobj.MainHeader.Rate=FramesPerSec*100;
	aviobj.MainHeader.Start=0;
	aviobj.MainHeader.Length=0;

	% Default stream header, from VFW Programmer's Guide
	aviobj.StreamHeader.fccType = 'vids';
	if (isunix),	aviobj.StreamHeader.fccHandler =  'DIB ';
	else			aviobj.StreamHeader.fccHandler =  'DIB ';
	end
	%else            aviobj.StreamHeader.fccHandler =  'iv32';   end
	aviobj.StreamHeader.Flags = 0;
	aviobj.StreamHeader.Reserved = 0;
	aviobj.StreamHeader.InitialFrames = 0;
	aviobj.StreamHeader.Scale = 100;
	aviobj.StreamHeader.Rate = FramesPerSec*100;
	aviobj.StreamHeader.Start = 0;
	aviobj.StreamHeader.Length = 0;
	aviobj.StreamHeader.SuggestedBufferSize = 0;
	aviobj.StreamHeader.Quality = Quality;
	aviobj.StreamHeader.SampleSize =  0;

	% Default data header
	aviobj.Bitmapheader.biSize = 40;
	aviobj.Bitmapheader.biWidth =0;
	aviobj.Bitmapheader.biHeight =0;
	aviobj.Bitmapheader.biPlanes =1;
	aviobj.Bitmapheader.biBitCount =0;
	if (isunix),	aviobj.Bitmapheader.biCompression ='DIB ';
	else			aviobj.Bitmapheader.biCompression ='DIB ';
	end
	%else            aviobj.Bitmapheader.biCompression ='iv32';      end
	aviobj.Bitmapheader.biSizeImage = 0;
	aviobj.Bitmapheader.biXPelsPerMeter =0;
	aviobj.Bitmapheader.biYPelsPerMeter =0;
	aviobj.Bitmapheader.biClrUsed =0;
	aviobj.Bitmapheader.biClrImportant =0;
	aviobj.Bitmapheader.Colormap = [];

	[path,name,ext] = fileparts(filename);
	if isempty(ext),    filename = strcat(filename,'.avi');     end

	aviobj.Filename = filename;

	if ispc
		aviobj = class(aviobj,'abifile');
		% Take care of any parameters set at the command line
		%aviobj = set(aviobj,varargin{:});
		aviobj.FileHandle = avi('open',filename);
	end
