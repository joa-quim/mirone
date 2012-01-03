function varargout = tfw_funs(opt,varargin)
% Functions to read, write or inquire if exist+read .tfw type files
% File extension is used to search for .tfw .jgw .pgw or .gfw
% for registration files to, respectively, .tif|tiff .jpg|jpeg .png .gif 

%	Copyright (c) 2004-2012 by J. Luis
%
% 	This program is part of Mirone and is free software; you can redistribute
% 	it and/or modify it under the terms of the GNU Lesser General Public
% 	License as published by the Free Software Foundation; either
% 	version 2.1 of the License, or any later version.
% 
% 	This program is distributed in the hope that it will be useful,
% 	but WITHOUT ANY WARRANTY; without even the implied warranty of
% 	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% 	Lesser General Public License for more details.
%
%	Contact info: w3.ualg.pt/~jluis/mirone
% --------------------------------------------------------------------

    switch opt(1)
        case 'w'
            write_tfw(varargin{:})
        case 'r'        % Currently, this option is not called directly
            [varargout{1} varargout{2}] = read_tfw(varargin{:});
        case 'i'
            [varargout{1} varargout{2}] = inquire_tfw(varargin{:});
    end

% -------------------------------------------------------------------------
function write_tfw(handles,fname)
% HANDLES is the Mirone handles.
% FNAME (optional) is the name of the .tfw file

	if (nargin == 1)			% Propose a name based on the image's name
		name = get(handles.figure1,'Name');
		[pato, fname, EXT] = fileparts(name);
		EXT = strtok(EXT);		% Remove the "@ ??%" part
		if (strcmpi(EXT,'jpg')),		fname = [fname '.jgw'];
		elseif (strcmpi(EXT,'png'))		fname = [fname '.pgw'];
		elseif (strcmpi(EXT,'gif'))		fname = [fname '.gfw'];
		else							fname = [fname '.tfw'];
		end

		if (~isempty(pato)),	handles.work_dir = pato;	end
		str1 = fname;
		[FileName,PathName] = put_or_get_file(handles,str1,'Select TFW File name','put');
		if isequal(FileName,0),		return,		end
		fname = [PathName FileName];
	end
	
	% Write the .tfw file
	str{1} = num2str(handles.head(8),'%.10f');
	str{2} = '0.0';
	str{3} = '0.0';
	str{4} = num2str(-handles.head(9),'%.10f');
	str{5} = num2str(handles.head(1),'%.10f');
	str{6} = num2str(handles.head(4),'%.10f');
	fid = fopen(fname,'wt');
	for (i=1:6),    fprintf(fid,'%s\n',str{i});     end
	fclose(fid);

% ---------------------------------------------------------------------    
function [head,msg] = read_tfw(fname)
% read a FNAME .tfw type file (or its relatives extensions) and return a partially filled head(9)
    
	msg = [];   head = zeros(1,9);
	fid = fopen(fname,'r');
	if (fid < 0)
		msg = 'Error opening file';		return
	end

	str = fread(fid,'*char');
	fclose(fid);
	fw = strread(str','%f');
	if (numel(fw) ~= 6)
		msg = 'Wrong number of lines in file (must be 6).';		return
	end

	head(8) = fw(1);		head(9) = -fw(4);
	head(1) = fw(5);		head(4) = fw(6);

% ---------------------------------------------------------------------    
function [head,msg] = inquire_tfw(imgSize,pato,name,ext)
% Give at least two inputs. imgSize is a [1x2] two element vector with image's [size(I,1) size(I,2)]
% There are no error tests

    head = [];      msg = [];
    
    if (nargin == 2)		% Second arg is in fact the full file name
		[pato,name,ext] = fileparts(pato);
    end

	fw_ext = [];
	switch lower(ext(2:3))		% ext(1) = '.'
		case 'ti',		fw_ext = '.tfw';		% .tif or .tiff
		case 'jp',		fw_ext = '.jgw';		% .jpg or .jpeg
		case 'pn',		fw_ext = '.pgw';		% .png
		case 'gi',		fw_ext = '.gfw';		% .gif
	end

	fs = filesep;
	if (isempty(fw_ext) || ~exist([pato fs name fw_ext],'file'))			% world file?
		if (~exist([pato fs name '.wld'],'file')),		return,		end		% Definitely no. bye bye
		fw_ext = '.wld';
	end
    
	[head0,msg] = read_tfw([pato fs name fw_ext]);   % Remember that HEAD is not complete
	if (~isempty(msg)),     return;     end         % STOP here, an error occured while reading file

	% OK, if we reach here we have to compute the remaining HEAD(0) elements
	n = imgSize(2);							% columns
	m = imgSize(1);							% rows
	head0(2) = head0(1) + (n-1)*head0(8);	% x_max
	head0(3) = head0(4) - (m-1)*head0(9);	% y_min
	head0(5) = 255;							% Since this applies to images even if false shouldn't be dramatic

	head = head0;
