function varargout = tfw_funs(opt,varargin)
    % Functions to read, write or inquire if exist+read .tfw type files
    % File extension is used to search for .tfw .jgw .pgw or .gfw
    % for registration files to, respectively, .tif|tiff .jpg|jpeg .png .gif 

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

	if (nargin == 1)    % Propose a name based on the image's name
        name = get(handles.figure1,'Name');
        [pato, fname, EXT] = fileparts(name);
        fname = [fname '.tfw'];
        if (~isempty(pato)),	cd(pato);   end
        [FileName,PathName] = uiputfile(fname,'Select TFW File name');
		pause(0.01)
		cd(handles.home_dir);       % allways go home
		if isequal(FileName,0);     return;     end
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
    % read a FNAME .tfw type file and return a partially filled head(9)
    
    msg = [];   head = zeros(1,9);
    fid = fopen(fname,'r');
    if (fid < 0)
        msg = 'Error opening file';
        return
    end
    
    str = fread(fid,'*char');
    fclose(fid);
    str = strread(str','%s','delimiter','\n');
    if (numel(str) ~= 6)
        msg = 'Wrong number of lines in file (must be 6).';
        return
    end

    head(8) = str{1};
    head(9) = -str{4};
    head(1) = str{5};
    head(4) = str{6};

% ---------------------------------------------------------------------    
function [head,msg] = inquire_tfw(head,pato,name,ext)
    % Give at least two inputs. HEAD is the non referenced image header
    % There are no error tests

    head = [];      msg = [];
    
    if (nargin == 2)    % Second arg is in fact the full file name
        [pato,name,ext] = fileparts(pato);
    end

    fw_ext = [];
    switch lower(ext(1:2))
        case 'ti',      fw_ext = '.tfw';       % .tif or .tiff
        case 'jp',      fw_ext = '.jgw';       % .jpg or .jpeg
        case 'pn',      fw_ext = '.pgw';       % .png
        case 'gi',      fw_ext = '.gfw';       % .gif
    end

    fs = filesep;
    if (isempty(fw_ext) || exist([pato fs name fw_ext]) == 0)   % no world file. bye bye
        return
    end
    
    [head0,msg] = read_tfw([pato fs name fw_ext]);   % Remember that HEAD is not complete
    if (~isempty(msg)),     return;     end         % STOP here, an error occured while reading file
    
    % OK, if we reach here we have to compute the remaining HEAD(0) elements
    n = round(diff(head(1:2)) / head(8) + 1);       % columns
    m = round(diff(head(3:4)) / head(9) + 1);       % rows
    head0(2) = head0(1) + (n-1)*head(8);            % x_max
    head0(3) = head0(4) - (m-1)*head(9);            % y_min
    head0(5) = 255;                                 % Since this applies to images even if false shouldn't be dramatic
    
    head = head0;