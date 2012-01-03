function img = screen_capture( h, name, format )
% This function starts a chain of calls that will finally arrive to the
% "hardcopy" built in function that does a screen capture.
[fp,fn,fe]=fileparts(name);
format = fe(2:end); % Do not want the '.'
[ops,dev,ext] = printtables;
i = strmatch( format, ext );
img = print_j( h, name, ['-d' dev{i}] );

% --------------------------------------------------------------------
function img = print_j( varargin )
%PRINT Print figure or model. 
pj = printjob;
%Check the input arguments and flesh out settings of PrintJob
pj = inputcheck_j( pj, varargin{:} );  
%Objects to print have their handles in a cell-array of vectors.
pj = render_j(pj, pj.Handles{1});
img = pj.Return;

% --------------------------------------------------------------------
function [pj, devices, options] = inputcheck_j( pj, varargin )
% Method to validate input arguments to PRINT.

	%Get cell arrays of driver specific information
	[options, devices, extensions, classes, colorDevs, destinations ] = printtables;

	for i = 1 : length(varargin)
		cur_arg = varargin{i};
		if (~ischar(cur_arg))
			pj.Handles = [ pj.Handles LocalCheckHandles(cur_arg) ];
		elseif (cur_arg(1) ~= '-')
			pj.FileName = cur_arg;
		else
			[ pj, devIndex ] = LocalCheckDevice( pj, cur_arg, devices );            
			pj.DriverExt = extensions{ devIndex };
			pj.DriverClass = classes{ devIndex };
			pj.DriverColor = strcmp( 'C', colorDevs{ devIndex } );
			pj.DriverColorSet = 1;
			pj.DriverExport = strcmp( 'X', destinations{ devIndex } );
		end
	end

% --------------------------------------------------------------------
function h = LocalCheckHandles( cur_arg )
% Checks that input matrix is full of handles to Figures and/or models.
	if ~iscell( cur_arg )
		cur_arg = {cur_arg};
	end
	h = cur_arg;

% --------------------------------------------------------------------
function [ pj, devIndex ] = LocalCheckDevice( pj, cur_arg, devices )
%LocalCheckDevice Verify device given is supported, and only one is given.
%    device proper starts after '-d', if only '-d'
%    we will later echo out possible choices

if ( size(cur_arg, 2) > 2 )    
    %Is there one unique match?
    devIndex = strmatch( cur_arg(3:end), devices, 'exact' );
    if length(devIndex) == 1
        pj.Driver = cur_arg(3:end);
    else
        %Is there one partial match, i.e. -dtiffn[ocompression]
        devIndex = strmatch( cur_arg(3:end), devices );
        if length( devIndex ) == 1
            %Save the full name
            pj.Driver = devices{devIndex};
        elseif length( devIndex ) > 1
            error( ['Device option ''' cur_arg '''is not unique'] )
        else
            % A special case, -djpegnn, where nn == quality level
            if strncmp( cur_arg, '-djpeg', 6 )
                if isempty( str2double(cur_arg(7:end)) )
                    error( 'JPEG quality level in device name must be numeric.' );
                end
                %We want to keep quality level in device name.
                pj.Driver = cur_arg(3:end); 
                devIndex = strmatch('jpeg',devices);
            else
                error(['Illegal device option, ''' cur_arg ''', specified.']);
            end
        end
    end
else
    devIndex = 0;
end

function pj = render_j( pj, h )
% Method to draw a model or Figure on current page of a print job.
% Figure or model is drawn in the output format specified by the device option to PRINT.

try
    pj.Error = 0;  %So caller knows there was an error
    %Make argument cell array for calling HARDCOPY
    inputargs = LocalPrintJob2Arguments( pj, h );
                       
	pj.Return = hardcopy( inputargs{:} );
catch
	pj.Error = 1;
    errordlg('An error occured while doing the screen capture','Error');
end

%------------------------------------------------------------------------------------
function imwriteArgs = LocalCreateImwriteArgs( pj )
%LOCALCREATEIMWRITEARGS Create a cell-array of input arguments for IMWRITE

imwriteArgs = {};
%We will have extra arguments for when we call IMWRITE.
if strcmp(pj.DriverClass, 'IM' )
    if strncmp( pj.Driver, 'jpeg', 4 )
        %Already checked that it is in acceptable format.
        imwriteArgs{end+1} = 'Quality';
        imwriteArgs{end+1} = sscanf(pj.Driver,'jpeg%d');    
        if isempty( imwriteArgs{end} )
            %Default quality level
            imwriteArgs{end} = 75;
        end
    end
end

%------------------------------------------------------------------------------------
function inputargs = LocalPrintJob2Arguments( pj, h )
%LOCALPRINTJOB2ARGUMENTS Make Cell-array of input arguments for old hardcopy from PrintJob.

inputargs{1} = h;   inputargs{2} = pj.FileName;
%If asking internal driver to create a zbuffer image to call
%IMWRITE with, set driver argument accordingly.
if strcmp( pj.DriverClass, 'IM')
  inputargs{3} = LocalGetImageDriver(h, pj.Renderer);
else
  inputargs{3} = [ '-d' pj.Driver];
end

if pj.DPI ~= -1,    inputargs{end+1} = ['-r' num2str(pj.DPI)];  end

%------------------------------------------------------------------------------------
function imageDriver = LocalGetImageDriver(h, renderer)
if isempty(renderer)
  renderer = get(h, 'Renderer');
end
if strcmp( renderer, 'painters' )
  renderer = 'zbuffer';
end
imageDriver = ['-d' renderer];
