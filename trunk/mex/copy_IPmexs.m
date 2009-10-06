% This script will copy the MEXs from the IPT that are needed by Mirone

% Do private do IPT
fs = filesep;
pato_in  = [matlabroot fs 'toolbox' fs 'images' fs 'images' fs 'private'];

% Where to copy the MEXs. It must be MIRONEROOT/lib_mex
pato_out = '/home/bagside/programs/mirone_dev/trunk/lib_mex';

% The suffix for MEXs in this OS. You may need to addapt
if (ispc)
	suffix = '.dll';
	suffix = '.mexw32';
elseif (strncmp(computer,'MAC',3))
	suffix = '.maci';
else
	suffix = '.mexglx';
end

ipts = {
	'applylutc'
	'bwlabel1'
	'bwlabel2'
	'cq'
	'ditherc'
	'intlutc'
	'grayto8'
	'grayto16'
	'imhistc'
	'inv_lwm'
	'ordf'
	'parityscan'
	'bwboundariesmex'
	'bwlabelnmex'
	'grayxform'
	'imfilter_mex'
	'imlincombc'
	'imreconstructmex'
	'inv_piecewiselinear'
	'morphmex'
	'resampsep'};

for k=1:numel(ipts)
	copyfile([pato_in fs ipts{k} suffix], pato_out, 'f')
end

% In iptutils of IPT
%iptcheckinput
copyfile([matlabroot fs 'toolbox' fs 'images' fs 'iptutils' fs 'iptcheckinput' suffix], pato_out, 'f')
