% This script will copy the MEXs from the IPT that are needed by Mirone
% Please revise for your OS (for instance, if it is a 64 bits suffixs need to be changed accordingly)

% Do private do IPT
fs = filesep;
pato_in  = [matlabroot fs 'toolbox' fs 'images' fs 'images' fs 'private'];

% Where to copy the MEXs. It must be MIRONEROOT/lib_mex
pato_out = '/home/bagside/programs/mirone_dev/trunk/lib_mex';

% The following 3 still need to be copied (either depend on IPP or not enough info to rebuild)
ipts = { 'imfilter_mex' 'imlincombc' 'inv_piecewiselinear'};

% On newer MATLAB instalations you may want to also copy the ones below as well (they might be faster)
%	'applylutc'
%	'bwlabel1'
%	'bwlabel2'
%	'cq'
%	'ditherc'
%	'intlutc'
%	'grayto8'
%	'grayto16'
%	'imhistc'
%	'inv_lwm'
%	'ordf'
%	'parityscan'
%	'bwboundariesmex'
%	'bwlabelnmex'
%	'grayxform'
%	'imreconstructmex'
%	'morphmex'
%	'resampsep'


% The suffix for MEXs in this OS. You may need to addapt
if (ispc)
	suffix = '.dll';
	suffix = '.mexw32';
elseif (strncmp(computer,'MAC',3))
	suffix = '.mexmaci';
else
	suffix = '.mexglx';
end


for k=1:numel(ipts)
	copyfile([pato_in fs ipts{k} suffix], pato_out, 'f')
end

% In iptutils of IPT
%iptcheckinput
copyfile([matlabroot fs 'toolbox' fs 'images' fs 'iptutils' fs 'iptcheckinput' suffix], pato_out, 'f')
