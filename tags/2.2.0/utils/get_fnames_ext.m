function [files,comp_files,comp_ext] = get_fnames_ext(pato, ext)
% Get the list of all files with extention "EXT" sitting in the "PATO" dir
%
% EXT may be either a char or a cell array. In the first case, only files with extension EXT
% will be returned (that is;  COMP_FILES & COMP_EXT are empty)
% On the second case, extra values of EXT will will be searched as well (that is; files with
% extension *.EXT{1}.EXT{2:length(EXT)}.
% FILES is a cell arrays of chars with the names that have extension EXT.
% COMP_FILES is a cell arrays of chars with the names that had extension EXT{2, or 3, or 4, etc...}.
% NOTE: the last extension is removed. E.G if file was lixo.dat.zip, it will become lixo.dat
% COMP_EXT is a cell arrays of chars with the extensions corresponding to COMP_FILES.
% An example is the search for files terminating in *.dat or *.dat.zip (EXT = {'dat' 'zip'})

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

	comp_files =[];     comp_ext = [];
	if (~(strcmp(pato(end),'\') || strcmp(pato(end),'/')))
		pato(end+1) = filesep;
	else
		pato(end) = filesep;
	end
	if (iscell(ext))
		ext1 = ext{1};
		for (k=2:length(ext))
			ext2{k-1} = ext{k};
		end
	else
		ext1 = ext;    ext2 = [];
	end

	tmp = dir([pato filesep '*.' ext1]);
	files = {tmp(:).name}';

	if (~isempty(ext2))         % That is, if we have one or more compression types (e.g. 'zip' 'gz')
		comp_files = [];    comp_ext = [];
		for (k=1:length(ext2))  % Loop over compression types
			tmp = dir([pato filesep '*.' ext1 '.' ext2{k}]);
			tmp = {tmp(:).name}';
			tmp1 = [];
			for m=1:length(tmp) % Loop over compressed files
				[PATH,FNAME,EXT] = fileparts(tmp{m});
				tmp{m} = [PATH,FNAME];
				tmp1{m} = EXT;  % Save File last extension as well
			end
			comp_files = [comp_files; tmp];
			comp_ext = [comp_ext; tmp1'];
		end
	end
