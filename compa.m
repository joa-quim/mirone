function compa(fname, varargin)
% Para criar so o Mirone em C e nao compilar, fazer
%       compa('exe','-c')

% $Id: compa.m 7787 2016-02-09 16:07:48Z j $

	if (nargin == 0)
		fname = 'exe';	varargin = {'-c'};
	end

	aqui = pwd;
	if (~strcmp(fname,'exe'))
		cd c:\SVN\mironeWC\mex
		clear mex

		% Detect which matlab version is beeing used.
		versao = version;
		if (str2double(versao(end-5:end-2)) >= 2007),   MEX_EXT = '.mexw32';
		else                                            MEX_EXT = '.dll';
		end
		computas = computer;
		if (MEX_EXT(end) == '2' && computas(end) == '4')
			MEX_EXT = '.mexw64';
		end
		
		try
			make_mexs(fname, varargin{:})
			if (strmatch(fname, {'all' 'gmt' 'gdal' 'simple'}) ),	fname = '*';	end
			dos(['move ' fname MEX_EXT ' ..\lib_mex']);
			cd(aqui)
		catch
			cd(aqui);	disp('O compa deu merda')
			disp(lasterr)
		end
	else
		try
			cd([aqui '\utils'])
			dos('copy zoom_j.m zoom_j_bak.m');
			dos('C:\j\bin\sed.exe 1,$s/scribefiglisten_j/scribefiglisten/g zoom_j_bak.m > zoom_j.m');
			cd(aqui)
			
			pato_m = [matlabroot '\toolbox\matlab\graph2d\'];
			cd(pato_m)
			dos('ren scribefiglisten.m scribefiglisten_.m');
			dos('ren scribefiglisten_j.m scribefiglisten.m');	cd(aqui)
			
			addpath([aqui '\lib_fig'])
			
			if (isempty(varargin))
				mcc -B sgl -O all mirone
			else
				disp('Traduz so e nao compa')
				mcc -c -B sgl -O all mirone
				disp('Agora zipa')
				!zip -m -q mirone_src *.c *.h
			end
			%mcc -c -B sgl -O all lolo
			
			rmpath([aqui '\lib_fig'])
			dos(['copy ' aqui '\utils\zoom_j_bak.m ' aqui '\utils\zoom_j.m']);
			cd(pato_m)
			dos('ren scribefiglisten.m scribefiglisten_j.m');
			dos('ren scribefiglisten_.m scribefiglisten.m');	cd(aqui)
		catch
			cd(pato_m)
			dos('ren scribefiglisten.m scribefiglisten_j.m');
			dos('ren scribefiglisten_.m scribefiglisten.m');	cd(aqui)
			
			disp('O compa deu merda')
            disp(lasterr)
		end
	end
