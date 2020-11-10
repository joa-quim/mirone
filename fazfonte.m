function fazfonte(nome, dir_opt)
% Converte to C code one or all of the figure files
%
% DIR_OPT use when compiling something not in the "src_figs" dir

	if (~strcmp(nome, 'todos'))
		if (nargin == 1)
			fig = [pwd '\src_figs\' nome '.m'];
		else
			fig = [pwd '\' dir_opt '\' nome '.m'];
		end
		mcc('-c', '-x', '-O', 'all', fig)
		disp('Agora zipa')
		dos(['zip -m -q ' nome ' *.c *.h']);
	else
		lst = dir([pwd '/src_figs/*.m']);
		for (k = 1:numel(lst))
			[pato,nome_] = fileparts(lst(k).name);
			mcc('-c', '-x', '-O', 'all', lst(k).name)
			disp(['Zipa o ' nome_])
			dos(['zip -m -q ' nome_ ' *.c *.h']);
		end
	end

