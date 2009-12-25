function out = deal_opts(opt, opt2)
% Function to provide fine (manual) control on some Mirone modules.
%
% The mechanism of this function relies on the contents of the "OPTcontrol.txt" control file.
% That file has keewords that trigger the usage of optional (specialized) features in Mirone.
%
% OPT directs to the option of interest (growing number)
% OPT2 is a uimenu handle for cases where the OPT case needs one

	MIRONE_DIRS = getappdata(0,'MIRONE_DIRS');
	opt_file = [MIRONE_DIRS.home_dir filesep 'data' filesep 'OPTcontrol.txt'];
	out = [];

	if ( ~exist(opt_file, 'file') == 2 ),	return,		end
	fid = fopen(opt_file, 'r');
	c = (fread(fid,'*char'))';      fclose(fid);
	lines = strread(c,'%s','delimiter','\n');   clear c fid;
	m = numel(lines);
	c = true(1, m);

	for (k = 1:m)
		if (~strncmp(lines{k},'MIR_',4)),	c(k) = false;	end
	end
	lines(c) = [];		% Delete non-keyword (comment) lines

	switch lower(opt)
		case 'gmtedit'		% Used when one wants to select what to plot in GMTEDIT slots
			for (k = 1:numel(lines))
				if (strncmp(lines{k}(5:end),'GMTEDIT',7))
					t = strtok(lines{k}(13:end));
					if ( strcmp(t(1:2), '-V') )		% Here we only check for a -V... and do not check for errors
						out = t;
					end
					break
				end
			end
		case 'mgg'			% When we want to find tracks inside calling rectangle
			for (k = 1:numel(lines))
				if (strncmp(lines{k}(5:end),'MGG',3))
					t = strtok(lines{k}(9:end));
					if (t == '1'),	out = true;		end
					break
				end
			end
			if (out),	custom_rect(opt2),	end		% Create a uimenu associated to a rectangle
	end
