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
