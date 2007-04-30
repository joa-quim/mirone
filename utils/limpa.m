function limpa(handles)

mata = true;
fid = fopen([handles.path_data 'dourada'],'r');

if (fid > 0)
    mata = false;
    fclose(fid);
end

hTelhometro = findobj(handles.figure1,'Label','entry_tl');
hFloodFil = findobj(handles.figure1,'Label','entry_sh');
hVitrinite = findobj(handles.figure1,'Label','entry_vtr');
hTsunamovie = findobj(handles.figure1,'Label','entry_tsm');

if (mata)
    delete([hTelhometro hFloodFil hVitrinite hTsunamovie])
else
    set(hTelhometro,'Callback','telhometro(gcf)','Label','Telhometro')
    set(hFloodFil,'Callback','mirone(''Transfer_CB'',gcbo,[],guidata(gcbo),''Shape'')','Label','Shape detector')
    set(hVitrinite,'Callback','vitrinite','Label','Vitrinite')
    set(hTsunamovie,'Label','Make tsunami movies')
    uimenu('Parent',hTsunamovie,'Callback','umDmovie(guidata(gcbo))','Label','1d movie')
    uimenu('Parent',hTsunamovie,'Callback','tsunamovie(guidata(gcbo))','Label','2d movie')
end

setappdata(handles.figure1,'esDourada',~mata)   % Save goldness info
