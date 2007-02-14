function writekml(handles,fname,Z)
% Write a GoogleEarth kml file that will allow displaing the current image in GE

    % Check that we are dealing with geog coordinates
    if (~handles.geog)
        errordlg('Haven''t you noticed that GoogleEarth works only with geographical coordinates?','Error')
        return
    end
    
    [PATH,FNAME,EXT] = fileparts( fname);
    fname_kml = [PATH filesep FNAME '.kml'];% Make sure it has the .kml extension
    fname_img = [PATH filesep FNAME];       % Extension will be added later

    % Get the image at its own resolution
    img = get(handles.hImg,'CData');
    flipa = false;
    if (strcmp(get(handles.axes1,'YDir'),'normal'))
        img = flipdim(img,1);               % Ghrrrrrrr
        flipa = true;
    end

    if (~isempty(handles.have_nans))        % We need transparency here. Note that nans imply that image derives from grid
        fname_img = [fname_img '.png'];     % And we'll use png
        cmap = get(handles.figure1,'Colormap');
        if (ndims(img) == 2)                % Indexed image
            imwrite(img,cmap,fname_img,'Transparency',0)
        else                                % Truecolor
            %imwrite(img,fname_img,'Transparenc',cmap(1,:));    % BUGGED
            m = size(Z,1);      n = size(Z,2);
            ind = isnan(Z);
            alfa = alloc_mex(m,n,'uint8');
            alfa(~ind) = 255;       clear ind;
            if (flipa)      alfa = flipdim(alfa,1);     end
            imwrite(img,fname_img,'Alpha',alfa);
        end
    else                                    % Eventual original image transparency was already lost
        if (ndims(img) == 2)                % Indexed image
            fname_img = [fname_img '.png']; % And we still use png because img is indexed
            imwrite(img,get(handles.figure1,'Colormap'),fname_img);
        else
            fname_img = [fname_img '.jpg']; % Now use jpeg because it allows compression
            imwrite(img,fname_img,'Quality',100);
        end
    end
    
    %fname_kml = 'lixo.kml';
    fid = fopen(fname_kml,'w');
    
    fprintf(fid,'%s\n','<?xml version="1.0" encoding="UTF-8"?>');
    fprintf(fid,'%s\n','<kml xmlns="http://earth.google.com/kml/2.1">');
        fprintf(fid,'\t%s\n','<GroundOverlay>');
            fprintf(fid,'\t\t%s%s%s\n','<name>',fname_kml,'</name>');
            fprintf(fid,'\t\t%s\n','<Icon>');
               fprintf(fid,'\t\t\t%s%s%s\n','<href>',fname_img,'</href>');
            fprintf(fid,'\t\t%s\n','</Icon>');
            fprintf(fid,'\t\t%s\n','<altitudeMode>clampToGround</altitudeMode>');
            fprintf(fid,'\t\t%s\n','<LatLonBox>');
                fprintf(fid,'\t\t\t%s%s%s\n','<north>',sprintf('%.6f',handles.head(4)),'</north>');
                fprintf(fid,'\t\t\t%s%s%s\n','<south>',sprintf('%.6f',handles.head(3)),'</south>');
                fprintf(fid,'\t\t\t%s%s%s\n','<east>',sprintf('%.6f',handles.head(2)),'</east>');
                fprintf(fid,'\t\t\t%s%s%s\n','<west>',sprintf('%.6f',handles.head(1)),'</west>');
            fprintf(fid,'\t\t%s\n','</LatLonBox>');
        fprintf(fid,'\t%s\n','</GroundOverlay>');
    fprintf(fid,'%s','</kml>');
    
    fclose(fid);