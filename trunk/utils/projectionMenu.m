% --------------------------- PROJECTIONS MENU ------------------------------------
function projectionMenu(hFig)
    % Creates the Projection Menu from data read from the 'GMTproj_def.txt' file 

    fs = filesep;
    fid = fopen([cd fs 'data' fs 'GMTproj_def.txt']);
    if (fid < 0),   return;     end

	c = (fread(fid,'*char'))';      fclose(fid);
	menus = strread(c,'%s','delimiter','\n');   clear c fid;
    if (numel(menus{end}) < 5)     % Empty last line
        menus(end) = [];
    end
    
    % Parse the file contents
	m = length(menus);
    c = false(m,1);
    mainMenu = cell(m,1);
    subMenu = cell(m,1);
    projStr = cell(m,1);
    for (k=1:m)
        [t,r] = strtok(menus{k});
        if (t(1) == '#'),  c(k) = true;  continue;   end
        mainMenu{k} = t;
        [t,r] = strtok(r);
        if (numel(t) == 1),     subMenu{k} = [];    % We don't have a submenu
        else                    subMenu{k} = t;
        end
        if (strcmp(mainMenu{k},'None')),    projStr{k} = [];
        else                                projStr{k} = r;
        end
    end
    mainMenu(c) = [];    subMenu(c) = [];    projStr(c) = [];   % Remove comments lines
	m = length(mainMenu);       % Update counting

    % Detect repeated main menus which occur when we have subMenus
    primos = true(m,1);
    for (k=2:m),    primos(k) = ~strcmp(mainMenu{k}, mainMenu{k-1});    end
    
    % See if ...
    projGMT = cell(m,3);        % 3 for -J -C -T
    for (k=1:m)
        if (isempty(projStr{k})),   continue;   end      % The 'None' line
        [t,r] = strtok(projStr{k});
        projGMT{k,1} = [t '/1'];    % Append scale
        if (~isempty(r))            % Either a -C<...> or -T
            projGMT{k,2} = r;
            [t,r] = strtok(r);
            if (~isempty(t))
                projGMT{k,3} = t;   % The other (of -C<...> -T above)
            end
        else
            projGMT{k,2} = '-C';    % A must have
        end
    end
    
    % We are ready -- Create the "Projections" menu
    hProj = uimenu('Parent',hFig,'Label','Projections','Tag','Projections');

    hMain = zeros(1,m);    hSec  = zeros(1,m);
    for (k=1:m)
        if (isempty(subMenu{k}))    % No subMenu, we have than a direct projection setting
            hMain(k) = uimenu('Parent',hProj,'Label',mainMenu{k},'Call',{@setPRJ,hFig,k,projGMT});
        else
            if (primos(k))          % Parent of a Submenu. No proj string settings
                hMain(k) = uimenu('Parent',hProj,'Label',mainMenu{k});
                hUI = hMain(k);     % Make a copy to use on SubMenus
                hSec(k) = uimenu('Parent',hUI,'Label',subMenu{k},'Call',{@setPRJ,hFig,k,projGMT});
            else                    % Child of a Parent with a submenu. We have a proj string to set
                hSec(k) = uimenu('Parent',hUI,'Label',subMenu{k},'Call',{@setPRJ,hFig,k,projGMT});
            end
        end
    end
    
    setappdata(hFig,'ProjGMT','')       % Start by none
    
    ind = (hMain == 0);     hMain(ind) = [];
    ind = (hSec == 0);      hSec(ind) = [];
    
    projList =[hMain hSec];
    setappdata(hFig,'ProjList',projList)

% -------------------------------------------------------------------------------
function setPRJ(obj,nikles,hFig,k,projGMT)
    % Set the projection string in Figure's appdata and a checkmark on the selected projection
    projList = getappdata(hFig,'ProjList');
    handles = guidata(hFig);
    unchk = setxor(obj,projList);
    set(obj,'checked','on');    set(unchk,'checked','off')
    if (strcmp(get(obj,'Label'),'None'))
        setappdata(hFig,'ProjGMT','')
        handles.is_projected = 0;
    else
        prj = projGMT(k,:);
        % If don't have the third element, remove it because no empties in mapproject
        if (isempty(prj{end})),     prj(end) = [];  end
        setappdata(hFig,'ProjGMT',prj)
        handles.is_projected = 1;
    end
    guidata(hFig, handles)
    setAxesDefCoordIn(handles)    % Sets the value of the axes uicontextmenu that selects whether project or not
    