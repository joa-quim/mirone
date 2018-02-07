function reconstruct_plates(hLine)
% Reconstruct the base image at the time of the HLINE isochron

%	Copyright (c) 2004-2018 by J. Luis
%
% 	This program is part of Mirone and is free software; you can redistribute
% 	it and/or modify it under the terms of the GNU Lesser General Public
% 	License as published by the Free Software Foundation; either
% 	version 2.1 of the License, or any later version.
% 
% 	This program is distributed in the hope that it will be useful,
% 	but WITHOUT ANY WARRANTY; without even the implied warranty of
% 	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
% 	Lesser General Public License for more details.
%
%	Contact info: w3.ualg.pt/~jluis/mirone
% --------------------------------------------------------------------

% $Id$

	handles = guidata(hLine);

	hLines = get_polygon(handles.figure1, 'multi');				% Get the line handles
	hLines = setxor(hLines, hLine);
	pause(0.1)

	hAllIsocs = findobj(get(hLine,'Parent'),'Tag', get(hLine,'Tag'));		% Fish all isocs, but only the isocs

	% Find the conjugate anomaly pairs.
	indConj = zeros(numel(hLines),1);		poles = zeros(numel(hLines), 3);
	for (k = 1:numel(hLines))
		[ind, pole, msg] = get_conjugates(hAllIsocs, hLines(k));
		if (~isempty(msg)),		errordlg(msg, 'ERROR'),	return,		end
		indConj(k) = ind;
		poles(k, :) = [pole.lon pole.lat pole.ang];
	end

	[plates_fixed, plates_mobile] = build_plate_polyg(handles, hAllIsocs, hLines, indConj, poles);
	name = strtok(getappdata(hLines(1),'LineInfo'));		% The anomaly name

	for (k = 1:numel(hLines))
		h = mirone('DrawLine_CB', handles, 'data', plates_fixed{k}(:,1),  plates_fixed{k}(:,2));
		set(h, 'Tag', 'PlateFix'),		setappdata(h, 'PlateFix', k)
		ui = get(h, 'UIContextMenu');
		uimenu(ui, 'Label', 'RECONSTRUCT PLATES','Sep','on', 'Callback', {@do_reconst, poles, name});

		h = mirone('DrawLine_CB', handles, 'data', plates_mobile{k}(:,1), plates_mobile{k}(:,2));
		set(h, 'Tag', 'PlateMov'),		setappdata(h, 'PlateMov', k)
		ui = get(h, 'UIContextMenu');
		uimenu(ui, 'Label', 'RECONSTRUCT PLATES', 'Sep','on', 'Callback', {@do_reconst, poles, name});
	end

% -----------------------------------------------------------------------------------------------------------------
function do_reconst(obj, evt, poles, name)
% Callback function that does the reconstruction work
	handles = guidata(obj);
	set(handles.figure1,'pointer','watch')

	% We do this again to allow for the possibility that polygons have been edited
	[plates_fixed, plates_mobile, globalBB, opt_A] = get_plate_polyg(handles, poles);

	n_rows = round((globalBB(4) - globalBB(3)) / handles.head(9) + 1 - handles.head(7));
	n_cols = round((globalBB(2) - globalBB(1)) / handles.head(8) + 1 - handles.head(7));
	Zg = alloc_mex(n_rows, n_cols, 'single', NaN);
	sG.head = [globalBB 0 0 0 handles.head(8:9)];
	nPlates = numel(plates_fixed);

	use_aguenta = false;					% Less than 20 tracks don't use aguentabar
 	if (nPlates > 1),	aguentabar(0,'title','Reconstruction Work'),	use_aguenta = true;		end

	for (k = 1:nPlates)
		[X,Y,Z,head] = mirone('ImageCrop_CB', handles, plates_fixed{k}, 'CropaGrid_pure-');
		s.head = head;		s.X = X;	s.Y = Y;	s.Z = Z;			% A wrap for transplants()
		sG.Z = Zg;
		Zg = transplants([], 'grid_brute', true, sG, s);				% Transplant this fixed plate to final grid

		[X,Y,Z,head] = mirone('ImageCrop_CB', handles, plates_mobile{k}, 'CropaGrid_pure-');
		G = gmt('wrapgrid', Z, head);
		Gr = gmtmex(sprintf('grdrotater %s -N -nn -E%g/%g/%g -F', opt_A{k}, poles(k,1), poles(k,2), -poles(k,3)), G, plates_mobile{k});
		s.head = [Gr.range 0 Gr.inc];	s.X = Gr.x;		s.Y = Gr.y;		s.Z = Gr.z;
		sG.Z = Zg;
		Zg = transplants([], 'grid_brute', true, sG, s);				% Transplant this rotated plate to final grid
		if (use_aguenta),	aguentabar(k/nPlates);		end
	end
	set(handles.figure1,'pointer','arrow'),		clear Gr Z s

	G = gmt('wrapgrid', Zg, sG.head);
	G.title = ['Reconstruction at anomaly ', name];
	h = mirone(G);

	% Overlay the plates polygons
	handMirNew = guidata(h);
	for (k = 1:nPlates)
		mirone('DrawLine_CB', handMirNew, 'data', plates_fixed{k}(:,1), plates_fixed{k}(:,2))
	end	
	for (k = 1:nPlates)
		[rlon,rlat] = rot_euler(plates_mobile{k}(:,1), plates_mobile{k}(:,2), poles(k,1), poles(k,2), -poles(k,3), -1);
		mirone('DrawLine_CB', handMirNew, 'data', rlon, rlat)
	end	

% -----------------------------------------------------------------------------------------------------------------
function [ind, pole, msg] = get_conjugates(hAllIsocs, hLine)
% Find the conjugate plate of HLINE. Return its index in HALLISOCS and the finite pole stored in HLINE appdata.

	ind = 0;	msg = '';
	lineInfo = getappdata(hLine, 'LineInfo');	% Ex: 9 NORTH AMERICA/IBERIA FIN"138.2 62 6.132 27.47"
	isoc = strtok(lineInfo);
	[P1, P2] = pole2neighbor('get_plate_pair', [], [], [], [], hLine);

	for (n = 1:numel(hAllIsocs))		% Loop over all Isochrons
		this_lineInfo  = getappdata(hAllIsocs(n),'LineInfo');
		this_isoc = strtok(this_lineInfo);
		if (strcmpi(this_isoc, isoc))	% OK, same isochron so now search for the same plate pair
			if (~isempty(strfind(this_lineInfo, P1)) && ~isempty(strfind(this_lineInfo, P2)))	% Then we found a target
				% but we still don't know the plate pair order
				[PP1, PP2] = pole2neighbor('get_plate_pair', [], [], [], [], this_lineInfo);
				if (strcmp(P1, PP2) && strcmp(P2, PP1))		% Found the conjugate pair
					ind = n;
					break
				end
			end
		end
	end
	if (~ind)
		msg = sprintf('Could not find the conjugate anomaly -> %s  %s/%s  Bye', isoc, PP2, PP1);
		return
	end
	pole = pole2neighbor('parse_finite_pole', [], [], [], [], lineInfo);
	if (isempty(pole))
		msg = 'The selected anomaly has no FINite rotation info. Bye';
	end

% -----------------------------------------------------------------------------------------------------------------
function [plates_fixed, plates_mobile] = build_plate_polyg(handles, hAllIsocs, hLines, indConj, poles)
% Construct the polygons of each pseudo-plate from each isochron. Do this crudly from an isoc outward rotation
% Comute also the global limits (in BB) of the final grid holding the fix and rotated plates.
% OPT_A hold the limits (-R) of the rotated plates as sub-grids of the global grid.

	nPlates = numel(hLines);
	plates_fixed = cell(nPlates,1);		plates_mobile = cell(nPlates,1);
	for (k = 1:nPlates)
		x = get(hLines(k), 'XData');			y = get(hLines(k), 'YData');					% A fixed isoc
		[x, y] = pseudo_plate_polyg(x, y, poles(k, :), 'left');
		plates_fixed{k} = [x y];
		x = get(hAllIsocs(indConj(k)), 'XData');	y = get(hAllIsocs(indConj(k)), 'YData');	% A moving isoc
		[x, y] = pseudo_plate_polyg(x, y, poles(k, :), 'right');
		plates_mobile{k} = [x y];
	end

% -----------------------------------------------------------------------------------------------------------------
function [plates_fixed, plates_mobile, BB, opt_A] = get_plate_polyg(handles, poles)
% ...
	hFix = findobj(handles.axes1, 'Tag', 'PlateFix');
	hMov = findobj(handles.axes1, 'Tag', 'PlateMov');
	if (numel(hFix) ~= numel(hMov))
		errordlg('Happy? You screw it by deleting plates polygons. Bye Bye!', 'Error')
		plates_fixed = [];	plates_mobile = [];		BB = [];	opt_A = '';
		return
	end

	nPlates = numel(hFix);
	orderFix = zeros(nPlates,1);	orderMov = zeros(nPlates,1);
	for (k = 1:nPlates)
		orderFix(k) = getappdata(hFix(k), 'PlateFix');
		orderMov(k) = getappdata(hMov(k), 'PlateMov');
	end

	[lix, ind] = sort(orderFix);
	hFix = hFix(ind);
	[lix, ind] = sort(orderMov);
	hMov = hMov(ind);

	plates_fixed = cell(nPlates,1);		plates_mobile = cell(nPlates,1);
	for (k = 1:nPlates)
		x = get(hFix(k), 'XData');		y = get(hFix(k), 'YData');
		plates_fixed{k} = [x(:) y(:)];
		x = get(hMov(k), 'XData');		y = get(hMov(k), 'YData');
		plates_mobile{k} = [x(:) y(:)];
	end

	% Now compute the global BB
	BB = [Inf -Inf Inf -Inf];
	for (k = 1:nPlates)
		BB(1) = min(BB(1), min(plates_fixed{k}(:,1)));		BB(2) = max(BB(2), max(plates_fixed{k}(:,1)));
		BB(3) = min(BB(3), min(plates_fixed{k}(:,2)));		BB(4) = max(BB(4), max(plates_fixed{k}(:,2)));
	end

	% For the moving plates, we have to rotate their individual BBs first.
	movBB = zeros(nPlates, 4);		thisBB = zeros(nPlates, 4);
	for (k = 1:nPlates)
		thisBB(k,1:2:3) = min(plates_mobile{k});	thisBB(k,2:2:4) = max(plates_mobile{k});
		[rlon,rlat] = rot_euler([thisBB(k,1) thisBB(k,1) thisBB(k,2) thisBB(k,2)]', [thisBB(k,3) thisBB(k,4) thisBB(k,4) thisBB(k,3)]', ...
		                        poles(k,1), poles(k,2), -poles(k,3), -1);
		movBB(k,1) = min(rlon);				movBB(k,2) = max(rlon);
		movBB(k,3) = min(rlat);				movBB(k,4) = max(rlat);
		BB(1) = min(BB(1), movBB(k,1));		BB(2) = max(BB(2), movBB(k,2));
		BB(3) = min(BB(3), movBB(k,3));		BB(4) = max(BB(4), movBB(k,4));
	end

	headOuter = [BB 0 0 0 handles.head(8:9)];
	BB = aux_funs('adjust_inner_BB', handles.head, headOuter);		% Adjust the global BB to fit nicely in base grid
	headOuter = [BB 0 0 0 handles.head(8:9)];		% To use in the next loop

	% We also compute the grdrotater opt_A's because it's here that we have the needed info.
	opt_A = cell(nPlates,1);
	for (k = 1:nPlates)
		headInner = [movBB(k,:) 0 0 0 handles.head(8:9)];
		thisBB = aux_funs('adjust_inner_BB', headOuter, headInner);
		opt_A{k} = sprintf('-A%.12g/%.12g/%.12g/%.12g', thisBB);
	end

% -----------------------------------------------------------------------------------------------------------------
function [xP, yP] = pseudo_plate_polyg(x, y, pole, side)
% ... The base and top should be flow-lines?
	x = x(:);	y = y(:);
	if (side(1) == 'l')			% Left
		ang = -4;
	else
		ang = 4;
	end

	[rlon,rlat] = rot_euler(x, y, pole(1), pole(2), ang, -1);
	D = gmtmex('gmtsimplify -T5m -fg', [rlon rlat]);		% DP simplify it. Usefull if we want line-edit it.
	rlon = D.data(:,1);		rlat = D.data(:,2);
	xP = [x; rlon(end:-1:1); x(1)+100*eps];			% Adding the 'eps' to prevent that polygons are plotted as patches
	yP = [y; rlat(end:-1:1); y(1)];
