function reconstruct_plates(hLine)
% Reconstruct the base image at the time of the HLINE isochron

%	Copyright (c) 2004-2019 by J. Luis
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

% $Id: reconstruct_plates.m 11406 2019-01-15 20:33:04Z j $

	handles = guidata(hLine);
	name = strtok(getappdata(hLine,'LineInfo'));			% The anomaly name

	is_polygon = aux_funs('isclosed', hLine, 0.5);
	if (is_polygon)
		hLines = hLine;						% For polygons only one at each time
% 	x = get(hLine, 'XData');	y = get(hLine, 'YData');
% 	if (x(1) == x(end) && y(1) == y(end))
% 		hLines = hLine;						% For polygons only one at each time
% 		is_polygon = true;
	else
		hLines = get_polygon(handles.figure1, 'multi');		% Get other line handles
		hLines = setxor(hLines, hLine);
		if (isempty(hLines)),	hLines = hLine;		end		% One selection only, which was the hLine
% 		is_polygon = false;
	end

	hAllIsocs = findobj(get(hLine,'Parent'),'Tag', get(hLine,'Tag'));	% Fish all isocs, but only the isocs

	% Find the conjugate anomaly pairs.
	indConj = zeros(numel(hLines),1);	poles = zeros(numel(hLines), 3);
	for (k = 1:numel(hLines))				% Loop over fixed plates
		if (is_polygon)
			[ind, p] = get_conjugates_polyg(hAllIsocs, hLines(k));
			if (k == 1),	poles = cell(numel(hLines),1);		indConj = cell(numel(hLines),1);	end
			poles{k} = p;	indConj{k} = ind;
		else
			[ind, pole] = get_conjugates(hAllIsocs, hLines(k));
			indConj(k)  = ind;
			poles(k,:) = [pole.lon pole.lat pole.ang];
		end
	end

	if (is_polygon)
		hPF = zeros(10,1);		hPM = zeros(10,1);		% Allocate in excess
		n = 1;
		for (k = 1:numel(indConj))
			for (m = 1:numel(indConj{k}))
				hPF(n) = hLines(k);
				hPM(n) = hAllIsocs(indConj{k}(m));
				n = n + 1;
			end
		end
		poles = cat(1, poles{:});			% We need an array in which the rows correspond to the line handles in hPM
		hPF(hPF == 0) = [];		hPM(hPM == 0) = [];
		hF = unique(hPF);
		for (k = 1:numel(hF))				% Loop over unique fixed plates (normally should be only 1)
			setappdata(hF(k), 'PlateFix', k)
			ui = get(hF(k), 'UIContextMenu');
			uimenu(ui, 'Label', 'RECONSTRUCT PLATES','Sep','on', 'Callback', {@do_reconst, poles, name, hPF, hPM});
		end
		for (k = 1:numel(hPM))
			setappdata(hPM(k), 'PlateMov', k)
			ui = get(hPM(k), 'UIContextMenu');
			uimenu(ui, 'Label', 'RECONSTRUCT PLATES','Sep','on', 'Callback', {@do_reconst, poles, name, hPF, hPM});
		end
	else
		[plates_fixed, plates_mobile] = build_plate_polyg(handles, hAllIsocs, hLines, indConj, poles);
		poles(:,3) = -poles(:,3);		% Because in the automatic polygons mode the poles were fetch from Fix plate
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
	end

% -----------------------------------------------------------------------------------------------------------------
function do_reconst(obj, evt, poles, name, hPlateFix, hPlateMov)
% Callback function that does the reconstruction work
	if (nargin == 4),	hPlateFix = [];		hPlateMov = [];		end
	handles = guidata(obj);
	set(handles.figure1,'pointer','watch'),		pause(0.01)

	% Do this again to allow for the possibility that polygons have been edited
	[plates_fixed, plates_mobile, globalBB, opt_A] = get_plate_polyg(handles, poles, hPlateFix, hPlateMov);

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
		Gr = gmtmex(sprintf('grdrotater %s -N -nn -E%g/%g/%g -F', opt_A{k}, poles(k,1), poles(k,2), poles(k,3)), G, plates_mobile{k});
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
		[rlon,rlat] = rot_euler(plates_mobile{k}(:,1), plates_mobile{k}(:,2), poles(k,1), poles(k,2), poles(k,3), -1);
		mirone('DrawLine_CB', handMirNew, 'data', rlon, rlat)
	end	

% -----------------------------------------------------------------------------------------------------------------
function [ind, pole] = get_conjugates(hAllIsocs, hLine)
% Find the conjugate plate of HLINE. Return its index in HALLISOCS and the finite pole stored in HLINE appdata.

	ind = 0;
	lineInfo = getappdata(hLine, 'LineInfo');	% Ex: 9 NORTH AMERICA/IBERIA FIN"138.2 62 6.132 27.47"
	pole = pole2neighbor('parse_finite_pole', [], [], [], [], lineInfo);
	if (isempty(pole)),		errordlg('The selected anomaly has no FINite rotation info. Bye'),	end
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
		errordlg(sprintf('Could not find the conjugate anomaly -> %s  %s/%s  Bye', isoc, PP2, PP1))
	end

% -----------------------------------------------------------------------------------------------------------------
function [ind, poles] = get_conjugates_polyg(hAllIsocs, hLine)
% Find the conjugate plate of HLINE. Return its index in HALLISOCS and the finite pole stored in HLINE appdata.

	lineInfo = getappdata(hLine, 'LineInfo');   % Ex: 9 NORTH AMERICA/IBERIA FIN"138.2 62 6.132 27.47"
	isoc = strtok(lineInfo);
	[P1, P2] = pole2neighbor('get_plate_pair', [], [], [], [], hLine);

	ind = zeros(10,1);		poles = zeros(10, 3);		k = 1;
	for (n = 1:numel(hAllIsocs))		% Loop over all Isochrons
		this_lineInfo  = getappdata(hAllIsocs(n),'LineInfo');
		this_isoc = strtok(this_lineInfo);
		if (strcmpi(this_isoc, isoc))	% OK, same isochron so now search for plates that connect to HLINE (P1 == PP2)
			if (~isempty(strfind(this_lineInfo, P1)))	% Potential target. But only want the P1 == PP2 combination
				[PP1, PP2] = pole2neighbor('get_plate_pair', [], [], [], [], this_lineInfo);
				if (strcmp(P1, PP2))	% Found one target
					p = pole2neighbor('parse_finite_pole', [], [], [], [], this_lineInfo);
					if (isempty(p)),	errordlg('The conjugated polygon has no FINite rotation info. Bye'),	end
					closed = aux_funs('isclosed', hAllIsocs(n), 0.5);
					if (~closed),	errordlg('This conjugated polygon is NOT closed as it must. Bye'),	end
					poles(k, :) = [p.lon p.lat p.ang];
					ind(k) = n;		k = k + 1;
				end
			end
		end
	end
	poles(ind == 0, :) = [];
	ind(ind == 0) = [];			% Remove empties
	if (isempty(ind)),	errordlg(sprintf('Could not find the paired polygons to %s. Bye', isoc)),	end

% -----------------------------------------------------------------------------------------------------------------
function [plates_fixed, plates_mobile] = build_plate_polyg(handles, hAllIsocs, hLines, indConj, poles)
% Construct the polygons of each pseudo-plate from each isochron. Do this crudly from an isoc outward rotation
% Compute also the global limits (in BB) of the final grid holding the fix and rotated plates.

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
function [plates_fixed, plates_mobile, BB, opt_A] = get_plate_polyg(handles, poles, hPlateFix, hPlateMov)
% BB    Global BB
% OPT_A holds the limits (-R) of the rotated plates as sub-grids of the global grid.

	if (~isempty(hPlateFix))		% Case of given polygons
		hFix = hPlateFix;		hMov = hPlateMov;
		nPlates = numel(hFix);
	else
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
	end

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
		                        poles(k,1), poles(k,2), poles(k,3), -1);
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
		ang = -8;
	else
		ang = 8;
	end

	[rlon,rlat] = rot_euler(x, y, pole(1), pole(2), ang, -1);
	D = gmtmex('gmtsimplify -T5m -fg', [rlon rlat]);		% DP simplify it. Usefull if we want line-edit it.
	rlon = D.data(:,1);		rlat = D.data(:,2);
	xP = [x; rlon(end:-1:1); x(1)+100*eps];			% Adding the 'eps' to prevent that polygons are plotted as patches
	yP = [y; rlat(end:-1:1); y(1)];
