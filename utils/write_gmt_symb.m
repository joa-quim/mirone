function write_gmt_symb(handles)
% Save vector contents of the Mirone fig as a GMT custom symbol

%	Copyright (c) 2004-2016 by J. Luis
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

% $Id: write_gmt_symb.m 7742 2016-01-08 18:09:32Z j $

	ALLlineHand = findobj(get(handles.axes1,'Child'),'Type','line');
	ALLpatchHand = findobj(get(handles.axes1,'Child'),'Type','patch');

	[FileName,PathName,handles] = put_or_get_file(handles, ...
		{'*.def','GMT symbol def file (*.def)'; '*.*', 'All Files (*.*)'},'Select output .def name','put','.def');
	if isequal(FileName,0),		return,		end
	fname = [PathName FileName];

	n_patch = numel(ALLpatchHand);
	n_pline = numel(ALLlineHand);

	% ---- Compute the scale and offset needed to bring all elements to the [-0.5 0.5] interval --------
	scale = 1;		shift_x = 0;	shift_y = 0;
	if (~isempty(ALLlineHand))      % 
		min_x = zeros(n_pline,1);		max_x = zeros(n_pline,1);
		min_y = zeros(n_pline,1);		max_y = zeros(n_pline,1);
		for (k = 1:n_pline)
			min_x(k) = min(get(ALLlineHand(k),'XData'));	max_x(k) = max(get(ALLlineHand(k),'XData'));
			min_y(k) = min(get(ALLlineHand(k),'YData'));	max_y(k) = max(get(ALLlineHand(k),'YData'));
		end
		shift_x = min(min_x);			shift_y = min(min_y);
		dx = max(max_x) - min(min_x);	dy = max(max_y) - min(min_y);
		d = max(dx, dy);
	end
	if (~isempty(ALLpatchHand))
		min_x = zeros(n_patch,1);		max_x = zeros(n_patch,1);
		min_y = zeros(n_patch,1);		max_y = zeros(n_patch,1);
		for (k = 1:n_patch)
			min_x(k) = min(get(ALLpatchHand(k),'XData'));	max_x(k) = max(get(ALLpatchHand(k),'XData'));
			min_y(k) = min(get(ALLpatchHand(k),'YData'));	max_y(k) = max(get(ALLpatchHand(k),'YData'));
		end
		shift_x = min(shift_x, min(min_x));		shift_y = min(shift_y, min(min_y));
		dx = max(max_x) - min(min_x);	dy = max(max_y) - min(min_y);
		d = max(d, max(dx, dy));
	end
	scale = 1 / d;
	% ---------------------------------------------------------------------------------------------------

	fid = fopen(fname,'wt');
	fprintf(fid,'#\n# Definition File created by Mirone Tech\n#\n');

	% ------------- Search for closed polygons ----------------------------
	if (~isempty(ALLpatchHand))
		ALLpatchHand = ALLpatchHand(end:-1:1);	% Don't know if this always good but respects stack order
		xx = get(ALLpatchHand,'XData');     yy = get(ALLpatchHand,'YData');
		LineStyle = get(ALLpatchHand,'LineStyle');
		LineWidth = get(ALLpatchHand,'LineWidth');
		if (iscell(LineWidth)),     LineWidth = cat(1,LineWidth{:});     end
		EdgeColor = get(ALLpatchHand,'EdgeColor');
		if (iscell(EdgeColor)),     EdgeColor = cat(1,EdgeColor{:});     end
		FillColor = get(ALLpatchHand,'FaceColor');
		if (iscell(FillColor))
			resp = strmatch('none',char(FillColor{:}));
			if (isempty(resp))
				FillColor = cat(1,FillColor{:});
			else
				for (i=1:length(resp))                  % Signal down that this is a non colored polygon
					FillColor{resp(i)} = [-1 -1 -1];    % FDS it worked. I wonder why FillColor{resp} = repmat([-1 -1 -1],length(resp),1); DOESN'T
				end
				FillColor = cat(1,FillColor{:});
			end
		else                % We have only one patch
			xx = num2cell(xx,1);   yy = num2cell(yy,1);   % Make it a cell for reducing the head-hakes
			resp = strmatch('none',FillColor);
			if (~isempty(resp))
				FillColor = [-1 -1 -1];                 % Signal down that this is a non colored polygon
			end
		end

		for (i = 1:n_patch)
			cor_edge = round(EdgeColor(i,1:3) * 255);
			cor_edge = sprintf('%d/%d/%d', cor_edge(1), cor_edge(2), cor_edge(3));
			cor_fill = round(FillColor(i,1:3) * 255);
			if (cor_fill(1) >= 0)       % Color filled polygon
				cor_fill = sprintf('%d/%d/%d', cor_fill(1), cor_fill(2), cor_fill(3));
				mlt_comm = [' -G' cor_fill ' -W' num2str(LineWidth(i)) 'p,' cor_edge];
			else                        % No filling color
				mlt_comm = [' -W' num2str(LineWidth(i)) 'p,' cor_edge];
			end

			x = (xx{i}(:)' - shift_x) * scale - 0.5;	y = (yy{i}(:)' - shift_y) * scale - 0.5;	% Shift-n-scale
			fprintf(fid,'%f\t%f\tM\t%s\n', x(1), y(1), mlt_comm);
			fprintf(fid,'%f\t%f\tD\n', [x(2:end); y(2:end)]);
			if (i < n_pline)
				fprintf(fid,'\n#\n');
			end
		end
	end

	% ------------- Search for lines or polylines ----------------------------
	if (~isempty(ALLlineHand))      % 
		LineWidth = get(ALLlineHand,'LineWidth');
		if (iscell(LineWidth)),     LineWidth = cat(1,LineWidth{:});     end
		LineColor = get(ALLlineHand,'Color');
		if (iscell(LineColor)),     LineColor = cat(1,LineColor{:});     end
		for (k = 1:n_pline)
			if (strcmp(get(ALLlineHand(k),'Tag'),'circleCart'))
				xyRad = getappdata(ALLlineHand(k),'LonLatRad');
				fprintf(fid,'%f\t%f\t%.4f\tc\n', (xyRad(1) - shift_x) * scale - 0.5, ...
						(xyRad(2) - shift_y) * scale - 0.5, xyRad(3)*2 * scale);
			else
				lcor = round(LineColor(k,1:3) * 255);
				lcor = sprintf('%d/%d/%d', lcor(1), lcor(2), lcor(3));
				mlt_comm = [' -W' num2str(LineWidth(k)) 'p,' lcor];
				x = get(ALLlineHand(k),'XData');		y = get(ALLlineHand(k),'YData');
				x = (x(:)' - shift_x) * scale - 0.5;	y = (y(:)' - shift_y) * scale - 0.5;	% Shift-n-scale
 				fprintf(fid,'%f\t%f\tM\t%s\n', x(1), y(1), mlt_comm);
 				fprintf(fid,'%f\t%f\tD\n', [x(2:end); y(2:end)]);
				%fprintf(fid,'%f\t%f\n', [x; y]);
				if ( (x(end) ~= x(1)) || (y(end) ~= y(1)) )		% An open polyline
					fprintf(fid,'S\n');
				end
				if (k < n_pline)
					fprintf(fid,'\n#\n');
				end
			end
		end	
	end
		
	fclose(fid);
