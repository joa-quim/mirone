function write_gmt_symb(handles)
% ...

	ALLlineHand = findobj(get(handles.axes1,'Child'),'Type','line');
	ALLpatchHand = findobj(get(handles.axes1,'Child'),'Type','patch');

	[FileName,PathName,handles] = put_or_get_file(handles, ...
		{'*.def','GMT symbol def file (*.def)'; '*.*', 'All Files (*.*)'},'Select output .def name','put','.def');
	if isequal(FileName,0),		return,		end
	fname = [PathName FileName];

	fid = fopen(fname,'wt');
	fprintf(fid,'#\n# Definition File created by Mirone Tech\n#\n');

	% ------------- Search for closed polygons ----------------------------
	if (~isempty(ALLpatchHand))
		ALLpatchHand = ALLpatchHand(end:-1:1);	% Don't know if this always good but respects stack order
		xx = get(ALLpatchHand,'XData');     yy = get(ALLpatchHand,'YData');
		n_patch = numel(ALLpatchHand);
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
			cor_edge = [num2str(cor_edge(1)) '/' num2str(cor_edge(2)) '/' num2str(cor_edge(3))];
			cor_fill = round(FillColor(i,1:3) * 255);
			if (cor_fill(1) >= 0)       % Color filled polygon
				cor_fill = [num2str(cor_fill(1)) '/' num2str(cor_fill(2)) '/' num2str(cor_fill(3))];
				mlt_comm = [' -G' cor_fill ' -W' num2str(LineWidth(i)) 'p,' cor_edge];
			else                        % No filling color
				mlt_comm = [' -W' num2str(LineWidth(i)) 'p,' cor_edge];
			end

			x = xx{i}(:)';			y = yy{i}(:)';			% Force row vecs
			fprintf(fid,'%f\t%f\tM\t%s\n', x(1), y(1), mlt_comm);
			fprintf(fid,'%f\t%f\tD\n', [x(2:end); y(2:end)]);
			fprintf(fid,'\n#\n');
		end
	end

	% ------------- Search for lines or polylines ----------------------------
	if (~isempty(ALLlineHand))      % 
		for (k = 1:numel(ALLlineHand))
			if (strcmp(get(ALLlineHand(k),'Tag'),'circleCart'))
				xyRad = getappdata(ALLlineHand(k),'LonLatRad');
				fprintf(fid,'%f\t%f\t%.4f\tc\n', xyRad(1), xyRad(2), xyRad(3)*2);
			else
				x = get(ALLlineHand(k),'XData');		y = get(ALLlineHand(k),'YData');
				x = x(:)';			y = y(:)';			% Force row vecs
				fprintf(fid,'%f\t%f\tM\n', x(1), y(1));
				fprintf(fid,'%f\t%f\tD\n', [x(2:end); y(2:end)]);
				if ( (x(end) ~= x(1)) || (y(end) ~= y(1)) )		% An open polyline
					fprintf(fid,'S\n');
				end
				fprintf(fid,'\n#\n');
			end
		end	
	end
		
	fclose(fid);