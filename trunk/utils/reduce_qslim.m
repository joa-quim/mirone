function [v, f] = reduce_qslim(varargin)
%REDUCE_QSLIM  Reduce number of triangle faces.
%   REDUCE_QSLIM(..., R) reduces the number of faces while trying to preserve the
%   overall shape of the surface.  If R is less than or equal to 1, R is interpreted
%   as a fraction of the original faces; for example, if R is 0.2, 20% of the faces
%   will be kept. If R is greater than 1, then R is the target number of faces.
%
%	[V, F] = REDUCE_QSLIM(handles,...) constructs the faces and vertices from Mirone
%	handles structure.
%
%   [V, F] = REDUCE_QSLIM(F, V, R) uses faces and vertices in arrays F and V.
%
%	In the above, if R is not provided assumes a reduction of .5.
%
%   REDUCE_QSLIM(...,'verbose') prints progress messages to the command window.  

%	Copyright (c) 2004-2012 by J. Luis
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

	v = [];		f = [];
	[faces, verts, reduction, verbose, got_handle, outFname, formato, tipo, base, thickness, msg] = ...
		parseargs(nargin,varargin);
	if (~isempty(msg)),		errordlg(msg, 'Error'),		return,		end

	if (got_handle)
		handles = varargin{1};
		[X,Y,Z,head,nRow,nCol] = load_grd(handles);
		if (isempty(Z)),	return,		end

		if (~isempty(base))
			faces = zeros((nRow-1) * (nCol-1) * 2 + (4*(nCol-1) + 4*(nRow-1) + 2), 3, 'int32');
		elseif (~isempty(thickness))
			faces = zeros((nRow-1) * (nCol-1) * 4 + (4*(nCol-1) + 4*(nRow-1)), 3, 'int32');
		else
			faces = zeros((nRow-1) * (nCol-1) * 2, 3, 'int32');
		end
		k = 1;		kk = 1;

		if (~handles.have_nans)			% Simpler case (don't waist time with ifnan tests)

			for (n = 1:nCol-1)
				j = (n - 1);
				for (m = 1:nRow-1)
					k1 = k + 1 + j;		k2 = k+nRow + j;
					faces(kk,1) = k+j;	faces(kk,2) = k1;		faces(kk,3) = k2;	% NW face
					kk = kk + 1;
					faces(kk,1) = k1;	faces(kk,2) = k1+nRow;	faces(kk,3) = k2;	% SE face
					kk = kk + 1;		k = k + 1;
				end
			end

			nFacets = kk - 1;

			[Xf, Yf] = meshgrid(single(X), single(Y));
			verts = [Xf(:) Yf(:) Z(:)];
		else

			for (n = 1:nCol-1)
				j = (n - 1);
				nanCol = isnan(Z(:,n));
				gotNaNs = any(nanCol); 
				for (m = 1:nRow-1)
					if (gotNaNs && nanCol(m)),		continue,		end
					k1 = k + 1 + j;			k2 = k+nRow + j;
					faces(kk,1) = k+j;		faces(kk,2) = k1;		faces(kk,3) = k2;		% NW face
					kk = kk + 1;
					faces(kk,1) = k1;		faces(kk,2) = k1+nRow;	faces(kk,3) = k2;		% NW face
					kk = kk + 1;			k = k + 1;
				end
				if (gotNaNs && nanCol(nRow))	% Check last row
					faces(kk-1,:) = 0;			% Start a kk-1 because meanwhile kk had advanced 1
					faces(kk-2,:) = 0;
				end
			end

			% We now need to see what's in last col which were not probed above
			nanCol = isnan(Z(:,nCol));
			if (any(nanCol))
				ind = find(nanCol);				ind = 2*ind(:);
				ind = [ind ind-1 ind-2];		% May have out of bounds values
				if (ind(1,3) == 0),				ind(1,3) = 1;			end			% First and last indices
				if (ind(end,3) == 2*nRow),		ind(end,:) = 2*nRow;	end			% need a special check
				ind = unique(ind(:));			% We could easily have several repeated faces
				ind = ind + nRow * (nCol - 2) * 2;
				faces(ind,:) = [];				% Remove them
			end

			% Now remove the faces which have zero indices
			ind = (faces(:,1) == 0);		faces(ind,:) = [];
			nFacets = size(faces,1);
		end
	end

	if (got_handle && ~isempty(base))
		[faces, verts] = closeSurf(faces, verts, X(:), Y(:), nFacets, base, 'base');
	elseif (got_handle && ~isempty(thickness))
		[faces, verts] = closeSurf(faces, verts, X(:), Y(:), nFacets, thickness, 'thick');
	end
	
	numFaces = size(faces,1);
	if (reduction > 0 && reduction <= 1)
		reduction = round(numFaces*reduction);
	end
	reduction = min(numFaces+1, reduction);

	% CALL MEX TO DO THE HARD WORK
	if (reduction ~= numFaces)			% Tests showed that if (nin == nout), nothing is done
		[v, f] = reducep_s(verts, faces, reduction, verbose);
		%[f, v] = reducepatch(double(faces), double(verts), reduction, 'fast', 'verbose');
	else
		v = verts;		f = faces;
	end

	if (~isempty(outFname))
		if (strcmp(formato,'nc'))
			if (got_handle),	geog = handles.geog;
			else				geog = 0;
			end
			saveFV_nc(outFname, f, v, geog)
		elseif (strcmp(formato,'stl') || strcmp(formato,'xyz'))
			write_mex(v, f, outFname, formato, tipo)
		else
			disp('Saving formats other than NC, STL or XYZ are not yet programmed')
		end
	end

% ---------------------------------------------------------------------------------------------
function [faces, verts, r, verbose, got_handle, outFname, formato, tipo, base, thickness, msg] = ...
	parseargs(nin, vargin)

	faces = [];		verts = [];		verbose = 0;	r = [];			got_handle = false;
	msg = [];		outFname = [];	tipo = 'asc';	formato = 'nc';	thickness = [];	base = [];

	% Count number of char argins
	nPV = 0;		nArgNoChar = 0;		n = 1;
	while (~ischar(vargin{n}) && n <= nin)
		nArgNoChar = nArgNoChar + 1;	n = n + 1;
	end
	if (nArgNoChar),	nPV = nin - nArgNoChar;		end

	for (j = (nin - nPV + 1):2:nin)
		if (strcmp(vargin{j}, 'name'))
			outFname = vargin{j+1};
		elseif (strcmp(vargin{j}, 'format'))
			if (strcmpi(vargin{j+1},'nc'))		% Confirm that rquired format is possible
				formato = 'nc';
			elseif (strcmpi(vargin{j+1},'stl'))
				formato = 'stl';
			elseif (strcmpi(vargin{j+1},'xyz'))
				formato = 'xyz';
			else
				msg = 'Unknown output format request';
				return
			end
		elseif (strcmp(vargin{j}, 'type'))
			if (strcmpi(vargin{j+1},'binary'))
				tipo = 'bin';
			end
		elseif (strcmp(vargin{j}, 'close'))		% Close the surface with a const base (or top) surface
			base = vargin{j+1};					% Flat base. Need to test that this value makes sense
		elseif (strcmp(vargin{j}, 'thick'))		% Close the surface withe a constant thickness layer
			thickness = vargin{j+1};
		elseif (strncmp(vargin{j}, 'ver', 3))	% Verbose (we don't care of next j, we take as a YES)
			verbose = 1;
		end
	end
	nin = nArgNoChar;

	if ( nin == 1 && isa(vargin{1}, 'struct') )		% reduce_qslim(handles)
		r = 0.5;		got_handle = true;
	elseif ( (nin == 2) && isa(vargin{1}, 'struct') && (numel(vargin{2}) == 1) )
		r = vargin{2};	got_handle = true;
	elseif ( (nin >= 2) && (size(vargin{1},2) == 3) && (size(vargin{2},2) == 3) )
		faces = vargin{1};
		verts = vargin{2};
		if (nin == 3 && (numel(vargin{3}) == 1))
			r = vargin{3};
		end
	else
		msg = 'reduce_qslim:  Wrong number of input arguments.';
	end
	if (isempty(r)),	r = 0.5;	end

% ----------------------------------------------------------------------------------------------
function [faces, verts] = closeSurf(faces, verts, X, Y, nFilled, base, modo)
% Close a surface either with flat base (MODO == 'base') or constant thickness
%
%	NFILLED is the number of facets already filled in main (top surface)
%	MODO is either	'base'	-> close surface with a flat base
%					'thick'	-> close surface making it a constant thickness layer

	kk = nFilled + 1;
	nCol = numel(X);		nRow = numel(Y);
	mn = nRow * nCol;		k = 1;

	if (strcmp(modo, 'base'))
		newVerts = [X(:) repmat(Y(1),nCol,1) repmat(base,nCol,1); ...				% To S face
					repmat(X(end),nRow,1) Y(:) repmat(base,nRow,1); ...				% To E face
					X(end:-1:1) repmat(Y(end),nCol,1) repmat(base,nCol,1); ...		% To N face
					repmat(X(1),nRow,1) Y(end:-1:1) repmat(base,nRow,1)];			% To W face

		verts = [verts; newVerts];

		for (n = 1:nCol-1)			% Fill the South wall
			faces(kk,1) = (n-1)*nRow+1;	faces(kk,2) = n*nRow+1;		faces(kk,3) = mn + k;		kk = kk + 1;
			faces(kk,1) = mn + k;		faces(kk,2) = n*nRow+1;		faces(kk,3) = mn + k + 1;	kk = kk + 1;	
			k = k + 1;
		end

		mn_ = (nCol-1)*nRow;	k = k + 1;
		for (n = 1:nRow-1)			% Fill the East wall
			faces(kk,1) = mn_ + n;		faces(kk,2) = mn_ + n + 1;	faces(kk,3) = mn + k;		kk = kk + 1;
			faces(kk,1) = mn + k;		faces(kk,2) = mn_ + n + 1;	faces(kk,3) = mn + k + 1;	kk = kk + 1;
			k = k + 1;			
		end

		k = k + 1;
		for (n = nCol:-1:2)			% Fill the North wall
			faces(kk,1) = n*nRow;		faces(kk,2) = (n-1)*nRow;	faces(kk,3) = mn + k;		kk = kk + 1;
			faces(kk,1) = mn + k;		faces(kk,2) = (n-1)*nRow;	faces(kk,3) = mn + k + 1;	kk = kk + 1;
			k = k + 1;			
		end

		k = k + 1;
		for (n = nRow:-1:2)			% Fill the West wall
			faces(kk,1) = n;			faces(kk,2) = n - 1;		faces(kk,3) = mn + k;		kk = kk + 1;
			faces(kk,1) =  mn + k;		faces(kk,2) = n - 1;		faces(kk,3) = mn + k + 1;	kk = kk + 1;
			k = k + 1;			
		end

		% And now close the base (2 more facets)
		faces(kk,1) = mn + 1;				faces(kk,2) = mn + nCol;	faces(kk,3) = mn + 2*nCol + nRow;
		kk = kk + 1;
		faces(kk,1) = mn + 2*nCol + nRow;	faces(kk,2) = mn + nCol;	faces(kk,3) = mn + nCol + nRow;

	elseif (strcmp(modo, 'thick'))		% Make a constant thickness layer

		newVerts = verts;		newVerts(:,3) = newVerts(:,3) + base;
		verts = [verts; newVerts];

		for (n = 1:nCol-1)			% Fill the South wall
			m1 = (n-1)*nRow+1;		m2 = n*nRow+1;
			faces(kk,1) = m1;		faces(kk,2) = m2;		faces(kk,3) = m1 + mn;		kk = kk + 1;
			faces(kk,1) = m1 + mn;	faces(kk,2) = m2;		faces(kk,3) = m2 + mn;		kk = kk + 1;
		end

		mn_ = (nCol-1)*nRow;
		for (n = 1:nRow-1)			% Fill the East wall
			m = mn_ + n;
			faces(kk,1) = m;		faces(kk,2) = m + 1;	faces(kk,3) = m + mn;		kk = kk + 1;
			faces(kk,1) = m + mn;	faces(kk,2) = m + 1;	faces(kk,3) = m + 1 + mn;	kk = kk + 1;			
		end

		for (n = nCol:-1:2)			% Fill the North wall
			m1 = n*nRow;			m2 = (n-1)*nRow;
			faces(kk,1) = m1;		faces(kk,2) = m2;		faces(kk,3) = m1 + mn;		kk = kk + 1;
			faces(kk,1) = m1 + mn;	faces(kk,2) = m2;		faces(kk,3) = m2 + mn;		kk = kk + 1;			
		end

		for (n = nRow:-1:2)			% Fill the West wall
			faces(kk,1) = n;		faces(kk,2) = n - 1;	faces(kk,3) = n + mn;		kk = kk + 1;
			faces(kk,1) = n + mn;	faces(kk,2) = n - 1;	faces(kk,3) = n - 1 + mn;	kk = kk + 1;
		end

		% Make the lower surface equal to top
		tmp = faces(1:nFilled,:);
		cvlib_mex('addS', tmp, mn);				% This works with ALL ML versions, and compiled
		faces(kk:end,:) = tmp;

	end


% ----------------------------------------------------------------------------------------------
function saveFV_nc(fname, faces, verts, is_geog)
% Save the FACES, VERTS matrix as 1D vectors in netCDF format.

	nc_funs('create_empty', fname);

	% ---------------------------- Write the dimensions --------------------------------
	nc_funs('add_dimension', fname, 'dimvertices', size(verts,1) )
	nc_funs('add_dimension', fname, 'dimfaces', size(faces,1) )

	% ---------------------------- Write Globals attributes -----------------------------
	nc_funs('attput', fname, -1, 'title', 'A TIN made after a regular grid');
	str = sprintf(['Set of triangles reduced with reduce_qslim.\nTo recover the triangles use the connection\n' ...
		'id stored in the TriVertInd_A,B,C variables.']);
	nc_funs('attput', fname, -1, 'Description', str);
	%nc_funs('attput', fname, -1, 'spatial_ref', spatial_ref);
	lm = [min(verts) max(verts)];
	limits = [lm(1) lm(4) lm(2) lm(5) lm(3) lm(6)];
	nc_funs('attput', fname, -1, 'BoundingBox', limits);

	% ------------------------------ Write the variables --------------------------------
	if (is_geog)
		long_name = {'Longitude' 'Latitude'};	units = {'degrees_east' 'degrees_north'};
		Xname = 'lon';				Yname = 'lat';
	else
		long_name = {'X' 'Y'};		units = {'xunits' 'yunits'};
		Xname = 'X';				Yname = 'Y';
	end

	Nctype = 5;			% NC_FLOAT
	write_var(fname, Xname, Nctype, 'dimvertices', long_name{1}, units{1})
	nc_funs('varput', fname, Xname, verts(:,1));
	write_var(fname, Yname, Nctype, 'dimvertices', long_name{1}, units{1})
	nc_funs('varput', fname, Yname, verts(:,2));
	write_var(fname, 'z', Nctype, 'dimvertices', 'z', 'yunits')
	nc_funs('varput', fname, 'z', verts(:,3));

	Nctype = 4;			% NC_INT
	write_var(fname, 'TriVertInd_A', Nctype, 'dimfaces', 'indices','numbers')
	nc_funs('varput', fname, 'TriVertInd_A', int32(faces(:,1)));
	write_var(fname, 'TriVertInd_B', Nctype, 'dimfaces', 'indices','numbers')
	nc_funs('varput', fname, 'TriVertInd_B', int32(faces(:,2)));
	write_var(fname, 'TriVertInd_C', Nctype, 'dimfaces', 'indices','numbers')
	nc_funs('varput', fname, 'TriVertInd_C', int32(faces(:,3)));


% ----------------------------------------------------------------------------------------------
function write_var(fname, name, tipo, dim, long_name, units, actual_range, comment, fillValue, missing_value, scale_factor)	
% NARGIN == 3 & NARGIN == 6 are special handled cases
	if (nargin == 3)
		dim = [];	long_name = [];	units = [];	actual_range = [];	comment = [];	fillValue = [];	missing_value = [];	scale_factor = [];
	elseif (nargin == 6)
		actual_range = [];	comment = [];	fillValue = [];	missing_value = [];	scale_factor = [];
	end
	varstruct.Name = name;
	varstruct.Nctype = tipo;
	if (~isempty(dim)),		varstruct.Dimension = {dim};	end
	nc_funs('addvar', fname, varstruct)
	if (~isempty(long_name)),		nc_funs('attput', fname, varstruct.Name, 'long_name', long_name ),	end
	if (~isempty(units)),			nc_funs('attput', fname, varstruct.Name, 'units', units),			end
	if (~isempty(comment)),			nc_funs('attput', fname, varstruct.Name, 'comment', comment),		end
	if (~isempty(actual_range)),	nc_funs('attput', fname, varstruct.Name, 'actual_range', actual_range),	end
	if (~isempty(fillValue)),		nc_funs('attput', fname, varstruct.Name, '_FillValue', fillValue),		end
	if (~isempty(missing_value)),	nc_funs('attput', fname, varstruct.Name, 'missing_value', missing_value),	end
	if (~isempty(scale_factor)),	nc_funs('attput', fname, varstruct.Name, 'scale_factor', scale_factor),	end
