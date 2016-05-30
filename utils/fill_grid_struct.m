function G = fill_grid_struct(Z, head)
% Fill the Grid struct used in gmtmex

% $Id$

	if (~isa(head, 'double')),	head = double(head);	end
	G.ProjectionRefPROJ4 = '';
	G.ProjectionRefWKT = '';	
	G.hdr = head;
	G.range = head(1:6);
	G.inc = head(8:9);
	G.dim = [size(Z,1) size(Z,2)];
	G.n_rows = G.dim(1);
	G.n_columns = G.dim(2);
	G.MinMax = head(5:6);
	G.NoDataValue = NaN;
	G.registration = head(7);
	G.title = '';
	G.remark = '';
	G.command = '';
	G.DataType = 'float32';
	G.LayerCount = 1;
	G.x = linspace(head(1), head(2), G.n_columns);
	G.y = linspace(head(3), head(4), G.n_rows);
	G.z = Z;
	G.x_units = '';
	G.y_units = '';
	G.z_units = '';	
