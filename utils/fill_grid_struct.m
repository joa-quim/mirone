function G = fill_grid_struct(Z, head)
% Fill the Grid struct used in gmtmex

% $Id$

	if (~isa(head, 'double')),	head = double(head);	end
	G.projection_ref_proj4 = '';
	G.projection_ref_wkt = '';	
	G.hdr = head;
	G.range = head(1:6);
	G.inc = head(8:9);
	G.no_data_value = NaN;
	G.registration = head(7);
	G.title = '';
	G.remark = '';
	G.command = '';
	G.datatype = 'float32';
	G.x = linspace(head(1), head(2), size(Z,2));
	G.y = linspace(head(3), head(4), size(Z,1));
	G.z = Z;
	G.x_unit = '';
	G.y_unit = '';
	G.z_unit = '';	
