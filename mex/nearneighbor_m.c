/*--------------------------------------------------------------------
 *	$Id: nearneighbor.c,v 1.19 2004/04/25 09:10:45 pwessel Exp $
 *
 *	Copyright (c) 1991-2004 by P. Wessel and W. H. F. Smith
 *	See COPYING file for copying and redistribution conditions.
 *
 *	This program is free software; you can redistribute it and/or modify
 *	it under the terms of the GNU General Public License as published by
 *	the Free Software Foundation; version 2 of the License.
 *
 *	This program is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *	GNU General Public License for more details.
 *
 *	Contact info: gmt.soest.hawaii.edu
 *--------------------------------------------------------------------*/
/*
 * Based on a specified grid size, nearneighbor reads an xyz file and
 * determines the nearest points to each node in sectors.  The default
 * looks for the nearest point for each quadrant.  The points must also
 * be within a maximum search-radius from the node.  For the nodes that
 * have a full set of nearest neighbors, a weighted average value is
 * computed.  New feature is full support for boundary conditions so
 * that geographic or periodic conditions are explicitly dealt with
 * in the sense that a data point may wrap around to serve as a
 * constraint on the other side of the periodic boundary.
 *
 * Author:	Paul Wessel
 * Date:	14-JUL-2000
 * Version:	4
 *	 
 *		04/06/06 J Luis, Updated to compile with version 4.1.3
 *		14/10/06 J Luis, Now includes the memory leak solving solution
 */
 
#include "gmt.h"
#include "mex.h"

struct NODE {	/* Structure with point id and distance pairs for all sectors */
	float *distance;	/* Distance of nearest datapoint to this node per sector */
	int *datum;		/* Point id of this data point */
};

struct POINT {	/* Structure with input data constraints */
	float x, y, z, w;
};

struct NODE *add_new_node(int n);
void assign_node (struct NODE **node, int n_sector, int sector, double distance, int id);

BOOLEAN GMTisLoaded = FALSE;	/* Used to know wether GMT stuff is already in memory or not */

/* --------------------------------------------------------------------------- */
/* Matlab Gateway routine */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	int i, j, k, ij, i0, j0, *di, dj, n_sectors = 4, sector, n, n_alloc = 5 * GMT_CHUNK, n_fields, nx_2;
	int n_set, n_almost, n_none, n_files = 0, n_args, fno, one_or_zero, n_expected_fields, pad[4], distance_flag = 0;
	int max_di, actual_max_di, ii, jj, x_wrap, y_wrap, ix, iy;

	BOOLEAN go, error = FALSE, done = FALSE, nofile = TRUE;
	BOOLEAN set_empty = FALSE, weighted = FALSE, wrap_180, replicate_x, replicate_y;

	double radius = 0.0, weight, weight_sum, grd_sum, *x0, *y0, dx, dy, delta, distance, factor;
	double *in, *shrink, km_pr_deg, x_left, x_right, y_top, y_bottom, offset, xinc2, yinc2, idx, idy;
	double half_y_width, y_width, half_x_width, x_width, three_over_radius, *info;
	float empty = 0.0, *grd, *z_4, *pdata;
	int	argc = 0, n_arg_no_char = 0, n_cols, nx, ny;
	char	line[BUFSIZ], **argv, buffer[BUFSIZ], *p;

	FILE *fp = NULL;

	struct GRD_HEADER header;
	struct GMT_EDGEINFO edgeinfo;
	struct NODE **grid_node;
	struct POINT  *point;
	
	argc = nrhs;
	for (i = 0; i < nrhs; i++) {		/* Check input to find how many arguments are of type char */
		if(!mxIsChar(prhs[i])) {
			argc--;
			n_arg_no_char++;	/* Number of arguments that have a type other than char */
		}
	}
	argc++;			/* to account for the program's name to be inserted in argv[0] */

	/* get the length of the input string */
	argv = (char **)mxCalloc(argc, sizeof(char *));
	argv[0] = "nearneighbor_m";
	for (i = 1; i < argc; i++) {
		argv[i] = (char *)mxArrayToString(prhs[i+n_arg_no_char-1]);
	}
	
	if (!GMTisLoaded) {
		argc = GMT_begin (argc, argv);
		GMTisLoaded = TRUE;
	}
	else
		argc = GMT_short_begin (argc, argv);

	GMT_boundcond_init (&edgeinfo);
	GMT_grd_init (&header, argc, argv, FALSE);

	pad[0] = pad[1] = pad[2] = pad[3] = 0;
	
	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {
				/* Common parameters */
			
				case 'R':
				case 'f':
				case '\0':
					error += GMT_get_common_args (argv[i], &header.x_min, &header.x_max, &header.y_min, &header.y_max);
					break;
				case 'H':
					gmtdefs.n_header_recs = atoi (&argv[i][2]);
					break;
				case ':':
					gmtdefs.xy_toggle[0] = TRUE;
					break;
				
				/* Supplemental parameters */
				
				case 'I':
					GMT_getinc (&argv[i][2], &header.x_inc, &header.y_inc);
					break;
				case 'E':
					if (!argv[i][2]) {
						mexPrintf ("%s: GMT SYNTAX ERROR -E option:  Must specify value or NaN\n", GMT_program);
						error++;
					}
					else
						empty = (argv[i][2] == 'N' || argv[i][2] == 'n') ? GMT_f_NaN : (float)atof (&argv[i][2]);	
					set_empty = TRUE;
					break;
				case 'F':
					header.node_offset = TRUE;
					break;
				case 'L':
					if (argv[i][2]) {
						error += GMT_boundcond_parse (&edgeinfo, &argv[i][2]);
						if (edgeinfo.gn) {
							GMT_io.in_col_type[0] = GMT_io.out_col_type[0] = GMT_IS_LON;
							GMT_io.in_col_type[1] = GMT_io.out_col_type[1] = GMT_IS_LAT;
						}
					}
					else {
						GMT_io.in_col_type[0] = GMT_io.out_col_type[0] = GMT_IS_LON;
						GMT_io.in_col_type[1] = GMT_io.out_col_type[1] = GMT_IS_LAT;
					}
					break;
				case 'N':
					n_sectors = atoi (&argv[i][2]);
					break;
				case 'S':
					GMT_getinc (&argv[i][2], &radius, &radius);
					if (argv[i][strlen(argv[i])-1] == 'k') distance_flag = 1;
					if (argv[i][strlen(argv[i])-1] == 'K') distance_flag = 2;
					break;
				case 'W':
					weighted = TRUE;
					break;
				default:
					error = TRUE;
					GMT_default_error (argv[i][1]);
					break;
			}
		}
		else {
			if (n_arg_no_char == 0) {
				if ((fp = fopen(argv[i], "r")) == NULL) {
					mexPrintf ("nearneighbor_m: cannot open input data file %s\n", argv[i]);
					mexErrMsgTxt("\n");
				}
			}
		}
	}
	
	if (argc == 1 || GMT_give_synopsis_and_exit) {
		mexPrintf ("nearneighbor %s - A \"Nearest neighbor\" gridding algorithm\n\n", GMT_VERSION);
		mexPrintf("usage: nearneighbor [xyzfile(s)] -I<dx>[m|c][/<dy>[m|c]]\n");
		mexPrintf("\t-N<sectors> -R<west/east/south/north> -S<radius>[m|c|k|K] [-E<empty>] [-F]\n");
		mexPrintf("\t[-H ] [-L<flags>] [-V ] [-W] [-:] [-bi[s][<n>]] [-f[i|o]<colinfo>]\n\n");
		mexPrintf("\t-I sets the grid spacing for the grid.  Append m for minutes, c for seconds.\n");
		mexPrintf("\t-N sets number of sectors. Default is quadrant search [4].\n");
		mexPrintf("\t-S sets search radius in -R, -I units; append m or c for minutes or seconds.\n");
		mexPrintf("\t   Append k for km (implies -R,-I in degrees), use flat Earth approximation.\n");
		mexPrintf("\t   Append K for km (implies -R,-I in degrees), use great circle distances.\n");
		mexPrintf("\n\tOPTIONS:\n");
		mexPrintf("\t-E value to use for empty nodes [Default is NaN].\n");
		mexPrintf("\t-F Force pixel registration [Default is gridline registration].\n");
		mexPrintf("\t-L sets boundary conditions.  <flags> can be either\n");
		mexPrintf("\t   g for geographic boundary conditions, or one or both of\n");
		mexPrintf("\t   x for periodic boundary conditions on x\n");
		mexPrintf("\t   y for periodic boundary conditions on y\n");
		mexPrintf("\t-W input file has observation weights in 4th column.\n");
		mexPrintf("\t   Default is 3 (or 4 if -W is set) columns\n");
		
		return;
	}

	if (!project_info.region_supplied) {
		mexPrintf ("%s: GMT SYNTAX ERROR:  Must specify -R option\n", GMT_program);
		error++;
	}
	if (n_sectors <= 0) {
		mexPrintf ("%s: GMT SYNTAX ERROR -N option:  Must specify a positive number of sectors\n", GMT_program);
		error++;
	}
	if (radius <= 0.0) {
		mexPrintf ("%s: GMT SYNTAX ERROR -S option:  Must specify a positive search radius\n", GMT_program);
		error++;
	}
	if (header.x_inc <= 0.0 || header.y_inc <= 0.0) {
		mexPrintf ("%s: GMT SYNTAX ERROR -I option.  Must specify positive increment(s)\n", GMT_program);
		error++;
	}	
	if (GMT_io.binary[GMT_IN] && gmtdefs.io_header[GMT_IN]) {
		mexPrintf ("%s: GMT SYNTAX ERROR.  Binary input data cannot have header -H\n", GMT_program);
		error++;
	}
	if (error) mexErrMsgTxt("\n");

	/*GMT_grd_RI_verify (&header, 1);*/		/* IF (IVAN == TRUE)  ==> Matlab = BOOM */

	n_expected_fields = 3 + weighted;
	
        if (n_files > 0)
        	nofile = FALSE;
        else
        	n_files = 1;
        n_args = (argc > 1) ? argc : 2;
        
	if (header.node_offset) {
		one_or_zero = 0;
		offset = 0.0;
		xinc2 = 0.5 * header.x_inc;
		yinc2 = 0.5 * header.y_inc;
	}
	else {
		one_or_zero = 1;
		offset = 0.5;
		xinc2 = yinc2 = 0.0;
	}
	
	idx = 1.0 / header.x_inc;
	idy = 1.0 / header.y_inc;

	header.nx = irint ( (header.x_max - header.x_min) * idx) + one_or_zero;
	header.ny = irint ( (header.y_max - header.y_min) * idy) + one_or_zero;
	
	GMT_boundcond_param_prep (&header, &edgeinfo);

	grid_node = (struct NODE **) GMT_memory (VNULL, (size_t)(header.nx * header.ny), sizeof (struct NODE *), GMT_program);
	point = (struct POINT *) GMT_memory (VNULL, (size_t)n_alloc, sizeof (struct POINT), GMT_program);
	
	di = (int *) GMT_memory (VNULL, (size_t)header.ny, sizeof (int), GMT_program);
	shrink = (double *) GMT_memory (VNULL, (size_t)header.ny, sizeof (double), GMT_program);

	x0 = (double *) GMT_memory (VNULL, (size_t)header.nx, sizeof (double), GMT_program);
	y0 = (double *) GMT_memory (VNULL, (size_t)header.ny, sizeof (double), GMT_program);
	for (i = 0; i < header.nx; i++) x0[i] = header.x_min + i * header.x_inc + xinc2;
	for (j = 0; j < header.ny; j++) y0[j] = header.y_max - j * header.y_inc - yinc2;
	if (distance_flag) {	/* Input data is geographical */
		km_pr_deg = 0.001 * 2.0 * M_PI * gmtdefs.ref_ellipsoid[gmtdefs.ellipsoid].eq_radius / 360.0;
		max_di = (int) (ceil (header.nx / 2.0) + 0.1);
		actual_max_di = 0;
		for (j = 0; j < header.ny; j++) {
			shrink[j] = cosd (y0[j]);
			di[j] = (fabs (y0[j]) == 90.0) ? max_di : (int)(ceil (radius / (km_pr_deg * header.x_inc * shrink[j])) + 0.1);
			if (di[j] > max_di) di[j] = max_di;
			if (di[j] > actual_max_di) actual_max_di = di[j];
		}
		dj = (int) (ceil (radius / (km_pr_deg * header.y_inc)) + 0.1);
	}
	else {	/* Plain Cartesian data */
		max_di = (int) (ceil (radius * idx) + 0.1);
		for (j = 0; j < header.ny; j++) di[j] = max_di;
		dj = (int) (ceil (radius * idy) + 0.1);
		actual_max_di = max_di;
	}
	factor = n_sectors / (2.0 * M_PI);

	x_left = header.x_min - actual_max_di * header.x_inc;	x_right = header.x_max + actual_max_di * header.x_inc;
	y_top = header.y_max + dj * header.y_inc;	y_bottom = header.y_min - dj * header.y_inc;
	x_width = header.x_max - header.x_min;		y_width = header.y_max - header.y_min;
	half_x_width = 0.5 * x_width;			half_y_width = 0.5 * y_width;
	nx_2 = edgeinfo.nxp / 2;
	n = 0;
	replicate_x = (edgeinfo.nxp && !header.node_offset);	/* Gridline registration has duplicate column */
	replicate_y = (edgeinfo.nyp && !header.node_offset);	/* Gridline registration has duplicate row */
	x_wrap = header.nx - 1;			/* Add to node index to go to right column */
	y_wrap = (header.ny - 1) * header.nx;	/* Add to node index to go to bottom row */
	
	in = (double *) calloc ((size_t)(4), sizeof(double));

	/* Use here the old toggle recipe because things inside GMT_get_common_args are behaving odly */
	if (gmtdefs.xy_toggle[0]) {
		ix = 1;		iy = 0;
	}
	else {
		ix = 0;		iy = 1;
	}
	for (i = 0; i < gmtdefs.n_header_recs; i++) fgets (line, 1024, fp);

	while (fgets (line, 1024, fp)) {
		if (n_cols == 0) {	/* First time, allocate # of columns */
			strcpy (buffer, line);
			p = (char *)strtok (buffer, " \t\n");
			while (p) {	/* Count # of fields */
				n_cols++;
				p = (char *)strtok ((char *)NULL, " \t\n");
			}
			/* Now we know # of columns */
			if (!weighted && n_cols < 3) {
				mexPrintf("NEARNEIGHBOR ERROR: input file must have at least 3 columns");
				mexErrMsgTxt("\n");
			}
			else if (weighted && n_cols < 4) {
				mexPrintf("NEARNEIGHBOR ERROR: input file must have at least 4 columns");
				mexErrMsgTxt("\n");
			}
		}
		p = (char *)strtok (line, " \t\n");
		jj = 0;
		while (p && jj < n_expected_fields) {
			sscanf (p, "%lf", &in[jj]);
			jj++;
			p = (char *)strtok ((char *)NULL, " \t\n");
		}
		if (jj != n_expected_fields) {
			mexPrintf ("Expected %d but found %d fields in record # %d\n", n_expected_fields, jj, k);
			continue;
		}
	
		if (in[ix] < x_left || in[ix] > x_right) continue;
		if (in[iy] < y_bottom || in[iy] > y_top) continue;
		
		point[n].x = (float)in[ix];
		point[n].y = (float)in[iy];
		point[n].z = (float)in[2];
		if (weighted) point[n].w = (float)in[3];
	
		/* Find indices of the node closest to this data point */

		i0 = GMT_x_to_i(in[0], header.x_min, idx, offset, header.nx);
		j0 = GMT_y_to_j(in[1], header.y_min, idy, offset, header.ny);

		/* Loop over all nodes within radius of this node */

		for (j = j0 - dj; j <= (j0 + dj); j++) {
			jj = j;
			if (GMT_y_out_of_bounds (&jj, &header, &edgeinfo, &wrap_180)) continue;	/* Outside y-range */

			for (i = i0 - di[jj]; i <= (i0 + di[jj]); i++) {
				ii = i;
				if (GMT_x_out_of_bounds (&ii, &header, &edgeinfo, wrap_180)) continue;	/* Outside x-range */ 

				/* Here, (ii,jj) is index of a node (k) inside the grid */

				k = jj * header.nx + ii;
				dx = in[0] - x0[ii];	dy = in[1] - y0[jj];

				/* Check for wrap-around in x or y.  This should only occur if the
				   search radius is larger than 1/2 the grid width/height so that
				   the shortest distance is going through the periodic boundary.
				   For longitudes the dx obviously cannot exceed 180 (half_x_width)
				   since we could then go the other direction instead.  */
				if (edgeinfo.nxp && fabs (dx) > half_x_width) dx -= copysign (x_width, dx);
				if (edgeinfo.nyp && fabs (dy) > half_y_width) dy -= copysign (y_width, dy);

				switch (distance_flag) {	/* Take different action depending on how we want distances calculated */
					case 0:		/* Cartesian distance */
						distance = hypot (dx, dy);
						break;
					case 1:		/* Flat Earth Approximation */
						distance = km_pr_deg * hypot (dx * shrink[jj], dy);
						break;
					case 2:		/* Full spherical calculation */
						distance = km_pr_deg * GMT_great_circle_dist (x0[ii], y0[jj], in[0], in[1]);
						break;
					default:
						break;
				}

				if (distance > radius) continue;	/* Data constraint is too far from this node */

				/* OK, this point should constrain this node.  Calculate which sector and assign the value */

				sector = ((int)((d_atan2 (dy, dx) + M_PI) * factor)) % n_sectors;
				assign_node (&grid_node[k], n_sectors, sector, distance, n);

				/* With periodic, gridline-registered grids there are duplicate rows and/or columns
				   so we may have to assign the point to more than one node.  The next section deals
				   with this situation.  */

				if (replicate_x) {	/* Must check if we have to replicate a column */
					if (ii == 0) 	/* Must replicate left to right column */
						assign_node (&grid_node[k+x_wrap], n_sectors, sector, distance, n);
					else if (ii == edgeinfo.nxp)	/* Must replicate right to left column */
						assign_node (&grid_node[k-x_wrap], n_sectors, sector, distance, n);
				}
				if (replicate_y) {	/* Must check if we have to replicate a row */
					if (jj == 0)	/* Must replicate top to bottom row */
						assign_node (&grid_node[k+y_wrap], n_sectors, sector, distance, n);
					else if (jj == edgeinfo.nyp)	/* Must replicate bottom to top row */
						assign_node (&grid_node[k-y_wrap], n_sectors, sector, distance, n);
				}
			}
		}
		n++;
		if (n == n_alloc) {
			n_alloc += GMT_CHUNK;
			point = (struct POINT *) GMT_memory ((void *)point, (size_t)n_alloc, sizeof (struct POINT), GMT_program);
		}
	}
				
	point = (struct POINT *) GMT_memory ((void *)point, (size_t)n, sizeof (struct POINT), GMT_program);
	grd = (float *) GMT_memory (VNULL, (size_t)(header.nx * header.ny), sizeof (float), GMT_program);

	/* Compute weighted averages based on the nearest neighbors */
	
	n_set = n_almost = n_none = 0;

	if (!set_empty) empty = GMT_f_NaN;
	three_over_radius = 3.0 / radius;
	
	for (j = ij = 0; j < header.ny; j++) {
		for (i = 0; i < header.nx; i++, ij++) {

			if (!grid_node[ij]) {	/* No nearest neighbors, set to empty and goto next node */
				n_none++;
				grd[ij] = empty;
				continue;
			}

			for (k = 0, go = TRUE; go && k < n_sectors; k++) if (grid_node[ij]->datum[k] < 0) go = FALSE;
			if (!go) { 	/* Not full set of neighbors in all sectors, set to empty and goto next node */
				n_almost++;
				grd[ij] = empty;
				continue;
			}

			/* OK, here we have enough data and need to calculate the weighted value */

			n_set++;
			weight_sum = grd_sum = 0.0;	/* Initialize sums */
			for (k = 0; k < n_sectors; k++) {
				delta = three_over_radius * grid_node[ij]->distance[k];
				weight = 1.0 / (1.0 + delta * delta);	/* This is distance weight */
				if (weighted) weight *= point[grid_node[ij]->datum[k]].w;	/* This is observation weight */
				grd_sum += weight * point[grid_node[ij]->datum[k]].z;
				weight_sum += weight;
			}
			grd[ij] = (float)(grd_sum / weight_sum);
		}
	}
	
	/* Transpose from gmt grd orientation to Matlab orientation */
	nx = header.nx;		ny = header.ny;
	z_4 = mxCalloc (nx*ny, sizeof (float));
	for (i = 0; i < ny; i++) for (j = 0; j < nx; j++) z_4[j*ny+ny-i-1] = grd[i*nx+j];

	plhs[0] = mxCreateNumericMatrix (ny,nx,mxSINGLE_CLASS,mxREAL);
	pdata = mxGetData(plhs[0]);
	memcpy(pdata, z_4, ny*nx * 4);
	mxFree(z_4);

	if (nlhs == 2) {	/* User also wants the header */
		plhs[1] = mxCreateDoubleMatrix (1, 9, mxREAL);
		info = mxGetPr (plhs[1]);
		info[0] = header.x_min;
		info[1] = header.x_max;
		info[2] = header.y_min;
		info[3] = header.y_max;
		info[4] = header.z_min;
		info[5] = header.z_max;
		info[6] = header.node_offset;
		info[7] = header.x_inc;
		info[8] = header.y_inc;
	}

	GMT_free ((void *)grd);
	GMT_free ((void *)point);
	GMT_free ((void *)grid_node);
	GMT_free ((void *)shrink);
	GMT_free ((void *)di);
	GMT_free ((void *)x0);
	GMT_free ((void *)y0);
	
	GMT_end_for_mex (argc, argv);
}

struct NODE *add_new_node(int n)
{
	struct NODE *new;
	
	new = (struct NODE *) GMT_memory (VNULL, (size_t)1, sizeof (struct NODE), GMT_program);
	new->distance = (float *) GMT_memory (VNULL, (size_t)n, sizeof (float), GMT_program);
	new->datum = (int *) GMT_memory (VNULL, (size_t)n, sizeof (int), GMT_program);
	while (n > 0) new->datum[--n] = -1;
	
	return (new);
}

void assign_node (struct NODE **node, int n_sector, int sector, double distance, int id)
{
	/* Allocates node space if not already used and updates the value if closer to node */

	if (!(*node)) *node = add_new_node (n_sector);
	if ((*node)->datum[sector] == -1 || (*node)->distance[sector] > distance) {
		(*node)->distance[sector] = (float)distance;
		(*node)->datum[sector] = id;
	}
}
