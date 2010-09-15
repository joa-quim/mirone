/*--------------------------------------------------------------------
 *	$Id: grdproject.c,v 1.13 2004/09/15 00:22:39 pwessel Exp $
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
 * grdproject reads a geographical grdfile and evaluates the grid at new grid positions
 * specified by the map projection and new dx/dy values using a weighted average of all
 * points within the search radius. Optionally, grdproject may perform the inverse
 * transformation, going from rectangular coordinates to geographical.
 *
 * Author:	Paul Wessel
 * Date:	15-JUL-2000
 * Ver:		4
 *
 *
 *--------------------------------------------------------------------*/
/* Mexified version of grdproject
 * Author:	J. Luis
 * Date: 	18-Sep-2004
 *
 * Usage
 * Zout = grdproject_m(Zin,head,'params');
 *
 * where	Zin is the array containing the input file
 *		head is the header descriptor in the format of grdread/grdwrite mexs
 *		and params are any of the parameters used in the gmt grdproject (e.g.: '-J...','-R...')
 *
 * IMPORTANT NOTE. The data type of Zin is preserved in Zout. That means you can send Zin
 * as a double, single, Int32, Int16, Uint16 or Uint8 and receive Zout in one of those types
 *
 *	 
 *		04/06/06 J Luis, Updated to compile with version 4.1.3
 *		14/10/06 J Luis, Now includes the memory leak solving solution
 *		02/07/07 J Luis, Fixed Bug on the "head[3]" return value
 *		21/02/08 J Luis, Patched version to compile with 4.1.2
 *		11/03/08 J Luis, Memory leak protection seams to make it hang on second call
 */

#include "gmt.h"
#include "mex.h"

/* int GMTisLoaded = FALSE;	/* Used to know wether GMT stuff is already in memory or not */
float *grd_out;

/* --------------------------------------------------------------------------- */
/* Matlab Gateway routine */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	int i, j, k, i_pad, dpi, nx, ny, nm, unit = 0;
	
	int error = FALSE, inverse = FALSE, n_set = FALSE, set_n = FALSE, one_to_one = FALSE;
	int d_set = FALSE, e_set = FALSE, m_set = FALSE, map_center = FALSE, offset, toggle_offset = FALSE;
	int bilinear = FALSE, shift_xy = FALSE;
	int is_double = FALSE, is_single = FALSE, is_int32 = FALSE, is_int16 = FALSE;
	int is_uint16 = FALSE, is_uint8 = FALSE;
	
	double w, e, s, n, x_inc = 0.0, y_inc = 0.0, max_radius = 0.0, one_or_zero;
	double xmin, xmax, ymin, ymax, inch_to_unit, unit_to_inch, fwd_scale, inv_scale;
	double false_easting = 0.0, false_northing = 0.0;
	
	char unit_name[80], scale_unit_name[80];
	
	int	argc = 0, n_arg_no_char = 0, n_used = 0, i2, nx_in, ny_in, mx, nc_h, nr_h, *i_4, *o_i4, *pdata_i4;
	short int *i_2, *o_i2, *pdata_i2;
	unsigned short int *ui_2, *o_ui2, *pdata_ui2;
	char	**argv;
	unsigned char *ui_1, *o_ui1, *pdata_ui1;
	float	*grd_in, *z_4, *pdata_s, *o_s;
	double	*pdata, *pdata_d, *z_8, *head, *o_d, major_axis, flat, head_o[9];

	struct GRD_HEADER g_head, r_head;
	struct GMT_EDGEINFO edgeinfo;

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
	argv[0] = "grdproject";
	for (i = 1; i < argc; i++) {
		argv[i] = (char *)mxArrayToString(prhs[i+n_arg_no_char-1]);
	}
	
	/*if (!GMTisLoaded) {
		argc = GMT_begin (argc, argv);
		GMTisLoaded = TRUE;
	}
	else
		argc = GMT_short_begin (argc, argv);*/
	argc = GMT_begin (argc, argv);
	
	w = e = s = n = 0.0;
	nx = ny = 0;

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {
				/* Common parameters */
			
				case 'J':
				case 'R':
				case '\0':
					error += GMT_get_common_args (argv[i], &w, &e, &s, &n);
					break;
				
				/* Supplemental parameters */
				
				case 'C':
					map_center = TRUE;
					if (argv[i][2]) {	/* Also gave shifts */
	 					n = sscanf (&argv[i][2], "%lf/%lf", &false_easting, &false_northing);
						if (n != 2) {
							mexPrintf ("%s: GMT SYNTAX ERROR.  Expected -C[<false_easting>/<false_northing>]\n", "grdproject");
							error++;
						}
						shift_xy = TRUE;
					}
					break;
				case 'D':
					GMT_getinc (&argv[i][2], &x_inc, &y_inc);
					d_set = TRUE;
					break;
				case 'E':
					dpi = atoi (&argv[i][2]);
					e_set = TRUE;
					break;
				case 'A':
					one_to_one = TRUE;
					unit = GMT_check_scalingopt ('A',argv[i][2], scale_unit_name);
					break;
				case 'F':
					toggle_offset = TRUE;
					break;
				case 'I':
					inverse = TRUE;
					break;
				case 'M':	/* Directly specify units */
					GMT_set_measure_unit (argv[i][2]);
					m_set = TRUE;
					break;
				case 'N':
					sscanf (&argv[i][2], "%d/%d", &nx, &ny);
					if (ny == 0) ny = nx;
					n_set = TRUE;
					break;
				default:
					error = TRUE;
					break;
			}
		}
	}
	
	if ((d_set + e_set + n_set) == 0) n_set = set_n = TRUE;

	if (argc == 1 || error) {
		mexPrintf("grdproject %s - Project geographical grid to/from rectangular grid\n\n", GMT_VERSION);
		mexPrintf("usage: [out,hdr] = grdproject_m(grdfile, grdhead, '-J<parameters>', '-R<west/east/south/north>'\n");
		mexPrintf("\t'[-A[k|m|n|i|c|p]]', '[-C]', '[-D<dx[m|c]>[/dy[m|c]]]', '[-E<dpi>]', '[-F]', '[-I]',\n");
		mexPrintf("\t'[-Mc|i|m]', '[-N<nx/ny>]',)\n\n");
		
		mexPrintf("\tgrdfile is data set to be transformed\n");
		mexPrintf("\n\tOPTIONS:\n");
		mexPrintf("\t-A force projected values to be in actual meters [Default uses the given map scale]\n");
		mexPrintf("\t   Specify another unit by appending k (km), m (miles), n (nautical miles), i (inch), c (cm), or p (points)\n");
		mexPrintf("\t-C coordinates relative to projection center [Default is relative to lower left corner]\n");
		mexPrintf("\t   Optionally append dx/dy to add (or subtract if -I) (i.e., false easting & northing) [0/0]\n");
		mexPrintf("\t-D sets the grid spacing for the new grid\n");
		mexPrintf("\t-E sets dpi for output grid\n");
		mexPrintf("\t-F toggle between pixel and grid registration  [Default is same as input]\n");
		mexPrintf("\t-I Inverse transformation from rectangular to geographical\n");
		mexPrintf("\t-M Temporarily reset MEASURE_UNIT to be c (cm), i (inch), m (meter), or p (point)\n");
		mexPrintf("\t   Cannot be used if -A is set.\n");
		mexPrintf("\t-N sets the number of nodes for the new grid\n");
		mexPrintf("\t   Only one of -D, -E, and -N can be specified!\n");
		mexPrintf("\t   If none are specified, nx,ny of the input file is used\n");
		
		return;
	}
	
	if (!project_info.region_supplied) {
		mexPrintf("%s: GMT SYNTAX ERROR:  Must specify -R option\n", GMT_program);
		error++;
	}
	if ((m_set + one_to_one) == 2) {
		mexPrintf("%s: GMT SYNTAX ERROR:  Can specify only one of -A and -M\n", GMT_program);
		error++;
	}
	if ((d_set + e_set + n_set) != 1) {
		mexPrintf("%s: GMT SYNTAX ERROR:  Must specify only one of -D, -E, or -N\n", GMT_program);
		error++;
	}
	if (d_set && (x_inc <= 0.0 || y_inc <= 0.0)) {
		mexPrintf("%s: GMT SYNTAX ERROR -D option.  Must specify positive increment(s)\n", GMT_program);
		error++;
	}
	if (n_set && !set_n && (nx <= 0 || ny <= 0)) {
		mexPrintf("%s: GMT SYNTAX ERROR -N option.  Must specify positive integers\n", GMT_program);
		error++;
	}
	if (e_set && dpi <= 0) {
		mexPrintf("%s: GMT SYNTAX ERROR -E option.  Must specify positive dpi\n", GMT_program);
		error++;
	}
	
	if (error) return;

	if (nlhs < 2)
		mexErrMsgTxt("GRDPROJECT ERROR: Must provide two outputs.\n");

	/* Find out in which data type was given the input array */
	if (mxIsDouble(prhs[0])) {
		z_8 = mxGetPr(prhs[0]);
		is_double = TRUE;
	}
	else if (mxIsSingle(prhs[0])) {
		z_4 = mxGetData(prhs[0]);
		is_single = TRUE;
	}
	else if (mxIsInt32(prhs[0])) {
		i_4 = mxGetData(prhs[0]);
		is_int32 = TRUE;
	}
	else if (mxIsInt16(prhs[0])) {
		i_2 = mxGetData(prhs[0]);
		is_int16 = TRUE;
	}
	else if (mxIsUint16(prhs[0])) {
		ui_2 = mxGetData(prhs[0]);
		is_uint16 = TRUE;
	}
	else if (mxIsUint8(prhs[0])) {
		ui_1 = mxGetData(prhs[0]);
		is_uint8 = TRUE;
	}
	else {
		mexPrintf("GRDPROJECT ERROR: Unknown input data type.\n");
		mexErrMsgTxt("Valid types are:double, single, In32, In16, UInt16 and Uint8.\n");
	}

	nx_in = mxGetN (prhs[0]);
	ny_in = mxGetM (prhs[0]);
	if (!mxIsNumeric(prhs[0]) || ny_in < 2 || nx_in < 2)
		mexErrMsgTxt("GRDPROJECT ERROR: First argument must contain a decent array\n");

	nc_h = mxGetN (prhs[1]);
	nr_h = mxGetM (prhs[1]);
	if (!mxIsNumeric(prhs[1]) || nr_h > 1 || nc_h < 9)
		mexErrMsgTxt("GRDPROJECT ERROR: Second argument must contain a valid header of the input array.\n");

	head  = mxGetPr(prhs[1]);		/* Get header info */

	GMT_init_scales (unit, &fwd_scale, &inv_scale, &inch_to_unit, &unit_to_inch, unit_name);

	GMT_map_setup (w, e, s, n);
	
	xmin = (map_center) ? project_info.xmin - project_info.x0 : project_info.xmin;
	xmax = (map_center) ? project_info.xmax - project_info.x0 : project_info.xmax;
	ymin = (map_center) ? project_info.ymin - project_info.y0 : project_info.ymin;
	ymax = (map_center) ? project_info.ymax - project_info.y0 : project_info.ymax;
	if (one_to_one) {	/* Convert to chosen units */
		strncpy (unit_name, scale_unit_name, 80);
		xmin /= project_info.x_scale;
		xmax /= project_info.x_scale;
		ymin /= project_info.y_scale;
		ymax /= project_info.y_scale;
		if (unit) {	/* Change the 1:1 unit used */
			xmin *= fwd_scale;
			xmax *= fwd_scale;
			ymin *= fwd_scale;
			ymax *= fwd_scale;
		}
	}
	else {	/* Convert inches to chosen MEASURE */
		xmin *= inch_to_unit;
		xmax *= inch_to_unit;
		ymin *= inch_to_unit;
		ymax *= inch_to_unit;
	}
	if (shift_xy) {
		xmin += false_easting;
		xmax += false_easting;
		ymin += false_northing;
		ymax += false_northing;
	}

	GMT_grd_init (&r_head, argc, argv, FALSE);
	GMT_grd_init (&g_head, argc, argv, FALSE);

	if (inverse) {
		r_head.x_min = head[0];		r_head.x_max = head[1];
		r_head.y_min = head[2];		r_head.y_max = head[3];
		r_head.z_min = head[4];		r_head.z_max = head[5];
		r_head.x_inc = head[7];		r_head.y_inc = head[8];
		r_head.nx = nx_in;		r_head.ny = ny_in;
		r_head.node_offset = irint(head[6]);
	}
	else {
		g_head.x_min = head[0];		g_head.x_max = head[1];
		g_head.y_min = head[2];		g_head.y_max = head[3];
		g_head.z_min = head[4];		g_head.z_max = head[5];
		g_head.x_inc = head[7];		g_head.y_inc = head[8];
		g_head.nx = nx_in;		g_head.ny = ny_in;
		g_head.node_offset = irint(head[6]);
	}

	mx = nx_in + 4;		i_pad = 2;
	GMT_pad[0] = GMT_pad[1] = GMT_pad[2] = GMT_pad[3] = 2;
	nm = (nx_in + 4) * (ny_in + 4);

	grd_in = mxMalloc (nm * sizeof (float));
	/* Transpose from Matlab orientation to gmt grd orientation */
	if (is_double) {
		for (i = 0, i2 = ny_in - 1; i < ny_in; i++, i2--) {
			k = i2 * mx + i_pad * (mx + 1); 
			for (j = 0; j < nx_in; j++) grd_in[j + k] = (float)z_8[j*ny_in+i];
		}
	}
	else if (is_single) {
		for (i = 0, i2 = ny_in - 1; i < ny_in; i++, i2--) {
			k = i2 * mx + i_pad * (mx + 1); 
			for (j = 0; j < nx_in; j++) grd_in[j + k] = z_4[j*ny_in+i];
		}
	}
	else if (is_int32) {
		for (i = 0, i2 = ny_in - 1; i < ny_in; i++, i2--) {
			k = i2 * mx + i_pad * (mx + 1); 
			for (j = 0; j < nx_in; j++) grd_in[j + k] = (float)i_4[j*ny_in+i];
		}
	}
	else if (is_int16) {
		for (i = 0, i2 = ny_in - 1; i < ny_in; i++, i2--) {
			k = i2 * mx + i_pad * (mx + 1); 
			for (j = 0; j < nx_in; j++) grd_in[j + k] = (float)i_2[j*ny_in+i];
		}
	}
	else if (is_uint16) {
		for (i = 0, i2 = ny_in - 1; i < ny_in; i++, i2--) {
			k = i2 * mx + i_pad * (mx + 1); 
			for (j = 0; j < nx_in; j++) grd_in[j + k] = (float)ui_2[j*ny_in+i];
		}
	}
	else if (is_uint8) {
		for (i = 0, i2 = ny_in - 1; i < ny_in; i++, i2--) {
			k = i2 * mx + i_pad * (mx + 1); 
			for (j = 0; j < nx_in; j++) grd_in[j + k] = (float)ui_1[j*ny_in+i];
		}
	}

	if (inverse) {	/* Transforming from rectangular projection to geographical */
	
		g_head.x_min = w;	g_head.x_max = e;	g_head.y_min = s;	g_head.y_max = n;
	
		GMT_boundcond_init (&edgeinfo);

		offset = r_head.node_offset;		/* Same as input */
		if (toggle_offset) offset = !offset;	/* Toggle */
		one_or_zero = (offset) ? 0 : 1;
		if (set_n) {
			nx = r_head.nx;
			ny = r_head.ny;
		}
		GMT_grdproject_init (&g_head, x_inc, y_inc, nx, ny, dpi, offset);
		grd_out = (float *) GMT_memory (VNULL, (size_t)(g_head.nx * g_head.ny), sizeof (float), GMT_program);

		/* Modify input rect header if -A, -C, -M have been set */
		
		if (shift_xy) {
			r_head.x_min -= false_easting;
			r_head.x_max -= false_easting;
			r_head.y_min -= false_northing;
			r_head.y_max -= false_northing;
		}
		if (one_to_one) {	/* Convert from 1:1 scale */
			if (unit) {	/* Undo the 1:1 unit used */
				r_head.x_min *= inv_scale;
				r_head.x_max *= inv_scale;
				r_head.y_min *= inv_scale;
				r_head.y_max *= inv_scale;
				max_radius *= inv_scale;
			}
			r_head.x_min *= project_info.x_scale;
			r_head.x_max *= project_info.x_scale;
			r_head.y_min *= project_info.y_scale;
			r_head.y_max *= project_info.y_scale;
			max_radius *= project_info.x_scale;
		}
		else if (gmtdefs.measure_unit != GMT_INCH) {	/* Convert from inch to whatever */
			r_head.x_min *= unit_to_inch;
			r_head.x_max *= unit_to_inch;
			r_head.y_min *= unit_to_inch;
			r_head.y_max *= unit_to_inch;
			max_radius *= unit_to_inch;
		}
		if (map_center) {	/* Then correct so lower left corner is (0,0) */
			r_head.x_min += project_info.x0;
			r_head.x_max += project_info.x0;
			r_head.y_min += project_info.y0;
			r_head.y_max += project_info.y0;
		}
		r_head.x_inc = (r_head.x_max - r_head.x_min) / (r_head.nx - one_or_zero);
		r_head.y_inc = (r_head.y_max - r_head.y_min) / (r_head.ny - one_or_zero);
		
		GMT_grd_project (grd_in, &r_head, grd_out, &g_head, &edgeinfo, TRUE, BCR_BICUBIC, 0.5, TRUE);
	}
	else {	/* Forward projection from geographical to rectangular grid */
	
		GMT_boundcond_init (&edgeinfo);

		r_head.x_min = project_info.xmin;	r_head.x_max = project_info.xmax;
		r_head.y_min = project_info.ymin;	r_head.y_max = project_info.ymax;
		if (one_to_one) {	/* Convert from 1:1 scale */
			if (unit) {	/* Undo the 1:1 unit used */
				x_inc *= inv_scale;
				y_inc *= inv_scale;
				max_radius *= inv_scale;
			}
			x_inc *= project_info.x_scale;
			y_inc *= project_info.y_scale;
			max_radius *= project_info.x_scale;
		}
		else if (gmtdefs.measure_unit != GMT_INCH) {	/* Convert from inch to whatever */
			x_inc *= unit_to_inch;
			y_inc *= unit_to_inch;
			max_radius *= unit_to_inch;
		}
		if (set_n) {
			nx = g_head.nx;
			ny = g_head.ny;
		}
		
		offset = g_head.node_offset;		/* Same as input */
		if (toggle_offset) offset = !offset;	/* Toggle */
		one_or_zero = (offset) ? 0 : 1;

		GMT_grdproject_init (&r_head, x_inc, y_inc, nx, ny, dpi, offset);
		grd_out = (float *) GMT_memory (VNULL, (size_t)(r_head.nx * r_head.ny), sizeof (float), GMT_program);
		GMT_grd_project (grd_in, &g_head, grd_out, &r_head, &edgeinfo, TRUE, BCR_BICUBIC, 0.5, FALSE);
		
		/* Modify output rect header if -A, -C, -M have been set */

		if (map_center) {	/* Change origin from lower left to projection center */
			r_head.x_min -= project_info.x0;
			r_head.x_max -= project_info.x0;
			r_head.y_min -= project_info.y0;
			r_head.y_max -= project_info.y0;
		}
		if (one_to_one) {	/* Convert to 1:1 scale */
			r_head.x_min /= project_info.x_scale;
			r_head.x_max /= project_info.x_scale;
			r_head.y_min /= project_info.y_scale;
			r_head.y_max /= project_info.y_scale;
			if (unit) {	/* Change the 1:1 unit used */
				r_head.x_min *= fwd_scale;
				r_head.x_max *= fwd_scale;
				r_head.y_min *= fwd_scale;
				r_head.y_max *= fwd_scale;
			}
		}
		else if (gmtdefs.measure_unit != GMT_INCH) {	/* Convert from inch to whatever */
			r_head.x_min *= unit_to_inch;
			r_head.x_max *= unit_to_inch;
			r_head.y_min *= unit_to_inch;
			r_head.y_max *= unit_to_inch;
		}
		if (shift_xy) {
			r_head.x_min += false_easting;
			r_head.x_max += false_easting;
			r_head.y_min += false_northing;
			r_head.y_max += false_northing;
		}
		r_head.x_inc = (r_head.x_max - r_head.x_min) / (r_head.nx - one_or_zero);
		r_head.y_inc = (r_head.y_max - r_head.y_min) / (r_head.ny - one_or_zero);

		/* rect xy values are here in GMT projected units chosen by user */
	}
	
	/* Transpose from gmt grd orientation to Matlab orientation */
	if (inverse) {			/* nx & ny maight have been changed by the -E<dpi> option */
		nx = g_head.nx;		ny = g_head.ny;
	}
	else {
		nx = r_head.nx;		ny = r_head.ny;
	}
	
	mxFree(grd_in);

	if (is_double) {
		plhs[0] = mxCreateDoubleMatrix (ny,nx, mxREAL);
		o_d = (double *)mxGetData(plhs[0]);
		for (i = 0; i < ny; i++) 
			for (j = 0; j < nx; j++) o_d[j*ny+ny-i-1] = (double)grd_out[i*nx+j];
	}
	else if (is_single) {
		plhs[0] = mxCreateNumericMatrix (ny,nx,mxSINGLE_CLASS,mxREAL);
		o_s = (float *)mxGetData(plhs[0]);
		for (i = 0; i < ny; i++) 
			for (j = 0; j < nx; j++) o_s[j*ny+ny-i-1] = grd_out[i*nx+j];
	}
	else if (is_int32) {
		plhs[0] = mxCreateNumericMatrix (ny,nx,mxINT32_CLASS,mxREAL);
		o_i4 = (int *)mxGetData(plhs[0]);
		for (i = 0; i < ny; i++) 
			for (j = 0; j < nx; j++) o_i4[j*ny+ny-i-1] = irint(grd_out[i*nx+j]);
	}
	else if (is_int16) {
		plhs[0] = mxCreateNumericMatrix (ny,nx,mxINT16_CLASS,mxREAL);
		o_i2 = (short int *)mxGetData(plhs[0]);
		for (i = 0; i < ny; i++) 
			for (j = 0; j < nx; j++) o_i2[j*ny+ny-i-1] = (short int)irint(grd_out[i*nx+j]);
	}
	else if (is_uint16) {
		plhs[0] = mxCreateNumericMatrix (ny,nx,mxUINT16_CLASS,mxREAL);
		o_ui2 = (unsigned short int *)mxGetData(plhs[0]);
		for (i = 0; i < ny; i++) 
			for (j = 0; j < nx; j++) o_ui2[j*ny+ny-i-1] = (unsigned short int)irint(grd_out[i*nx+j]);
	}
	else if (is_uint8) {
		plhs[0] = mxCreateNumericMatrix (ny,nx,mxUINT8_CLASS ,mxREAL);
		o_ui1 = (unsigned char *)mxGetData(plhs[0]);
		for (i = 0; i < ny; i++) 
			for (j = 0; j < nx; j++) o_ui1[j*ny+ny-i-1] = (unsigned char)grd_out[i*nx+j];
	}

	GMT_free ((void *)grd_out);
	GMT_end (argc, argv);

	if (inverse) {
		head_o[0] = g_head.x_min;		head_o[1] = g_head.x_max;
		head_o[2] = g_head.y_min;		head_o[3] = g_head.y_max;
		head_o[4] = 0;				head_o[5] = 0;
		head_o[6] = offset;
		head_o[7] = g_head.x_inc;		head_o[8] = g_head.y_inc;
	}
	else {
		head_o[0] = r_head.x_min;		head_o[1] = r_head.x_max;
		head_o[2] = r_head.y_min;		head_o[3] = r_head.y_max;
		head_o[4] = 0;				head_o[5] = 0;
		head_o[6] = offset;
		head_o[7] = r_head.x_inc;		head_o[8] = r_head.y_inc;
	}
	plhs[1] = mxCreateDoubleMatrix (1,9, mxREAL);
	pdata = mxGetPr(plhs[1]);
	memcpy(pdata, head_o, 8*9);
}
