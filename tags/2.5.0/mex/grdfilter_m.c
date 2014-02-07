/*--------------------------------------------------------------------
 *	$Id$
 *	grdfilter.c,v 1.18 2004/08/14 00:38:19
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
 * grdfilter.c  reads a grdfile and creates filtered grd file
 *
 * user selects boxcar, gaussian, cosine arch, median, or mode filter
 * user selects distance scaling as appropriate for map work, etc.
 *
 * Author:  W.H.F. Smith
 * Date: 	16 Aug 88
 *
 * Modified:	27 Sep 88 by W.H.F. Smith, to use new median routine.
 *
 * Updated:	PW: 13-Jun-1991 to v2.0
 *		PW: 13-Jul-1992 to actually do v2 i/o.
 *		PW: 15-Jul-1992 to handle arbitrary new -R -I
 *		PW: 03-Jan-1995 to offer mode-filter (from filter1d)
 *		WS: 21-May-1998 to make warnings silent unless -V on.
 *		PW: 03-Jun-1998 upgrade to GMT 3.1
 *		PW: 02-Jun-1999 upgrade to GMT 3.3 + added SINCOS option
 *		PW: 18-Oct-1999 Use sincos directly
 *		PW: 18-JUN-2000 3.3.5
 * Version:	4
 *		PW: 27-MAY-2004 Added extreme values filter options l, L, u, U
*/
 /*--------------------------------------------------------------------
 * Mexified version of grdfilter
 * Author:	J. Luis
 * Date: 	2 Nov 2004
 *
 * Usage
 * [Zout,head] = grdfilter_m(Z,head,'params');
 * 	OR
 * Zout = grdfilter_m(Z,head,'params');
 *
 * where	Z is the array containing the input file 
 *		head is the header descriptor in the format of grdread/grdwrite mexs
 *		and params are any of the parameters used in the gmt grdfilter (e.g.: '-F...','-R...')
 *	 
 *		04/06/06 J Luis, Updated to compile with version 4.1.3
 *		14/10/06 J Luis, Now includes the memory leak solving solution
 */
#define debuga	0
 
#include "gmt.h"
#include "mex.h"

double	get_wt();
void	set_weight_matrix(int nx_f, int ny_f, double y_0, double north, double south, double dx, double dy, double f_wid, int f_flag, int d_flag, double x_off, double y_off, int fast);
void DEBUGA(int n);

int	*i_origin;
float	*input, *output;
double	*weight, *work_array, *x_shift;
double	DEG2KM;

char *filter_name[9] = {
	"Boxcar",
	"Cosine Arch",
	"Gaussian",
	"Median",
	"Mode",
	"Lower",
	"Lower+",
	"Upper",
	"Upper-"
};

GMT_LONG GMT_mode_selection = 0, GMT_n_multiples = 0;

/* int GMTisLoaded = FALSE;	/* Used to know wether GMT stuff is already in memory or not */

/* --------------------------------------------------------------------------- */
/* Matlab Gateway routine */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	GMT_LONG	nx_out, ny_out, nx_fil, ny_fil, n_in_median, n_nan = 0;
	GMT_LONG	x_half_width, y_half_width, j_origin, i_out, j_out, i2, nx, ny, k1, k2;
	GMT_LONG	i_in, j_in, ii, jj, i, j, ij_in, ij_out, ij_wt, effort_level;
	GMT_LONG	distance_flag, filter_type, one_or_zero = 1, nr_h, nc_h;
	int	argc = 0, n_arg_no_char = 0, *i_4, *o_i4, *pdata_i4;
	char	**argv, c;
	short int *i_2, *o_i2, *pdata_i2;
	unsigned short int *ui_2, *o_ui2, *pdata_ui2;
	unsigned char *ui_1, *o_ui1, *pdata_ui1;
	
	int	error, new_range, new_increment, fast_way, shift = FALSE, slow, toggle = FALSE, pole_trouble = FALSE;
	int is_double = FALSE, is_single = FALSE, is_int32 = FALSE, is_int16 = FALSE;
	int is_uint16 = FALSE, is_uint8 = FALSE;
	
	float	*z_4, *pdata_s, *o_s;
	double	west_new, east_new, south_new, north_new, dx_new, dy_new, offset;
	double	filter_width, x_scale, y_scale, x_width, y_width;
	double	x_out, y_out, wt_sum, value, last_median, this_median, xincnew2, yincnew2;
	double	xincold2, yincold2, y_shift, x_fix = 0.0, y_fix = 0.0;
	struct	GRD_HEADER h, test_h;
	double	*pdata, *pdata_d, *z_8, *head, *o_d, head_o[9];

	DEG2KM = 2.0 * M_PI * gmtdefs.ref_ellipsoid[GMT_N_ELLIPSOIDS-1].eq_radius / 360.0 * 0.001;	/* GRS-80 sphere degree->m  */

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
	argv[0] = "grdfilter";
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
	
	error = new_range = new_increment = FALSE;
	filter_width = dx_new = dy_new = west_new = east_new = 0.0;
	filter_type = distance_flag = -1;
	
	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {
				/* Common parameters */
			
				case 'R':
				case '\0':
					error += GMT_get_common_args (argv[i], &west_new, &east_new, &south_new, &north_new);
					break;
				
				case 'V':
					gmtdefs.verbose = 1;
					break;
				case 'F':
					switch (argv[i][2]) {
						case 'b':
							filter_type = 0;
							break;
						case 'c':
							filter_type = 1;
							break;
						case 'g':
							filter_type = 2;
							break;
						case 'm':
							filter_type = 3;
							break;
						case 'p':
							filter_type = 4;
							c = argv[i][strlen(argv[i]-1)];
							if (c == '-') GMT_mode_selection = -1;
							if (c == '+') GMT_mode_selection = +1;
							break;
						case 'l':
							filter_type = 5;
							break;
						case 'L':
							filter_type = 6;
							break;
						case 'u':
							filter_type = 7;
							break;
						case 'U':
							filter_type = 8;
							break;
						default:
							error = TRUE;
							break;
					}
					filter_width = atof(&argv[i][3]);
					break;
				case 'D':
					distance_flag = atoi(&argv[i][2]);
					break;
				case 'I':
					GMT_getinc (&argv[i][2], &dx_new, &dy_new);
					new_increment = TRUE;
					break;
				case 'N':	/* Pixel registration OBSOLETE but BACKWARD COMPATIBLE */
					one_or_zero = 0;
					break;
				case 'T':	/* Toggle registration */
					toggle = TRUE;
					break;
				default:
					error = TRUE;
					/*GMT_default_error (argv[i][1]);*/
					break;
			}
		}
	}
	
	if (argc == 1 || GMT_give_synopsis_and_exit) {
		mexPrintf("grdfilter - Filter a 2-D netCDF grdfile in the Time domain\n\n");
		mexPrintf("usage: [Zout,[head]] = grdfilter_m(Z, '-D<distance_flag>', '-F<type><filter_width>[<mode>]',\n");
		mexPrintf("\t '[-I<xinc>[m|c][/<yinc>[m|c]]]', '[-R<west/east/south/north>]', '[-T]', '[-V]')\n");
		
		mexPrintf("\tDistance flag determines how grid (x,y) maps into distance units of filter width as follows:\n");
		mexPrintf("\t   -D0 grid x,y same units as <filter_width>, cartesian Distances.\n");
		mexPrintf("\t   -D1 grid x,y in degrees, <filter_width> in km, cartesian Distances.\n");
		mexPrintf("\t   -D2 grid x,y in degrees, <filter_width> in km, x_scaled by cos(middle y), cartesian Distances.\n");
		mexPrintf("\t   These first three options are faster; they allow weight matrix to be computed only once.\n");
		mexPrintf("\t   Next two options are slower; weights must be recomputed for each scan line.\n");
		mexPrintf("\t   -D3 grid x,y in degrees, <filter_width> in km, x_scale varies as cos(y), cartesian Distances.\n");
		mexPrintf("\t   -D4 grid x,y in degrees, <filter_width> in km, spherical Distances.\n");
		mexPrintf("\t-F sets the filter type and full diameter (6 sigma) filter-width.  Choose between\n");
		mexPrintf("\t   convolution-type filters which differ in how weights are assigned and geospatial\n");
		mexPrintf("\t   filters that seek to return a representative value.\n");
		mexPrintf("\t   Convolution filters:\n");
		mexPrintf("\t     b: Boxcar : a simple averaging of all points inside filter radius.\n");
		mexPrintf("\t     c: Cosine arch : a weighted averaging with cosine arc weights\n");
		mexPrintf("\t     g: Gaussian : weighted averaging with Gaussian weights.\n");
		mexPrintf("\t   Geospatial filters:\n");
		mexPrintf("\t     l: Lower : return minimum of all points.\n");
		mexPrintf("\t     L: Lower+ : return minimum of all +ve points.\n");
		mexPrintf("\t     m: Median : return the median value of all points.\n");
		mexPrintf("\t     p: Maximum likelihood probability estimator : return mode of all points.\n");
		mexPrintf("\t        By default, we return the average if more than one mode is found.\n");
		mexPrintf("\t        Append - or + to the width to instead return the smallest or largest mode.\n");
		mexPrintf("\t     u: Upper : return maximum of all points.\n");
		mexPrintf("\t     U: Upper- : return maximum of all -ve points.\n");
		mexPrintf("\n\tOPTIONS:\n");
		mexPrintf("\t-I for new Increment of output grid; enter xinc, optionally xinc/yinc.\n");
		mexPrintf("\t   Default is yinc = xinc.  Append an m [or c] to xinc or yinc to indicate minutes [or seconds];\n");
		mexPrintf("\t   The new xinc and yinc should be divisible by the old ones (new lattice is subset of old).\n");
		mexPrintf("\t-T Toggles between grid and pixel registration for output grid [Default is same as input registration]\n");
		mexPrintf("\t-R for new Range of output grid; enter <WESN> (xmin, xmax, ymin, ymax) separated by slashes.\n");
		mexErrMsgTxt("\n");
	}

	if (distance_flag < 0 || distance_flag > 4) {
		mexPrintf ("%s: GMT SYNTAX ERROR -D option:  Choose from the range 0-4\n", GMT_program);
		error++;
	}
	if (filter_type < 0 || filter_width <= 0.0) {
		mexPrintf ("%s: GMT SYNTAX ERROR -F option:  Correct syntax:\n", GMT_program);
		mexPrintf ("\t-FX<width>, with X one of bcgmp, width is filter fullwidth\n");
		error++;
	}
	if (toggle && one_or_zero == 0) {	/* Both -N and -T set, not good */
		mexPrintf ("%s: GMT SYNTAX ERROR -T option:  Not allowed with obsolete -N option\n", GMT_program);
		error++;
	}
	if (error) mexErrMsgTxt("\n");

	if (nlhs == 0)
		mexErrMsgTxt("GRDSAMPLE ERROR: Must provide an output.\n");

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
		mexPrintf("GRDFILTER ERROR: Unknown input data type.\n");
		mexErrMsgTxt("Valid types are:double, single, In32, In16, UInt16 and Uint8.\n");
	}

	nx = mxGetN (prhs[0]);
	ny = mxGetM (prhs[0]);
	if (!mxIsNumeric(prhs[0]) || ny < 2 || nx < 2)
		mexErrMsgTxt("First argument must contain a decent array\n");

	nc_h = mxGetN (prhs[1]);
	nr_h = mxGetM (prhs[1]);
	if (!mxIsNumeric(prhs[1]) || nr_h > 1 || nc_h < 9)
		mexErrMsgTxt("Second argument must contain a valid header of the input array\n");

	head = mxGetPr(prhs[1]);		/* Get header info */
	h.x_min = head[0];	h.x_max = head[1];
	h.y_min = head[2];	h.y_max = head[3];
	h.z_min = head[4];	h.z_max = head[5];
	h.x_inc = head[7];	h.y_inc = head[8];
	h.nx = nx;		h.ny = ny;
	h.node_offset = irint(head[6]);

	if (project_info.region_supplied) new_range = TRUE;

	GMT_grd_init (&h, argc, argv, TRUE);	/* Update command history only */

	if (toggle) {	/* Make output grid of the opposite registration */
		one_or_zero = (h.node_offset) ? 1 : 0;
	}
	else
		one_or_zero = (h.node_offset) ? 0 : 1;
	
	input = (float *)mxMalloc ((size_t)(h.nx * h.ny) * sizeof (float));

	/* Transpose from Matlab orientation to gmt grd orientation */
	if (is_double) {
		for (i = 0, i2 = ny - 1; i < ny; i++, i2--) 
			for (j = 0; j < nx; j++) input[i2*nx + j] = (float)z_8[j*ny+i];
	}
	else if (is_single) {
		for (i = 0, i2 = ny - 1; i < ny; i++, i2--) 
			for (j = 0; j < nx; j++) input[i2*nx + j] = z_4[j*ny+i];
	}
	else if (is_int32) {
		for (i = 0, i2 = ny - 1; i < ny; i++, i2--) 
			for (j = 0; j < nx; j++) input[i2*nx + j] = (float)i_4[j*ny+i];
	}
	else if (is_int16) {
		for (i = 0, i2 = ny - 1; i < ny; i++, i2--) 
			for (j = 0; j < nx; j++) input[i2*nx + j] = (float)i_2[j*ny+i];
	}
	else if (is_uint16) {
		for (i = 0, i2 = ny - 1; i < ny; i++, i2--) 
			for (j = 0; j < nx; j++) input[i2*nx + j] = (float)ui_2[j*ny+i];
	}
	else if (is_uint8) {
		for (i = 0, i2 = ny - 1; i < ny; i++, i2--) 
			for (j = 0; j < nx; j++) input[i2*nx + j] = (float)ui_1[j*ny+i];
	}

	last_median = 0.5 * (h.z_min + h.z_max);

	/* Check range of output area and set i,j offsets, etc.  */

	if (!new_range) {
		west_new = h.x_min;
		east_new = h.x_max;
		south_new = h.y_min;
		north_new = h.y_max;
	}
	if (!new_increment) {
		dx_new = h.x_inc;
		dy_new = h.y_inc;
	}

	if (west_new < h.x_min) error = TRUE;
	if (east_new > h.x_max) error = TRUE;
	if (south_new < h.y_min) error = TRUE;
	if (north_new > h.y_max) error = TRUE;
	if (dx_new <= 0.0) error = TRUE;
	if (dy_new <= 0.0) error = TRUE;

	if (error) {
		mexPrintf("%s: New WESN incompatible with old.\n", GMT_program);
		mexErrMsgTxt("\n");
	}

	/* Make sure output grid is kosher */

	test_h.x_min = west_new;	test_h.x_max = east_new;	test_h.x_inc = dx_new;
	test_h.y_min = south_new;	test_h.y_max = north_new;	test_h.y_inc = dy_new;
	/*GMT_grd_RI_verify (&grd_b, 1);*/	/* IF (IVAN == TRUE)  ==> Matlab = BOOM */

	/* We can save time by computing a weight matrix once [or once pr scanline] only
	   if new grid spacing is multiple of old spacing */

	fast_way = (fabs (fmod (dx_new / h.x_inc, 1.0)) < GMT_SMALL && fabs (fmod (dy_new / h.y_inc, 1.0)) < GMT_SMALL);

	if (!fast_way && gmtdefs.verbose) {
		mexPrintf ("%s: Warning - Your output grid spacing is such that filter-weights must\n", GMT_program);
		mexPrintf ("be recomputed for every output node, so expect this run to be slow.  Calculations\n");
		mexPrintf ("can be speeded up significantly if output grid spacing is chosen to be a multiple\n");
		mexPrintf ("of the input grid spacing.  If the odd output grid is necessary, consider using\n");
		mexPrintf ("a \'fast\' grid for filtering and then resample onto your desired grid with grdsample.\n");
	}

	nx_out = one_or_zero + irint ( (east_new - west_new) / dx_new);
	ny_out = one_or_zero + irint ( (north_new - south_new) / dy_new);

	output = (float *)mxMalloc ((size_t)(nx_out*ny_out) * sizeof (float));
	i_origin = (int *)mxMalloc ((size_t)nx_out * sizeof (float));
	if (!fast_way) x_shift = (double *) mxMalloc ((size_t)nx_out * sizeof(double));

	xincnew2 = (one_or_zero) ? 0.0 : 0.5 * dx_new;
	yincnew2 = (one_or_zero) ? 0.0 : 0.5 * dy_new;
	/* xincold2 = (h.node_offset) ? 0.0 : 0.5 * h.x_inc;
	yincold2 = (h.node_offset) ? 0.0 : 0.5 * h.y_inc; */
	xincold2 = (h.node_offset) ? 0.5 * h.x_inc : 0.0;
	yincold2 = (h.node_offset) ? 0.5 * h.y_inc : 0.0;
	offset = (h.node_offset) ? 0.0 : 0.5;
	
	if (fast_way && h.node_offset == one_or_zero) {	/* multiple grid but one is pix, other is grid */
		x_fix = 0.5 * h.x_inc;
		y_fix = 0.5 * h.y_inc;
		shift = (x_fix != 0.0 || y_fix != 0.0);
	}
	
	/* Set up weight matrix and i,j range to test  */

	x_scale = y_scale = 1.0;
	if (distance_flag > 0) {
		x_scale *= DEG2KM;
		y_scale *= DEG2KM;
	}
	if (distance_flag == 2)
		x_scale *= cos (D2R * 0.5 * (north_new + south_new) );
	else if (distance_flag > 2) {
		if ( fabs(south_new) > north_new ) {
			x_scale *= cos (D2R * south_new);
			pole_trouble = (fabs (south_new + 90.0) < GMT_CONV_LIMIT);
		}
		else {
			x_scale *= cos (D2R * north_new);
			pole_trouble = (fabs (north_new - 90.0) < GMT_CONV_LIMIT);
		}
	}
	x_width = filter_width / (h.x_inc * x_scale);
	y_width = filter_width / (h.y_inc * y_scale);
	y_half_width = (int) (ceil(y_width) / 2.0);
	x_half_width = (int) (ceil(x_width) / 2.0);
	nx_fil = 2 * x_half_width + 1;
	ny_fil = 2 * y_half_width + 1;
	if (pole_trouble || nx_fil > h.nx) {	/* Safety valve when x_scale -> 0.0 */
		nx_fil = h.nx + 1;
		x_half_width = h.nx / 2;
	}
	if (ny_fil > h.ny) {
		ny_fil = h.ny;
		y_half_width = h.ny / 2;
	}
	weight = (double *)mxMalloc ((size_t)(nx_fil*ny_fil) * sizeof (double));

	slow = (filter_type >= 3);	/* Will require sorting or comparisons */
	
	if (slow) work_array = (double *) mxMalloc ((size_t)(nx_fil*ny_fil) * sizeof(double));

	/* Compute nearest xoutput i-indices and shifts once */
	
	for (i_out = 0; i_out < nx_out; i_out++) {
		x_out = west_new + i_out * dx_new + xincnew2;
		i_origin[i_out] = (int)floor(((x_out - h.x_min) / h.x_inc) + offset);
		if (!fast_way) x_shift[i_out] = x_out - (h.x_min + i_origin[i_out] * h.x_inc + xincold2);
	}
	
	/* Now we can do the filtering  */

	/* Determine how much effort to compute weights:
		1 = Compute weights once for entire grid
		2 = Compute weights once per scanline
		3 = Compute weights for every output point [slow]
	*/
	
	if (fast_way && distance_flag <= 2)
		effort_level = 1;
	else if (fast_way && distance_flag > 2)
		effort_level = 2;
	else
		effort_level = 3;
		
	if (effort_level == 1) {	/* Only need this once */
		y_out = north_new - yincnew2;
		set_weight_matrix (nx_fil, ny_fil, y_out, north_new, south_new, h.x_inc, h.y_inc, filter_width, filter_type, distance_flag, x_fix, y_fix, shift);
	}
	
	for (j_out = 0; j_out < ny_out; j_out++) {
	
		y_out = north_new - j_out * dy_new - yincnew2;
		j_origin = (int)floor(((h.y_max - y_out) / h.y_inc) + offset);
		if (effort_level == 2)
			set_weight_matrix (nx_fil, ny_fil, y_out, north_new, south_new, h.x_inc, h.y_inc, filter_width, filter_type, distance_flag, x_fix, y_fix, shift);
		if (!fast_way) y_shift = y_out - (h.y_max - j_origin * h.y_inc - yincold2);
		
		for (i_out = 0; i_out < nx_out; i_out++) {
		
			if (effort_level == 3)
				set_weight_matrix (nx_fil, ny_fil, y_out, north_new, south_new, h.x_inc, h.y_inc, filter_width, filter_type, distance_flag, x_shift[i_out], y_shift, fast_way);
			wt_sum = value = 0.0;
			n_in_median = 0;
			ij_out = j_out * nx_out + i_out;
			
			for (ii = -x_half_width; ii <= x_half_width; ii++) {
				i_in = i_origin[i_out] + ii;
				if ( (i_in < 0) || (i_in >= h.nx) ) continue;

				for (jj = -y_half_width; jj <= y_half_width; jj++) {
					j_in = j_origin + jj;
					if ( (j_in < 0) || (j_in >= h.ny) ) continue;
										
					ij_wt = (jj + y_half_width) * nx_fil + ii + x_half_width;
					if (weight[ij_wt] < 0.0) continue;
					
					ij_in = j_in*h.nx + i_in;
					if (GMT_is_fnan (input[ij_in])) continue;

					/* Get here when point is usable  */
					if (slow) {
						work_array[n_in_median] = input[ij_in];
						n_in_median++;
					}
					else {
						value += input[ij_in] * weight[ij_wt];
						wt_sum += weight[ij_wt];
					}
				}
			}
			
			/* Now we have done the convolution and we can get the value  */
			
			if (slow) {
				if (n_in_median) {
					switch (filter_type) {
						case 3:	/* Median */
							GMT_median (work_array, n_in_median, h.z_min, h.z_max, last_median, &this_median);
							last_median = this_median;
							break;
						case 4:	/* Mode */
							GMT_mode (work_array, n_in_median, n_in_median/2, TRUE, GMT_mode_selection, &GMT_n_multiples, &this_median);
							break;
						case 5:	/* Lowest of all */
							this_median = GMT_extreme (work_array, n_in_median, DBL_MAX, 0, -1);
							break;
						case 6:	/* Lowest of positive values */
							this_median = GMT_extreme (work_array, n_in_median, 0.0, +1, -1);
							break;
						case 7:	/* Upper of all values */
							this_median = GMT_extreme (work_array, n_in_median, -DBL_MAX, 0, +1);
							break;
						case 8:	/* Upper of negative values */
							this_median = GMT_extreme (work_array, n_in_median, 0.0, -1, +1);
							break;
					}
					output[ij_out] = (float)this_median;
				}
				else {
					output[ij_out] = GMT_f_NaN;
					n_nan++;
				}
			}
			else {
				if (wt_sum == 0.0) {	/* Assign value = GMT_f_NaN */
					n_nan++;
					output[ij_out] = GMT_f_NaN;
				}
				else
					output[ij_out] = (float)(value / wt_sum);
			}
		}
	}
	
	/* At last, that's it!  Output: */

	h.nx = nx_out;
	h.ny = ny_out;
	h.x_min = west_new;
	h.x_max = east_new;
	h.y_min = south_new;
	h.y_max = north_new;
	h.x_inc = dx_new;
	h.y_inc = dy_new;
	h.node_offset = !one_or_zero;
	
	if (nlhs == 2) {
		head_o[0] = h.x_min;		head_o[1] = h.x_max;
		head_o[2] = h.y_min;		head_o[3] = h.y_max;
		head_o[4] = 0;			head_o[5] = 0;
		head_o[6] = h.node_offset;
		head_o[7] = h.x_inc;		head_o[8] = h.y_inc;
		plhs[1] = mxCreateDoubleMatrix (1,9, mxREAL);
		pdata = mxGetPr(plhs[1]);
		memcpy(pdata, head_o, 8*9);
	}

	mxFree((void *) input);
	mxFree((void *) i_origin);
	mxFree ((void *) weight);
	if (slow) mxFree((void *) work_array);
	if (!fast_way) mxFree(x_shift);
	
	GMT_end (argc, argv);

	/* Transpose from gmt grd orientation to Matlab orientation */
	/* Because we need to do the transposition and also a type conversion, we need a extra array */

	nx = h.nx;		ny = h.ny;
	if (is_double) {
		plhs[0] = mxCreateDoubleMatrix (ny,nx, mxREAL);
		pdata_d = mxGetPr(plhs[0]);
		for (i = 0; i < ny; i++) {
			k1 = ny - i - 1;	k2 = i * nx;
			for (j = 0; j < nx; j++) pdata_d[j*ny+k1] = (double)output[k2+j];
		}
	}
	else if (is_single) {
		plhs[0] = mxCreateNumericMatrix (ny,nx,mxSINGLE_CLASS,mxREAL);
		pdata_s = (float *)mxGetData(plhs[0]);
		for (i = 0; i < ny; i++) {
			k1 = ny - i - 1;	k2 = i * nx;
			for (j = 0; j < nx; j++) pdata_s[j*ny+k1] = output[k2+j];
		}
	}
	else if (is_int32) {
		plhs[0] = mxCreateNumericMatrix (ny,nx,mxINT32_CLASS,mxREAL);
		pdata_i4 = (int *)mxGetData(plhs[0]);
		for (i = 0; i < ny; i++) {
			k1 = ny - i - 1;	k2 = i * nx;
			for (j = 0; j < nx; j++) pdata_i4[j*ny+k1] = irint(output[k2+j]);
		}
	}
	else if (is_int16) {
		plhs[0] = mxCreateNumericMatrix (ny,nx,mxINT16_CLASS,mxREAL);
		pdata_i2 = (short int *)mxGetData(plhs[0]);
		for (i = 0; i < ny; i++) {
			k1 = ny - i - 1;	k2 = i * nx;
			for (j = 0; j < nx; j++) pdata_i2[j*ny+k1] = (short int)irint(output[k2+j]);
		}
	}
	else if (is_uint16) {
		plhs[0] = mxCreateNumericMatrix (ny,nx,mxUINT16_CLASS,mxREAL);
		pdata_ui2 = (unsigned short int *)mxGetData(plhs[0]);
		for (i = 0; i < ny; i++) {
			k1 = ny - i - 1;	k2 = i * nx;
			for (j = 0; j < nx; j++) pdata_ui2[j*ny+k1] = (unsigned short int)irint(output[k2+j]);
		}
	}
	else if (is_uint8) {
		plhs[0] = mxCreateNumericMatrix (ny,nx,mxUINT8_CLASS ,mxREAL);
		pdata_ui1 = (unsigned char *)mxGetData(plhs[0]);
		for (i = 0; i < ny; i++) {
			k1 = ny - i - 1;	k2 = i * nx;
			for (j = 0; j < nx; j++) pdata_ui1[j*ny+k1] = (unsigned char)irint(output[k2+j]);
		}
	}

	mxFree((void *) output);
}

void	set_weight_matrix (int nx_f, int ny_f, double y_0, double north, double south, double dx, double dy, double f_wid, int f_flag, int d_flag, double x_off, double y_off, int fast) {

	/* Last two gives offset between output node and 'origin' input node for this window (0,0 for integral grids) */
	/* TRUE when input/output grids are offset by integer values in dx/dy */
            
	int	i, j, ij, i_half, j_half;
	double	x_scl, y_scl, f_half, r_f_half, sigma, sig_2;
	double	y1, y2, theta, x, y, r, s_y1, c_y1, s_y2, c_y2;

	/* Set Scales  */

	y_scl = (d_flag) ? DEG2KM : 1.0;
	if (d_flag < 2)
		x_scl = y_scl;
	else if (d_flag == 2)
		x_scl = DEG2KM * cos (0.5 * D2R * (north + south));
	else
		x_scl = DEG2KM * cos (D2R * y_0);

	/* Get radius, weight, etc.  */

	i_half = nx_f / 2;
	j_half = ny_f / 2;
	f_half = 0.5 * f_wid;
	r_f_half = 1.0 / f_half;
	sigma = f_wid / 6.0;
	sig_2 = -0.5 / (sigma * sigma);
	for (i = -i_half; i <= i_half; i++) {
		for (j = -j_half; j <= j_half; j++) {
			ij = (j + j_half) * nx_f + i + i_half;
			if (fast && i == 0)
				r = (j == 0) ? 0.0 : j * y_scl * dy;
			else if (fast && j == 0)
				r = i * x_scl * dx;
			else if (d_flag < 4) {
				x = x_scl * (i * dx - x_off);
				y = y_scl * (j * dy - y_off);
				r = hypot (x, y);
			}
			else {
				theta = D2R * (i * dx - x_off);
				y1 = D2R * (90.0 - y_0);
				y2 = D2R * (90.0 - (y_0 + (j * dy - y_off)) );
				sincos (y1, &s_y1, &c_y1);
				sincos (y2, &s_y2, &c_y2);
				r = d_acos (c_y1 * c_y2 + s_y1 * s_y2 * cos(theta)) * DEG2KM * R2D;
			}
			/* Now we know r in f_wid units  */
			
			if (r > f_half) {
				weight[ij] = -1.0;
				continue;
			}
			else if (f_flag >= 3) {
				weight[ij] = 1.0;
				continue;
			}
			else {
				if (f_flag == 0)
					weight[ij] = 1.0;
				else if (f_flag == 1)
					weight[ij] = 1.0 + cos (M_PI * r * r_f_half);
				else
					weight[ij] = exp (r * r * sig_2);
			}
		}
	}
}

void DEBUGA(int n) {
#if debuga
	mexPrintf("Merda %d\n",n);
#endif
}
