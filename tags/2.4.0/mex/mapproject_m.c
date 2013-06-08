/*--------------------------------------------------------------------
 *	$Id: mapproject.c,v 1.41 2004/09/20 20:48:49 pwessel Exp $
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
 * mapproject reads a pair of coordinates [+ optional data fields] from
 * standard input or file(s) and transforms the coordinates according to the
 * map projection selected. See man gmt-system for projections currently supported.
 *
 * The default is to expect longitude, latitude, [and optional datavalues],
 * and return x, y, [ and optional datavalues], but if the -I option is used,
 * the reverse is true.  Specifying -C means that the origin of projected coordinates
 * should be set to origin of projection  [Default origin is lower left corner of "map"].
 * If your data is lat lon instead of lon lat [Default], use -: to toggle x/y -> y/x.
 * Note that only unprojected values are affected by the -: switch.  True x,y values are
 * always printed out as x,y.  Option -G allows calculation of distances along track or
 * to a fixed point, while -L calculates shortest distances to a line.
 *
 *
 * Author:	Paul Wessel
 * Date:	1-MAR-1990
 * Version:	4
 *
 */
/*--------------------------------------------------------------------
 * Mexified version of mapproject
 * Author:	J. Luis
 * Date: 	18-Sep-2004
 *
 * Usage
 * Zout = mapproject_m([lon lat ...],'params');
 *
 * where	[lon lat] is the mxn table with m lines and x,y,... columns
 *		and params are any of the parameters used in the gmt mapproject (e.g.: '-J...','-R...')
 *
 * NOTE. Tested only with coordinate conversions. The other options probably wont work (and crash Matlab)
 *
 *	 
 *		04/06/06 J Luis, Updated to compile with version 4.1.3
 *		14/10/06 J Luis, Now includes the memory leak solving solution
 *		18/02/08 J Luis, T_from & T_to not global (otherwise they remember previous values)
 */ 

#include "gmt.h"
#include "mex.h"

void lon_range_adjust (int range, double *lon);

/* int GMTisLoaded = FALSE;	/* Used to know wether GMT stuff is already in memory or not */

/* --------------------------------------------------------------------------- */
/* Matlab Gateway routine */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	GMT_LONG	i, j, k, n = 0, n_files = 0, unit = 0, n_slash, g_report = 0;
	GMT_LONG	n_fields, distance = 0, proj_type = 0, save[2], two, line_mode = 2;
	
	int error = FALSE, inverse = FALSE, suppress = FALSE, one_to_one = FALSE, ECEF_conv = FALSE;
	int map_center = FALSE, nofile = TRUE, done = FALSE, first = TRUE, datum_conv = FALSE;
	int back_az = FALSE, d_set = FALSE, line_start = TRUE, do_az = FALSE;
	int geodetic_calc = FALSE, do_geocentric, do_line_dist = FALSE, greenwich = FALSE;
	int datum_conv_only = FALSE, double_whammy = FALSE, shift_xy = FALSE, T_heights = FALSE;
	
	double west = 0.0, east = 0.0, south = 0.0, north = 0.0, d, s, x0, y0;
	double *in, *out, fwd_scale, inv_scale, xtmp, ytmp;
	double x_in_min, x_in_max, y_in_min, y_in_max, inch_to_unit, unit_to_inch;
	double x_out_min, x_out_max, y_out_min, y_out_max, u_scale, d_scale;
	double xnear, ynear, false_easting = 0.0, false_northing = 0.0;
	
	char line_file[BUFSIZ], unit_name[80], scale_unit_name[80];
	char txt_a[32], txt_b[32], c;
	char from[GMT_LONG_TEXT], to[GMT_LONG_TEXT];

	struct GMT_DATUM datum;	/* Contains a, f, xyz[3] */
	struct GMT_DATUM T_from;
	struct GMT_DATUM T_to;

	int	argc = 0, n_arg_no_char = 0, n_pts, range = -1, n_comp_here = 1;
	char	**argv;
	double	*in_m, *out_m, *pdata;
	mxArray *cell_array_ptr;

	struct GMT_TABLE *xyline;

	PFD distance_func;
	PFI near_a_line;
	PFD azimuth_func;
	
	out = (double *)NULL;
	in = (double *)NULL;

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
	argv[0] = "mapproject";
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

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {
		
				/* Common parameters */
			
				case 'R':
				case 'J':
				case 'f':
#ifdef GMT_MINOR_VERSION
				case 'm':
				case 'M':
#endif
				case '\0':
					error += GMT_parse_common_options (argv[i], &west, &east, &south, &north);
					break;
				
				/* Supplemental parameters */

				case 'A':
					do_az = TRUE;
	 				n = sscanf (&argv[i][2], "%c%[^/]/%s", &c, txt_a, txt_b);
					if (n < 3) {
						mexPrintf("%s: GMT SYNTAX ERROR.  Expected -Ab|B|f|F<lon0>/<lat0>\n", "mapproject");
						error++;
					}
					else {
						do_geocentric = back_az = FALSE;
						switch (c) {
							case 'B':
								do_geocentric = TRUE;
							case 'b':
								back_az = TRUE;
								break;
							case 'F':
								do_geocentric = TRUE;
							case 'f':
								break;
							default:
								mexPrintf("%s: GMT SYNTAX ERROR.  Expected -Ab|B|f|F<lon0>/<lat0>\n", "mapproject");
								error++;
								break;
						}
						error += GMT_verify_expectations (GMT_io.in_col_type[0], GMT_scanf_arg (txt_a, GMT_io.in_col_type[0], &x0), txt_a);
						error += GMT_verify_expectations (GMT_io.in_col_type[1], GMT_scanf_arg (txt_b, GMT_io.in_col_type[1], &y0), txt_b);
					}
					break;
				case 'C':
					map_center = TRUE;
					if (argv[i][2]) {	/* Also gave shifts */
	 					n = sscanf (&argv[i][2], "%lf/%lf", &false_easting, &false_northing);
						if (n != 2) {
							mexPrintf ("%s: GMT SYNTAX ERROR.  Expected -C[<false_easting>/<false_northing>]\n", "mapproject");
							error++;
						}
						shift_xy = TRUE;
					}
					break;
				case 'D':
					d_set = TRUE;
					GMT_set_measure_unit (argv[i][2]);
					break;
				case 'd':
					if (argv[i][2] == '0') range = 0;
					if (argv[i][2] == '1') range = 1;
					break;
				case 'E':
					ECEF_conv = TRUE;
					if (GMT_set_datum (&argv[i][2], &datum) == -1) error++;
					GMT_ECEF_init (&datum);
					break;
				case 'F':
					one_to_one = TRUE;
					unit = GMT_check_scalingopt ('F',argv[i][2], scale_unit_name);
					break;
				case 'G':
					for (n_slash = 0, k = 2; argv[i][k]; k++) if (argv[i][k] == '/') n_slash++;
					if (n_slash == 2 || n_slash == 1) {	/* Got -Glon0/lat0[/units] */
						distance = 1;
	 					n = sscanf (&argv[i][2], "%[^/]/%[^/]/%c", txt_a, txt_b, &c);
						if (n < 2) {
							mexPrintf("%s: GMT SYNTAX ERROR.  Expected -G<lon0>/<lat0>[/e|E|k|K|m|M|n|N|c|C|d|D]\n", "mapproject");
							error++;
						}
						error += GMT_verify_expectations (GMT_io.in_col_type[0], GMT_scanf_arg (txt_a, GMT_io.in_col_type[0], &x0), txt_a);
						error += GMT_verify_expectations (GMT_io.in_col_type[1], GMT_scanf_arg (txt_b, GMT_io.in_col_type[1], &y0), txt_b);
						if (n_slash == 2)
							error += GMT_get_dist_scale (c, &d_scale, &proj_type, &distance_func);
						else
							error += GMT_get_dist_scale ('\0', &d_scale, &proj_type, &distance_func);
					}
					else {				/* Got -G[units] */
						distance = 2;
						error += GMT_get_dist_scale (argv[i][2], &d_scale, &proj_type, &distance_func);
					}
					break;
				case 'I':
					inverse = TRUE;
					break;
				case 'L':
	 				strcpy (line_file, &argv[i][2]);
					k = strlen (line_file) - 1;
					if (line_file[k] == '+') {	/* Flag to get point number instead of coordinates at nearest point on line */
						line_mode = 3;
						line_file[k] = '\0';
						k--;
					}
					k--;
					if (line_file[k] == '/' && strchr ("ekmndcC", line_file[k+1])) {
						error += GMT_get_dist_scale (line_file[k+1], &d_scale, &proj_type, &distance_func);
						line_file[k] = 0;
					}
					else
						error += GMT_get_dist_scale ('\0', &d_scale, &proj_type, &distance_func);
					do_line_dist = TRUE;
					break;
				case 'Q':
					if (argv[i][2] == 'e') g_report += 1;
					if (argv[i][2] == 'd') g_report += 2;
					if (argv[i][2] == '\0') g_report = 3;
					break;
#ifndef GMT_MINOR_VERSION
				case 'M':               /* Multiple line segments */
					GMT_multisegment (&argv[i][2]);
					break;*/
#endif
				case 'S':
					suppress = TRUE;
					break;
				case 'T':
					datum_conv = TRUE;
					k = 2;
					if (argv[i][k] == 'h') {	/* We will process lon, lat, height data */
						k = 3;
						T_heights = TRUE;	/* If FALSE we set height = 0 */
					}

					if (strchr (&argv[i][k], '/')) {	/* Gave from/to */
						sscanf (&argv[i][k], "%[^/]/%s", from, to);
					}
					else {	/* to not given, set to - which means WGS-84 */
						strcpy (to, "-");
						strcpy (from, &argv[i][k]);
					}
					if (GMT_set_datum (to, &T_to) == -1 || GMT_set_datum (from, &T_from) == -1) {
						error++;
						fprintf (stderr, "%s: GMT SYNTAX ERROR -T: Usage -T[h]<from>[/<to>]\n", GMT_program);
					}
					break;
				default:
					error = TRUE;
					mexPrintf("%s: GMT SYNTAX ERROR.\n", "mapproject");
					break;
			}
		}
	}
	
	if (argc == 1 || error) {
		mexPrintf("mapproject - Forward and Inverse map transformations and geodesy\n\n");
		mexPrintf("usage: mapproject <infiles> -J<parameters> -R<west/east/south/north> [-C]\n");
		mexPrintf("\t[-Ab|B|f|F<lon0/lat0>] [-Dc|i|m|p] [d-1|0|1] [-E[<datum>]] [-F[k|m|n|i|c|p]] [-G[<lon0/lat0>/][<unit>]\n");
		mexPrintf("\t[-H[<nrec>]] [-I] [-L<line.xy>[/<unit>] [-M[<flag>]] [-S] [-T[h]<from>[/<to>]\n");
		mexPrintf("\t[-V[l]] [-:] [-bi[s][<n>]] [-bo[s][<n>]] [-f[i|o]<colinfo>]\n\n");
		
		mexPrintf("\tinfiles (in ASCII or binary) has 2 or more columns.  If no file(s) is given, standard input is read.\n");
		mexPrintf("\n\tOPTIONS:\n");
		mexPrintf("\t-A Calculate azimuths from specified point to each input data point with -Af.\n");
		mexPrintf("\t   Use -Ab to calculate backazimuths from data to specified point.\n");
		mexPrintf("\t   Upper case B or F gives azimuths of geodesics using current ellipsoid.\n");
		mexPrintf("\t-C returns x/y relative to projection center [Default is relative to lower left corner]\n");
		mexPrintf ("\t   Optionally append dx/dy to add (or subtract if -I) (i.e., false easting & northing) [0/0]\n");
		mexPrintf ("\t   Units are plot units unless -F is set in which case the unit is meters.\n");
		mexPrintf("\t-D Temporarily reset MEASURE_UNIT to be c (cm), i (inch), m (meter), or p (point)\n");
		mexPrintf("\t   Cannot be used if -F is set.\n");
		mexPrintf("\t-d Specify longitude range by appending 0 (0 <= lon < 360), 1 (-360 < lon <= 0) [default is -180 < lon < +180]\n");
		mexPrintf("\t-E Convert (lon, lat, h) to Earth Centered Earth Fixed (ECEF) coordinates [-I for inverse].\n");
		mexPrintf("\t   Specify <datum> using datum ID (see -Qd or man page) or as <ellipsoid>:<dx,dy,dz>\n");
		mexPrintf("\t   where <ellipsoid> may be ellipsoid ID (see -Qe or man page) or <major,inv_flattening>.\n");
		mexPrintf("\t   If <datum> = - or not given we assume WGS-84.\n");
		mexPrintf("\t-F force projected values to be in actual meters [Default uses the given plot scale]\n");
		mexPrintf("\t   Specify unit by appending k (km), m (miles), n (nautical miles), i (inch), c (cm), or p (points)\n");
		mexPrintf("\t-G Calculate distances to specified point OR culumlative distances along track (if point not given).\n");
		mexPrintf("\t   Specify unit as m(e)ter, (k)m, (m)ile, (n)autical mile, (d)egree, or (c)artesian in user units.\n");
		mexPrintf("\t   Unit C means Cartesian distances after first projecting the input coordinates (-R, -J).\n");
		mexPrintf("\t   Units E, K, M, N, D mean geodesic distance using current ellipsoid [lower case is spherical].\n");
		mexPrintf("\t   Default is meters on spherical earth.\n");
		mexPrintf("\t-I means Inverse, i.e., get lon/lat from x/y input. [Default is lon/lat -> x/y]\n");
		mexPrintf("\t-L Calculate minimum distances to specified line(s) in the file <line.xy>.\n");
		mexPrintf("\t   Specify unit as m(e)ter, (k)m, (m)ile, (n)autical mile, (d)egree, or (c)artesian in user units.\n");
		mexPrintf("\t   Unit C means Cartesian distances after first projecting the input coordinates (-R, -J).\n");
		mexPrintf("\t   Calculations uses spherical approximations.  Default unit is meters.\n");
		mexPrintf("\t-Q list projection parameters and exit.  For subsets [Default is all] use\n");
		mexPrintf("\t   -Qe shows ellipsoid parameters\n");
		mexPrintf("\t   -Qd shows datum parameters\n");
		mexPrintf("\t-S means Suppress points outside region\n");
		mexPrintf("\t-T means coordinate transformation from datum <from> to datum <to>.\n");
		mexPrintf("\t   Prepend h if input data are lon, lat, height [Default sets height = 0].\n");
		mexPrintf("\t   Specify datums using datum ID (see -Qd or man page) or as <ellipsoid>:<dx,dy,dz>\n");
		mexPrintf("\t   where <ellipsoid> may be ellipsoid ID (see -Qe or man page) or <major,inv_flattening>.\n");
		mexPrintf("\t   <from> = - means WGS-84.  If /<to> is not given we assume WGS-84.\n");
		return;
	}
	
	geodetic_calc = (distance || do_az || do_line_dist);
	
	if (datum_conv && (distance + ECEF_conv + do_line_dist) > 0) {	/* No good... */
		mexPrintf("%s: GMT SYNTAX ERROR: -T cannot work with -E, -G or -L\n", "mapproject");
		error++;
	}
	if (geodetic_calc && inverse) {	/* No good... */
		mexPrintf("%s: GMT SYNTAX ERROR: -A, -G, and -L cannot work with -I\n", "mapproject");
		error++;
	}
	if (project_info.projection < 0 && (distance || do_line_dist) && proj_type == 2) {	/* Must have -J */
		mexPrintf("%s: GMT SYNTAX ERROR:  Must specify -J option with selected form of -G or -L\n", "mapproject");
		error++;
	}
	if (!project_info.region_supplied && !(geodetic_calc || datum_conv || ECEF_conv)) {
		mexPrintf("%s: GMT SYNTAX ERROR:  Must specify -R option\n", "mapproject");
		error++;
	}
 	if ((d_set + one_to_one) == 2) {
		mexPrintf("%s: GMT SYNTAX ERROR:  Can specify only one of -D and -F\n", "mapproject");
		error++;
	}
	if (GMT_io.binary[GMT_IN] && GMT_io.ncol[GMT_IN] == 0) GMT_io.ncol[GMT_IN] = 2;
        if (((datum_conv && GMT_datum.h_given) || ECEF_conv) && GMT_io.binary[GMT_IN] && GMT_io.ncol[GMT_IN] < 3) {
                mexPrintf("%s: GMT SYNTAX ERROR.  For -E or -T, binary input data (-bi) must have at least 3 columns\n", "mapproject");
		error++;
	}

	if (error) mexErrMsgTxt ("\n");

	if (nlhs == 0)
		mexErrMsgTxt("MAPPROJECT ERROR: Must provide an output.\n");

	/* Check that first argument contains at least a mx2 table */
	n_pts = mxGetM (prhs[0]);
	n_fields = mxGetN(prhs[0]);
	if (!mxIsNumeric(prhs[0]) || (n_fields < 2))
		mexErrMsgTxt("MAPPROJECT ERROR: first argument must contain a mx2 table with the x,y positions to convert.\n");

	/* Read the input points and convert them to double */
	if (mxIsDouble(prhs[0]))
		in_m = mxGetPr(prhs[0]);
	else if (mxIsSingle(prhs[0]))
		in_m = mxGetData(prhs[0]);

	if (datum_conv) GMT_datum_init (&T_from, &T_to, T_heights);

	if (datum_conv && project_info.projection && project_info.region_supplied) {	/* Do datum shift & project coordinates */
		double_whammy = TRUE;
		if (inverse) {	/* Need to set the ellipsoid to that of the old datum */
			if (GMT_datum.from.ellipsoid_id < 0) {
				gmtdefs.ellipsoid = GMT_N_ELLIPSOIDS - 1;
				gmtdefs.ref_ellipsoid[i].eq_radius = GMT_datum.from.a;
				gmtdefs.ref_ellipsoid[i].flattening = GMT_datum.from.f;
			}
			else
				gmtdefs.ellipsoid = GMT_datum.from.ellipsoid_id;
		}
		else {	/* Need to set the ellipsoid to that of the new datum */
			if (GMT_datum.to.ellipsoid_id < 0) {
				gmtdefs.ellipsoid = GMT_N_ELLIPSOIDS - 1;
				gmtdefs.ref_ellipsoid[i].eq_radius = GMT_datum.to.a;
				gmtdefs.ref_ellipsoid[i].flattening = GMT_datum.to.f;
			}
			else
				gmtdefs.ellipsoid = GMT_datum.to.ellipsoid_id;
		}
	}
	else
		datum_conv_only = datum_conv;

	GMT_init_scales (unit, &fwd_scale, &inv_scale, &inch_to_unit, &unit_to_inch, unit_name);

	if (distance) {	/* save output format in case -J changes it */
		save[0] = GMT_io.out_col_type[0];
		save[1] = GMT_io.out_col_type[1];
	}
	u_scale = (inverse) ? inv_scale : fwd_scale;

	if (project_info.projection < 0) {	/* Supply dummy linear proj */
		GMT_parse_J_option ("x1d");
		if (!project_info.region_supplied) {
			west = 0.0;	east = 360.0;
			south = -90.0;	north = 90.0;
		}
	}

	if (west == east && south == north)
		mexErrMsgTxt ("GMT Fatal Error: No region selected - Aborts!\n");

	GMT_map_setup (west, east, south, north);
	
	azimuth_func = (GMT_IS_SPHERICAL || !do_geocentric) ? GMT_az_backaz_sphere : GMT_az_backaz_geodesic;
	 
	if (distance && proj_type < 2) {	/* Ensure we use the selected output coordinates */
		GMT_io.out_col_type[0] = save[0];
		GMT_io.out_col_type[1] = save[1];
	}

	if (GMT_io.in_col_type[0] & GMT_IS_GEO && proj_type == 0) {	/* Geographic data */
		GMT_distance_func = (PFD) GMT_great_circle_dist;
		near_a_line = (PFI) GMT_near_a_line_spherical;
		greenwich = (west < 0.0 && east > 0.0);
	}
	else {
		GMT_distance_func = (PFD) GMT_cartesian_dist;
		near_a_line = (PFI) GMT_near_a_line_cartesian;
	}
	if (do_line_dist) {
		GMT_import_table ((void *)line_file, GMT_IS_FILE, &xyline, 0.0, greenwich, FALSE, FALSE);
		if (proj_type == 2) {	/* Must convert the line points first */
			for (i = 0; i < xyline->n_segments; i++) {
				for (j = 0; j < xyline->segment[i]->n_rows; j++) {
					GMT_geo_to_xy (xyline->segment[i]->coord[GMT_X][j], xyline->segment[i]->coord[GMT_Y][j], &xtmp, &ytmp);
					xyline->segment[i]->coord[GMT_X][j] = xtmp;
					xyline->segment[i]->coord[GMT_Y][j] = ytmp;
				}
			}
		}
	}

	/* Now we are ready to take on some input values */
	
	if ((GMT_IS_MAPPING || ECEF_conv) && inverse) {
		GMT_io.out_col_type[0] = GMT_IS_LON;	GMT_io.out_col_type[1] = GMT_IS_LAT;	/* Inverse projection expects x,y and gives lon, lat */
		GMT_io.in_col_type[0] = GMT_io.in_col_type[1] = GMT_IS_FLOAT;
	}
	if (datum_conv_only) {	/* Both in and out are geographic */
		GMT_io.in_col_type[0] = GMT_io.out_col_type[0] = GMT_IS_LON;
		GMT_io.in_col_type[1] = GMT_io.out_col_type[1] = GMT_IS_LAT;
		GMT_io.in_col_type[2] = GMT_io.out_col_type[2] = GMT_IS_FLOAT;
	}

	x_in_min = y_in_min = x_out_min = y_out_min = DBL_MAX;
	x_in_max = y_in_max = x_out_max = y_out_max = -DBL_MAX;
	
	two = (ECEF_conv || (datum_conv && GMT_datum.h_given)) ? 3 : 2;	/* # of output points from conversion */
	
	if (shift_xy && one_to_one) {	/* Use same units in -C and -F */
		false_easting *= u_scale;
		false_northing *= u_scale;		
	}
	if (distance == 2 && proj_type == 2) {	/* Must project the fixed point here */
		GMT_geo_to_xy (x0, y0, &xtmp, &ytmp);
		if (map_center) {	/* Change origin from lower left to projection center */
			xtmp -= project_info.x0;
			ytmp -= project_info.y0;
		}
		if (one_to_one) {	/* Convert to 1:1 scale */
			xtmp /= project_info.x_scale;
			ytmp /= project_info.y_scale;
			if (unit) {
				xtmp *= u_scale;
				ytmp *= u_scale;
			}
		}
		else if (gmtdefs.measure_unit != GMT_INCH) {	/* Convert from inch to whatever */
			xtmp *= inch_to_unit;
			ytmp *= inch_to_unit;
		}
		if (shift_xy) {
			xtmp += false_easting;
			ytmp += false_northing;
		}
		x0 = xtmp;
		y0 = ytmp;
	}

	n_comp_here = (do_line_dist) ? 3 : 1;	/* Number of vars computed here and to be output */

	if ((out_m = mxCalloc(n_pts * (n_fields+n_comp_here), sizeof (double))) == 0)
		mexErrMsgTxt("MAPPROJECT ERROR: Could not allocate memory\n");

	for (i = 0; i < n_pts; i++) {
		
		while ( (GMT_is_dnan(in_m[i]) || GMT_is_dnan(in_m[i+n_pts])) ) {
			for (j = 0; j < n_fields; j++) out_m[j*n_pts+i] = in_m[j*n_pts+i];
			i++;
		}

		if (inverse) {		/* Do inverse transformation */

			if (!out) out = mxCalloc (n_fields, sizeof (double));
			if (!in) in = mxCalloc (n_fields, sizeof (double));
			for (k = 0; k < n_fields; k++)		/* copy the ith line to a temporary var*/
				in[k] = in_m[i+k*n_pts];
			if (shift_xy) {
				in[0] -= false_easting;
				in[1] -= false_northing;
			}
			if (ECEF_conv) {
				GMT_ECEF_inverse (in, out);
			}
			else {
				if (one_to_one) {	/* Convert from 1:1 scale */
					if (unit) {
						in[0] *= u_scale;
						in[1] *= u_scale;
					}
					in[0] *= project_info.x_scale;
					in[1] *= project_info.y_scale;
				}
				else if (gmtdefs.measure_unit != GMT_INCH) {	/* Convert from whatever to inch */
					in[0] *= unit_to_inch;
					in[1] *= unit_to_inch;
				}
				if (map_center) {	/* Then correct so lower left corner is (0,0) */
					in[0] += project_info.x0;
					in[1] += project_info.y0;
				}
				GMT_xy_to_geo (&out[0], &out[1], in[0], in[1]);
			}
			if (double_whammy) {	/* Now apply datum shift */
				in[0] = out[0];
				in[1] = out[1];
				GMT_conv_datum (in, out);
			}

			if (suppress && GMT_map_outside (out[0], out[1])) continue;

			/* Simply copy other columns and output */
			lon_range_adjust (range, &out[0]);
			for (k = 0; k < two; k++) out_m[k*n_pts+i] = out[k];		/* first copy the projected columns */
			for (k = two; k < n_fields; k++) out_m[k*n_pts+i] = in_m[k*n_pts+i];	/* now copy whatever follows it */
		}
		else {		/* Do forward transformation */

			if (geodetic_calc) { 
				while ( (GMT_is_dnan(in_m[i]) || GMT_is_dnan(in_m[i+n_pts])) ) {
					for (j = 0; j < n_fields; j++) out_m[j*n_pts+i] = in_m[j*n_pts+i];
					i++;
					line_start = TRUE;
				}
			}

			if (!out) out = mxCalloc (n_fields+n_comp_here, sizeof (double));
			if (!in) in = mxCalloc (n_fields, sizeof (double));
			for (k = 0; k < n_fields; k++)		/* copy the ith line to a temporary var*/
				in[k] = in_m[i+k*n_pts];

			if (suppress && GMT_map_outside (in[0], in[1])) continue;
				
			if (datum_conv_only) {
				GMT_conv_datum (in, out);
			}
			else if (ECEF_conv) {
				GMT_ECEF_forward (in, out);
			}
			else if (do_line_dist && proj_type != 2) {	/* Do nothing below as we are not using "out" */
			}
			else {
				if (double_whammy) {	/* Apply datum shift first */
					GMT_conv_datum (in, out);
					in[0] = out[0];
					in[1] = out[1];
				}
				GMT_geo_to_xy (in[0], in[1], &out[0], &out[1]);
				if (map_center) {	/* Change origin from lower left to projection center */
					out[0] -= project_info.x0;
					out[1] -= project_info.y0;
				}
				if (one_to_one) {	/* Convert to 1:1 scale */
					out[0] /= project_info.x_scale;
					out[1] /= project_info.y_scale;
					if (unit) {
						out[0] *= u_scale;
						out[1] *= u_scale;
					}
				}
				else if (gmtdefs.measure_unit != GMT_INCH) {	/* Convert from inch to whatever */
					out[0] *= inch_to_unit;
					out[1] *= inch_to_unit;
				}
				if (shift_xy) {
					out[0] += false_easting;
					out[1] += false_northing;
				}
			}

			if (geodetic_calc) {	/* Get either distances or azimuths */
				if (distance) {	/* Cumulative distances along track */
					if (distance == 2 && line_start)
						s = d = 0.0;
					else if (proj_type == 2)	/* Calculate Cartesian distances using projected units */
						s = hypot (x0 - out[0], y0 - out[1]);
					else if (proj_type == 1)	/* Plain Cartesian distances using input points */
						s = hypot (x0 - in[0], y0 - in[1]);
					else				/* Great circle distances */
						s = d_scale * (*distance_func) (x0, y0, in[0], in[1]);
					if (distance == 2) {
						line_start = FALSE;
						d += s;
						if (proj_type == 2) {	/* Calculate distances using projected units */
							x0 = out[0];
							y0 = out[1];
						}
						else {
							x0 = in[0];
							y0 = in[1];
						}
					}
					else
						d = s;
				}
				else if (do_line_dist) {	/* Compute closest distance to line */
					if (proj_type == 2)	/* Using projected coordinates */
						(void *) near_a_line (out[0], out[1], xyline, line_mode, &d, &xnear, &ynear);
					else			/* Using input coordinates */
						(void *) near_a_line (in[0], in[1], xyline, line_mode, &d, &xnear, &ynear);
					d *= d_scale;
				}
				else {	/* Azimuths */
					d = (*azimuth_func) (x0, y0, in[0], in[1], back_az);
				}

				/* Simply copy other columns and output */
				for (k = 0; k < n_fields; k++) out_m[k*n_pts+i] = in[k];
				out_m[n_fields*n_pts+i] = d;
				if (do_line_dist) { 
					out_m[(n_fields+1)*n_pts+i] = xnear;
					out_m[(n_fields+2)*n_pts+i] = ynear;
				}
			}
			else {
				/* Simply copy other columns and output */
				for (k = 0; k < two; k++) out_m[k*n_pts+i] = out[k];		/* first copy the projected columns */
				for (k = two; k < n_fields; k++) out_m[k*n_pts+i] = in_m[k*n_pts+i];	/* now copy whatever follows it */
			}
		}
	}
	
	mxFree(in);
	mxFree(out);
	GMT_end (argc, argv);

	if (geodetic_calc)
		k = n_fields + n_comp_here;
	else
		k = n_fields;
	plhs[0] = mxCreateDoubleMatrix (n_pts,k, mxREAL);
	pdata = mxGetPr(plhs[0]);
	memcpy(pdata, out_m, n_pts*k*8);
	mxFree(out_m);
	mxFree(argv);
}

void lon_range_adjust (int range, double *lon) {
	switch (range) {	/* Adjust to the desired range */
		case 0:		/* Make 0 <= lon < 360 */
			while ((*lon) < 0.0) (*lon) += 360.0;
			while ((*lon) >= 360.0) (*lon) -= 360.0;
			break;
		case 1:		/* Make -360 < lon <= 0 */
			while ((*lon) <= -360.0) (*lon) += 360.0;
			while ((*lon) > 0) (*lon) -= 360.0;
			break;
		default:	/* Make -180 < lon < +180 */
			while ((*lon) < -180.0) (*lon) += 360.0;
			while ((*lon) > 180.0) (*lon) -= 360.0;
			break;
	}
}
