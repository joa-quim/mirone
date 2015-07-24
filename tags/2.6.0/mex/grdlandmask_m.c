/*--------------------------------------------------------------------
 *	$Id: grdlandmask.c,v 1.42 2008/01/23 03:22:48 guru Exp $
 *
 *	The Mirone/GMT-system:	@(#)grdlandmask.c
 *
 *	Copyright (c) 1991-2008 by P. Wessel and W. H. F. Smith
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
 * grdlandmask defines a grid based on region and xinc/yinc values,
 * reads a shoreline data base, and sets the grid nodes inside, on the
 * boundary, and outside of the polygons to the user-defined values
 * <in>, <on>, and <out>.  These may be any number, including NaN.
 *
 * Author:	P. Wessel
 * Date:	23-Sep-1994
 * Version:	3.0
 * Modified:	24-JUN-1998, for GMT 3.1
 *		18-AUG-1999, for GMT 3.3.2
 *		13-JUL-2000, for GMT 3.3.5
 * Version:	4
 *
 * Mexified by	Joaquim Luis
 * Date:	26-MAR-2008
 *		18/01/10 	J Luis, Add option to call aguentabar.dll via mexEvalString
 */
 
#include "gmt.h"
#include "mex.h"

#define GRDLANDMASK_N_CLASSES	(GMT_MAX_GSHHS_LEVEL + 1)	/* Number of bands separated by the levels */

struct GRDLANDMASK_CTRL {	/* All control options for this program (except common args) */
	/* ctive is TRUE if the option has been activated */
	struct A {	/* -A<min_area>[/<min_level>/<max_level>] */
		int active;
#ifdef GMT_MINOR_VERSION
		struct GMT_SHORE_SELECT info;
#else
		int low;	/* Lowest hierarchical level to use [0] */
		int high;	/* Highest hierarchical level to use [4] */
		double area;	/* Area of smallest geographical feature to include [0] */
#endif
	} A;
	struct D {	/* -D<resolution> */
		int active;
		char set;	/* One of f, h, i, l, c */
	} D;
	struct e {	/* -e For use in compiled Mirone */
		int active;
	} e;
	struct F {	/* -F */
		int active;
	} F;
	struct G {	/* -G<maskfile> */
		int active;
		char *file;
	} G;
	struct I {	/* -Idx[/dy] */
		int active;
		double xinc, yinc;
	} I;
	struct N {	/* -N<maskvalues>[o] */
		int active;
		int edge;	/* TRUE if edges are considere outside */
		double mask[GRDLANDMASK_N_CLASSES];	/* values for each level */
	} N;
};

/* int GMTisLoaded = FALSE;	/* Used to know wether GMT stuff is already in memory or not */

/* --------------------------------------------------------------------------- */
/* Matlab Gateway routine */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	int	error = FALSE, dry_wet_only = FALSE, greenwich = FALSE;
	int	temp_shift = FALSE, wrap, used_polygons;
	int	out_logic = TRUE, out_uint = FALSE, out_float = FALSE;

	char line[GMT_LONG_TEXT], ptr[BUFSIZ], cmd[24];
	char *shore_resolution[5] = {"full", "high", "intermediate", "low", "crude"};

	GMT_LONG i, j, k, ij, ii, bin, ind, nm, np, side, i_min, i_max, j_min, j_max, nx1, ny1, np_new, pos;
	int      base = 3, direction, is_inside = 1, err;

	unsigned char *data_8;

	float	*data_32;
	double	*x, *y, xmin, xmax, ymin, ymax, west_border, east_border, i_dx_inch, i_dy_inch;
	double	xinc2, yinc2, i_dx, i_dy, edge = 0.0, del_off, dummy, *hdr, *ptr_wb;

	struct GMT_SHORE c;
	struct GRD_HEADER header;
	struct GMT_GSHHS_POL *p;
	struct GRDLANDMASK_CTRL *Ctrl;

	void *New_Grdlandmask_Ctrl (), Free_Grdlandmask_Ctrl (struct GRDLANDMASK_CTRL *C);

	char **argv, res = 'l';
	int argc, dims[] = {0,0}, n;

	mxArray *rhs[3];	/* For the aguentabar */

	argc = nrhs;
	for (i = 0; i < argc; i++) {		/* Check input to be sure it is of type char. */
		if(!mxIsChar(prhs[i]))
			mexErrMsgTxt("Input must be of type char.");
	}
	argc++;			/* to account for the program's name to be inserted in argv[0] */

	/* get the length of the input string */
	argv=(char **)mxCalloc(argc, sizeof(char *));
	argv[0] = "grdlandmask_m";
	for (i = 1; i < argc; i++) {
		argv[i] = (char *)mxArrayToString(prhs[i-1]);
	}

	/*if (!GMTisLoaded) {
		argc = GMT_begin (argc, argv);
		GMTisLoaded = TRUE;
	}
	else
		argc = GMT_short_begin (argc, argv);*/
	argc = GMT_begin (argc, argv);

	/* Check and interpret the command line arguments */


	Ctrl = (struct GRDLANDMASK_CTRL *) New_Grdlandmask_Ctrl ();	/* Allocate and initialize defaults in a new control structure */

	GMT_grd_init (&header, argc, argv, FALSE);

	/* Check command line arguments */

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {

				/* Common parameters */

				case 'R':
				case 'V':
				case ':':
				case '\0':
					error += GMT_parse_common_options (argv[i], &header.x_min, &header.x_max, &header.y_min, &header.y_max);
					break;

				/* Supplemental parameters */

				case 'A':
					Ctrl->A.active = TRUE;
#ifdef GMT_MINOR_VERSION
					Ctrl->A.info.fraction = Ctrl->A.info.flag = Ctrl->A.info.low = 0;
					Ctrl->A.info.high = GMT_MAX_GSHHS_LEVEL;
					Ctrl->A.info.area = 0.;
					GMT_set_levels (&argv[i][2], &Ctrl->A.info);
#else
					j = sscanf (&argv[i][2], "%lf/%d/%d", &Ctrl->A.area, &Ctrl->A.low, &Ctrl->A.high);
					if (j == 1) Ctrl->A.low = 0, Ctrl->A.high = GMT_MAX_GSHHS_LEVEL;
#endif
					break;
				case 'D':
					Ctrl->D.active = TRUE;
					Ctrl->D.set = argv[i][2];
					break;
				case 'e':
					Ctrl->e.active = TRUE;
					break;
				case 'N':
					Ctrl->N.active = TRUE;
					strcpy (line, &argv[i][2]);
					if (line[strlen(line)-1] == 'o') { /* Edge is considered outside */
						Ctrl->N.edge = TRUE;
						line[strlen(line)-1] = 0;
					}
					j = pos = 0;
					while (j < 5 && (GMT_strtok (line, "/", &pos, ptr))) {
						Ctrl->N.mask[j] = (ptr[0] == 'N' || ptr[0] == 'n') ? GMT_f_NaN : (float)atof (ptr);
						j++;
					}
					if (!(j == 2 || j == 5))
						mexErrMsgTxt ("grdlandmask_m: SYNTAX ERROR -N option:  Specify 2 or 5 arguments\n");

					dry_wet_only = (j == 2);
					for (k = 0; k < j; k++) {	/* Pick up what the output data type will be */
						if (Ctrl->N.mask[k] != 0 && Ctrl->N.mask[k] != 1) {
							if (Ctrl->N.mask[k] < 0 || Ctrl->N.mask[k] > 255) {
								out_float = TRUE;
								out_logic = out_uint = FALSE;
								break;
							}
							else {
								out_uint = TRUE;
								out_logic = FALSE;
							}
						}
						if (mxIsNaN(Ctrl->N.mask[k])) {
							out_float = TRUE;
							out_logic = out_uint = FALSE;
							break;
						}
					}
					break;
				case 'F':
					Ctrl->F.active = TRUE;
					break;
				case 'I':
					if (GMT_getinc (&argv[i][2], &Ctrl->I.xinc, &Ctrl->I.yinc)) {
						GMT_inc_syntax ('I', 1);
						error = TRUE;
					}
					Ctrl->I.active = TRUE;
					break;
				default:
					error = TRUE;
					GMT_default_error (argv[i][1]);
					break;
			}
		}
	}

	if (argc == 1 || GMT_give_synopsis_and_exit) {
		mexPrintf ("grdlandmask_m - Create \"wet-dry\" mask grid file from shoreline data base\n\n");
		mexPrintf ("usage: [MASK,hdr,X,Y] = grdlandmask_m( %s, %s, \n", GMT_I_OPT, GMT_Rgeo_OPT);
		mexPrintf ("\t[-A<min_area>[/<min_level>/<max_level>]], [-D<resolution>], [-e] [-F], [-N<maskvalues>[o]])\n\n");

		if (GMT_give_synopsis_and_exit) return;

		mexPrintf ("\t-I sets the grid spacing (dx, dy) for the new grid\n");
		mexPrintf ("\n\tOPTIONS:\n");
		mexPrintf ("\t-A features smaller than <min_area> (in km^2) or of levels (0-4) outside the min-max levels\n");
		mexPrintf ("\t   will be skipped [0/4] (see pscoast for details)]\n");
		mexPrintf ("\t-D Choose one of the following resolutions:\n");
		mexPrintf ("\t   f - full resolution (may be very slow for large regions)\n");
		mexPrintf ("\t   h - high resolution (may be slow for large regions)\n");
		mexPrintf ("\t   i - intermediate resolution\n");
		mexPrintf ("\t   l - low resolution [Default]\n");
		mexPrintf ("\t   c - crude resolution, for tasks that need crude continent outlines only\n");
		mexPrintf ("\t-F Force pixel registration for output grid [Default is gridline orientation]\n");
		mexPrintf ("\t-N gives values to use if a node is outside or inside a feature.\n");
		mexPrintf ("\t   Append o to let feature boundary be considered outside [Default is inside].\n");
		mexPrintf ("\t   Specify this information using 1 of 2 formats:\n");
		mexPrintf ("\t   -N<wet>/<dry>.\n");
		mexPrintf ("\t   -N<ocean>/<land>/<lake>/<island>/<pond>.\n");
		mexPrintf ("\t   NaN is a valid entry. Default values are 0/1/0/1/0 (i.e., 0/1)\n\n");
		mexPrintf ("\t   By default (-N with 0 & 1) MASK is of type logical but that depends on the -N option\n");
		mexPrintf ("\t   If NaN or values outside the [0-255] range are transmitted, MASK is of type float.\n");
		mexPrintf ("\t   For values in the [0-255] range, MASK is of type uint8.\n");
		mexPrintf ("\t-e To be used from the Mirone stand-alone version.\n");
		return;
	}

	if (!project_info.region_supplied) {
		mexPrintf ("%s: GMT SYNTAX ERROR:  Must specify -R option\n", GMT_program);
		error++;
	}
	if (Ctrl->I.xinc <= 0.0 || Ctrl->I.yinc <= 0.0) {
		mexPrintf ("%s: GMT SYNTAX ERROR -I option.  Must specify positive increment(s)\n", GMT_program);
		error = TRUE;
	}

	if (error) return;

	header.x_inc = Ctrl->I.xinc;
	header.y_inc = Ctrl->I.yinc;
	header.node_offset = Ctrl->F.active;
	
	GMT_RI_prepare (&header);	/* Ensure -R -I consistency and set nx, ny */

	if (header.x_min < 0.0 && header.x_max < 0.0) {	/* Shift longitudes */
		temp_shift = TRUE;
		header.x_min += 360.0;
		header.x_max += 360.0;
	}

	if (GMT_grd_RI_verify (&header, 1))
		mexErrMsgTxt("GRDLANDMASK_M: failed in 'GMT_grd_RI_verify' test"); 

	base = GMT_set_resolution (&Ctrl->D.set, 'D');
	
	is_inside = (Ctrl->N.edge) ? 2 : 1;		/* Whether of not being exactly on an edge is outside */
	if (dry_wet_only) {
		Ctrl->N.mask[3] = Ctrl->N.mask[1];
		Ctrl->N.mask[2] = Ctrl->N.mask[4] = Ctrl->N.mask[0];
	}

#ifdef GMT_MINOR_VERSION
	if (GMT_init_shore(Ctrl->D.set, &c, header.x_min, header.x_max, header.y_min, header.y_max, &Ctrl->A.info))  {
#else
	if (GMT_init_shore(Ctrl->D.set, &c, header.x_min, header.x_max, header.y_min, header.y_max))  {
#endif
		mexPrintf ("%s: %s resolution shoreline data base not installed\n", GMT_program, shore_resolution[base]);
		mexErrMsgTxt ("");
	}

	sprintf (line, "%s\n", gmtdefs.d_format);

	i_dx = 1.0 / header.x_inc;
	i_dy = 1.0 / header.y_inc;
	del_off = (header.node_offset) ? 0.5 : 0.0;
	xinc2 = (header.node_offset) ? 0.5 * header.x_inc : 0.0;
	yinc2 = (header.node_offset) ? 0.5 * header.y_inc : 0.0;
	nm = header.nx * header.ny;

	/* The output data type depends on what was asked via -N option */
	if (out_logic) {
		plhs[0] = mxCreateNumericMatrix (header.ny,header.nx,mxLOGICAL_CLASS,mxREAL);
		data_8 = (unsigned char *)mxGetData(plhs[0]);
	}
	else if (out_uint) {
		plhs[0] = mxCreateNumericMatrix (header.ny,header.nx,mxUINT8_CLASS,mxREAL);
		data_8 = (unsigned char *)mxGetData(plhs[0]);
	}
	else {
		plhs[0] = mxCreateNumericMatrix (header.ny,header.nx,mxSINGLE_CLASS,mxREAL);
		data_32 = (float *)mxGetData(plhs[0]);
	}

	/* All data nodes are thus initialized to 0 */
	x = (double *) GMT_memory (VNULL, (size_t)header.nx, sizeof(double), GMT_program);
	y = (double *) GMT_memory (VNULL, (size_t)header.ny, sizeof(double), GMT_program);

	nx1 = header.nx - 1;	ny1 = header.ny - 1;

	if (header.x_min < 0.0 && header.x_max > 0.0) {	/* Must shift longitudes */
		greenwich = TRUE;
		edge = header.x_min;
	}

	GMT_parse_J_option ("x1d");	/* Fake linear projection */
	GMT_err_fail (GMT_map_setup (header.x_min, header.x_max, header.y_min, header.y_max), "");
	wrap = GMT_360_RANGE (header.x_max, header.x_min);
	
	if (!Ctrl->e.active && gmtdefs.verbose) {
		rhs[0] = mxCreateDoubleScalar(0.0);
		ptr_wb = mxGetPr(rhs[0]);
		rhs[1] = mxCreateString("title");
		rhs[2] = mxCreateString("Creating Mask ...");
		mexCallMATLAB(0,NULL,3,rhs,"aguentabar");
	}
	else if (Ctrl->e.active && gmtdefs.verbose)
		mexEvalString("aguentabar(0,\"title\",\"Creating Mask ...\")");

	/* Fill out gridnode coordinates and apply the implicit linear projection */

	for (i = 0; i < header.nx; i++) 
		GMT_geo_to_xy (GMT_i_to_x (i, header.x_min, header.x_max, header.x_inc, header.xy_off, header.nx), 0.0, &x[i], &dummy);
	for (j = 0; j < header.ny; j++) 
		GMT_geo_to_xy (0.0, GMT_j_to_y (j, header.y_min, header.y_max, header.y_inc, header.xy_off, header.ny), &dummy, &y[j]);
	i_dx_inch = 1.0 / fabs (x[1] - x[0]);
	i_dy_inch = 1.0 / fabs (y[1] - y[0]);

        west_border = floor (project_info.w / c.bsize) * c.bsize;
        east_border = ceil (project_info.e / c.bsize) * c.bsize;
	for (ind = 0; ind < c.nb; ind++) {	/* Loop over necessary bins only */

		bin = c.bins[ind];
		if (!Ctrl->e.active && gmtdefs.verbose) {
			*ptr_wb = (double)(ind+1) / c.nb;
			mexCallMATLAB(0,NULL,1,rhs,"aguentabar");
		}
		else if (Ctrl->e.active && gmtdefs.verbose) {
			sprintf(cmd, "aguentabar(%f)", (double)(ind+1) / c.nb);
			mexEvalString(cmd);
		}

#ifdef GMT_MINOR_VERSION
			if ((err = GMT_get_shore_bin (ind, &c))) {
#else
			if ((err = GMT_get_shore_bin (ind, &c,  Ctrl->A.area, Ctrl->A.low, Ctrl->A.high))) {
#endif
			mexPrintf ("%s: %s [%s resolution shoreline]\n", GMT_program, GMT_strerror(err), shore_resolution[base]);
			mexErrMsgTxt ("");
		}

		/* Use polygons, if any.  Go in both directions to cover both land and sea */

		used_polygons = FALSE;

		for (direction = -1; c.ns > 0 && direction < 2; direction += 2) {

			/* Assemble one or more segments into polygons */
#ifdef GMT_MINOR_VERSION
			np = GMT_assemble_shore (&c, direction, TRUE, greenwich, west_border, east_border, &p);
#else
			np = GMT_assemble_shore (&c, direction, Ctrl->A.low, TRUE, greenwich, west_border, east_border, &p);
#endif

			/* Get clipped polygons in x,y inches that can be processed */

			np_new = GMT_prep_polygons (&p, np, FALSE, 0.0, -1);

			for (k = 0; k < np_new; k++) {

				if (p[k].n == 0) continue;

				used_polygons = TRUE;	/* At least some points made it to here */

				/* Find min/max of polygon in inches */

				xmin = xmax = p[k].lon[0];
				ymin = ymax = p[k].lat[0];
				for (i = 1; i < p[k].n; i++) {
					if (p[k].lon[i] < xmin) xmin = p[k].lon[i];
					if (p[k].lon[i] > xmax) xmax = p[k].lon[i];
					if (p[k].lat[i] < ymin) ymin = p[k].lat[i];
					if (p[k].lat[i] > ymax) ymax = p[k].lat[i];
				}
				i_min = (int)MAX (0, ceil (xmin * i_dx_inch - del_off - GMT_CONV_LIMIT));
				if (i_min > nx1) i_min = 0;
				i_max = (int)MIN (nx1, floor (xmax * i_dx_inch - del_off + GMT_CONV_LIMIT));
				if (i_max <= 0 || i_max < i_min) i_max = nx1;
				j_min = (int)MAX (0, ceil ((project_info.ymax - ymax) * i_dy_inch - del_off - GMT_CONV_LIMIT));
				j_max = (int)MIN (ny1, floor ((project_info.ymax - ymin) * i_dy_inch - del_off + GMT_CONV_LIMIT));

				if (out_logic || out_uint) {
					for (i = i_min; i <= i_max; i++) {
						ii = header.ny * (i + 1) - 1;
						for (j = j_min; j <= j_max; j++) {
							if ((side = GMT_non_zero_winding (x[i], y[j], p[k].lon, p[k].lat, p[k].n)) < 
									is_inside) continue;	/* Outside */
							/* Here, point is inside, we must assign value */
							ij = ii - j;
							if (p[k].level > data_8[ij]) data_8[ij] = (unsigned char)p[k].level;
						}
					}
				}
				else {
					for (i = i_min; i <= i_max; i++) {
						ii = header.ny * (i + 1) - 1;
						for (j = j_min; j <= j_max; j++) {
							if ((side = GMT_non_zero_winding (x[i], y[j], p[k].lon, p[k].lat, p[k].n)) < 
									is_inside) continue;	/* Outside */
							/* Here, point is inside, we must assign value */
							ij = ii - j;
							if (p[k].level > data_32[ij]) data_32[ij] = (float)p[k].level;
						}
					}
				}
			}

			GMT_free_polygons (p, np_new);
			GMT_free ((void *)p);
		}

		if (!used_polygons) {	/* Lack of polygons or clipping etc resulted in no polygons after all, must deal with background */

			k = INT_MAX;	/* Initialize to outside range of levels (4 is highest) */
			/* Visit each of the 4 nodes, test if it is inside -R, and if so update lowest level found so far */

			if (!GMT_map_outside (c.lon_sw, c.lat_sw)) k = MIN (k, c.node_level[0]);			/* SW */
			if (!GMT_map_outside (c.lon_sw + c.bsize, c.lat_sw)) k = MIN (k, c.node_level[1]);		/* SE */
			if (!GMT_map_outside (c.lon_sw + c.bsize, c.lat_sw - c.bsize)) k = MIN (k, c.node_level[2]);	/* NE */
			if (!GMT_map_outside (c.lon_sw, c.lat_sw - c.bsize)) k = MIN (k, c.node_level[3]);		/* NW */

			/* If k is still INT_MAX we must assume this patch should have the min level of the bin */

			if (k == INT_MAX) k = MIN (MIN (c.node_level[0], c.node_level[1]) , MIN (c.node_level[2], c.node_level[3]));

			/* Determine nodes to initialize */

			j_min = (int)MAX (0, ceil ((header.y_max - c.lat_sw - c.bsize) * i_dy - del_off));
			j_max = (int)MIN (ny1, floor ((header.y_max - c.lat_sw) * i_dy - del_off));
			i_min = (int)ceil (fmod (c.lon_sw - header.x_min + 360.0, 360.0) * i_dx - del_off);
			i_max = (int)floor (fmod (c.lon_sw + c.bsize - header.x_min + 360.0, 360.0) * i_dx - del_off);
			if (wrap && i_max < i_min) i_max += nx1;
			if (out_logic || out_uint) {
				for (i = i_min; i <= i_max; i++) {
					ii = (wrap) ? i % nx1 : i;
					if (ii < 0 || ii > nx1) continue;
					ii = header.ny * (ii + 1) - 1;
					for (j = j_min; j <= j_max; j++) {
						/*ij = j * header.nx + ii;*/
						ij = ii - j;
						data_8[ij] = (unsigned char)k;
					}
				}
			}
			else {
				for (i = i_min; i <= i_max; i++) {
					ii = (wrap) ? i % nx1 : i;
					if (ii < 0 || ii > nx1) continue;
					ii = header.ny * (ii + 1) - 1;
					for (j = j_min; j <= j_max; j++) {
						ij = ii - j;
						data_32[ij] = (float)k;
					}
				}
			}
		}

		GMT_free_shore (&c);
	}

	GMT_shore_cleanup (&c);

	if (out_logic || out_uint) {
		for (ij = 0; ij < nm; ij++) {
			k = data_8[ij];
			data_8[ij] = (unsigned char)Ctrl->N.mask[k];
		}

		if (wrap && header.node_offset == 0) { /* Copy over values to the repeating right column */
			for (j = ij = 0; j < header.ny; j++, ij += header.nx) data_8[ij+nx1] = data_8[ij];
		}
	}
	else {
		for (ij = 0; ij < nm; ij++) {
			k = irint (data_32[ij]);
			data_32[ij] = (float)Ctrl->N.mask[k];
		}

		if (wrap && header.node_offset == 0) { /* Copy over values to the repeating right column */
			for (j = ij = 0; j < header.ny; j++, ij += header.nx) data_32[ij+nx1] = data_32[ij];
		}
	}
	
	if (temp_shift) {
		header.x_min -= 360.0;
		header.x_max -= 360.0;
	}

	if (nlhs >= 2) {
		plhs[1] = mxCreateDoubleMatrix(1, 9, mxREAL);
		hdr = mxGetPr(plhs[1]);
		hdr[0] = header.x_min;		hdr[1] = header.x_max;
		hdr[2] = header.y_min;		hdr[3] = header.y_max;
		hdr[4] = 0;			hdr[5] = 1;
		hdr[6] = header.node_offset;
		hdr[7] = header.x_inc;		hdr[8] = header.y_inc;
		if (nlhs >= 3) {
			double *xx;
			plhs[2] = mxCreateDoubleMatrix(1, header.nx, mxREAL);
			xx = mxGetPr(plhs[2]);
			for (i = 0; i <= header.nx; i++) xx[i] = header.x_min + (i + header.node_offset) * header.x_inc;
			if (nlhs == 4) {
				double *yy;
				plhs[3] = mxCreateDoubleMatrix(1, header.ny, mxREAL);
				yy = mxGetPr(plhs[3]);
				for (i = 0; i <= header.ny; i++) yy[i] = header.y_min + (i + header.node_offset) * header.y_inc;
			}
		}
	}

	GMT_free ((void *)x);
	GMT_free ((void *)y);

	Free_Grdlandmask_Ctrl (Ctrl);	/* Deallocate control structure */

	GMT_end (argc, argv);
}

void *New_Grdlandmask_Ctrl () {	/* Allocate and initialize a new control structure */
	struct GRDLANDMASK_CTRL *C;
	
	C = (struct GRDLANDMASK_CTRL *) GMT_memory (VNULL, 1, sizeof (struct GRDLANDMASK_CTRL), "New_Grdlandmask_Ctrl");
	
	/* Initialize values whose defaults are not 0/FALSE/NULL */
	
#ifdef GMT_MINOR_VERSION
	C->A.info.high = GMT_MAX_GSHHS_LEVEL;			/* Include all GSHHS levels */
#else
	C->A.high = GMT_MAX_GSHHS_LEVEL;				/* Include all GSHHS levels */
#endif
	C->D.set = 'l';							/* Low-resolution coastline data */
	memset ((void *)C->N.mask, 0, (size_t)(GRDLANDMASK_N_CLASSES * sizeof (float)));	/* Default "wet" value = 0 */
	C->N.mask[1] = C->N.mask[3] = 1.0;				/* Default for "dry" areas = 1 (inside) */
	
	return ((void *)C);
}

void Free_Grdlandmask_Ctrl (struct GRDLANDMASK_CTRL *C) {	/* Deallocate control structure */
	GMT_free ((void *)C);	
}
