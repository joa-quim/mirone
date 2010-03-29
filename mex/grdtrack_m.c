/*--------------------------------------------------------------------
 *	$Id: grdtrack.c,v 1.19 2004/07/16 17:07:02 pwessel Exp $
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
 * grdtrack reads a xyfile, opens the 2d binary gridded grdfile, 
 * and samples the dataset at the xy positions with a bilinear or bicubic
 * interpolant.  This new data is added to the input as an extra column
 * and printed to standard output.  In order to evaluate derivatives along
 * the edges of the grdfile region, we assume natural bicubic spline
 * boundary conditions (d2z/dn2 = 0, n being the normal to the edge;
 * d2z/dxdy = 0 in the corners).  Rectangles of size x_inc by y_inc are 
 * mapped to [0,1] x [0,1] by affine transformation, and the interpolation
 * done on the normalized rectangle.
 *
 * Author:	Walter H F Smith
 * Date:	23-SEP-1993
 * 
 * Based on the original grdtrack, which had this authorship/date/history:
 *
 * Author:	Paul Wessel
 * Date:	29-JUN-1988
 * Revised:	5-JAN-1990	PW: Updated to v.2.0
 *		4-AUG-1993	PW: Added -Q
 *		14-AUG-1998	PW: GMT 3.1
 *  Modified:	10 Jul 2000 3.3.5  by PW to allow plain -L to indicate geographic coordinates
 * Version:	4
 *
 *--------------------------------------------------------------------*/
/*
 * Mexified version of grdtrack
 * Author:	J. Luis
 * Date: 	17 Aug 2004
 *
 * Usage
 * out = grdtrack_m(Zin,head,xy,'options');
 *
 * where	Z is the array containing the input file
 *		head is the header descriptor in the format of grdread/grdwrite mexs
 *		xy is an array with at least two columns (first column x and second y)
 *		and options may be for example: '-Z'
 *
 * IMPORTANT NOTE. The data type of Zin can be any of: double, single, Int32, Int16,
 * Uint16 or Uint8 but "out" will allways be of type double
 *	 
 *		04/06/06 J Luis, Updated to compile with version 4.1.3
 *		14/10/06 J Luis, Now includes the memory leak solving solution
 *		03/03/08 J Luis, Adapted for 4.2.1
 *		12/07/08 J Luis, Made it standalone (no GMT lib dependency)
 *		13/09/08 J Luis, Added column major code (when Z input in singles)
 */

#include "mex.h"
#include <math.h>
#include <string.h>
#include <float.h>

#define GMT_SMALL		1.0e-4	/* Needed when results aren't exactly zero but close */
#define GMT_CONV_LIMIT	1.0e-8	/* Fairly tight convergence limit or "close to zero" limit */

#define	FALSE	0
#define	TRUE	1
#ifndef M_PI
#define M_PI	3.14159265358979323846
#endif

#define D2R (M_PI / 180.0)
#define R2D (180.0 / M_PI)
#define cosd(x) cos ((x) * D2R)
#define d_sqrt(x) ((x) < 0.0 ? 0.0 : sqrt (x))
#define CNULL	((char *)NULL)
#define Loc_copysign(x,y) ((y) < 0.0 ? -fabs(x) : fabs(x))

#ifndef M_SQRT2
#define	M_SQRT2		1.41421356237309504880
#endif

#ifndef MIN
#define MIN(x, y) (((x) < (y)) ? (x) : (y))	/* min and max value macros */
#endif
#ifndef MAX
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#endif

#ifndef rint
#define rint(x) (floor((x)+0.5))
#endif
#ifndef irint
#define irint(x) ((int)rint(x))
#endif

#define EQ_RAD 6371.0087714
#define M_PR_DEG (EQ_RAD * 1000 * M_PI / 180.0)

#define BCR_NEARNEIGHBOR	0
#define BCR_BILINEAR		1
#define BCR_BSPLINE		2
#define BCR_BICUBIC		3

#define original_GMT_code	0	/* Set it to 1 to always use original (rowmajor) code */

struct GRD_HEADER {
/* Do not change the first three items. They are copied verbatim to the native grid header */
	int nx;				/* Number of columns */
	int ny;				/* Number of rows */
	int node_offset;		/* 0 for node grids, 1 for pixel grids */
/* This section is flexible. It is not copied to any grid header */
	int type;			/* Grid format */
	char name[256];			/* Actual name of the file after any ?<varname> and =<stuff> has been removed */
	char varname[80];		/* NetCDF: variable name */
	int y_order;			/* NetCDF: 1 if S->N, -1 if N->S */
	int z_id;			/* NetCDF: id of z field */
	int ncid;			/* NetCDF: file ID */
	int t_index[3];			/* NetCDF: index of higher coordinates */
	double nan_value;		/* Missing value as stored in grid file */
	double xy_off;			/* 0.0 (node_offset == 0) or 0.5 ( == 1) */
/* The following elements should not be changed. They are copied verbatim to the native grid header */
	double x_min;			/* Minimum x coordinate */
	double x_max;			/* Maximum x coordinate */
	double y_min;			/* Minimum y coordinate */
	double y_max;			/* Maximum y coordinate */
	double z_min;			/* Minimum z value */
	double z_max;			/* Maximum z value */
	double x_inc;			/* x increment */
	double y_inc;			/* y increment */
	double z_scale_factor;		/* grd values must be multiplied by this */
	double z_add_offset;		/* After scaling, add this */
	char x_units[80];		/* units in x-direction */
	char y_units[80];		/* units in y-direction */
	char z_units[80];		/* grid value units */
	char title[80];			/* name of data set */
	char command[320];		/* name of generating command */
	char remark[160];		/* comments re this data set */
}; 

struct GMT_EDGEINFO {
	/* Description below is the final outcome after parse and verify */
	int	nxp;	/* if X periodic, nxp > 0 is the period in pixels  */
	int	nyp;	/* if Y periodic, nxp > 0 is the period in pixels  */
	int	gn;	/* TRUE if top    edge will be set as N pole  */
	int	gs;	/* TRUE if bottom edge will be set as S pole  */
};

struct GMT_BCR {	/* Used mostly in gmt_support.c */
	double	rx_inc;			/* 1.0 / grd.x_inc  */
	double	ry_inc;			/* 1.0 / grd.y_inc  */
	double	offset;			/* 0 or 0.5 for grid or pixel registration  */
	double	threshold;		/* sum of cardinals must >= threshold in bilinear; else NaN */
	int	interpolant;		/* Interpolation function used (0, 1, 2, 3) */
	int	n;			/* Width of the interpolation function */
	int	ioff;			/* Padding on west side of array  */
	int	joff;			/* Padding on north side of array  */
	int	mx;			/* Padded array dimension  */
	int	my;			/* Ditto  */
};

void GMT_boundcond_init (struct GMT_EDGEINFO *edgeinfo);
int GMT_boundcond_set (struct GRD_HEADER *h, struct GMT_EDGEINFO *edgeinfo, int *pad, float *a);
int GMT_boundcond_param_prep (struct GRD_HEADER *h, struct GMT_EDGEINFO *edgeinfo);
int GMT_boundcond_parse (struct GMT_EDGEINFO *edgeinfo, char *edgestring);
double GMT_get_bcr_z (struct GRD_HEADER *grd, double xx, double yy, float *data, struct GMT_EDGEINFO *edgeinfo, struct GMT_BCR *bcr, int row_maj);
void GMT_bcr_init (struct GRD_HEADER *grd, int *pad, int interpolant, double threshold, struct GMT_BCR *bcr);
double row_or_column_major (struct GMT_BCR *bcr, float *data, double *wx, double *wy, int ij, int row_maj);

int decode_R (char *item, double *w, double *e, double *s, double *n);
int check_region (double w, double e, double s, double n);
double ddmmss_to_degree (char *text);

/* --------------------------------------------------------------------------- */
/* Matlab Gateway routine */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	int i, j, nx, ny, n_read = 0, n_points = 0, one_or_zero;
	int n_output = 0, n_fields, n_pts, ii, jj, GMT_pad[4];
	
	int error = FALSE, suppress = FALSE, node = FALSE, z_only = FALSE;
	int is_double = FALSE, is_single = FALSE, is_int32 = FALSE, is_int16 = FALSE;
	int is_uint16 = FALSE, is_uint8 = FALSE, is_int8 = FALSE;
	int free_copy = TRUE, need_padding = FALSE, row_maj = FALSE;
	
	double value, west, east, south, north, threshold = 1.0, i_dx, i_dy, half, *in, *out;
	float *f;

	int	i2, argc = 0, nc_h, nr_h, mx, n_arg_no_char = 0, *i_4, interpolant = BCR_BICUBIC;
	short int *i_2;
	unsigned short int *ui_2;
	char	**argv, *i_1;
	unsigned char *ui_1;
	float	*z_4;
	double	*pdata_d, *z_8, *head;
	
	struct GRD_HEADER grd;
	struct GMT_EDGEINFO edgeinfo;
	struct GMT_BCR bcr;

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
	argv[0] = "grdtrack_m";
	for (i = 1; i < argc; i++)
		argv[i] = (char *)mxArrayToString(prhs[i+n_arg_no_char-1]);

	west = east = south = north = 0.0;
	
	GMT_boundcond_init (&edgeinfo);

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {
			
				case 'R':
					error += decode_R (argv[i], &west, &east, &south, &north);
					break;
				case 'L':
					if (argv[i][2]) {
						error += GMT_boundcond_parse (&edgeinfo, &argv[i][2]);
						/*if (edgeinfo.gn) {
							GMT_io.in_col_type[0] = GMT_io.out_col_type[0] = GMT_IS_LON;
							GMT_io.in_col_type[1] = GMT_io.out_col_type[1] = GMT_IS_LAT;
						}*/
					}
					/*else {
						GMT_io.in_col_type[0] = GMT_io.out_col_type[0] = GMT_IS_LON;
						GMT_io.in_col_type[1] = GMT_io.out_col_type[1] = GMT_IS_LAT;
					}*/
					break;
				case 'N':
					node = TRUE;
					break;
				case 'Q':
					interpolant = BCR_BILINEAR;
					threshold = (argv[i][2]) ? atof (&argv[i][2]) : 1.0;
					break;
				case 'S':
					suppress = TRUE;
					break;
				case 'Z':
					z_only = TRUE;
					break;
				default:
					error = TRUE;
					break;
			}
		}
	}
	
	if (argc == 1 || error) {
		mexPrintf ("grdtrack - Sampling of a 2-D gridded netCDF grdfile along 1-D trackline\n\n");
		mexPrintf ("usage: out = grdtrack_m(grd,head,xydata, ['-L<flag>'], ['-N']\n"); 
		mexPrintf ("\t['-Q[<value>]'], ['-R<west/east/south/north>[r]'] ['-S'] ['-Z'] ['-f[i|o]<colinfo>']\n");
				
		mexPrintf ("\t<xydata> is an multicolumn array with (lon,lat) in the first two columns\n");
		mexPrintf ("\n\tOPTIONS:\n");
		mexPrintf ("\t-L sets boundary conditions.  <flag> can be either\n");
		mexPrintf ("\t   g for geographic boundary conditions\n");
		mexPrintf ("\t   or one or both of\n");
		mexPrintf ("\t   x for periodic boundary conditions on x\n");
		mexPrintf ("\t   y for periodic boundary conditions on y\n");
		mexPrintf ("\t-N Report value at nearest node instead of interpolating\n");
		mexPrintf ("\t-Q Quick mode, use bilinear rather than bicubic interpolation.\n");
		mexPrintf ("\t   Optionally, append <value> in the 0 < value <= 1 range.\n");
		mexPrintf ("\t   [Default = 1 requires all 4 nodes to be non-NaN.], <value> = 0.5\n");
		mexPrintf ("\t   will interpolate about 1/2 way from a non-NaN to a NaN node, while\n");
		mexPrintf ("\t   0.1 will go about 90%% of the way, etc.\n");
		mexPrintf ("\t-R specifies a subregion [Default is old region]\n");
		mexPrintf ("\t-S Suppress output when result equals NaN\n");
		mexPrintf ("\t-Z only output z-values [Default gives all columns]\n");
		mexPrintf ("\n\tSECRET INFO:\n");
		mexPrintf ("\t   When input points are inside the outer skirt of 2 rows and columns of\n");
		mexPrintf ("\t   the 2-D grid we don't need to set boundary conditions and as\n");
		mexPrintf ("\t   such we can use the input array without further to C order\n");
		mexPrintf ("\t   conversion (to row major). This save a lot of memory and execution\n");
		mexPrintf ("\t   time. However, this possibility works only when input 2-D array is\n");
		mexPrintf ("\t   of tipe SINGLE (though the computations are all done in doubles).\n");
		mexPrintf ("\t   The other cases request using a temporary array of size (M+2)x(N+2)\n");
		return;
	}

	if (threshold <= 0.0 || threshold > 1.0) {
		mexPrintf ("GRDTRACK_M SYNTAX ERROR -Q:  threshold must be in <0,1] range\n");
		error++;
	}
	if (error) return;

	if (nlhs == 0) {
		mexPrintf("ERROR: Must provide an output.\n");
		return;
	}

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
	else if (mxIsInt8(prhs[0])) {
		i_1 = mxGetData(prhs[0]);
		is_int8 = TRUE;
	}
	else {
		mexPrintf("GRDTRACK ERROR: Unknown input data type.\n");
		mexErrMsgTxt("Valid types are:double, single, Int32, Int16, UInt16, UInt8 and Int8.\n");
	}

	nx = mxGetN (prhs[0]);
	ny = mxGetM (prhs[0]);
	if (!mxIsNumeric(prhs[0]) || ny < 2 || nx < 2)
		mexErrMsgTxt("First argument must contain a decent array\n");

	nc_h = mxGetN (prhs[1]);
	nr_h = mxGetM (prhs[1]);
	if (!mxIsNumeric(prhs[1]) || nr_h > 1 || nc_h < 9)
		mexErrMsgTxt("Second argument must contain a valid header of the input array.\n");
	
	head  = mxGetPr(prhs[1]);		/* Get header info */

	/* Check that thirth argument contains at least a mx2 table */
	n_pts = mxGetM (prhs[2]);
	n_fields = mxGetN(prhs[2]);
	if (!mxIsNumeric(prhs[2]) || (n_fields < 2))
		mexErrMsgTxt("GRDTRACK ERROR: thirth argument must contain the x,y positions where to interpolate.\n");
	if (z_only)
		n_fields = 0;

	/* Read the interpolation points and convert them to double */
	if (mxIsDouble(prhs[2]))
		in = mxGetPr(prhs[2]);
	else if (mxIsSingle(prhs[2]))
		in = mxGetData(prhs[2]);

	grd.x_min = head[0];	grd.x_max = head[1];
	grd.y_min = head[2];	grd.y_max = head[3];
	grd.z_min = head[4];	grd.z_max = head[5];
	grd.x_inc = head[7];	grd.y_inc = head[8];
	grd.nx = nx;			grd.ny = ny;
	grd.node_offset = irint(head[6]);
	mx = nx + 4;

	if (west == east) {	/* No subset asked for */
		west = grd.x_min;
		east = grd.x_max;
		south = grd.y_min;
		north = grd.y_max;
	}
	one_or_zero = (grd.node_offset) ? 0 : 1;
	half = (grd.node_offset) ? 0.5 : 0.0;
	nx = irint ( (east - west) / grd.x_inc) + one_or_zero;
	ny = irint ( (north - south) / grd.y_inc) + one_or_zero;
	i_dx = 1.0 / grd.x_inc;
	i_dy = 1.0 / grd.y_inc;

	if (!node) {
		/* If we don't have any point inside the two outer row/columns
		   there is no need to set boundary conditions plus all the extra
		   ovehead that it implies. So check it out here. */
		int n;
		double this_xmin, this_xmax, this_ymin, this_ymax;
		n = (interpolant == BCR_BILINEAR) ? 1 : 2;
		this_xmin = grd.x_min + n * grd.x_inc;
		this_xmax = grd.x_max - n * grd.x_inc;
		this_ymin = grd.y_min + n * grd.y_inc;
		this_ymax = grd.y_max - n * grd.y_inc;
		for (i = 0; i < n_pts; i++) {
			if (in[i] < this_xmin || in[i] > this_xmax) {
				need_padding = TRUE;
				break;
			}
			if (in[i+n_pts] < this_ymin || in[i+n_pts] > this_ymax) {
				need_padding = TRUE;
				break;
			}
		}
	}

#if original_GMT_code
	need_padding = TRUE;
#endif

	if (need_padding) row_maj = TRUE;	/* Here we have to use the old row major code */

	if (!need_padding) {		/* We can use the column major order of the Matlab array */

		if (!is_single)
			f = mxCalloc (nx * ny, sizeof (float));

		if (is_double) 
			for (j = 0; j < nx*ny; j++) f[j] = (float)z_8[j];

		else if (is_single) {
			f = z_4;
			free_copy = FALSE;	/* Signal that we shouldn't free f */
		}

		else if (is_int32)
			for (j = 0; j < nx*ny; j++) f[j] = (float)i_4[j];

		else if (is_int16)
			for (j = 0; j < nx*ny; j++) f[j] = (float)i_2[j];

		else if (is_uint16)
			for (j = 0; j < nx*ny; j++) f[j] = (float)ui_2[j];

		else if (is_uint8)
			for (j = 0; j < nx*ny; j++) f[j] = (float)ui_1[j];

		else if (is_int8)
			for (j = 0; j < nx*ny; j++) f[j] = (float)i_1[j];

		GMT_pad[0] = GMT_pad[1] = GMT_pad[2] = GMT_pad[3] = 0;
	}

	else {

		f = mxCalloc ((nx+4)*(ny+4), sizeof (float));

		/* Transpose from Matlab orientation to gmt grd orientation */
		if (is_double) {
			for (i = 0, i2 = ny - 1; i < ny; i++, i2--) {
				ii = (i2 + 2)*mx + 2;
				for (j = 0; j < nx; j++) f[ii + j] = (float)z_8[j*ny+i];
			}
		}
		else if (is_single) {
			for (i = 0, i2 = ny - 1; i < ny; i++, i2--) {
				ii = (i2 + 2)*mx + 2;
				for (j = 0; j < nx; j++) f[ii + j] = z_4[j*ny+i];
			}
		}
		else if (is_int32) {
			for (i = 0, i2 = ny - 1; i < ny; i++, i2--) {
				ii = (i2 + 2)*mx + 2;
				for (j = 0; j < nx; j++) f[ii + j] = (float)i_4[j*ny+i];
			}
		}
		else if (is_int16) {
			for (i = 0, i2 = ny - 1; i < ny; i++, i2--) {
			ii = (i2 + 2)*mx + 2;
				for (j = 0; j < nx; j++) f[ii + j] = (float)i_2[j*ny+i];
			}
		}
		else if (is_uint16) {
			for (i = 0, i2 = ny - 1; i < ny; i++, i2--) {
				ii = (i2 + 2)*mx + 2;
				for (j = 0; j < nx; j++) f[ii + j] = (float)ui_2[j*ny+i];
			}
		}
		else if (is_uint8) {
			for (i = 0, i2 = ny - 1; i < ny; i++, i2--) {
				ii = (i2 + 2)*mx + 2;
				for (j = 0; j < nx; j++) f[ii + j] = (float)ui_1[j*ny+i];
			}
		}
		else if (is_int8) {
			for (i = 0, i2 = ny - 1; i < ny; i++, i2--) {
				ii = (i2 + 2)*mx + 2;
				for (j = 0; j < nx; j++) f[ii + j] = (float)i_1[j*ny+i];
			}
		}

		GMT_pad[0] = GMT_pad[1] = GMT_pad[2] = GMT_pad[3] = 2;

		GMT_boundcond_param_prep (&grd, &edgeinfo);
	}
	
	/*project_info.w = west;	project_info.e = east;
	project_info.s = south;	project_info.n = north;*/
	
	/* Initialize bcr structure:  */

	GMT_bcr_init (&grd, GMT_pad, interpolant, threshold, &bcr);

	if (need_padding)
		/* Set boundary conditions  */
		GMT_boundcond_set (&grd, &edgeinfo, GMT_pad, f);
	
	if ((out = mxCalloc(n_pts * (n_fields+1), sizeof (double))) == 0)
		mexErrMsgTxt("GRDTRACK ERROR: Could not allocate memory\n");

	for (i = 0; i < n_pts; i++) {
		while ( (mxIsNaN(in[i]) || mxIsNaN(in[i+n_pts])) && !z_only) {
			for (j = 0; j < n_fields; j++) out[j*n_pts+i] = in[j*n_pts+i];
			out[j*n_pts+i] = mxGetNaN();
			i++;
		}

		/* If point is outside grd area, shift it using periodicity or skip if not periodic. */
		while ( (in[i+n_pts] < grd.y_min) && (edgeinfo.nyp > 0) ) in[i+n_pts] += (grd.y_inc * edgeinfo.nyp);
		if (in[i+n_pts] < grd.y_min) continue;

		while ( (in[i+n_pts] > grd.y_max) && (edgeinfo.nyp > 0) ) in[i+n_pts] -= (grd.y_inc * edgeinfo.nyp);
		if (in[i+n_pts] > grd.y_max) continue;

		while ( (in[i] < grd.x_min) && (edgeinfo.nxp > 0) ) in[i] += (grd.x_inc * edgeinfo.nxp);
		if (in[i] < grd.x_min) continue;

		while ( (in[i] > grd.x_max) && (edgeinfo.nxp > 0) ) in[i] -= (grd.x_inc * edgeinfo.nxp);
		if (in[i] > grd.x_max) continue;
		
		if (node) {
			ii = irint ((in[i] - grd.x_min) * i_dx - half) + one_or_zero;
			jj = irint ((grd.y_max - in[i+n_pts]) * i_dy - half) + one_or_zero;
			value = f[(jj+GMT_pad[3])*mx+ii+GMT_pad[0]];
		}
		else
			value = GMT_get_bcr_z(&grd, in[i], in[i+n_pts], f, &edgeinfo, &bcr, row_maj);

		if (suppress && mxIsNaN (value)) continue;

		if (z_only) {	/* Simply print out value */
			out[i] = value;
		}
		else {	/* Simply copy other columns, append value, and output */
			for (j = 0; j < n_fields; j++) out[j*n_pts+i] = in[j*n_pts+i];
			out[j*n_pts+i] = value;
		}
	}

	/*if (!(!need_padding && !is_single)) {
		mexPrintf("Merda vou Friar %d\t%d\n", need_padding, is_single);
		mxFree((void *)f);
	}*/
	if (free_copy)
		mxFree((void *)f);

	plhs[0] = mxCreateDoubleMatrix (n_pts,n_fields+1, mxREAL);
	pdata_d = mxGetPr(plhs[0]);
	memcpy(pdata_d, out, n_pts*(n_fields+1)*8);
	mxFree(out);
}


void GMT_bcr_init (struct GRD_HEADER *grd, int *pad, int interpolant, double threshold, struct GMT_BCR *bcr) {
	/* Initialize interpolant and threshold */
	bcr->interpolant = interpolant;
	bcr->threshold = threshold;
	if (interpolant == BCR_NEARNEIGHBOR)
		bcr->n = 1;
	else if (interpolant == BCR_BILINEAR)
		bcr->n = 2;
	else
		bcr->n = 4;

	/* Initialize ioff, joff, mx, my according to grd and pad:  */
	bcr->ioff = pad[0];
	bcr->joff = pad[3];
	bcr->mx = (int)(grd->nx + pad[0] + pad[1]);
	bcr->my = (int)(grd->ny + pad[2] + pad[3]);

	/* Initialize rx_inc, ry_inc, and offset:  */
	bcr->rx_inc = 1.0 / grd->x_inc;
	bcr->ry_inc = 1.0 / grd->y_inc;
	bcr->offset = (grd->node_offset) ? 0.5 : 0.0;
}

double GMT_get_bcr_z (struct GRD_HEADER *grd, double xx, double yy, float *data, struct GMT_EDGEINFO *edgeinfo, struct GMT_BCR *bcr, int row_maj) {
	/* Given xx, yy in user's grid file (in non-normalized units)
	   this routine returns the desired interpolated value (nearest-neighbor, bilinear
	   B-spline or bicubic) at xx, yy. */

	int i, j, ij;
	double	x, y, wx[4], wy[4], w, wp, wq;

	/* First check that xx,yy are not Nan - if so return NaN */
	
	if (mxIsNaN (xx) || mxIsNaN (yy)) return (mxGetNaN());
	
	/* First check if the xx and yy are within the grid.
	   16-Sep-2007: Added some slack (GMT_SMALL) here to avoid setting to NaN points
	   that are really on the edge but because of rounding errors are regarded outside.
	   Remember that we have padded the grid with 2 extra values, so this should not be
	   a problem. */

	if (xx < grd->x_min - GMT_SMALL || xx > grd->x_max + GMT_SMALL) return (mxGetNaN());
	if (yy < grd->y_min - GMT_SMALL || yy > grd->y_max + GMT_SMALL) return (mxGetNaN());

	/* Compute the normalized real indices (x,y) of the point (xx,yy) within the grid.
	   Note that the y axis points down from the upper left corner of the grid. */

	x = (xx - grd->x_min) * bcr->rx_inc - bcr->offset;
	y = (grd->y_max - yy) * bcr->ry_inc - bcr->offset;

	if (bcr->interpolant == BCR_NEARNEIGHBOR) {
		/* Find the indices (i,j) of the closest node. */
		i = irint(x);
		j = irint(y);
	}
	else {
		/* Find the indices (i,j) of the node to the upper left of that.
	   	   Because of padding, i and j can be on the edge. */
		i = (int)floor(x);
		j = (int)floor(y);

		/* Determine the offset of (x,y) with respect to (i,j). */
		x -= (double)i;
		y -= (double)j;

		/* For 4x4 interpolants, move over one more cell to the upper left corner */
		if (bcr->n == 4) { i--; j--; }
	}

	/* Normally, one would expect here a check on the value (i,j) to make sure that the
	   corners of the convolution kernel, (i,j) and (i+bcr->n-1,j+bcr->n-1), are both within
	   the padded grid. However, the check on (xx, yy) above, even with the slack, ensures
	   that the corner points are between (-2,-2) and (grd->nx+1,grd->ny+1), the corners
	   of the padding.

	if (i < -2 || j < -2 || i+bcr->n > grd_nx+2 || j+bcr->n > grd_ny+2) return (GMT_d_NaN);
	*/

	/* Build weights */

	switch (bcr->interpolant) {
	case BCR_NEARNEIGHBOR:
		wx[0] = wy[0] = 1.0;
		break;
	case BCR_BILINEAR:
		/* Simple 1-D linear weights */
		wx[0] = 1.0 - x;
		wx[1] = x;

		wy[0] = 1.0 - y;
		wy[1] = y;
		break;
	case BCR_BSPLINE:
		/* These are B-spline weights */
		wp = x * x;
		wq = wp * x;
		wx[1] = wq / 2 - wp + 2.0 / 3.0;
		wx[3] = wq / 6;
		w = 1.0 - x;
		wp = w * w;
		wq = wp * w;
		wx[2] = wq / 2 - wp + 2.0 / 3.0;
		wx[0] = wq / 6;

		wp = y * y;
		wq = wp * y;
		wy[1] = wq / 2 - wp + 2.0 / 3.0;
		wy[3] = wq / 6;
		w = 1.0 - y;
		wp = w * w;
		wq = wp * w;
		wy[2] = wq / 2 - wp + 2.0 / 3.0;
		wy[0] = wq / 6;
		break;
	default:
		/* These weights are based on the cubic convolution kernel, see for example
		   http://undergraduate.csse.uwa.edu.au/units/CITS4241/Handouts/Lecture04.html
		   These weights include a free parameter (a), which is set to -0.5 in this case.

		   In the absence of NaNs, the result of this is identical to the scheme introduced
		   by Walter Smith. The current implementation, however, is much less complex, faster,
		   allows NaNs to be skipped, and much more similar to the bilinear case.

		   Remko Scharroo, 10 Sep 2007.
		*/
		w = 1.0 - x;
		wp = w * x;
		wq = -0.5 * wp;
		wx[0] = wq * w;
		wx[3] = wq * x;
		wx[1] = 3 * wx[3] + w + wp;
		wx[2] = 3 * wx[0] + x + wp;

		w = 1.0 - y;
		wp = w * y;
		wq = -0.5 * wp;
		wy[0] = wq * w;
		wy[3] = wq * y;
		wy[1] = 3 * wy[3] + w + wp;
		wy[2] = 3 * wy[0] + y + wp;
		break;
	}

	/* Save the location of the upper left corner point of the convolution kernel */
	if (row_maj)		/* C order */
		ij = (j + bcr->joff) * bcr->mx + (i + bcr->ioff);
	else
		/* I'm not sure why we need that ... - j - 1. I think it is because in
		   original GMT code y has origin ij upper left corner and j = (int)floor(y); 
		   above. Whilst here Y origin his at lower left. So it might be that floor().
		   Anyway, I tested it against GMT grdtrack and it gave the same result.
		   So I believe this is ok. The NEARNEIGHBOR  case was not tested, but we don't
		   use it in Mirone. 		*/
		ij = (i + bcr->ioff) * bcr->my + (bcr->my - j - 1 + bcr->joff);

	return (row_or_column_major (bcr, data, wx, wy, ij, row_maj));
}

double row_or_column_major (struct GMT_BCR *bcr, float *data, double *wx, double *wy, int ij, int row_maj) {

	int i, j;
	double	x, y, retval, wsum, w;

	retval = wsum = 0.0;
	if (row_maj) {	/* C order */
		for (j = 0; j < bcr->n; j++) {
			for (i = 0; i < bcr->n; i++) {
				if (!mxIsNaN(data[ij+i])) {
					w = wx[i] * wy[j];
					retval += data[ij+i] * w;
					wsum += w;
				}
			}
			ij += bcr->mx;
		}
	}
	else {		/* Matlab (fortran) ordering */
		for (i = 0; i < bcr->n; i++) {
			for (j = 0; j < bcr->n; j++) {
				if (!mxIsNaN(data[ij-j])) {
					w = wx[i] * wy[j];
					retval += data[ij-j] * w;
					wsum += w;
				}
			}
			ij += bcr->my;
		}
	}

	return ( ((wsum + GMT_CONV_LIMIT - bcr->threshold) > 0.0) ? retval / wsum : mxGetNaN());
}

void GMT_boundcond_init (struct GMT_EDGEINFO *edgeinfo) {
	edgeinfo->nxp = 0;
	edgeinfo->nyp = 0;
	edgeinfo->gn = FALSE;
	edgeinfo->gs = FALSE;
	return;
}

int GMT_boundcond_parse (struct GMT_EDGEINFO *edgeinfo, char *edgestring) {
	/* Parse string beginning at argv[i][2] and load user's
		requests in edgeinfo->  Return success or failure.
		Requires that edgeinfo previously initialized to
		zero/FALSE stuff.  Expects g or (x and or y) is
		all that is in string.  */

	int	i, ier;

	i = 0;
	ier = FALSE;
	while (!ier && edgestring[i]) {
		switch (edgestring[i]) {
			case 'g':
			case 'G':
				edgeinfo->gn = TRUE;
				edgeinfo->gs = TRUE;
				break;
			case 'x':
			case 'X':
				edgeinfo->nxp = -1;
				break;
			case 'y':
			case 'Y':
				edgeinfo->nyp = -1;
				break;
			default:
				ier = TRUE;
				break;

		}
		i++;
	}

	if (ier) return (-1);

	return (0);
}

int GMT_boundcond_param_prep (struct GRD_HEADER *h, struct GMT_EDGEINFO *edgeinfo) {
	/* Called when edgeinfo holds user's choices.  Sets edgeinfo according to choices and h.  */

	double	xtest;

	/*if (edgeinfo->gn || GMT_grd_is_global(h)) {	I'm not using the is_global test */
	if (edgeinfo->gn) {
		/* User has requested geographical conditions.  */
		if ( (h->x_max - h->x_min) < (360.0 - GMT_SMALL * h->x_inc) ) {
			mexPrintf ("GRDTRACK_M Warning: x range too small; g boundary condition ignored.\n");
			edgeinfo->nxp = edgeinfo->nyp = 0;
			edgeinfo->gn  = edgeinfo->gs = FALSE;
			return (0);
		}
		xtest = fmod (180.0, h->x_inc) / h->x_inc;
		/* xtest should be within GMT_SMALL of zero or of one.  */
		if ( xtest > GMT_SMALL && xtest < (1.0 - GMT_SMALL) ) {
			/* Error.  We need it to divide into 180 so we can phase-shift at poles.  */
			mexPrintf ("GRDTRACK_M Warning: x_inc does not divide 180; g boundary condition ignored.\n");
			edgeinfo->nxp = edgeinfo->nyp = 0;
			edgeinfo->gn  = edgeinfo->gs = FALSE;
			return (0);
		}
		edgeinfo->nxp = irint(360.0/h->x_inc);
		edgeinfo->nyp = 0;
		edgeinfo->gn = ( (fabs(h->y_max - 90.0) ) < (GMT_SMALL * h->y_inc) );
		edgeinfo->gs = ( (fabs(h->y_min + 90.0) ) < (GMT_SMALL * h->y_inc) );
	}
	else {
		if (edgeinfo->nxp != 0) edgeinfo->nxp = (h->node_offset) ? h->nx : h->nx - 1;
		if (edgeinfo->nyp != 0) edgeinfo->nyp = (h->node_offset) ? h->ny : h->ny - 1;
	}
	return (0);
}

int GMT_boundcond_set (struct GRD_HEADER *h, struct GMT_EDGEINFO *edgeinfo, int *pad, float *a) {
	/* Set two rows of padding (pad[] can be larger) around data according
		to desired boundary condition info in edgeinfo.
		Returns -1 on problem, 0 on success.
		If either x or y is periodic, the padding is entirely set.
		However, if neither is true (this rules out geographical also)
		then all but three corner-most points in each corner are set.

		As written, not ready to use with "surface" for GMT v4, because
		assumes left/right is +/- 1 and down/up is +/- mx.  In "surface"
		the amount to move depends on the current mesh size, a parameter
		not used here.

		This is the revised, two-rows version (WHFS 6 May 1998).
	*/

	int	bok;	/* Counter used to test that things are OK  */
	int	mx;	/* Width of padded array; width as malloc'ed  */
	int	mxnyp;	/* distance to periodic constraint in j direction  */
	int	i, jmx;	/* Current i, j * mx  */
	int	nxp2;	/* 1/2 the xg period (180 degrees) in cells  */
	int	i180;	/* index to 180 degree phase shift  */
	int	iw, iwo1, iwo2, iwi1, ie, ieo1, ieo2, iei1;  /* see below  */
	int	jn, jno1, jno2, jni1, js, jso1, jso2, jsi1;  /* see below  */
	int	jno1k, jno2k, jso1k, jso2k, iwo1k, iwo2k, ieo1k, ieo2k;
	int	j1p, j2p;	/* j_o1 and j_o2 pole constraint rows  */


	/* Check pad  */
	bok = 0;
	for (i = 0; i < 4; i++) {
		if (pad[i] < 2) bok++;
	}
	if (bok > 0) {
		mexPrintf ("GMT BUG:  bad pad for GMT_boundcond_set.\n");
		return (-1);
	}

	/* Check minimum size:  */
	if (h->nx < 2 || h->ny < 2) {
		mexPrintf ("GMT ERROR:  GMT_boundcond_set requires nx,ny at least 2.\n");
		return (-1);
	}

	/* Initialize stuff:  */

	mx = h->nx + pad[0] + pad[1];
	nxp2 = edgeinfo->nxp / 2;	/* Used for 180 phase shift at poles  */

	iw = pad[0];		/* i for west-most data column */
	iwo1 = iw - 1;		/* 1st column outside west  */
	iwo2 = iwo1 - 1;	/* 2nd column outside west  */
	iwi1 = iw + 1;		/* 1st column  inside west  */

	ie = pad[0] + h->nx - 1;	/* i for east-most data column */
	ieo1 = ie + 1;		/* 1st column outside east  */
	ieo2 = ieo1 + 1;	/* 2nd column outside east  */
	iei1 = ie - 1;		/* 1st column  inside east  */

	jn = mx * pad[3];	/* j*mx for north-most data row  */
	jno1 = jn - mx;		/* 1st row outside north  */
	jno2 = jno1 - mx;	/* 2nd row outside north  */
	jni1 = jn + mx;		/* 1st row  inside north  */

	js = mx * (pad[3] + h->ny - 1);	/* j*mx for south-most data row  */
	jso1 = js + mx;		/* 1st row outside south  */
	jso2 = jso1 + mx;	/* 2nd row outside south  */
	jsi1 = js - mx;		/* 1st row  inside south  */

	mxnyp = mx * edgeinfo->nyp;

	jno1k = jno1 + mxnyp;	/* data rows periodic to boundary rows  */
	jno2k = jno2 + mxnyp;
	jso1k = jso1 - mxnyp;
	jso2k = jso2 - mxnyp;

	iwo1k = iwo1 + edgeinfo->nxp;	/* data cols periodic to bndry cols  */
	iwo2k = iwo2 + edgeinfo->nxp;
	ieo1k = ieo1 - edgeinfo->nxp;
	ieo2k = ieo2 - edgeinfo->nxp;

	/* Check poles for grid case.  It would be nice to have done this
		in GMT_boundcond_param_prep() but at that point the data
		array isn't passed into that routine, and may not have been
		read yet.  Also, as coded here, this bombs with error if
		the pole data is wrong.  But there could be an option to
		to change the condition to Natural in that case, with warning.  */

	if (h->node_offset == 0) {
		if (edgeinfo->gn) {
			bok = 0;
			if (mxIsNaN (a[jn + iw])) {
				for (i = iw+1; i <= ie; i++) {
					if (!mxIsNaN (a[jn + i])) bok++;
				}
			}
			else {
				for (i = iw+1; i <= ie; i++) {
					if (a[jn + i] != a[jn + iw]) bok++;
				}
			}
			if (bok > 0) mexPrintf ("GRDTRACK_M Warning: Inconsistent grid values at North pole.\n");
		}

		if (edgeinfo->gs) {
			bok = 0;
			if (mxIsNaN (a[js + iw])) {
				for (i = iw+1; i <= ie; i++)
					if (!mxIsNaN (a[js + i])) bok++;
			}
			else {
				for (i = iw+1; i <= ie; i++)
					if (a[js + i] != a[js + iw]) bok++;
			}
			if (bok > 0) mexPrintf ("GRDTRACK_M Warning: Inconsistent grid values at South pole.\n");
		}
	}

	/* Start with the case that x is not periodic, because in that
		case we also know that y cannot be polar.  */

	if (edgeinfo->nxp <= 0) {

		/* x is not periodic  */

		if (edgeinfo->nyp > 0) {

			/* y is periodic  */

			for (i = iw; i <= ie; i++) {
				a[jno1 + i] = a[jno1k + i];
				a[jno2 + i] = a[jno2k + i];
				a[jso1 + i] = a[jso1k + i];
				a[jso2 + i] = a[jso2k + i];
			}

			/* periodic Y rows copied.  Now do X naturals.
				This is easy since y's are done; no corner problems.
				Begin with Laplacian = 0, and include 1st outside rows
				in loop, since y's already loaded to 2nd outside.  */

			for (jmx = jno1; jmx <= jso1; jmx += mx) {
				a[jmx + iwo1] = (float)(4.0 * a[jmx + iw])
					- (a[jmx + iw + mx] + a[jmx + iw - mx] + a[jmx + iwi1]);
				a[jmx + ieo1] = (float)(4.0 * a[jmx + ie])
					- (a[jmx + ie + mx] + a[jmx + ie - mx] + a[jmx + iei1]);
			}

			/* Copy that result to 2nd outside row using periodicity.  */
			a[jno2 + iwo1] = a[jno2k + iwo1];
			a[jso2 + iwo1] = a[jso2k + iwo1];
			a[jno2 + ieo1] = a[jno2k + ieo1];
			a[jso2 + ieo1] = a[jso2k + ieo1];

			/* Now set d[laplacian]/dx = 0 on 2nd outside column.  Include
				1st outside rows in loop.  */
			for (jmx = jno1; jmx <= jso1; jmx += mx) {
				a[jmx + iwo2] = (a[jmx + iw - mx] + a[jmx + iw + mx] + a[jmx + iwi1])
					- (a[jmx + iwo1 - mx] + a[jmx + iwo1 + mx])
					+ (float)(5.0 * (a[jmx + iwo1] - a[jmx + iw]));

				a[jmx + ieo2] = (a[jmx + ie - mx] + a[jmx + ie + mx] + a[jmx + iei1])
					- (a[jmx + ieo1 - mx] + a[jmx + ieo1 + mx])
					+ (float)(5.0 * (a[jmx + ieo1] - a[jmx + ie]));
			}

			/* Now copy that result also, for complete periodicity's sake  */
			a[jno2 + iwo2] = a[jno2k + iwo2];
			a[jso2 + iwo2] = a[jso2k + iwo2];
			a[jno2 + ieo2] = a[jno2k + ieo2];
			a[jso2 + ieo2] = a[jso2k + ieo2];

			/* DONE with X not periodic, Y periodic case.  Fully loaded.  */

			return (0);
		}
		else {
			/* Here begins the X not periodic, Y not periodic case  */

			/* First, set corner points.  Need not merely Laplacian(f) = 0
				but explicitly that d2f/dx2 = 0 and d2f/dy2 = 0.
				Also set d2f/dxdy = 0.  Then can set remaining points.  */

	/* d2/dx2 */	a[jn + iwo1]   = (float)(2.0 * a[jn + iw] - a[jn + iwi1]);
	/* d2/dy2 */	a[jno1 + iw]   = (float)(2.0 * a[jn + iw] - a[jni1 + iw]);
	/* d2/dxdy */	a[jno1 + iwo1] = -(a[jno1 + iwi1] - a[jni1 + iwi1] + a[jni1 + iwo1]);


	/* d2/dx2 */	a[jn + ieo1]   = (float)(2.0 * a[jn + ie] - a[jn + iei1]);
	/* d2/dy2 */	a[jno1 + ie]   = (float)(2.0 * a[jn + ie] - a[jni1 + ie]);
	/* d2/dxdy */	a[jno1 + ieo1] = -(a[jno1 + iei1] - a[jni1 + iei1] + a[jni1 + ieo1]);

	/* d2/dx2 */	a[js + iwo1]   = (float)(2.0 * a[js + iw] - a[js + iwi1]);
	/* d2/dy2 */	a[jso1 + iw]   = (float)(2.0 * a[js + iw] - a[jsi1 + iw]);
	/* d2/dxdy */	a[jso1 + iwo1] = -(a[jso1 + iwi1] - a[jsi1 + iwi1] + a[jsi1 + iwo1]);

	/* d2/dx2 */	a[js + ieo1]   = (float)(2.0 * a[js + ie] - a[js + iei1]);
	/* d2/dy2 */	a[jso1 + ie]   = (float)(2.0 * a[js + ie] - a[jsi1 + ie]);
	/* d2/dxdy */	a[jso1 + ieo1] = -(a[jso1 + iei1] - a[jsi1 + iei1] + a[jsi1 + ieo1]);

			/* Now set Laplacian = 0 on interior edge points,
				skipping corners:  */
			for (i = iwi1; i <= iei1; i++) {
				a[jno1 + i] = (float)(4.0 * a[jn + i])
					- (a[jn + i - 1] + a[jn + i + 1] + a[jni1 + i]);

				a[jso1 + i] = (float)(4.0 * a[js + i])
					- (a[js + i - 1] + a[js + i + 1] + a[jsi1 + i]);
			}
			for (jmx = jni1; jmx <= jsi1; jmx += mx) {
				a[iwo1 + jmx] = (float)(4.0 * a[iw + jmx])
					- (a[iw + jmx + mx] + a[iw + jmx - mx] + a[iwi1 + jmx]);
				a[ieo1 + jmx] = (float)(4.0 * a[ie + jmx])
					- (a[ie + jmx + mx] + a[ie + jmx - mx] + a[iei1 + jmx]);
			}

			/* Now set d[Laplacian]/dn = 0 on all edge pts, including
				corners, since the points needed in this are now set.  */
			for (i = iw; i <= ie; i++) {
				a[jno2 + i] = a[jni1 + i]
					+ (float)(5.0 * (a[jno1 + i] - a[jn + i]))
					+ (a[jn + i - 1] - a[jno1 + i - 1])
					+ (a[jn + i + 1] - a[jno1 + i + 1]);
				a[jso2 + i] = a[jsi1 + i]
					+ (float)(5.0 * (a[jso1 + i] - a[js + i]))
					+ (a[js + i - 1] - a[jso1 + i - 1])
					+ (a[js + i + 1] - a[jso1 + i + 1]);
			}
			for (jmx = jn; jmx <= js; jmx += mx) {
				a[iwo2 + jmx] = a[iwi1 + jmx]
					+ (float)(5.0 * (a[iwo1 + jmx] - a[iw + jmx]))
					+ (a[iw + jmx - mx] - a[iwo1 + jmx - mx])
					+ (a[iw + jmx + mx] - a[iwo1 + jmx + mx]);
				a[ieo2 + jmx] = a[iei1 + jmx]
					+ (float)(5.0 * (a[ieo1 + jmx] - a[ie + jmx]))
					+ (a[ie + jmx - mx] - a[ieo1 + jmx - mx])
					+ (a[ie + jmx + mx] - a[ieo1 + jmx + mx]);
			}
			/* DONE with X not periodic, Y not periodic case.
				Loaded all but three cornermost points at each corner.  */

			return (0);
		}
		/* DONE with all X not periodic cases  */
	}
	else {
		/* X is periodic.  Load x cols first, then do Y cases.  */

		for (jmx = jn; jmx <= js; jmx += mx) {
			a[iwo1 + jmx] = a[iwo1k + jmx];
			a[iwo2 + jmx] = a[iwo2k + jmx];
			a[ieo1 + jmx] = a[ieo1k + jmx];
			a[ieo2 + jmx] = a[ieo2k + jmx];
		}

		if (edgeinfo->nyp > 0) {
			/* Y is periodic.  copy all, including boundary cols:  */
			for (i = iwo2; i <= ieo2; i++) {
				a[jno1 + i] = a[jno1k + i];
				a[jno2 + i] = a[jno2k + i];
				a[jso1 + i] = a[jso1k + i];
				a[jso2 + i] = a[jso2k + i];
			}
			/* DONE with X and Y both periodic.  Fully loaded.  */

			return (0);
		}

		/* Do north (top) boundary:  */

		if (edgeinfo->gn) {
			/* Y is at north pole.  Phase-shift all, incl. bndry cols. */
			if (h->node_offset) {
				j1p = jn;	/* constraint for jno1  */
				j2p = jni1;	/* constraint for jno2  */
			}
			else {
				j1p = jni1;	/* constraint for jno1  */
				j2p = jni1 + mx;	/* constraint for jno2  */
			}
			for (i = iwo2; i <= ieo2; i++) {
				i180 = pad[0] + ((i + nxp2)%edgeinfo->nxp);
				a[jno1 + i] = a[j1p + i180];
				a[jno2 + i] = a[j2p + i180];
			}
		}
		else {
			/* Y needs natural conditions.  x bndry cols periodic.
				First do Laplacian.  Start/end loop 1 col outside,
				then use periodicity to set 2nd col outside.  */

			for (i = iwo1; i <= ieo1; i++) {
				a[jno1 + i] = (float)(4.0 * a[jn + i])
					- (a[jn + i - 1] + a[jn + i + 1] + a[jni1 + i]);
			}
			a[jno1 + iwo2] = a[jno1 + iwo2 + edgeinfo->nxp];
			a[jno1 + ieo2] = a[jno1 + ieo2 - edgeinfo->nxp];


			/* Now set d[Laplacian]/dn = 0, start/end loop 1 col out,
				use periodicity to set 2nd out col after loop.  */

			for (i = iwo1; i <= ieo1; i++) {
				a[jno2 + i] = a[jni1 + i]
					+ (float)(5.0 * (a[jno1 + i] - a[jn + i]))
					+ (a[jn + i - 1] - a[jno1 + i - 1])
					+ (a[jn + i + 1] - a[jno1 + i + 1]);
			}
			a[jno2 + iwo2] = a[jno2 + iwo2 + edgeinfo->nxp];
			a[jno2 + ieo2] = a[jno2 + ieo2 - edgeinfo->nxp];

			/* End of X is periodic, north (top) is Natural.  */

		}

		/* Done with north (top) BC in X is periodic case.  Do south (bottom)  */

		if (edgeinfo->gs) {
			/* Y is at south pole.  Phase-shift all, incl. bndry cols. */
			if (h->node_offset) {
				j1p = js;	/* constraint for jso1  */
				j2p = jsi1;	/* constraint for jso2  */
			}
			else {
				j1p = jsi1;	/* constraint for jso1  */
				j2p = jsi1 - mx;	/* constraint for jso2  */
			}
			for (i = iwo2; i <= ieo2; i++) {
				i180 = pad[0] + ((i + nxp2)%edgeinfo->nxp);
				a[jso1 + i] = a[j1p + i180];
				a[jso2 + i] = a[j2p + i180];
			}
		}
		else {
			/* Y needs natural conditions.  x bndry cols periodic.
				First do Laplacian.  Start/end loop 1 col outside,
				then use periodicity to set 2nd col outside.  */

			for (i = iwo1; i <= ieo1; i++) {
				a[jso1 + i] = (float)(4.0 * a[js + i])
					- (a[js + i - 1] + a[js + i + 1] + a[jsi1 + i]);
			}
			a[jso1 + iwo2] = a[jso1 + iwo2 + edgeinfo->nxp];
			a[jso1 + ieo2] = a[jso1 + ieo2 - edgeinfo->nxp];


			/* Now set d[Laplacian]/dn = 0, start/end loop 1 col out,
				use periodicity to set 2nd out col after loop.  */

			for (i = iwo1; i <= ieo1; i++) {
				a[jso2 + i] = a[jsi1 + i]
					+ (float)(5.0 * (a[jso1 + i] - a[js + i]))
					+ (a[js + i - 1] - a[jso1 + i - 1])
					+ (a[js + i + 1] - a[jso1 + i + 1]);
			}
			a[jso2 + iwo2] = a[jso2 + iwo2 + edgeinfo->nxp];
			a[jso2 + ieo2] = a[jso2 + ieo2 - edgeinfo->nxp];

			/* End of X is periodic, south (bottom) is Natural.  */

		}

		/* Done with X is periodic cases.  */

		return (0);
	}
}


/* -------------------------------------------------------------------- */
int decode_R (char *item, double *w, double *e, double *s, double *n) {
	char *text, string[BUFSIZ];
	
	/* Minimalist code to decode option -R extracted from GMT_get_common_args */
	
	int i, error = 0;
	double *p[4];
	
	p[0] = w;	p[1] = e;	p[2] = s;	p[3] = n;
			
	i = 0;
	strcpy (string, &item[2]);
	text = strtok (string, "/");
	while (text) {
		*p[i] = ddmmss_to_degree (text);
		i++;
		text = strtok (CNULL, "/");
	}
	if (item[strlen(item)-1] == 'r')	/* Rectangular box given, but valid here */
		error++;
	if (i != 4 || check_region (*p[0], *p[1], *p[2], *p[3]))
		error++;
	w = p[0];	e = p[1];
	s = p[2];	n = p[3];
	return (error);
}

/* -------------------------------------------------------------------- */
int check_region (double w, double e, double s, double n) {
	/* If region is given then we must have w < e and s < n */
	return ((w >= e || s >= n));
}

/* -------------------------------------------------------------------- */
double ddmmss_to_degree (char *text) {
	int i, colons = 0, suffix;
	double degree, minute, degfrac, second;

	for (i = 0; text[i]; i++) if (text[i] == ':') colons++;
	suffix = (int)text[i-1];	/* Last character in string */
	if (colons == 2) {	/* dd:mm:ss format */
		sscanf (text, "%lf:%lf:%lf", &degree, &minute, &second);
		degfrac = degree + Loc_copysign (minute / 60.0 + second / 3600.0, degree);
	}
	else if (colons == 1) {	/* dd:mm format */
		sscanf (text, "%lf:%lf", &degree, &minute);
		degfrac = degree + Loc_copysign (minute / 60.0, degree);
	}
	else
		degfrac = atof (text);
	if (suffix == 'W' || suffix == 'w' || suffix == 'S' || suffix == 's') degfrac = -degfrac;	/* Sign was given implicitly */
	return (degfrac);
}
