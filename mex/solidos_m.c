/*--------------------------------------------------------------------
 *	$Id: solidos_m.c 3527 2012-04-25 15:24:25Z j $
 *
 *	Copyright (c) 2004-2012 by J. Luis
 *
 * 	This program is part of Mirone and is free software; you can redistribute
 * 	it and/or modify it under the terms of the GNU Lesser General Public
 * 	License as published by the Free Software Foundation; either
 * 	version 2.1 of the License, or any later version.
 * 
 * 	This program is distributed in the hope that it will be useful,
 * 	but WITHOUT ANY WARRANTY; without even the implied warranty of
 * 	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * 	Lesser General Public License for more details.
 *
 *	Contact info: w3.ualg.pt/~jluis/mirone
 *--------------------------------------------------------------------*/

/*
 * 
 * 	constroi corpos solidos
 * 
 * Author:	J. Luis
 * Date:	28 Outubro, 2000
 *		09-12-00 acrecentei um corpo gaussiano 
 * Revision	03-02-01 Bug na opcao piramide
 */

#include "mex.h"
#include "mwsize.h"
#include <float.h>
#include <string.h>

#define	CHUNK	2000
#define	FALSE	0
#define	TRUE	1
#define M_PI	3.14159265358979323846
#define D2R (M_PI / 180.0)
#define DPI (2 * M_PI)
#define PI2 (M_PI /2.)
#define R2D (180.0 / M_PI)
//#define ln(x)  (((x) <= 0.0) ? 0.0 : (double) log((double)x))

struct CIRC {
	double  x;
	double  y;
} *circ;

struct TRIANG {
	double  x[3];
	double  y[3];
	double  z[3];
} *tri;

int cilindro (double rad_c, double height_c, double z_c, double x0, double y0, int n_pts);
int paralelo (double a, double b, double c, double z_c, double x0, double y0);
int five_psoid (double a, double b, double c, double z_c, int n_pts, int n_slice, double x0, double y0, int cone, int piram, int sino, int hemi);
void euler (double phi, double teta, double psi, int n_pts, double x0, double y0, double z_c);
double	gaussian (double rad, double half_width);
double area (int n_pts);

float	n_sig = 2;	/* Number of sigmas which will determine the bell's base width */

/* --------------------------------------------------------------------------- */
/* Matlab Gateway routine */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	int     error = FALSE, cilinder = FALSE, rotate = FALSE, multi = FALSE;
	int     cone = FALSE, piram = FALSE, ellipsoid = FALSE, sphere = FALSE;
	int     hemi = FALSE, verbose = FALSE, paralelogram = FALSE;
	int     sino = FALSE;
	int     i, j, j1, n_body = 0;
	int     n_pts = 24;	/* divide a circle in 15 degrees intervals */
	int     n_slice = 5;	/* some bodies are constructed as pile of n_slice */
	int	    argc = 0, n_arg_no_char = 0;
	char	**argv, format[512], *out_file;
	double	x0 = 0, y0 = 0, rad_c, z_c, height_c, azim, inc, dec;
	double	a, b, c, h_bell, sigma_x, sigma_y, *ptr_d;

	azim = inc = dec = 0.;
 
	argc = nrhs;
	for (i = 0; i < nrhs; i++) {		/* Check input to find how many arguments are of type char */
		if (!mxIsChar(prhs[i])) {
			argc--;
			n_arg_no_char++;	/* Number of arguments that have a type other than char */
		}
	}
	argc++;			/* to account for the program's name to be inserted in argv[0] */

	/* get the length of the input string */
	argv = (char **)mxCalloc(argc, sizeof(char *));
	argv[0] = "solidos";
	for (i = 1; i < argc; i++)
		argv[i] = (char *)mxArrayToString(prhs[i+n_arg_no_char-1]);

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {
 		
				/* Common parameters */

				case 'B':
					j = sscanf(&argv[i][2], "%lf/%lf/%lf/%lf/%f/%d/%d", &h_bell, &sigma_x, &sigma_y, &z_c, &n_sig, &n_pts, &n_slice);
					if (j < 4) {
						mexPrintf ("%s: SYNTAX ERROR -G option: Wrong number of arguments %d\n", argv[0], j);
						error = TRUE;
					}	
					sino = TRUE;
					n_body++;
					break;
				case 'C':
					j = sscanf(&argv[i][2], "%lf/%lf/%lf/%d", &rad_c, &height_c, &z_c, &n_pts);
					if (j < 3) {
						mexPrintf ("%s: SYNTAX ERROR -C option: Wrong number of arguments %d\n", argv[0], j);
						error = TRUE;
					}	
					cilinder = TRUE;
					n_body++;
					break;
				case 'c':
					j = sscanf(&argv[i][2], "%lf/%lf/%lf/%lf/%d", &a, &b, &c, &z_c, &n_pts);
					if (j < 3) {
						mexPrintf ("%s: SYNTAX ERROR -c option: Wrong number of arguments %d\n", argv[0], j);
						error = TRUE;
					}	
					cone = TRUE;
					n_body++;
					break;
				case 'R':
					sscanf(&argv[i][2], "%lf/%lf/%lf", &azim, &inc, &dec);
					rotate = TRUE;
					break;
				case 'E':
					j = sscanf(&argv[i][2], "%lf/%lf/%lf/%lf/%d/%d", &a, &b, &c, &z_c, &n_pts, &n_slice);
					if (j < 4) {
						mexPrintf ("%s: SYNTAX ERROR -E option: Wrong number of arguments %d\n", argv[0], j);
						error = TRUE;
					}	
					ellipsoid = TRUE;
					n_body++;
					break;
				case 'F':
					out_file = &argv[i][2];
					break;
				case 'H':
					hemi = TRUE;
					break;
				case 'M':
					multi = TRUE;
					break;
				case 'P':
					j = sscanf(&argv[i][2], "%lf/%lf/%lf/%lf", &a, &b, &c, &z_c);
					if (j != 4) {
						mexPrintf ("%s: SYNTAX ERROR -P option: Wrong number of arguments %d\n", argv[0], j);
						error = TRUE;
					}	
					paralelogram = TRUE;
					n_body++;
					break;
				case 'p':
					j = sscanf(&argv[i][2], "%lf/%lf/%lf/%lf", &a, &b, &c, &z_c);
					if (j != 4) {
						mexPrintf ("%s: SYNTAX ERROR -p option: Wrong number of arguments %d\n", argv[0], j);
						error = TRUE;
					}	
					a /= 2.0;	b /= 2.0;
					n_pts = 4;
					piram = TRUE;
					n_body++;
					break;
				case 'S':
					j = sscanf(&argv[i][2], "%lf/%lf/%d/%d", &a, &z_c, &n_pts, &n_slice);
					if (j < 2) {
						mexPrintf ("%s: SYNTAX ERROR -S option: Wrong number of arguments %d\n", argv[0], j);
						error = TRUE;
					}	
					sphere = TRUE;
					n_body++;
					break;
				case 'V':
					verbose = TRUE;
					break;
				case 'X':
					sscanf(&argv[i][2], "%lf", &x0);
					break;
				case 'Y':
					sscanf(&argv[i][2], "%lf", &y0);
					break;
				default:
					error = TRUE;
					break;
			}
		}
	}

	if (nlhs == 0) {	/* Display usage */
 		mexPrintf("solidos - calculates a solid surface as a union of triangular facets.\n\n");
		mexPrintf("usage:  solidos [-E<semi_x/semi_y/semi_z/z_center[/npts/n_slice]>]\n");
		mexPrintf("  [-S<rad/z_center[/npts/n_slice]>] [-B<height/sx/sy/z0[/n_sig/npts/n_slice]>]\n");
		mexPrintf("  [-C<rad/height/z0[/npts]>] [-c<semi_x/semi_y/height/z0[/npts]>]\n");
		mexPrintf("  [-P<side_x/side_y/side_z/z0>] [-R<azim/inc/dec>]\n");
		mexPrintf("  [-p<side_x/side_y/height/z0>] [-M] [-H] [-X<x0>] [-Y<y0>]\n\n");
 		
		mexPrintf("\tBODIES:\n");
		mexPrintf("\t-B bell (gaussian) of height <height> with caracteristic standard\n");
		mexPrintf("\t   deviation <sx> and <sy>. The base width (located at depth <z0>) is\n");
		mexPrintf("\t   controled by the number of sigma-limits (<n_sig>) [Default = 2.0]\n");
		mexPrintf("\t-C cilinder of radius <rad> height <height> and base at depth <z0>\n");
		mexPrintf("\t-c cone of semi axes <semi_x/semi_y> (that means it can be elliptic)\n");
		mexPrintf("\t   height <height> and base at depth <z0>\n");
		mexPrintf("\t-E ellipsoid of semi axes <semi_x/semi_y/semi_z> and center\n");
		mexPrintf("\t   depth <z_center>\n");
		mexPrintf("\t-P paralelogram of sides <x/y/z> and base at depth <z0>\n");
		mexPrintf("\t-p piramid of sides <x/y> height <height> and base at depth <z0>\n");
		mexPrintf("\t-S sphere of radius <rad> and center at depth <z_center>\n");
		mexPrintf("\n\tCOMMON MEANING PARAMETERS:\n");
		mexPrintf("\t npts     -> number of points in which the circle is descretized\n");
		mexPrintf("\t             [Default = 24].\n");
		mexPrintf("\t n_slice  -> some bodies are constructed by a pile of slices.\n");
		mexPrintf("\t             Spheres and Ellipsoides are made by 2*n_slice.\n");
		mexPrintf("\t             Bells are made by n_slice. [Default = 5]\n");
		mexPrintf("\t z_center -> body's center of mass depth (applies to sphere\n");
		mexPrintf("\t             and ellipsoide) [positive up].\n");
		mexPrintf("\t z0       -> depth of body's base (applies to remanent bodies)\n");
		mexPrintf("\n\tOPTIONS:\n");
		mexPrintf("\t-M output triangles as multiple segments separated by > flag\n");
		mexPrintf("\t-R rotate body by azimuth inclination and declination\n");
		mexPrintf("\t   (follows the convention of Euler angles)\n");
		mexPrintf("\t-X -Y to shift origin of the body by <x0> and/or <y0>\n");
		return;
	}
	if (error) mexErrMsgTxt("");
	
	if (!n_body)
		mexErrMsgTxt ("SOLIDOS: Error. Must select one of -C -E -G -S -c -P options\n");

	else if (n_body > 1)
		mexErrMsgTxt ("SOLIDOS: Error. Choose only one of -C -E -G -S -c -P options\n");

	if (cilinder)
		n_pts = cilindro (rad_c, height_c, z_c, x0, y0, n_pts);
	else if (ellipsoid)
		n_pts = five_psoid (a, b, c, z_c, n_pts, n_slice, x0, y0, FALSE, FALSE, FALSE, hemi);
	else if (sphere)
		n_pts = five_psoid (a, a, a, z_c, n_pts, n_slice, x0, y0, FALSE, FALSE, FALSE, hemi);
	else if (cone)
		n_pts = five_psoid (a, b, c, z_c, n_pts, 1, x0, y0, cone, FALSE, FALSE, hemi);
	else if (piram)
		n_pts = five_psoid (a, b, c, z_c, n_pts, 1, x0, y0, FALSE, piram, FALSE, hemi);
	else if (sino)
		n_pts = five_psoid (sigma_x, sigma_y, h_bell, z_c, n_pts, n_slice, x0, y0, FALSE, FALSE, sino, hemi);
	else if (paralelogram)
		n_pts = paralelo (a, b, c, z_c, x0, y0);

	if (rotate)
		euler (azim, inc, dec, n_pts, x0, y0, z_c);

	if (!multi)
		sprintf (format,"%s\n\0", "%.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f");
	else
		sprintf (format,"%s\n\0", ">\n%.3f %.3f %.3f\n%.3f %.3f %.3f\n%.3f %.3f %.3f");

	if (nlhs == 0) {
		FILE *fp;
		fp = fopen(out_file, "wt");
		for (i = 0; i < n_pts; i++) {
			fprintf (fp, format, tri[i].x[0], tri[i].y[0], tri[i].z[0], tri[i].x[1], tri[i].y[1],
					tri[i].z[1], tri[i].x[2], tri[i].y[2], tri[i].z[2]);
		}
		fclose(fp);
	}
	else {
		plhs[0] = mxCreateDoubleMatrix(n_pts, 9, mxREAL);
		ptr_d = mxGetPr(plhs[0]);
		for (i = 0; i < n_pts; i++) {
			ptr_d[i          ] = tri[i].x[0];
			ptr_d[i +   n_pts] = tri[i].y[0];
			ptr_d[i + 2*n_pts] = tri[i].z[0];
			ptr_d[i + 3*n_pts] = tri[i].x[1];
			ptr_d[i + 4*n_pts] = tri[i].y[1];
			ptr_d[i + 5*n_pts] = tri[i].z[1];
			ptr_d[i + 6*n_pts] = tri[i].x[2];
			ptr_d[i + 7*n_pts] = tri[i].y[2];
			ptr_d[i + 8*n_pts] = tri[i].z[2];
		}
	}

	if (verbose) {
		mexPrintf ("%d triangles\n",n_pts);
		mexPrintf ("body area = %g\n", area (n_pts));
	}
	mxFree(tri);
}

int five_psoid (double a, double b, double c, double z_c, int n_pts, int n_slice, double x0, double y0, int cone, int piram, int sino, int hemi) {
	/*	Constructs either a sphere, ellipsoid, cone, piramid, or a bell
	/*	as a union of triagular facets. Returns number of triangles. */

	double z_top, z_bot, dfi, d_sli, ai0, ai1, bi0, bi1, ci, zi0, zi1;
	double d_tet, half_width_x, half_width_y, dx, dy, dz, rad_x, rad_y;
	int i, j, j1, k, l, m, m1, m2, n = 0, n_tri, first = TRUE;
	struct CIRC *ellipse[2];

	n_tri = (hemi) ? 2 * n_pts * n_slice: 2 * (n_pts * (n_slice*2 - 1));
	ellipse[0] = (struct CIRC *) mxCalloc ((size_t)(n_pts+1), sizeof(struct CIRC));
	ellipse[1] = (struct CIRC *) mxCalloc ((size_t)(n_pts+1), sizeof(struct CIRC));
	if ((tri = (struct TRIANG *) mxCalloc ((size_t)(n_tri), sizeof(struct TRIANG)) ) == NULL) {
		mexPrintf ("Fatal Error: %s (five_psoid) could not allocate memory, n = %d\n", "solidos", n_tri);
		mexErrMsgTxt("");
	}	

	dfi = (DPI/n_pts);	d_tet = (PI2/n_slice);
	d_sli = c / n_slice;
	z_top = z_c + c;	z_bot = z_c;
	half_width_x = 0.5 * a;	half_width_y = 0.5 * b;
	dx = n_sig * a / n_slice;	dy = n_sig * b / n_slice; 
	dz = c / n_slice;	/* repeated but ok */

	for (j = 0; j < n_slice; j++) {
		j1 = j + 1;
		if (cone || piram) {
			ai0 = j * a / n_slice;	bi0 = j * b / n_slice;
			zi0 = z_top - j * d_sli;
			ai1 = j1 * a / n_slice;	bi1 = j1 * b / n_slice;
			zi1 = z_top - j1 * d_sli;
		}
		else if (sino) { /* Bell shaped volume */
			rad_x = j * dx;		rad_y = j * dy;
			ai0 = rad_x;		bi0 = rad_y;
			zi0 = z_bot + c * gaussian (rad_x, half_width_x); 
			rad_x = j1 * dx;	rad_y = j1 * dy;
			ai1 = rad_x;		bi1 = rad_y;
			zi1 = z_bot + c * gaussian (rad_x, half_width_x);

/* A fatia constante. Preciso de descobrir como controlar o tamanho da base
			rad_x = a*sqrt(2*log(c/(c-j*dz)));
			rad_y = b*sqrt(2*log(c/(c-j*dz)));
			ai0 = rad_x;		bi0 = rad_y;
			zi0 = z_bot + c - j*dz;
			rad_x = a*sqrt(2*log(c/(c-j1*dz)));
			rad_y = b*sqrt(2*log(c/(c-j1*dz)));
			ai1 = rad_x;		bi1 = rad_y;
			zi1 = z_bot + c - j1*dz;
			if (j1 == n_slice) {ai1 = ai0;	bi1 = bi0;} */
		}
		else { /* Ellipsoide or Sphere */
			ai0 = a*cos(PI2-j*d_tet);	bi0 = b*cos(PI2-j*d_tet);
			zi0 = z_top - c * (1. - sqrt(1. - (ai0/a)*(ai0/a)));
			ai1 = a*cos(PI2-j1*d_tet);	bi1 = b*cos(PI2-j1*d_tet);
			zi1 = z_top - c * (1. - sqrt(1. - (ai1/a)*(ai1/a)));
		}
		for (i = 0; i < n_pts; i++) { /* compute slice j */
			ellipse[0][i].x = x0 + ai0 * cos (i*dfi);
			ellipse[0][i].y = y0 + bi0 * sin (i*dfi);
			ellipse[1][i].x = x0 + ai1 * cos (i*dfi);
			ellipse[1][i].y = y0 + bi1 * sin (i*dfi);
		}
		ellipse[0][n_pts].x = ellipse[0][0].x; /* close slice "contour" */
		ellipse[0][n_pts].y = ellipse[0][0].y;
		ellipse[1][n_pts].x = ellipse[1][0].x;
		ellipse[1][n_pts].y = ellipse[1][0].y;

		/* Calculates vertex of triangles in slice j */
		i = 0;
		if (first) {
			for (m = 0; m < n_pts ; m++) {
				tri[m].x[0] = ellipse[0][i].x;
				tri[m].y[0] = ellipse[0][i].y;
				tri[m].z[0] = zi0;
				tri[m].x[1] = ellipse[1][i+1].x;
				tri[m].y[1] = ellipse[1][i+1].y;
				tri[m].z[1] = zi1;
				tri[m].x[2] = ellipse[1][i].x;
				tri[m].y[2] = ellipse[1][i].y;
				tri[m].z[2] = zi1;
				i++;
			}
		}
		else {
			for (m = (j-1)*n_pts; m < j*n_pts ; m++) {
				/* First triangle */
				m1 = 2 * m + n_pts;	m2 = 2 * m + 1 + n_pts;
				tri[m1].x[0] = ellipse[0][i].x;
				tri[m1].y[0] = ellipse[0][i].y;
				tri[m1].z[0] = zi0;
				tri[m1].x[1] = ellipse[1][i+1].x;
				tri[m1].y[1] = ellipse[1][i+1].y;
				tri[m1].z[1] = zi1;
				tri[m1].x[2] = ellipse[1][i].x;
				tri[m1].y[2] = ellipse[1][i].y;
				tri[m1].z[2] = zi1;
				/* Second triangle */
				tri[m2].x[0] = ellipse[0][i].x;
				tri[m2].y[0] = ellipse[0][i].y;
				tri[m2].z[0] = zi0;
				tri[m2].x[1] = ellipse[0][i+1].x;
				tri[m2].y[1] = ellipse[0][i+1].y;
				tri[m2].z[1] = zi0;
				tri[m2].x[2] = ellipse[1][i+1].x;
				tri[m2].y[2] = ellipse[1][i+1].y;
				tri[m2].z[2] = zi1;
				i++;
			}
		}
		first = FALSE;
	}
	/* First half is ready. Now, either close it and return or 
	construct the other half by simetry */
	if (cone || piram || sino || hemi) { /* close the base and return */
		if (sino && fabs ((zi1 - z_c) / z_c) > 0.01) { 
			/* bell's last slice is 1% far from base, so we add a vertical wall */ 
			z_top = zi1;	/* update z_top value */
			for (n = 0; n < n_pts; n++) {
				/* First triangle */
				j = n_pts * (n_slice * 2 - 1) + 2 * n;
				tri[j].x[0] = ellipse[1][n].x;
				tri[j].y[0] = ellipse[1][n].y;
				tri[j].z[0] = z_top;
				tri[j].x[1] = ellipse[1][n+1].x;
				tri[j].y[1] = ellipse[1][n+1].y;
				tri[j].z[1] = z_top;
				tri[j].x[2] = ellipse[1][n].x;
				tri[j].y[2] = ellipse[1][n].y;
				tri[j].z[2] = z_bot;
				/* Second triangle */
				j1 = n_pts * (n_slice * 2 - 1) + 2 * n + 1;
				tri[j1].x[0] = ellipse[1][n+1].x;
				tri[j1].y[0] = ellipse[1][n+1].y;
				tri[j1].z[0] = z_top;
				tri[j1].x[1] = ellipse[1][n+1].x;
				tri[j1].y[1] = ellipse[1][n+1].y;
				tri[j1].z[1] = z_bot;
				tri[j1].x[2] = ellipse[1][n].x;
				tri[j1].y[2] = ellipse[1][n].y;
				tri[j1].z[2] = z_bot;
			}
		}
		else if (sino) /* slightly change base depth to force a closed volume */
			z_c = zi1;

		i = n_pts * (n_slice*2 - 1) + 2 * n;
		for (k = i, l = 0; k < i+n_pts; k++, l++) {
			tri[k].x[0] = x0;	tri[k].y[0] = y0;
			tri[k].z[0] = z_c;
			tri[k].x[1] = ellipse[1][l].x;
			tri[k].y[1] = ellipse[1][l].y;
			tri[k].z[1] = z_c;
			tri[k].x[2] = ellipse[1][l+1].x;
			tri[k].y[2] = ellipse[1][l+1].y;
			tri[k].z[2] = z_c;
		}
		free ((void *)ellipse[0]);
		free ((void *)ellipse[1]);
//		return (n_pts*n_slice*2 + n);
		return (k);
	}
	n_tri = n_pts * (n_slice*2 - 1);
	for (j = n_tri-1, i = n_tri; j >= 0; j--, i++) {
		tri[i].x[0] = tri[j].x[0];	tri[i].y[0] = tri[j].y[0];
		tri[i].z[0] = z_c + (z_c - tri[j].z[0]);
		tri[i].x[1] = tri[j].x[2];	tri[i].y[1] = tri[j].y[2];
		tri[i].z[1] = z_c + (z_c - tri[j].z[2]);
		tri[i].x[2] = tri[j].x[1];	tri[i].y[2] = tri[j].y[1];
		tri[i].z[2] = z_c + (z_c - tri[j].z[1]);
	}
	mxFree(ellipse[0]);
	mxFree(ellipse[1]);
	n_tri = 2 * (n_pts * (n_slice*2 - 1));
	return (n_tri);
}

double	gaussian (double rad, double half_width) {
	double y, gauss_ct;
	gauss_ct = -0.5 / (half_width * half_width);
	y = exp (rad * rad * gauss_ct);
	return (y);
}

int cilindro (double rad_c, double height_c, double z_c, double x0, double y0, int n_pts) {
	double z_top, z_bot, dfi;
	int i, j, j1, n_tri;

	n_tri = n_pts * 4;
	circ = (struct CIRC *) mxCalloc ((size_t)(n_pts+1), sizeof(struct CIRC));
	tri = (struct TRIANG *) mxCalloc ((size_t)(n_tri), sizeof(struct TRIANG));
	dfi = (DPI/n_pts);
	z_top = z_c + height_c;		z_bot = z_c ;

	for (i = 0; i < n_pts; i++) { /* compute circle */
		circ[i].x = x0 + rad_c * cos (i*dfi);
		circ[i].y = y0 + rad_c * sin (i*dfi);
	}
	circ[n_pts].x = circ[0].x;	circ[n_pts].y = circ[0].y;

	for (i = 0; i < n_pts; i++) { /* Calculates vertex of top circle */
		tri[i].x[0] = x0;	tri[i].y[0] = y0;	tri[i].z[0] = z_top;
		tri[i].x[1] = circ[i+1].x;	tri[i].y[1] = circ[i+1].y;
		tri[i].z[1] = z_top;
		tri[i].x[2] = circ[i].x;	tri[i].y[2] = circ[i].y;
		tri[i].z[2] = z_top;
	}

	for (i = 0; i < n_pts; i++) { /* Calculates vertex of side rectangle, where each one is decomposed into two triangles */
		/* First triangle */
		j = 2 * i + n_pts;
		tri[j].x[0] = circ[i].x;	tri[j].y[0] = circ[i].y;
		tri[j].z[0] = z_top;
		tri[j].x[1] = circ[i+1].x;	tri[j].y[1] = circ[i+1].y;
		tri[j].z[1] = z_top;
		tri[j].x[2] = circ[i].x;	tri[j].y[2] = circ[i].y;
		tri[j].z[2] = z_bot;
		/* Second triangle */
		j1 = 2 * i + n_pts + 1;
		tri[j1].x[0] = circ[i+1].x;	tri[j1].y[0] = circ[i+1].y;
		tri[j1].z[0] = z_top;
		tri[j1].x[1] = circ[i+1].x;	tri[j1].y[1] = circ[i+1].y;
		tri[j1].z[1] = z_bot;
		tri[j1].x[2] = circ[i].x;	tri[j1].y[2] = circ[i].y;
		tri[j1].z[2] = z_bot;
	}

	for (i = 0; i < n_pts; i++) { /* Now the vertex of bottom circle */
		j = i + 3 * n_pts;
		tri[j].x[0] = x0;	tri[j].y[0] = y0;	tri[j].z[0] = z_bot;
		tri[j].x[1] = circ[i].x;	tri[j].y[1] = circ[i].y;
		tri[j].z[1] = z_bot;
		tri[j].x[2] = circ[i+1].x;	tri[j].y[2] = circ[i+1].y;
		tri[j].z[2] = z_bot;
	}
	mxFree (circ);
	return (n_tri);
}

int paralelo (double a, double b, double c, double z_c, double x0, double y0) {
	double z_top, z_bot;
	int i, n_tri = 12;

	tri = (struct TRIANG *) calloc ((size_t) (n_tri), sizeof(struct TRIANG));
	z_top = z_c + c;
	z_bot = z_c ;

	/* vertex of top rectangle */
		/* first triangle */
	tri[0].x[0] = - a/2. + x0;	tri[0].y[0] = - b/2. + y0;
	tri[0].z[0] = z_top;
	tri[0].x[1] = - a/2. + x0;	tri[0].y[1] = b/2. + y0;
	tri[0].z[1] = z_top;
	tri[0].x[2] = a/2. + x0;	tri[0].y[2] = b/2. + y0;
	tri[0].z[2] = z_top;
		/* second triangle */
	tri[1].x[0] = - a/2. + x0;	tri[1].y[0] = - b/2. + y0;
	tri[1].z[0] = z_top;
	tri[1].x[1] = a/2. + x0;	tri[1].y[1] = b/2. + y0;
	tri[1].z[1] = z_top;
	tri[1].x[2] = a/2. + x0;	tri[1].y[2] = - b/2. + y0;
	tri[1].z[2] = z_top;
	/* vertex of est rectangle */
		/* first triangle */
	tri[2].x[0] = a/2. + x0;	tri[2].y[0] = - b/2. + y0;
	tri[2].z[0] = z_bot;
	tri[2].x[1] = a/2. + x0;	tri[2].y[1] = - b/2. + y0;
	tri[2].z[1] = z_top;
	tri[2].x[2] = a/2. + x0;	tri[2].y[2] = b/2. + y0;
	tri[2].z[2] = z_top;
		/* second triangle */
	tri[3].x[0] = a/2. + x0;	tri[3].y[0] = - b/2. + y0;
	tri[3].z[0] = z_bot;
	tri[3].x[1] = a/2. + x0;	tri[3].y[1] = b/2. + y0;
	tri[3].z[1] = z_top;
	tri[3].x[2] = a/2. + x0;	tri[3].y[2] = b/2. + y0;
	tri[3].z[2] = z_bot;
	/* vertex of north rectangle */
		/* first triangle */
	tri[4].x[0] = a/2. + x0;	tri[4].y[0] = b/2. + y0;
	tri[4].z[0] = z_bot;
	tri[4].x[1] = a/2. + x0;	tri[4].y[1] = b/2. + y0;
	tri[4].z[1] = z_top;
	tri[4].x[2] = - a/2. + x0;	tri[4].y[2] = b/2. + y0;
	tri[4].z[2] = z_top;
		/* second triangle */
	tri[5].x[0] = a/2. + x0;	tri[5].y[0] = b/2. + y0;
	tri[5].z[0] = z_bot;
	tri[5].x[1] = - a/2. + x0;	tri[5].y[1] = b/2. + y0;
	tri[5].z[1] = z_top;
	tri[5].x[2] = - a/2. + x0;	tri[5].y[2] = b/2. + y0;
	tri[5].z[2] = z_bot;
	/* vertex of west rectangle */
		/* first triangle */
	tri[6].x[0] = - a/2. + x0;	tri[6].y[0] = b/2. + y0;
	tri[6].z[0] = z_bot;
	tri[6].x[1] = - a/2. + x0;	tri[6].y[1] = b/2. + y0;
	tri[6].z[1] = z_top;
	tri[6].x[2] = - a/2. + x0;	tri[6].y[2] = - b/2. + y0;
	tri[6].z[2] = z_top;
		/* second triangle */
	tri[7].x[0] = - a/2. + x0;	tri[7].y[0] = b/2. + y0;
	tri[7].z[0] = z_bot;
	tri[7].x[1] = - a/2. + x0;	tri[7].y[1] = - b/2. + y0;
	tri[7].z[1] = z_top;
	tri[7].x[2] = - a/2. + x0;	tri[7].y[2] = - b/2. + y0;
	tri[7].z[2] = z_bot;
	/* vertex of south rectangle */
		/* first triangle */
	tri[8].x[0] = - a/2. + x0;	tri[8].y[0] = - b/2. + y0;
	tri[8].z[0] = z_bot;
	tri[8].x[1] = - a/2. + x0;	tri[8].y[1] = - b/2. + y0;
	tri[8].z[1] = z_top;
	tri[8].x[2] = a/2. + x0;	tri[8].y[2] = - b/2. + y0;
	tri[8].z[2] = z_top;
		/* second triangle */
	tri[9].x[0] = - a/2. + x0;	tri[9].y[0] = - b/2. + y0;
	tri[9].z[0] = z_bot;
	tri[9].x[1] = a/2. + x0;	tri[9].y[1] = - b/2. + y0;
	tri[9].z[1] = z_top;
	tri[9].x[2] = a/2. + x0;	tri[9].y[2] = - b/2. + y0;
	tri[9].z[2] = z_bot;
	/* vertex of bottom rectangle */
		/* first triangle */
	tri[10].x[0] = - a/2. + x0;	tri[10].y[0] = - b/2. + y0;
	tri[10].z[0] = z_bot;
	tri[10].x[1] = a/2. + x0;	tri[10].y[1] = b/2. + y0;
	tri[10].z[1] = z_bot;
	tri[10].x[2] = - a/2. + x0;	tri[10].y[2] = b/2. + y0;
	tri[10].z[2] = z_bot;
		/* second triangle */
	tri[11].x[0] = - a/2. + x0;	tri[11].y[0] = - b/2. + y0;
	tri[11].z[0] = z_bot;
	tri[11].x[1] = a/2. + x0;	tri[11].y[1] = - b/2. + y0;
	tri[11].z[1] = z_bot;
	tri[11].x[2] = a/2. + x0;	tri[11].y[2] = b/2. + y0;
	tri[11].z[2] = z_bot;

	return (12);
}

void euler (double phi, double teta, double psi, int n_pts, double x0, double y0, double z_c) {
	int i;
	double c_tet, s_tet, c_phi, s_phi, c_psi, s_psi, a[3][3];
	double	x1, y1, z1, x2, y2, z2, x3, y3, z3;

	c_phi = cos (phi*D2R);	s_phi = sin(phi*D2R);
	c_tet = cos (teta*D2R);	s_tet = sin(teta*D2R);
	c_psi = cos (psi*D2R);	s_psi = sin(psi*D2R);
	a[0][0] = c_psi*c_phi - c_tet*s_phi*s_psi;
	a[0][1] = c_psi*s_phi + c_tet*c_phi*s_psi;
	a[0][2] = s_psi*s_tet;
	a[1][0] = -s_psi*c_phi - c_tet*s_phi*c_psi;
	a[1][1] = -s_psi*s_phi + c_tet*c_phi*c_psi;
	a[1][2] = c_psi*s_tet;
	a[2][0] = s_tet*s_phi;
	a[2][1] = -s_tet*c_phi;
	a[2][2] = c_tet;

	for (i = 0; i < n_pts; i++) {
		x1 = tri[i].x[0] - x0;	y1 = tri[i].y[0] - y0;
		z1 = tri[i].z[0] - z_c;
		x2 = tri[i].x[1] - x0;	y2 = tri[i].y[1] - y0;
		z2 = tri[i].z[1] - z_c;
		x3 = tri[i].x[2] - x0;	y3 = tri[i].y[2] - y0;
		z3 = tri[i].z[2] - z_c;
		tri[i].x[0] = a[0][0]*x1 + a[0][1]*y1 + a[0][2]*z1 + x0;
		tri[i].y[0] = a[1][0]*x1 + a[1][1]*y1 + a[1][2]*z1 + y0;
		tri[i].z[0] = a[2][0]*x1 + a[2][1]*y1 + a[2][2]*z1 + z_c;
		tri[i].x[1] = a[0][0]*x2 + a[0][1]*y2 + a[0][2]*z2 + x0; 
		tri[i].y[1] = a[1][0]*x2 + a[1][1]*y2 + a[1][2]*z2 + y0; 
		tri[i].z[1] = a[2][0]*x2 + a[2][1]*y2 + a[2][2]*z2 + z_c; 
		tri[i].x[2] = a[0][0]*x3 + a[0][1]*y3 + a[0][2]*z3 + x0; 
		tri[i].y[2] = a[1][0]*x3 + a[1][1]*y3 + a[1][2]*z3 + y0; 
		tri[i].z[2] = a[2][0]*x3 + a[2][1]*y3 + a[2][2]*z3 + z_c; 
	}
}

double area (int n_pts) {
	/* Compute area of body */
	int i;
	double S, A, B, C, sup, dx0, dx1, dx2, dy0, dy1, dy2, dz0, dz1, dz2, x;

	sup = 0.;
	for (i = 0; i < n_pts; i++) {
		dx0 = tri[i].x[1] - tri[i].x[0];
		dx1 = tri[i].x[2] - tri[i].x[1];
		dx2 = tri[i].x[0] - tri[i].x[2];
		dy0 = tri[i].y[1] - tri[i].y[0];
		dy1 = tri[i].y[2] - tri[i].y[1];
		dy2 = tri[i].y[0] - tri[i].y[2];
		dz0 = tri[i].z[1] - tri[i].z[0];
		dz1 = tri[i].z[2] - tri[i].z[1];
		dz2 = tri[i].z[0] - tri[i].z[2];
		A = sqrt (dx0*dx0 + dy0*dy0 + dz0*dz0);
		B = sqrt (dx1*dx1 + dy1*dy1 + dz1*dz1);
		C = sqrt (dx2*dx2 + dy2*dy2 + dz2*dz2);
		S = (A + B + C) / 2.;
		sup += sqrt ( S*(S-A)*(S-B)*(S-C) );
	}
	return (sup);
}
