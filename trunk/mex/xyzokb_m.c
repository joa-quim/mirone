/*--------------------------------------------------------------------
 *	$Id:$
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
 * Author:	Joaquim Luis
 * Date:	06-OCT-2000
 * Revised:	07-FEB-2001	
 *
 *	Versao com -S e onde a okabe foi transformada em double
 *	Comecei a ler um fiche em stl 14-12-00
 *
 *	Muita cosmetica e passagem hypot->sqrt 05-01-01
 *	Resolvida(?) a questao da aceleracao para prismas de espessura constante
 *
 *	Mudei DBL_EPSILON para FLT_EPSILON prque senao dava merda 04-02-01
 *	Tirei uma origem orlat orlon  11-02-01
 *
 *	Corrected potential bug if -z and -P	25-5-01
 *
 *	Retirei o orlat e orlon (era estupido)	27-5-01
 *	Added a test for checking (and swaping if needed) the triangle order.	27-5-01
 *
 *	06-02-03
 *	Changed the way computations in geographical coordinates were donne.
 *	Given that the okabe routine requires coordinates in meters, the geographical
 *	coordinates were transformed to meters along paralels (relative to grenwich) and
 *	meridians (relative to equator). This produces a change in the body azimuth as
 *	comparing to coodinates in a map projection (e.g.UTM).
 *	To reduce this effect, I now remove the body's central longitude. However, for bodies
 *	with large longitude span the problem might still persist.
 *
 *	05-10-03
 *	Ta uma grande mixordia com a historia dos z_dir e se "bat" ou "topo". Nao me oriento.
 *	Pus tudo a z positivo up, por isso as opcoes -D e -T nao devem ser usadas
 *
 *	07-11-03
 *	Na modif anterior so tinha retirado o orlon. Agora retiro tambem o orlat.
 *	Nao sei se faz diferenca
 *
 *	16-11-03
 *	Adicionei um teste para nao gastar tempo em calculos se a mag do triang for nula
 *
 *	17-02-04
 *	No seguimento da mixordice voltei a por o z0 *= z_dir. Isto faz a coisa mais
 *	coerente. Assim, na opcao -Z o <level> deve obedecer ao sinal. Ou seja, -Z-2300
 *	significa que a base está a -2300 m de prof.
 */

#include "mex.h"
#include "mwsize.h"
#include <float.h>
#include <string.h>

#ifdef WIN32	/* Start of Windows setup */
//#pragma warning( disable : 4244 )	/* SHUT UP THAT SANITY KILLER WARNING "... possible loss of data" */
#endif

#define	FALSE	0
#define	TRUE	1
#ifndef M_PI
#define M_PI	3.14159265358979323846
#endif

#define D2R (M_PI / 180.0)
#define cosd(x) cos ((x) * D2R)
#define CNULL	((char *)NULL)
#define Loc_copysign(x,y) ((y) < 0.0 ? -fabs(x) : fabs(x))

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

struct  DATA    {
        double  x;
        double  y;
}       *data;

struct  BODY_DESC {
	int n_f;
	int *n_v;
	int *ind;
}       bd_desc;

struct  LOC_OR    {
        double  x;
        double  y;
        double  z;
}       *loc_or;

struct  TRIANG    {
        double  x;
        double  y;
        double  z;
}       *triang;

struct  VERT    {
        int  a;
        int  b;
        int  c;
}       *vert;

struct  TRI_CENTER {
        double  x;
        double  y;
        double  z;
}       *t_center;

struct  RAW    {
        double  t1[3];
        double  t2[3];
        double  t3[3];
}       *raw_mesh;

struct MAG_PARAM {
	double	rim[3];
}	*mag_param;

struct MAG_VAR {		/* Used when only the modulus of magnetization varies */
	double	rk[3];
}	*mag_var;

struct MAG_VAR2 {
	double	m;
	double	m_dip;
}	*mag_var2;

struct MAG_VAR3 {
	double	m;
	double	m_dec;
	double	m_dip;
}	*mag_var3;

struct MAG_VAR4 {
	double	t_dec;
	double	t_dip;
	double	m;
	double	m_dec;
	double	m_dip;
}	*mag_var4;


/* Old habits. Need to get rid of these globals */
double	xx[24], yy[24], zz[24], d_to_m, *mag_int;
double	c_tet, s_tet, c_phi, s_phi, central_long, central_lat;
int grav = TRUE, mag = FALSE, const_th = FALSE;
int m_var1 = FALSE, m_var2 = FALSE, m_var3 = FALSE, m_var4 = FALSE;

int read_xyz (FILE *fp_xyz, char *xyz_file, double z_dir, int m_var, int geo);
int read_t (FILE *fp_t, char *t_file);
int read_stl (FILE *fp_s, char *stl_file, double z_dir);
void set_center (int n_triang);
double okabe (double rho, double x, double y, double z, struct BODY_DESC bd_desc, int km, int pm, struct LOC_OR *loc_or);
double okb_grv (int n_vert, struct LOC_OR *loc_or);
double okb_mag (int n_vert, int km, int pm, struct LOC_OR *loc_or); //
double eq_30 (double c, double s, double x, double y, double z);
double eq_43 (double mz, double c, double tg, double auxil, double x, double y, double z);
void rot_17 (int n_vert, int top, struct LOC_OR *loc_or);
int facet_triangulate (int i, double z0, double dz, int bat);
int facet_raw (int i, int geo);
int check_triang_cw (int n, int type);
void parse_FV(float *v, mwSize *f, int nFacet, int nv, double z_dir, struct RAW *raw_mesh);
int	no_sys_mem (char *where, int n);

int GMT_getinc (char *line, double *dx, double *dy);
int decode_R (char *item, double *w, double *e, double *s, double *n);
int check_region (double w, double e, double s, double n);
double ddmmss_to_degree (char *text);
int GMT_getinc (char *line, double *dx, double *dy);
int GMT_getincn (char *line, double inc[], int n);
int GMT_strtok (const char *string, const char *sep, int *pos, char *token);

int GMT_inc_code[2] = {0, 0};

/* --------------------------------------------------------------------------- */
/* Matlab Gateway routine */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	int	error = FALSE, global = FALSE, geo = FALSE; 
	int	zero = FALSE, bat = TRUE, do_track = FALSE, min_thick = FALSE;
	int	inc_set = FALSE, do_grid = FALSE;
	int	m_var = FALSE, DO = TRUE, exact = TRUE, verbose = FALSE;
	int triangulate = FALSE, raw = FALSE, stl = FALSE;
	int	i, j, ii, nx, ny, k, kk, ij, ndata_r = 0;
	int	ndata_p, ndata_xyz, ndata_t, nx_p, ny_p, n_vert_max;
	int	z_th, n_triang, ndata_m, ndata_s, n_swap = 0, nFacet;
	int	argc = 0, n_arg_no_char = 0;
	int	do_FV = FALSE, nv, node_offset;
	int	km, pm;		/* index of current body facet (for mag only) */
	mwSize	*faces = NULL;
	float	one_100, *ptr_s = NULL, *vertices = NULL, *out_s;
	double	s_rad = 50000, s_rad2, z_dir = -1, z0 = 0.1, rho = 0.0, zobs = 0, dz = 0;
	double	w = 0.0, e = 0.0, s = 0.0, n = 0.0, t_mag, a, DX, DY;
	double	t_dec, t_dip, m_int, m_dec, m_dip, cc_t, cs_t, s_t, *ptr_d, *out_d;
	double	*x_obs = NULL, *y_obs = NULL, *z_obs = NULL, *x = NULL, *y = NULL, *cos_vec = NULL;
	double	west, east, south, north, x_inc, y_inc;
	char	*xyz_file = CNULL, *t_file = CNULL;
	char	line[512], *ptr, *raw_file, *mag_file, *stl_file;
	char	**argv, *type;
	mxArray	*mx_ptr;

	FILE	*fp = NULL, *fp_t = NULL, *fp_xyz = NULL, *fp_r = NULL, *fp_m = NULL;
	FILE	*fp_s = NULL;

	d_to_m = 2.0 * M_PI * 6371008.7714 / 360.0;

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
	argv[0] = "xyzokb";
	for (i = 1; i < argc; i++)
		argv[i] = (char *)mxArrayToString(prhs[i+n_arg_no_char-1]);

	west = east = south = north = 0.0;
	
	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {
			
				case 'R':
					error += decode_R (argv[i], &west, &east, &south, &north);
					break;
				case 'B':
					if ((sscanf(&argv[i][2], "%lf/%lf/%lf/%lf/%lf", &t_dec, &t_dip, &m_int, &m_dec, &m_dip)) != 5) {
						mexPrintf ("SYNTAX ERROR -B option: Can't dechiper values\n");
						error = TRUE;
					}
					mag = TRUE;
					grav = FALSE;
                                        break;
				case 'C':
					rho = atof (&argv[i][2]);
					rho *= 6.674e-6;
					break;
				case 'D':
					z_dir = 1;
					break;
				case 'I':
					GMT_getinc (&argv[i][2], &x_inc, &y_inc);
					inc_set = TRUE;
					break;
	 			case 'L':
					zobs = atof (&argv[i][2]);
					break;
	 			case 'M':
					geo = TRUE;
					break;
	 			case 'P':
					dz = atof (&argv[i][2]);
					const_th = TRUE;
					break;
	 			case 'S':
					s_rad = atof (&argv[i][2]);
					s_rad *= 1000.;
					exact = FALSE;
					break;
	 			case 'T':
					bat = FALSE;
					break;
				case 't': /* Selected input mesh format */ 
					switch (argv[i][2]) {
						case 'd':	/* Surface computed by triangulate */
							strcpy (line, &argv[i][3]);
							ptr = strtok (line, "/");
							j = 0;
							while (ptr) {
								switch (j) {
									case 0:
										xyz_file = ptr;
										break;
									case 1:
										t_file = ptr;
										break;
									case 2:
										m_var = TRUE;
										mag = TRUE;
										grav = FALSE;
										if (ptr[1] == '2') m_var2 = TRUE;
										else if (ptr[1] == '3') m_var3 = TRUE;
										else if (ptr[1] == '4') m_var4 = TRUE;
										else m_var1 = TRUE;
										break;
									default:
										break;
								}
								ptr = strtok (CNULL, "/");
								j++;
							}
							if (j != 2 && j != 3) {
								mexPrintf ("SYNTAX ERROR -t option: Must give names for data points and vertex files\n");
								error = TRUE;
							}
							triangulate = TRUE;
							break;
						case 'r':	/* Closed volume in RAW format */
	 						raw_file = &argv[i][3];
							raw = TRUE;
							break;
						case 's':	/* Closed volume in STL1format */
	 						stl_file = &argv[i][3];
							stl = TRUE;
							break;
						default:
							error = TRUE;
							break;
					}
					break;
	 			case 'V':
					verbose = TRUE;
					break;
	 			case 'z':
					dz = atof (&argv[i][2]);
					min_thick = TRUE;
					break;
				case 'Z':
					z0 = atof(&argv[i][2]);
					break;
	 			case '0':
					zero = TRUE;
					break;
				default:
					error = TRUE;
					break;
			}
		}
		else {
			error = TRUE;
			mexPrintf ("Unrecognized option %s\n", argv[i]);
		}
	}
	
	if (nlhs == 0) {
		mexPrintf("xyzokb - Computes the gravity/magnetic effect of a body by the method of Okabe\n\n");
		mexPrintf("usage: xyzokb [-C<density>] [-R<w/e/s/n>]\n");
		mexPrintf("\t[-Z<level>] [-T] [-L<z_observation>] [-S<radius>] [-D]\n");
		mexPrintf("\t[z<dz>] [-Bf_dec/f_dip/m_int/m_dec/m_dip]\n");
		mexPrintf("\t[-0] [-t[<d>xyz_file/vert_file[/m]]|<r>raw_file] [-P<thick>] [-V]\n");

		mexPrintf("\txyzfile is the file whose grav/mag efect is to be computed.\n");
		mexPrintf("\t-B sets parameters for computation of magnetic anomaly\n");
		mexPrintf("\t   f_dec/f_dip -> geomagnetic declination/inclination\n");
		mexPrintf("\t   m_int/m_dec/m_dip -> body magnetic intensity/declination/inclination\n");
		mexPrintf("\t-C sets body density in SI\n");
		mexPrintf("\t-I sets the grid spacing for the grid.  Append m for minutes, c for seconds\n");
		mexPrintf("\t-L sets level of observation [Default = 0]\n");
		mexPrintf("\t-M Map units TRUE; x,y in degrees, dist units in m.\n\t   [Default dist unit = x,y unit].\n");
		mexPrintf("\t-P give layer thickness in m [Default = 0 m]\n");
		mexPrintf("\t-S search radius in km\n");
		mexPrintf("\t-t Give either names of xyz[m] and vertex files or of a raw file defining a close surface.\n");
		mexPrintf("\t   In the first case append an d imediatly after -t and optionaly a /m after the vertex file name.\n");
		mexPrintf("\t   In the second case append an r imediatly after -t and before the raw file name.\n");
		mexPrintf("\t-Z z level of reference plane [Default = 0]\n");
		mexPrintf("\t-z don't compute prism efect if it's thickness < dz [Default = 0 m]\n");
		mexPrintf("\t-0 for some unknown reason some observation points far from poligon\n");
		mexPrintf("\t   which effect is beeing computed have a slightly negative gravity\n");
		mexPrintf("\t   contribution. This option turns it to zero\n");
		return;
	}

	if (nlhs == 0 || nlhs > 2)
		mexErrMsgTxt("XYZOKB: Error. Must select one or 2 outputs.");

	if (!mxIsStruct(prhs[0])) 
		mexErrMsgTxt("XYZOKB: Error. First argument must be a structure.");

	/* ------------------------------------------------------------------------------------------ */
	mx_ptr = mxGetField(prhs[0], 0, "type");
	if (mx_ptr == NULL)
		mexErrMsgTxt("XYZOKB: Input Error. Structure must have a member named 'type'.");
	else
		type = (char *)mxArrayToString(mx_ptr);

	if (!strcmp(type, "FV") || !strcmp(type, "fv")) {
		mx_ptr = mxGetField(prhs[0], 0, "faces");
		if (mx_ptr == NULL)
			mexErrMsgTxt("XYZOKB: Input Error. Structure must have a member named 'faces'.");
		if (mxGetN(mx_ptr) != 3)
			mexErrMsgTxt("XYZOKB: faces array must be a Mx3 array");
		if (!mxIsInt32( mx_ptr))
			mexErrMsgTxt("XYZOKB: The Faces array must be of type INT32");
	
		nFacet = mxGetM( mx_ptr );
		faces = (mwSize *)mxGetData( mx_ptr );

		mx_ptr = mxGetField(prhs[0], 0, "vertices");
		if (mx_ptr == NULL)
			mexErrMsgTxt("XYZOKB: Input Error. Structure must have a member named 'vertices'.");
		if (mxGetN(mx_ptr) != 3)
			mexErrMsgTxt("XYZOKB: faces array must be a Nx3 array");
		if (!mxIsSingle( mx_ptr))
			mexErrMsgTxt("XYZOKB: The Faces array must be of type SINGLE");

		nv = mxGetM( mx_ptr );
		vertices = (float *)mxGetData( mx_ptr );
		do_FV = TRUE;

		mx_ptr = mxGetField(prhs[0], 0, "hdr");
		if (mx_ptr != NULL) {
			if ( (mxGetM( mx_ptr ) * mxGetN( mx_ptr )) != 9)
				mexErrMsgTxt("XYZOKB: The 'hdr' member must contain a 9 elements header array");
			if (!mxIsDouble( mx_ptr))
				mexErrMsgTxt("XYZOKB: The 'hdr' array must be of type DOUBLE");

			ptr_d = mxGetPr( mx_ptr );
			west = ptr_d[0];	east = ptr_d[1];
			south = ptr_d[2];	north = ptr_d[3];
			x_inc = ptr_d[7];	y_inc = ptr_d[8];
			node_offset = (int)ptr_d[6];
			do_grid = TRUE;
		}

		mx_ptr = mxGetField(prhs[0], 0, "track");
		if (mx_ptr != NULL && !do_grid) {
			if (mxGetN( mx_ptr ) > 2)
				mexErrMsgTxt("XYZOKB: The XY (track) array must be a Mx2 array.");

			ndata_p = mxGetM( mx_ptr );
			if ((data = (struct DATA *) mxCalloc ((size_t) ndata_p, sizeof(struct DATA)) ) == NULL)
				mexErrMsgTxt("XYZOKB: Could not allocate memory to hold the track profile structure");

			if (mxIsDouble( mx_ptr)) {
				ptr_d = mxGetPr( mx_ptr );
				for (i = 0; i < ndata_p; i++) {
					data[i].x = ptr_d[i];
					data[i].y = ptr_d[i + ndata_p];
				}
			}
			else if (mxIsSingle( mx_ptr)) {
				ptr_s = (float *)mxGetData( mx_ptr );
				for (i = 0; i < ndata_p; i++) {
					data[i].x = ptr_s[i];
					data[i].y = ptr_s[i + ndata_p];
				}
			}
			else
				mexErrMsgTxt("XYZOKB: The XY (track) array must be of type DOUBLE or SINGLE");

			do_track = TRUE;
		}
		else if (mx_ptr != NULL && do_grid)
			mexPrintf ("XYZOKB: Warning: Given both grid and profile options. Ignoring profile.\n");
	}

	/* ------------------------------------------------------------------------------------------ */

	if (rho == 0.0 && !mag && !m_var4) {
		mexPrintf ("XYZOKB: ERROR:  Must specify either -Cdensity or -B<stuff>\n");
		error++;
	}
	if (const_th && min_thick) {
		mexPrintf ("XYZOKB: ERROR:  Cannot have -P and -z options\n");
		error++;
	}
	if (raw && !exact) {
		mexPrintf ("XYZOKB: Warning -tr overrides -S\n");
		exact = TRUE;
	}

	if (zero && mag) {
		mexPrintf ("XYZOKB: Warning -0 option will wipe out the negative part of,\n");
		mexPrintf ("magnetic anomaly so don't be surprized with the end result\n");
	}

	if (!do_track && !do_grid)
		do_grid = (inc_set && x_inc > 0.0 && y_inc > 0.0);

	if (!do_track && !do_grid && (!inc_set || x_inc == 0 || y_inc == 0)) {
		mexPrintf ("XYZOKB: ERROR. Must specify -R & -I or select what to compute via input swtructure.\n", argv[0]);
		error++;
	}

	if (triangulate) {
		if ((fp_xyz = fopen (xyz_file, "r")) == NULL) {
			mexPrintf ("XYZOKB: Cannot open file %s\n", xyz_file);
			error++;
		}
		if ((fp_t = fopen (t_file, "r")) == NULL) {
			mexPrintf ("XYZOKB: Cannot open file %s\n", t_file);
			error++;
		}
	}
	else if (do_FV) {

		if ((raw_mesh = (struct RAW *) mxCalloc ((size_t) nFacet, sizeof(struct RAW)) ) == NULL) 
			mexErrMsgTxt("XYZOKB: Could not allocate memory to hold the triangles structure");

		parse_FV(vertices, faces, nFacet, nv, z_dir, raw_mesh);
		//n_swap = check_triang_cw (nFacet, 1);
	}
	else
		error++;

	if (error) 
		mexErrMsgTxt("XYZOKB: Cannot continue");

/* ---- Read files section ---------------------------------------------------- */

	if (!geo) d_to_m = 1.;
	z0 *= z_dir;
	if (triangulate) { /* Read triangle file output from triangulate */
		ndata_xyz = read_xyz (fp_xyz, xyz_file, z_dir, m_var, geo);
		/* read vertex file */
		ndata_t = read_t (fp_t, t_file);
		t_center = (struct TRI_CENTER *) mxCalloc ((size_t) ndata_t, sizeof(struct TRI_CENTER));
		/* compute aproximate center of each triangle */
		n_swap = check_triang_cw (ndata_t, 0);
		set_center (ndata_t);
	}

	if (verbose && n_swap > 0)
		mexPrintf ("WARNING: %d triangles had ccw order\n", n_swap);

/* ---------------------------------------------------------------------------- */

	nx = irint ((east - west) / x_inc) + !node_offset;
	ny = irint ((north - south) / y_inc) + !node_offset;
	nx_p = (do_grid) ? nx : ndata_p;
	ny_p = (do_grid) ? ny : ndata_p;

	/*  Allocate the output array */
	if (do_grid) {
		x = (double *) mxCalloc ((mwSize) nx, sizeof (double));
		y = (double *) mxCalloc ((mwSize) ny, sizeof (double));
		plhs[0] = mxCreateNumericMatrix (ny,nx,mxSINGLE_CLASS,mxREAL);
		out_s   = (float *)mxGetData(plhs[0]);
	}
	else {
		plhs[0] = mxCreateDoubleMatrix(ndata_p, 3, mxREAL);
		out_d   = mxGetPr(plhs[0]);
		for (k = 0; k < ndata_p; k++) {
			out_d[k] = data[k].x;
			out_d[k+ndata_p] = data[k].y;
		}
	}

	x_obs = (double *) mxCalloc ((mwSize) nx_p, sizeof (double));
	y_obs = (double *) mxCalloc ((mwSize) ny_p, sizeof (double));
	/*z_obs = (double *) mxCalloc ((mwSize) nx_p, sizeof (double));

	for (i = 0; i < nx_p; i++) z_obs[i] = zobs;*/

	if (do_track) { /* Need to compute observation coords only once. */
		for (i = 0; i < ndata_p; i++) {
			x_obs[i] = (geo) ? (data[i].x-central_long)*d_to_m*cos(data[i].y*D2R): data[i].x;
			y_obs[i] = (geo) ? -(data[i].y-central_lat)*d_to_m: data[i].y; /* - because y positive 'south' */
		}
	}

	/* Build observation point vectors */
	if (do_grid) {
		double x_inc2 = 0.0, y_inc2 = 0.0;
		if (!node_offset) {		/* Pixel reg */
			x_inc2 = x_inc / 2;
			y_inc2 = y_inc / 2;
		}
		for (i = 0; i < nx; i++)
			x[i] = (i == (nx-1)) ?  east - x_inc2 :   west + i * x_inc + x_inc2;
		for (j = 0; j < ny; j++)
			y[j] = (j == (ny-1)) ? -south - y_inc2 : -(north - j * y_inc - y_inc2);
	}
	if (triangulate) {
		n_triang = ndata_t;

		bd_desc.n_f = 5;		/* Number of prism facets */
		bd_desc.n_v = (int *)mxCalloc ((size_t) (bd_desc.n_f), sizeof (int));
		bd_desc.n_v[0] = 3;	bd_desc.n_v[1] = 3;
		bd_desc.n_v[2] = 4;	bd_desc.n_v[3] = 4;
		bd_desc.n_v[4] = 4;
		bd_desc.ind = (int *)mxCalloc ((size_t) (bd_desc.n_v[0] + bd_desc.n_v[1] +
			bd_desc.n_v[2] + bd_desc.n_v[3] + bd_desc.n_v[4]), sizeof (int));
		bd_desc.ind[0] = 0;	bd_desc.ind[1] = 1; 	bd_desc.ind[2] = 2;	/* top triang */
		bd_desc.ind[3] = 3;	bd_desc.ind[4] = 5; 	bd_desc.ind[5] = 4;	/* bot triang */
		bd_desc.ind[6] = 1;	bd_desc.ind[7] = 4; 	bd_desc.ind[8] = 5;	bd_desc.ind[9] = 2;
		bd_desc.ind[10] = 0;	bd_desc.ind[11] = 3;	bd_desc.ind[12] = 4;	bd_desc.ind[13] = 1;
	 	bd_desc.ind[14] = 0;	bd_desc.ind[15] = 2;	bd_desc.ind[16] = 5;	bd_desc.ind[17] = 3;
	}
	else if (raw || do_FV) {
		n_triang = (raw) ? ndata_r : nFacet;
		bd_desc.n_f = 1;
		bd_desc.n_v = (int *)mxCalloc ((size_t) (bd_desc.n_f), sizeof (int));
		bd_desc.n_v[0] = 3;
		bd_desc.ind = (int *)mxCalloc ((size_t) (bd_desc.n_v[0]), sizeof (int));
		bd_desc.ind[0] = 0;	bd_desc.ind[1] = 1; 	bd_desc.ind[2] = 2;
	}
	else
		mexPrintf ("%s It shouldn't pass here\n", argv[0]);

	/* Allocate a structure that will be used inside okabe().
	   We do it here to avoid thousands of alloc/free that would result if done in okabe() */ 
	n_vert_max = bd_desc.n_v[0];
	for (i = 1; i < bd_desc.n_f; i++) {
		n_vert_max = MAX(bd_desc.n_v[i], n_vert_max);
	}
	loc_or = (struct LOC_OR *) mxCalloc ((size_t) (n_vert_max+1), sizeof(struct LOC_OR));


	if (mag) { /* 1e2 is a factor to obtain nT from magnetization in A/m */
		cc_t = cos(m_dip*D2R)*cos((m_dec - 90.)*D2R);
		cs_t = cos(m_dip*D2R)*sin((m_dec - 90.)*D2R);
		s_t = sin(m_dip*D2R);
		if (!m_var4) {		/* In all the other cases the field parameters are constatnt */
			mag_param = (struct MAG_PARAM *)mxCalloc ((size_t) 1, sizeof(struct MAG_PARAM));
			mag_param[0].rim[0] = 1e2*cos(t_dip*D2R) * cos((t_dec - 90.)*D2R);
			mag_param[0].rim[1] = 1e2*cos(t_dip*D2R) * sin((t_dec - 90.)*D2R);
			mag_param[0].rim[2] = 1e2*sin(t_dip*D2R);
		}
		if (!m_var) { /* Case of constant magnetization */
			mag_var = (struct MAG_VAR *) mxCalloc ((size_t) 1, sizeof(struct MAG_VAR));
			mag_var[0].rk[0] = m_int * cc_t;
			mag_var[0].rk[1] = m_int * cs_t;
			mag_var[0].rk[2] = m_int * s_t;
		}
		else { /* The triangles have a non-constant magnetization */
			if ((mag_var = (struct MAG_VAR *)mxCalloc ((size_t) n_triang, sizeof(struct MAG_VAR))) == NULL)
				{no_sys_mem("xyzokb", n_triang);}
			if (m_var1) {		/* Only the mag intensity changes, Dec & Dip are the same as the Field */
				for (i = 0; i < n_triang; i++) {
					t_mag = (mag_int[vert[i].a] + mag_int[vert[i].b] + mag_int[vert[i].c])/3.;
					mag_var[i].rk[0] = t_mag * cc_t;
					mag_var[i].rk[1] = t_mag * cs_t; 
					mag_var[i].rk[2] = t_mag * s_t; 
				}
			}
			else if (m_var2) {	/* Both mag intensity & dip varies. Dec is Zero (axial dipole) */
				for (i = 0; i < n_triang; i++) {
					t_mag = (mag_var2[vert[i].a].m + mag_var2[vert[i].b].m + mag_var2[vert[i].c].m)/3.;
					t_dip = (mag_var2[vert[i].a].m_dip + mag_var2[vert[i].b].m_dip + mag_var2[vert[i].c].m_dip)/3.;
					mag_var[i].rk[0] = 0.;
					mag_var[i].rk[1] = -t_mag * cos(t_dip*D2R); 
					mag_var[i].rk[2] = t_mag * sin(t_dip*D2R); 
				}
			}
			else if (m_var3) { 	/* Both mag intensity, mag_dec & mag_dip varies. */
				for (i = 0; i < n_triang; i++) {
					t_mag = (mag_var3[vert[i].a].m + mag_var3[vert[i].b].m + mag_var3[vert[i].c].m)/3.;
					t_dec = (mag_var3[vert[i].a].m_dec + mag_var3[vert[i].b].m_dec + mag_var3[vert[i].c].m_dec)/3.;
					t_dip = (mag_var3[vert[i].a].m_dip + mag_var3[vert[i].b].m_dip + mag_var3[vert[i].c].m_dip)/3.;
					mag_var[i].rk[0] = t_mag * cos(t_dip*D2R) * cos((t_dec - 90)*D2R);
					mag_var[i].rk[1] = t_mag * cos(t_dip*D2R) * sin((t_dec - 90)*D2R); 
					mag_var[i].rk[2] = t_mag * sin(t_dip*D2R); 
				}
			}
			else {			/* Everything varies. */ 
				if ((mag_param = (struct MAG_PARAM *)mxCalloc ((size_t) n_triang, sizeof(struct MAG_PARAM))) == NULL)
					{no_sys_mem("xyzokb", n_triang);}
				for (i = 0; i < n_triang; i++) {
					t_dec = (mag_var4[vert[i].a].t_dec + mag_var4[vert[i].b].t_dec + mag_var4[vert[i].c].t_dec)/3.;
					t_dip = (mag_var4[vert[i].a].t_dip + mag_var4[vert[i].b].t_dip + mag_var4[vert[i].c].t_dip)/3.;
					mag_param[i].rim[0] = 1e2*cos(t_dip*D2R) * cos((t_dec - 90.)*D2R);
					mag_param[i].rim[1] = 1e2*cos(t_dip*D2R) * sin((t_dec - 90.)*D2R);
					mag_param[i].rim[2] = 1e2*sin(t_dip*D2R);
					t_mag = (mag_var4[vert[i].a].m + mag_var4[vert[i].b].m + mag_var4[vert[i].c].m)/3.;
					t_dec = (mag_var4[vert[i].a].m_dec + mag_var4[vert[i].b].m_dec + mag_var4[vert[i].c].m_dec)/3.;
					t_dip = (mag_var4[vert[i].a].m_dip + mag_var4[vert[i].b].m_dip + mag_var4[vert[i].c].m_dip)/3.;
					mag_var[i].rk[0] = t_mag * cos(t_dip*D2R) * cos((t_dec - 90)*D2R);
					mag_var[i].rk[1] = t_mag * cos(t_dip*D2R) * sin((t_dec - 90)*D2R); 
					mag_var[i].rk[2] = t_mag * sin(t_dip*D2R); 
				}
			}
		}
	}

/* ---------------> Now start computing <------------------------------------- */
	one_100 = (float)(n_triang / 100.);
	s_rad2 = s_rad*s_rad;
	if (do_grid) {		/* Compute the cos(lat) vector only once */
		cos_vec = (double *) mxCalloc ((mwSize) ny, sizeof (double));
		for (i = 0; i < ny; i++)
			cos_vec[i] = (geo) ? cos(y[i]*D2R): 1;
	}
	for (i = j = 0; i < n_triang; i++) {
		if (verbose && i > j*one_100) {
			mexPrintf ("computed %.2d%s of %d prisms\r", j, "%", n_triang);
			j++;
		}
		km = (m_var)  ? i : 0;
		pm = (m_var4) ? i : 0;

		/* Don't wast time with zero mag triangles */
		if (mag && m_var && mag_var[i].rk[0] == 0 && mag_var[i].rk[1] == 0 && mag_var[i].rk[2] == 0)
			continue;
		if (triangulate)
			z_th = facet_triangulate (i, z0, dz, bat);
		else if (do_FV || raw)
			z_th = facet_raw (i, geo);
		if (z_th) {
			if (do_grid) {
				for (ii = ij = 0; ii < nx; ii++) {
					for (k = 0; k < ny; k++) {
						/* In the next line +central_lat because y was already *= -1 */
						y_obs[k]  = (geo) ? (y[k]  + central_lat) * d_to_m : y[k];
						x_obs[ii] = (geo) ? (x[ii] - central_long)*d_to_m*cos_vec[k] : x[ii]; 
						if (exact)
							DO = TRUE;
						else {
							DX = t_center[i].x - x_obs[ii];
							DY = t_center[i].y - y_obs[k];
							DO = (DX*DX + DY*DY) < s_rad2; 
						}
						if (DO) {
							a = okabe (rho, x_obs[ii], y_obs[k], zobs, bd_desc, km, pm, loc_or);
							if (zero && a < 0) continue;
							if (!mxIsNaN(a))
								out_s[ij] += (float)a;
						}
						ij++;
					}
				}
			}
			else { 		/* do Track */
				for (kk = 0; kk < ndata_p; kk++){
					if (exact)
						DO = TRUE;
					else {
						DX = t_center[i].x - x_obs[kk];
						DY = t_center[i].y - y_obs[kk];
						DO = (DX*DX + DY*DY) < s_rad2; 
					}
					if (DO) {
						a = okabe (rho, x_obs[kk], y_obs[kk], zobs, bd_desc, km, pm, loc_or);
						if (zero && a < 0) continue;
						if (!mxIsNaN(a))
							out_d[kk + 2*ndata_p] += a;
					}
				}
			}
		}
	}

	if (do_grid && nlhs == 2) {	/* Create the output header vector */
		plhs[1] = mxCreateDoubleMatrix(1, 9, mxREAL);
		ptr_d = mxGetPr( plhs[1] );
		ptr_d[0] = west;	ptr_d[1] = east;
		ptr_d[2] = south;	ptr_d[3] = north;
		ptr_d[7] = x_inc;	ptr_d[8] = y_inc;
		ptr_d[6] = node_offset;
	}

	if (x) mxFree ((void *)x);
	if (y) mxFree ((void *)y);
	if (z_obs) mxFree ((void *)z_obs);
	mxFree ((void *)x_obs);
	mxFree ((void *)y_obs);
	mxFree ((void *)data);
	mxFree ((void *)triang);
	mxFree ((void *)raw_mesh);
	mxFree ((void *)t_center);
	mxFree ((void *)vert);
	mxFree ((void *)mag_param);
	mxFree ((void *)mag_var);
	mxFree ((void *)bd_desc.n_v);
	mxFree ((void *)bd_desc.ind);
	if (cos_vec) mxFree ((void *)cos_vec);
	if (m_var1) mxFree ((void *)mag_int);
	if (m_var2) mxFree ((void *)mag_var2);
	if (m_var3) mxFree ((void *)mag_var3);
	if (m_var4) mxFree ((void *)mag_var4);
}
/* -------------------------------------------------------------------------*/
int read_xyz (FILE *fp_xyz, char *xyz_file, double z_dir, int m_var, int geo) {
	/* read xyz[m] file with point data coordinates */
	return (0);
}

/* -----------------------------------------------------------------*/
int read_t (FILE *fp_t, char *t_file) {
	/* read file with vertex indexes of triangles */
	return (0);
}

/* -----------------------------------------------------------------*/
void parse_FV(float *v, mwSize *f, int nFacet, int nv, double z_dir, struct RAW *raw_mesh) {
	int	i, i1, i2, i3;

	for (i = 0; i < nFacet; i++) {
		i1 = f[i];		i2 = f[i + nFacet];		i3 = f[i + 2*nFacet];
		raw_mesh[i].t1[0] = (double)v[i1-1];
		raw_mesh[i].t1[1] = (double)-v[i1 + nv-1];
		raw_mesh[i].t1[2] = (double)v[i1 + 2*nv-1] * z_dir;

		raw_mesh[i].t2[0] = (double)v[i2-1];
		raw_mesh[i].t2[1] = (double)-v[i2 + nv-1];
		raw_mesh[i].t2[2] = (double)v[i2 + 2*nv-1] * z_dir;

		raw_mesh[i].t3[0] = (double)v[i3-1];
		raw_mesh[i].t3[1] = (double)-v[i3 + nv-1];
		raw_mesh[i].t3[2] = (double)v[i3 + 2*nv-1] * z_dir;
	}
}

/* -----------------------------------------------------------------*/
int read_stl (FILE *fp_s, char *stl_file, double z_dir) {
	/* read a file with triagles in the stl format and returns nb of triangles */
	int n_alloc, ndata_s, chunk = 50000;
	double in[3];
	char line[512], text[100], ver_txt[100];

	n_alloc = chunk;
	ndata_s = 0;
	if ((raw_mesh = (struct RAW *) mxCalloc ((size_t) n_alloc, sizeof(struct RAW)) ) == NULL) 
		{no_sys_mem("xyzokb --> (read_stl)", n_alloc);}
	
	while (fgets (line, 512, fp_s)) { 
		sscanf (line, "%s", &text);
		if (strcmp (text, "outer") == 0) {
			fgets (line, 512, fp_s); /* get first vertex */
			if(sscanf (line, "%s %lg %lg %lg", &ver_txt, &in[0], &in[1], &in[2]) !=4)
				mexPrintf ("ERROR deciphering triangle %d of %s\n", ndata_s+1, stl_file);
			raw_mesh[ndata_s].t1[0] = in[0];
			raw_mesh[ndata_s].t1[1] = -in[1];
			raw_mesh[ndata_s].t1[2] = in[2] * z_dir;
			fgets (line, 512, fp_s); /* get second vertex */
			if(sscanf (line, "%s %lg %lg %lg", &ver_txt, &in[0], &in[1], &in[2]) !=4)
				mexPrintf ("ERROR deciphering triangle %d of %s\n", ndata_s+1, stl_file);
			raw_mesh[ndata_s].t3[0] = in[0];
			raw_mesh[ndata_s].t3[1] = -in[1];
			raw_mesh[ndata_s].t3[2] = in[2] * z_dir;
			fgets (line, 512, fp_s); /* get third vertex */
			if(sscanf (line, "%s %lg %lg %lg", &ver_txt, &in[0], &in[1], &in[2]) !=4)
				mexPrintf ("ERROR deciphering triangle %d of %s\n", ndata_s+1, stl_file);
			raw_mesh[ndata_s].t2[0] = in[0];
			raw_mesh[ndata_s].t2[1] = -in[1];
			raw_mesh[ndata_s].t2[2] = in[2] * z_dir;
			ndata_s++;
              		if (ndata_s == n_alloc) { /* with bad luck we have a flaw here */
				n_alloc += chunk;
				if ((raw_mesh = (struct RAW *) realloc ((void *)raw_mesh, (size_t)(n_alloc * sizeof(struct RAW)))) == NULL) 
					{no_sys_mem("xyzokb --> (read_stl)", n_alloc);}
			}
		}
		else
			continue;
	}
	fclose(fp_s);
	return (ndata_s);
}

/* -----------------------------------------------------------------*/
int facet_triangulate (int i, double z0, double dz, int bat) {
	/* Sets coodinates for the facet whose effect is beeing calculated */
	double x_a, x_b, x_c, y_a, y_b, y_c, z_a, z_b, z_c;

	x_a = triang[vert[i].a].x;	x_b = triang[vert[i].b].x;	x_c = triang[vert[i].c].x;
	y_a = triang[vert[i].a].y;	y_b = triang[vert[i].b].y;	y_c = triang[vert[i].c].y;
	z_a = triang[vert[i].a].z;	z_b = triang[vert[i].b].z;	z_c = triang[vert[i].c].z;
	/* top triang */
	xx[0] = x_a;	yy[0] = y_a;
	xx[1] = x_b;	yy[1] = y_b;
	xx[2] = x_c;	yy[2] = y_c;
	/* bot triang */
	xx[3] = xx[0];	yy[3] = yy[0];
	xx[4] = xx[1];	yy[4] = yy[1];
	xx[5] = xx[2];	yy[5] = yy[2];

	xx[6] = xx[1];	yy[6] = yy[1];
	xx[7] = xx[4];	yy[7] = yy[4];
	xx[8] = xx[5];	yy[8] = yy[5];
	xx[9] = xx[2];	yy[9] = yy[2];

	xx[10] = xx[1];	yy[10] = yy[1];
	xx[11] = xx[0];	yy[11] = yy[0];
	xx[12] = xx[3];	yy[12] = yy[3];
	xx[13] = xx[4];	yy[13] = yy[4];

	xx[14] = xx[0];	yy[14] = yy[0];
	xx[15] = xx[2];	yy[15] = yy[2];
	xx[16] = xx[5];	yy[16] = yy[5];
	xx[17] = xx[3];	yy[17] = yy[3];

	if (const_th) { /* Layer of constant thickness */
		zz[0] = z_a;		zz[1] = z_b;		zz[2] = z_c;
		zz[3] = z_a + dz;	zz[4] = z_b + dz;	zz[5] = z_c + dz;
		zz[6] = zz[1];		zz[7] = zz[4];		zz[8] = zz[5];
		zz[9] = zz[5];		zz[10] = zz[1];		zz[11] = zz[0];
		zz[12] = zz[3];		zz[13] = zz[4];		zz[14] = zz[0];
		zz[15] = zz[2];		zz[16] = zz[5];		zz[17] = zz[3];

		return (1);
	}
	if (bat) { /* Triangle mesh defines a bathymetric surface (TA MIXORDADO (== NOS DOIS CASOS)) */
		zz[0] = z_a;	zz[1] = z_b;	zz[2] = z_c;
		zz[3] = z0;	zz[4] = z0;	zz[5] = z0;
		/*zz[0] = z0;	zz[1] = z0;	zz[2] = z0;
		zz[3] = z_a;	zz[4] = z_b;	zz[5] = z_c; */
		if (fabs(zz[0] - zz[3]) > dz || fabs(zz[1] - zz[4]) > dz || fabs(zz[2] - zz[5]) > dz) return 1;
		else return (0);
	}
	else { /* Triangle mesh defines a topographic surface */
		zz[0] = z_a;	zz[1] = z_b;	zz[2] = z_c;
		zz[3] = z0;	zz[4] = z0;	zz[5] = z0;
		if (fabs(zz[0] - zz[3]) > dz || fabs(zz[1] - zz[4]) > dz || fabs(zz[2] - zz[5]) > dz) return 1;
		else return (0);
	}
}

/* -----------------------------------------------------------------*/
int facet_raw (int i, int geo) {
	/* Sets coodinates for the facet in the RAW format */
	double cos_a, cos_b, cos_c, x_a, x_b, x_c, y_a, y_b, y_c, z_a, z_b, z_c;

	x_a = raw_mesh[i].t1[0];   x_b = raw_mesh[i].t2[0];   x_c = raw_mesh[i].t3[0];
	y_a = raw_mesh[i].t1[1];   y_b = raw_mesh[i].t2[1];   y_c = raw_mesh[i].t3[1];
	z_a = raw_mesh[i].t1[2];   z_b = raw_mesh[i].t2[2];   z_c = raw_mesh[i].t3[2];
	if (geo) {
		cos_a = cos(y_a*D2R);	cos_b = cos(y_b*D2R);	cos_c = cos(y_c*D2R);
	}
	else {cos_a = cos_b = cos_c = 1.;}
	xx[0] = x_a*d_to_m*cos_a;	yy[0] = y_a*d_to_m;
	xx[1] = x_b*d_to_m*cos_b;	yy[1] = y_b*d_to_m;
	xx[2] = x_c*d_to_m*cos_c;	yy[2] = y_c*d_to_m;
	zz[0] = z_a;	zz[1] = z_b;	zz[2] = z_c;
	return 1; /* Allways return 1 */
}

int	no_sys_mem (char *where, int n) {	
		mexPrintf ("Fatal Error: %s could not allocate memory, n = %d\n", where, n);
		return (-1);
}

/* --------------------------------------------------------------------- */
double okabe (double rho, double x_o, double y_o, double z_o, struct BODY_DESC bd_desc, int km, int pm, struct LOC_OR *loc_or) {
	double okb = 0.0, tot = 0.;
	int i, n_vert, l, k, cnt_v = 0;
	int top = TRUE;

/* x_o, y_o, z_o are the coordinates of the observation point
 * rho is the body density times G constant
 * km is an: index of current body facet (if they have different mags); or 0 if mag=const
 * bd_desc is a structure containing the body's description. It contains the following members
 * n_f -> number of facets (int)
 * n_v -> number of vertex of each facet (pointer)
 * ind -> index describing the vertex order of each facet. These index must
 *	  describe the facet in a clock-wise order when viewed from outside.
	
/*  _________________________________________________________________ 
    |                                                               | 
    |  Reference : Okabe, M., Analytical expressions for gravity    |
    |     anomalies due to polyhedral bodies and translation into   |
    |     magnetic anomalies, Geophysics, 44, (1979), p 730-741.    |
    |_______________________________________________________________|
    _____________________________________________________________________
    |                                                                   |
    |  Ifac decrit le corps (ATTN : Integer*2) :                        |
    |   - Il y a Nff facettes qui sont decrites dans le sens des        |
    |          aiguilles d'une montre si on regarde le corps de         |
    |          l'exterieur. Mxsomf = Max de sommets / face              |
    |   - Le premier nombre indique le nombre de factettes. Suivent     |
    |       alors des groupes de nombres dont le premier de chaque      |
    |       groupe est le nombre de sommets de la facette, suivi par    |
    |       les indices (pointeurs) des sommets (rang dans Xx,Yy,Zz)    |
    |       correspondant a la facette.                                 |
    |                                                                   |
    |  Par exemple pour un cube                _________________        |
    |  (Nff=6; 4 sommets /face)              /|         X (Nord)        |
    |  [Ifac] = { 6,  4, 1,2,3,4,           / |                         |
    |                 4, 2,6,7,3,          /  |     1 ________ 2        |
    |                 4, 4,3,7,8,         /   |      /       /|         |
    |                 4, 5,1,4,8,      Y /    |     /       / |         |
    |                 4, 1,5,6,2,      (Est)  |  4 /_______/3 |         |
    |                 4, 5,8,7,6 }            |    |       |  |         |
    |                                         |    |       | / 6        |
    |                                       Z |    |       |/           |
    |                                         V    |_______/            |
    |                                             8         7           |
    |___________________________________________________________________|
    |                                                                   |
    |  X,Y ET Z sont les tableaux des coordonness des pts de mesure     |
    |___________________________________________________________________| */

	for (i = 0; i < bd_desc.n_f; i++) {	/* Loop over facets */
		n_vert = bd_desc.n_v[i];	/* Number of vertices of each face */
		if (n_vert < 3) 
			fprintf (stderr, "Warning: facet with less than 3 vertex\n");
		for (l = 0; l < n_vert; l++) {
			k = bd_desc.ind[l+cnt_v];
			loc_or[l].x = xx[k] - x_o;
			loc_or[l].y = yy[k] - y_o;
			loc_or[l].z = zz[k] - z_o;
		}
		rot_17 (n_vert, top, loc_or); /* rotate coords by eq (17) of okb */
		okb += (grav) ? okb_grv (n_vert, loc_or): okb_mag (n_vert, km, pm, loc_or);
		cnt_v += n_vert;
	}
	tot = (grav) ? okb * rho: okb;
	return (tot);
}

/* ---------------------------------------------------------------------- */
void rot_17 (int n_vert, int top, struct LOC_OR *loc_or) {
	/* Rotates coordinates by teta and phi acording to equation (17) of Okabe */
	/* store the result in external structure loc_or and angles c_tet s_tet c_phi s_phi */
	double xi, xj, xk, yi, yj, yk, zi, zj, zk, v, x, y, z;
	double r, r2, r_3d, Sxy, Szx, Syz;
	int i = 0, j, k, l;

	loc_or[n_vert].x = loc_or[0].x;		loc_or[n_vert].y = loc_or[0].y;	
	loc_or[n_vert].z = loc_or[0].z;		/* Last point = first point */

	if (top) { /* Currently, this is always true */
		j = i + 1;	k = i + 2;
		xi = loc_or[i].x;	xj = loc_or[j].x;	xk = loc_or[k].x;
		yi = loc_or[i].y;	yj = loc_or[j].y;	yk = loc_or[k].y;
		zi = loc_or[i].z;	zj = loc_or[j].z;	zk = loc_or[k].z;
		Sxy = xi * (yj - yk) + xj * (yk - yi) + xk * (yi - yj);
		Syz = yi * (zj - zk) + yj * (zk - zi) + yk * (zi - zj);
		Szx = zi * (xj - xk) + zj * (xk - xi) + zk * (xi - xj);
		r2 = Syz * Syz + Szx * Szx;	r = sqrt(r2);
		r_3d = sqrt(r2 + Sxy * Sxy);
		c_phi = - Sxy / r_3d;
		s_phi = r / r_3d;

		if (Szx == 0.0 && Syz == 0.0) { c_tet = 1.0;	s_tet = 0.0;}
		else { c_tet = - Syz / r;	s_tet = - Szx / r;}
		}
	else { /* Don't need to recompute angles, only do this */
		c_tet *= -1;	s_tet *= -1;	c_phi *= -1;
	}

	for (l = 0; l < n_vert + 1; l++) {
		x = loc_or[l].x;	y = loc_or[l].y;	z = loc_or[l].z;
		v = x * c_tet + y * s_tet;
		loc_or[l].x = v * c_phi - z * s_phi;
		loc_or[l].y = y * c_tet - x * s_tet;
		loc_or[l].z = v * s_phi + z * c_phi;
	}
}

/* ---------------------------------------------------------------------- */
double okb_grv (int n_vert, struct LOC_OR *loc_or) {
/*  Computes the gravity anomaly due to a facet. */
 
	double x1, x2, y1, y2, z2, z1, dx, dy, r, r_1, c_psi, s_psi;
	double grv = 0.0, grv_p;
	int l;

	if (fabs(c_phi) < FLT_EPSILON) return 0.0;
	for (l = 0; l < n_vert; l++) {
		x1 = loc_or[l].x;	x2 = loc_or[l+1].x;
		y1 = loc_or[l].y;	y2 = loc_or[l+1].y;
		dx = x2 - x1;	dy = y2 - y1;
		r = sqrt(dx*dx + dy*dy);
		r_1 = 1. / r;
		if (r > FLT_EPSILON) {
			c_psi = dx * r_1;	s_psi = dy * r_1;
			z2 = loc_or[l+1].z;	z1 = loc_or[l].z;
			grv_p = eq_30(c_psi, s_psi, x2, y2, z2) - eq_30(c_psi, s_psi, x1, y1, z1);
		}
		else
			grv_p = 0.0;
		grv += grv_p;
	}
	return (grv * c_phi);
}

/* ---------------------------------------------------------------------- */
double eq_30 (double c, double s, double x, double y, double z) {
	double r, Ji = 0.0, log_arg;

	r = sqrt(x * x + y * y + z * z);
	if (r > FLT_EPSILON) {
		if (fabs(z) > FLT_EPSILON && fabs(c) > FLT_EPSILON)
			Ji = -2. * z * atan ((x * c + (s + 1) * (y + r)) / (z * c));
		log_arg = x * c + y * s + r;
		if (log_arg > FLT_EPSILON)
			Ji += (x * s - y * c) * log(log_arg);
	}
	return Ji;
}

/* ---------------------------------------------------------------------- */
double okb_mag (int n_vert, int km, int pm, struct LOC_OR *loc_or) {
/*  Computes the total magnetic anomaly due to a facet. */
 
	double qsi1, qsi2, eta1, eta2, z2, z1, dx, dy, kx, ky, kz, v, r, c_psi, s_psi;
	double ano = 0.0, ano_p, mag_fac, xi, xi1, yi, yi1, mx, my, mz, r_1, tg_psi, auxil;
	int i;

	/*mag_fac = s_phi * (mag_param.rim[0] * c_tet + mag_param.rim[1] * s_tet) + mag_param.rim[2] * c_phi;*/
	mag_fac = s_phi * (mag_param[pm].rim[0] * c_tet + mag_param[pm].rim[1] * s_tet) + mag_param[pm].rim[2] * c_phi;

	if (fabs(mag_fac) < FLT_EPSILON) return 0.0;

	kx = mag_var[km].rk[0];	ky = mag_var[km].rk[1];	kz = mag_var[km].rk[2];
	v = kx * c_tet + ky * s_tet;
	mx = v * c_phi - kz * s_phi;	my = ky * c_tet - kx * s_tet;
	mz = v * s_phi + kz * c_phi;

	for (i = 0; i < n_vert; i++) {
		xi = loc_or[i].x;	xi1 = loc_or[i+1].x;
		yi = loc_or[i].y;	yi1 = loc_or[i+1].y;
		dx = xi1 - xi;	dy = yi1 - yi;
		r = sqrt(dx*dx + dy*dy);
		r_1 = 1. / r;
		if (r > FLT_EPSILON) {
			c_psi = dx * r_1;	s_psi = dy * r_1;
			tg_psi = dy / dx;
			auxil = my * c_psi - mx * s_psi;
			qsi1 = yi * s_psi + xi * c_psi;	qsi2 = yi1 * s_psi + xi1 * c_psi;
			eta1 = yi * c_psi - xi * s_psi;	eta2 = yi1 * c_psi - xi1 * s_psi;
			z1 = loc_or[i].z;	z2 = loc_or[i+1].z;
			ano_p = eq_43(mz, c_psi, tg_psi, auxil, qsi2, eta2, z2) - eq_43(mz, c_psi, tg_psi, auxil, qsi1, eta1, z1);
		}
		else
			ano_p = 0.0;
		ano += ano_p;
	}
	return (ano * mag_fac);
}

/* ---------------------------------------------------------------------- */
double eq_43 (double mz, double c, double tg, double auxil, double x, double y, double z) {
	double r, ez, Li = 0.0, tmp;

	ez = y * y + z * z;
	r = sqrt(x * x + ez);

	if (r > FLT_EPSILON) {
		if (fabs(z) > FLT_EPSILON && fabs(c) > FLT_EPSILON)
			Li = mz * atan((ez * tg - x * y) / (z * r));
		else
			Li = 0.0;
		tmp = x + r;
		if (tmp <= 0.)
			Li -= log(r - x) * auxil;
		else
			Li += log(tmp) * auxil;
	}
	return Li;
}
/* ---------------------------------------------------------------------- */
void set_center (int n_triang) {
	/* Calculates triangle center by an aproximate (iterative) formula */
	int i, j, k = 5;
	double x, y, z, xa[6], ya[6], xb[6], yb[6], xc[6], yc[6];

	for (i = 0; i < n_triang; i++) {
		xa[0] = (triang[vert[i].b].x + triang[vert[i].c].x) / 2.;
		ya[0] = (triang[vert[i].b].y + triang[vert[i].c].y) / 2.;
		xb[0] = (triang[vert[i].c].x + triang[vert[i].a].x) / 2.;
		yb[0] = (triang[vert[i].c].y + triang[vert[i].a].y) / 2.;
		xc[0] = (triang[vert[i].a].x + triang[vert[i].b].x) / 2.;
		yc[0] = (triang[vert[i].a].y + triang[vert[i].b].y) / 2.;
		for (j = 1; j <= k; j++) {
			xa[j] = (xb[j-1] + xc[j-1]) / 2.;
			ya[j] = (yb[j-1] + yc[j-1]) / 2.;
			xb[j] = (xc[j-1] + xa[j-1]) / 2.;
			yb[j] = (yc[j-1] + ya[j-1]) / 2.;
			xc[j] = (xa[j-1] + xb[j-1]) / 2.;
			yc[j] = (ya[j-1] + yb[j-1]) / 2.;
		}
		x = (xa[k]+xb[k]+xc[k])/3.;
		y = (ya[k]+yb[k]+yc[k])/3.;
		z = (triang[vert[i].a].z+triang[vert[i].b].z+triang[vert[i].c].z)/3.; 
		t_center[i].x = x;
		t_center[i].y = y;
		t_center[i].z = z;
	}
}

int check_triang_cw (int n, int type) {
	/* Checks that triangles are given in the correct clock-wise order.
	If not swap them. This is a tricky issue. In the case of "classic" 
	trihedron (x positive right; y positive "north" and z positive up),
	positive determinants signify counter clockwise order. However, in
	geomagnetic reference (x positive right; y positive "south" and z 
	positive down (OK, I know it's not exactly like this but instead 
	x->north; y->east; z->down)), counter clockwise order follows if
	determinant is negative. */

	int i, n_swaped = 0, tmp;
	double x1, x2, x3, y1, y2, y3, determ, d_tmp[3];

	for (i = 0; i < n; i++) {
		if (type == 0) { /* triangulate */
			x1 = triang[vert[i].a].x;	 y1 = triang[vert[i].a].y;
			x2 = triang[vert[i].b].x;	 y2 = triang[vert[i].b].y;
			x3 = triang[vert[i].c].x;	 y3 = triang[vert[i].c].y;
		}
		else if (type == 1) { /* raw */
			x1 = raw_mesh[i].t1[0];		y1 = raw_mesh[i].t1[1];
			x2 = raw_mesh[i].t2[0];		y2 = raw_mesh[i].t2[1];
			x3 = raw_mesh[i].t3[0];		y3 = raw_mesh[i].t3[1];
		}

		determ = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1);

		if (determ < 0.0) { /* counter clockwise triangle -> swap vertex order */
			if (type == 0) {
				tmp = vert[i].b;
				vert[i].b = vert[i].c;
				vert[i].c = tmp;
				n_swaped++;
			}
			else if (type == 1) {
				d_tmp[0] = raw_mesh[i].t2[0]; 
				d_tmp[1] = raw_mesh[i].t2[1]; 
				d_tmp[2] = raw_mesh[i].t2[2]; 
				raw_mesh[i].t2[0] = raw_mesh[i].t3[0]; 
				raw_mesh[i].t2[1] = raw_mesh[i].t3[1]; 
				raw_mesh[i].t2[2] = raw_mesh[i].t3[2]; 
				raw_mesh[i].t3[0] = d_tmp[0]; 
				raw_mesh[i].t3[1] = d_tmp[1]; 
				raw_mesh[i].t3[2] = d_tmp[2]; 
				n_swaped++;
			}
		}
	}
	return (n_swaped);
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

/* -------------------------------------------------------------------- */
#define GMT_INC_IS_M		1
#define GMT_INC_IS_KM		2
#define GMT_INC_IS_MILES	4
#define GMT_INC_IS_NMILES	8
#define GMT_INC_IS_NNODES	16
#define GMT_INC_IS_EXACT	32
#define GMT_INC_UNITS		15

#define GMT_DEG2MIN_F	60.0
#define GMT_MIN2DEG	(1.0 / GMT_DEG2MIN_F)
#define GMT_DEG2SEC_F	3600.0
#define GMT_SEC2DEG	(1.0 / GMT_DEG2SEC_F)

int GMT_getinc (char *line, double *dx, double *dy) {
	/* Special case of getincn use where n is two. */

	int n;
	double inc[2];

	/* Syntax: -I<xinc>[m|c|e|i|k|n|+|=][/<yinc>][m|c|e|i|k|n|+|=]
	 * Units: m = minutes
	 *	  c = seconds
	 *	  e = meter [Convert to degrees]
	 *	  i = miles [Convert to degrees]
	 *	  k = km [Convert to degrees]
	 *	  n = nautical miles [Convert to degrees]
	 * Flags: = = Adjust -R to fit exact -I [Default modifies -I to fit -R]
	 *	  + = incs are actually nx/ny - convert to get xinc/yinc
	 */

	n = GMT_getincn (line, inc, 2);
	*dx = inc[0] ; *dy = inc[1];
	if (n == 1) {	/* Must copy y info from x */
		*dy = *dx;
		GMT_inc_code[1] = GMT_inc_code[0];	/* Use exact inc codes for both x and y */
	}

	if (GMT_inc_code[0] & GMT_INC_IS_NNODES && GMT_inc_code[0] & GMT_INC_UNITS) {
		mexPrintf ("GMTMBGRID: ERROR: number of x nodes cannot have units\n");
		return (EXIT_FAILURE);
	}
	if (GMT_inc_code[1] & GMT_INC_IS_NNODES && GMT_inc_code[1] & GMT_INC_UNITS) {
		mexPrintf ("GMTMBGRID: ERROR: number of y nodes cannot have units\n");
		return (EXIT_FAILURE);
	}
	return (0);
}


int GMT_getincn (char *line, double inc[], int n) {
	int last, i, pos;
	char p[BUFSIZ];
	double scale = 1.0;

	/* Dechipers dx/dy/dz/dw/du/dv/... increment strings with n items */

	memset ((void *)inc, 0, (size_t)(n * sizeof (double)));

	i = pos = GMT_inc_code[0] = GMT_inc_code[1] = 0;

	while (i < n && (GMT_strtok (line, "/", &pos, p))) {
		last = strlen (p) - 1;
		if (p[last] == '=') {	/* Let -I override -R */
			p[last] = 0;
			if (i < 2) GMT_inc_code[i] |= GMT_INC_IS_EXACT;
			last--;
		}
		else if (p[last] == '+' || p[last] == '!') {	/* Number of nodes given, determine inc from domain (! added since documentation mentioned this once... */
			p[last] = 0;
			if (i < 2) GMT_inc_code[i] |= GMT_INC_IS_NNODES;
			last--;
		}
		switch (p[last]) {
			case 'm':
			case 'M':	/* Gave arc minutes */
				p[last] = 0;
				scale = GMT_MIN2DEG;
				break;
			case 'c':
			case 'C':	/* Gave arc seconds */
				p[last] = 0;
				scale = GMT_SEC2DEG;
				break;
			case 'e':
			case 'E':	/* Gave meters along mid latitude */
				p[last] = 0;
				if (i < 2) GMT_inc_code[i] |= GMT_INC_IS_M;
				break;
			case 'K':	/* Gave km along mid latitude */
			case 'k':
				p[last] = 0;
				if (i < 2) GMT_inc_code[i] |= GMT_INC_IS_KM;
				break;
			case 'I':	/* Gave miles along mid latitude */
			case 'i':
				p[last] = 0;
				if (i < 2) GMT_inc_code[i] |= GMT_INC_IS_MILES;
				break;
			case 'N':	/* Gave nautical miles along mid latitude */
			case 'n':
				p[last] = 0;
				if (i < 2) GMT_inc_code[i] |= GMT_INC_IS_NMILES;
				break;
			default:	/* No special flags or units */
				scale = 1.0;
				break;
		}
		if ( (sscanf(p, "%lf", &inc[i])) != 1) {
			mexPrintf ("GMTMBGRID ERROR: Unable to decode %s as a floating point number\n", p);
			return (EXIT_FAILURE);
		}
		inc[i] *= scale;
		i++;	/* Goto next increment */
	}

	return (i);	/* Returns the number of increments found */
}

int GMT_strtok (const char *string, const char *sep, int *pos, char *token) {
	/* Reentrant replacement for strtok that uses no static variables.
	 * Breaks string into tokens separated by one of more separator
	 * characters (in sep).  Set *pos to 0 before first call.  Unlike
	 * strtok, always pass the original string as first argument.
	 * Returns 1 if it finds a token and 0 if no more tokens left.
	 * pos is updated and token is returned.  char *token must point
	 * to memory of length >= strlen (string).
	 */

	int i, j, string_len;

	string_len = strlen (string);

	/* Wind up *pos to first non-separating character: */
	while (string[*pos] && strchr (sep, (int)string[*pos])) (*pos)++;

	token[0] = 0;	/* Initialize token to NULL in case we are at end */

	if (*pos >= string_len || string_len == 0) return 0;	/* Got NULL string or no more string left to search */

	/* Search for next non-separating character */
	i = *pos; j = 0;
	while (string[i] && !strchr (sep, (int)string[i])) token[j++] = string[i++];
	token[j] = 0;	/* Add terminating \0 */

	/* Wind up *pos to next non-separating character */
	while (string[i] && strchr (sep, (int)string[i])) i++;
	*pos = i;

	return 1;
}
