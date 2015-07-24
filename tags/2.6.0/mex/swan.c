/*    Copyright 1992, Mader Consulting Co., 1049 Kamehame Dr.,
 *    Honolulu, HI 96825, 808-396-9855

 *     dminx is longitude increment in minutes
 *     dminy is latitude increment in minutes
 *     cf is coriolis parameter
 *     cc is bottom stress of de chezy
 *     sfx is x wind stress
 *     sfy is y wind stress
 *     ROUGH is the roughness coefficient for FLOODING - Seldom Used
 *        surface roughness = actual height/ideal height
 *        1.0 for no roughness
 *     polar=0 for no polar, 1 for north lat, -1 for south lat
 *     polar=2 for Mercator
 *     cumint is cycle interval between points on cumulative plots
 *     ROUGH is the roughness coefficient for FLOODING - seldom used
 *     ind is 1 for piston, 2 for continuative , 3 for reflective
 *     in order of left,bottom,right,top of rectangle
 *     iopt is = 0 for extrapolated continuative boundary
 *     iopt is = 1 for linear continuative boundary
 *
 *
 *	Translated to C & mexified (+ some bug corrections) By
 *	Joaquim Luis - 2005

		20-10-2007 	(Killed gotos)
				removed write_grd_ascii() 		Blheght!!
				replaced several ij(i,j) occurences by its precomputed value
		22-10-2007	Added option -R
				Moved time_h += dt; to the end of the loop so that output naming start a time = 0
		04-11-2007	Fixed momentum ouput bug (negatives were being set to zero)
				Water heights < 5e-3 = 0
				Killed many gotos in bndy_()
		03-01-2008	Added output to ANUGA and MOST formats (netCDF)
		05-03-2008	Create empty (NaNs) global attribs to hold fault parameters in ANUGA format
		12-05-2008	SWW stage comes out with/without land as controlled by -L option
 		20-01-2010 	Add option to call aguentabar.dll via mexEvalString
				64 bits ready.
 *
 *	version WITH waitbar
 */

#include "mex.h"
#include "mwsize.h"
#include <netcdf.h>
#include <float.h>
#include <math.h>
#include <string.h>

#define	FALSE	0
#define	TRUE	1
#ifndef M_PI
#define M_PI	3.14159265358979323846
#endif
#define D2R		M_PI / 180.
#define LF		0x0A
#define CR		0x0D

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

/* Safe math macros that check arguments */

#define d_sqrt(x) ((x) < 0.0 ? 0.0 : sqrt (x))
#define d_acos(x) (fabs (x) >= 1.0 ? ((x) < 0.0 ? M_PI : 0.0) : acos (x))
#define d_asin(x) (fabs (x) >= 1.0 ? copysign (M_PI_2, (x)) : asin (x))

#define ij(i,j) ((i) + (j)*ip2)
#define ijs(i,j,n) ((i) + (j)*n)
#define ijc(i,j) ((i) + (j)*n_ptmar)

struct srf_header {		/* Surfer file header structure */
	char id[4];		/* ASCII Binary identifier (DSAA/DSBB) */
	short int nx;		/* Number of columns */
	short int ny;		/* Number of rows */
	double x_min;		/* Minimum x coordinate */
	double x_max;		/* Maximum x coordinate */
	double y_min;		/* Minimum y coordinate */
	double y_max;		/* Maximum y coordinate */
	double z_min;		/* Minimum z value */
	double z_max;		/* Maximum z value */
};

void no_sys_mem (char *where, mwSize n);
int count_col (char *line);
int read_grd_info_ascii (char *file, struct srf_header *hdr);
int read_header_bin (FILE *fp, struct srf_header *hdr);
int write_grd_bin(char *name, double x_min, double y_min, double x_inc, double y_inc, mwSize i_start, 
		mwSize j_start, mwSize i_end, mwSize j_end, mwSize nX, float *work);
int read_grd_ascii (char *file, struct srf_header *hdr, double *work);
int read_grd_bin (char *file, struct srf_header *hdr, double *work);
int read_params(char *file);
int read_maregs(char *file);
int count_n_maregs(char *file);
int uvh_(double *dep, double *r, double *rn, double *u, double *un, double *v, double *vn, double *h, 
		double *hn, double *dxp, double *cca, double m_per_deg);
int bndy_(double *dep, double *r, double *rn, double *u, double *un, double *v, double *vn, double *h, 
		double *hn, double *dxp);
void max_z (double *zm, double *h_bak);
void change (double *h, double *h_bak);
int decode_R (char *item, double *w, double *e, double *s, double *n);
int check_region (double w, double e, double s, double n);
double ddmmss_to_degree (char *text);
void err_trap(mwSize status);
void write_most_slice(mwSize *ncid_most, mwSize *ids_most, mwSize i_start, mwSize j_start, mwSize i_end, mwSize j_end,
		mwSize nX, float *work, double *h, double *dep, double *u, double *v, float *tmp,
		size_t *start, size_t *count);
int open_most_nc (char *basename, char *name_var, mwSize *ids, mwSize nx, mwSize ny, 
		double dtx, double dty, double xMinOut, double yMinOut);
int open_anuga_sww (char *fname_sww, mwSize *ids, mwSize i_start, mwSize j_start, mwSize i_end, mwSize j_end, 
		mwSize nX, double dtx, double dty, double *dep, double xMinOut, double yMinOut, 
		float z_min, float z_max);
void write_anuga_slice(mwSize ncid, mwSize z_id, mwSize i_start, mwSize j_start, mwSize i_end, mwSize j_end, mwSize nX, 
		float *work, double *h, double *dep, double *u, double *v, float *tmp,
		size_t *start, size_t *count, float *slice_range, mwSize idx, mwSize with_land);

mwSize	ip, jp, polar, indl, indb, indr, indt, iopt;
mwSize	ip1, jp1, ip2, jp2, grn, cumint, *lcum_p = NULL;
double	dx, dy, dt, cf, cc, sfx, sfy, time_h, rough;
double	pistal, pistbl, pistab, pistbb, pistar, pistbr, pistat, pistbt;
double	dangx, dangy, *anglt = NULL;
mwSize first_in_uvh = TRUE;
struct	srf_header hdr_b;
struct	srf_header hdr_f;

/* --------------------------------------------------------------------------- */
/* Matlab Gateway routine */

void mexFunction(mwSize nlhs, mxArray *plhs[], mwSize nrhs, const mxArray *prhs[]) {

	float	*work, dz, *ptr_mov_32, *mov_32, *time_p = NULL;
	float	work_min = FLT_MAX, work_max = -FLT_MAX;
	double	*inicial, *r, *rn, *u, *un, *v, *vn, *h, *h_bak, *hn, *zm;
	double	*dxp, *xpp, *cca, x, y, small = 1e-6, m_per_deg = 111317.1;
	double	x_min, y_min, dminx, dminy, dtx, dty, *head, *tmp;
	double	*dep, *dep1, *cum_p = NULL, cang, angltt, dumb;
	double	x_inc, y_inc, x_tmp, y_tmp;		/* Used in the maregs positiojn test */
	double	*ptr_wb; 		/* Pointer to be used in the aguentabar */
	double	dfXmin = 0.0, dfYmin = 0.0, dfXmax = 0.0, dfYmax = 0.0, xMinOut, yMinOut;
	double	time_jump = 0, time0, time_for_anuga, prc;
	mwSize	i, j, i2, ij, k, ix, jy, ncl, lcum, n_mareg, n_ptmar, cycle;
	mwSize	i_start, j_start, i_end, j_end;
	mwSize	w_bin = TRUE, cumpt = FALSE, error = FALSE;
	mwSize	first_anuga_time = TRUE, out_sww = FALSE, out_most = FALSE;
	mwSize	r_bin_b, r_bin_f, surf_level = TRUE, max_level = FALSE, water_depth = FALSE;
	mwSize	argc, n_arg_no_char = 0, nx, ny, dims[3], n_frames = 0;
	mwSize	ncid, ncid_most[3], z_id = -1, ids[13], ids_ha[6], ids_ua[6], ids_va[6], ids_most[3];
	mwSize	n_of_cycles = 1010;	/* Numero de ciclos a calcular */
	char	*bathy = NULL;		/* Name pointer for bathymetry file */
	char 	*params = NULL;		/* Name pointer for parameters file */
	char 	*hcum = NULL;		/* Name pointer for cumulative hight file */
	char 	*maregs = NULL;		/* Name pointer for maregraph positions file */
	char 	*fname_sww = NULL;	/* Name pointer for Anuga's .sww file */
	char 	*basename_most = NULL;	/* Name pointer for basename of MOST .nc files */
	char	*fonte = NULL;		/* Name pointer for tsunami source file */
	char	stem[80], prenome[128], cmd[16];
	char	**argv;
	unsigned char	*ptr_mov_8, *mov_8, *mov_8_tmp;
	mwSize	params_in_input = FALSE, bat_in_input = FALSE, source_in_input = FALSE;
	mwSize	write_grids = FALSE, movie = FALSE, movie_char = FALSE, movie_float = FALSE;
	mwSize	maregs_in_input = FALSE, out_velocity =	FALSE, out_momentum = FALSE, got_R = FALSE;
	mwSize	with_land = FALSE, IamCompiled = FALSE;
	mxArray *rhs[3];
	size_t	start0 = 0, count0 = 1, start1_A[2] = {0,0}, count1_A[2], start1_M[3] = {0,0,0}, count1_M[3];
	float	stage_range[2], xmom_range[2], ymom_range[2], *tmp_slice;
	FILE	*fp;

	movie_char = TRUE;	/* temporary */

	bathy = "bathy.grd";
	params = "swan.par";
	hcum = "maregs.dat";
	fonte = "fonte.grd";
	maregs = "mareg.xy";

	argc = nrhs;
	for (i = 0; i < nrhs; i++) {		/* Check input to find how many arguments are of type char */
		if(!mxIsChar(prhs[i])) {
			argc--;
			n_arg_no_char++;	/* Number of arguments that have a type other than char */
		}
	}

	if (n_arg_no_char >= 5) {		/* Not all combinations are tested */
		/* Bathymetry */
		dep1 = mxGetPr(prhs[0]);
		nx = mxGetN (prhs[0]);
		ny = mxGetM (prhs[0]);
		head  = mxGetPr(prhs[1]);	/* Get bathymetry header info */
		hdr_b.x_min = head[0];		hdr_b.x_max = head[1];
		hdr_b.y_min = head[2];		hdr_b.y_max = head[3];
		hdr_b.z_min = head[4];		hdr_b.z_max = head[5];
		hdr_b.nx = nx;			hdr_b.ny = ny;
		bat_in_input = TRUE;
		/* Tsunami source */
		inicial = mxGetPr(prhs[2]);
		nx = mxGetN (prhs[2]);
		ny = mxGetM (prhs[2]);
		head  = mxGetPr(prhs[3]);		/* Get bathymetry header info */
		hdr_f.x_min = head[0];		hdr_f.x_max = head[1];
		hdr_f.y_min = head[2];		hdr_f.y_max = head[3];
		hdr_f.z_min = head[4];		hdr_f.z_max = head[5];
		hdr_f.nx = nx;			hdr_f.ny = ny;
		source_in_input = TRUE;
		/* Now test if bat & source are compatible. If they are not, we go out right here. */
		if (hdr_f.nx != hdr_b.nx || hdr_f.ny != hdr_b.ny) {
			mexPrintf ("Bathymetry & Source grids have different nx and/or ny\n"); 
			mexPrintf ("%d %d %d %d\n", hdr_f.nx, hdr_b.nx, hdr_f.ny, hdr_b.ny); 
			return;
		}
		if ( fabs(hdr_f.x_min - hdr_b.x_min) > small || fabs(hdr_f.x_max - hdr_b.x_max) > small ||
			fabs(hdr_f.y_min - hdr_b.y_min) > small || fabs(hdr_f.y_max - hdr_b.y_max) > small ) {
			mexPrintf ("Bathymetry & Source grids do not cover the same region\n"); 
			mexPrintf ("%lf %lf %lf %lf\n", hdr_f.x_min, hdr_b.x_min, hdr_f.x_max, hdr_b.x_max); 
			mexPrintf ("%lf %lf %lf %lf\n", hdr_f.y_min, hdr_b.y_min, hdr_f.y_max, hdr_b.y_max); 
			return;
		}

		/* Now "read" the parameters */
		tmp = mxGetPr(prhs[4]);
		i = mxGetN (prhs[4]);
		if (i != 22) {
			mexPrintf("Params input argument has a wrong (=%d) number of arguments (should be 22)\n", i);
			return;
		}
		dt = tmp[0];		grn = (mwSize)tmp[1];	cf = tmp[2];
		cc = tmp[3];		sfx = tmp[4];		sfy = tmp[5];
		polar = (mwSize)tmp[6];	rough = tmp[7];		cumint = (mwSize)tmp[8];
		pistal = tmp[9];	pistbl = tmp[10];	pistab = tmp[11];
		pistbb = tmp[12];	pistar = tmp[13];	pistbr = tmp[14];
		pistat = tmp[15];	pistbt = tmp[16];	indl = (mwSize)tmp[17];
		indb = (mwSize)tmp[18];	indr = (mwSize)tmp[19];	indt = (mwSize)tmp[20];
		iopt = (mwSize)tmp[21];
		params_in_input = TRUE;
	}
	if(n_arg_no_char == 6) {		/* A maregraph vector was given as the sixth argument*/
		tmp = mxGetPr(prhs[5]);
		n_mareg = mxGetM(prhs[5]);
		dx = head[7];		dy = head[8];
		lcum_p = (mwSize *) mxCalloc ((size_t)(n_mareg), sizeof(mwSize));
		for (i = 0; i < n_mareg; i++) {
			x = tmp[i];		y = tmp[i+n_mareg];	/* Matlab vectors are stored by columns */
			lcum_p[i] = (irint((y - hdr_b.y_min) / dy) ) * hdr_b.nx + irint((x - hdr_b.x_min) / dx);
		}
		maregs_in_input = TRUE;
		cumpt = TRUE;
	}

	/* get the length of the input string */
	argv = (char **)mxCalloc(argc, sizeof(char *));
	for (i = 0; i < argc; i++)
		argv[i] = (char *)mxArrayToString(prhs[i+n_arg_no_char]);


	for (i = 0; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {
				case 'e':
					IamCompiled = TRUE;
					break;
				case 'f':	/* Movie */
					movie = TRUE;
					movie_char = FALSE;
					movie_float = TRUE;
					break;
				case 'm':	/* Movie */
					movie = TRUE;
					break;
				case 'n':	/* Write MOST files (*.nc) */
					basename_most  = &argv[i][2];
					out_most = TRUE;
					break;
				case 'A':	/* Name for the Anuga .sww netCDF file */
					fname_sww  = &argv[i][2];
					out_sww = TRUE;
					break;
				case 'B':	/* File with batymetry */
					bathy  = &argv[i][2];
					break;
				case 'D':
					water_depth = TRUE;
					surf_level = FALSE;
					break;
				case 'F':	/* File with the source grid */
					fonte  = &argv[i][2];
					break;
				case 'G':	/* Write grids at grn (see params) intervals */
					write_grids = TRUE;
					strcpy (stem, &argv[i][2]);
					break;
				case 'J':
					sscanf (&argv[i][2], "%lf", &time_jump);
					break;
				case 'L':	/* Output land nodes in SWW file */ 
					with_land = TRUE;
					break;
				case 'M':
					max_level = TRUE;
					surf_level = FALSE;
					break;
				case 'N':	/* Number of cycles to compute */
					n_of_cycles = atoi(&argv[i][2]);
					break;
				case 'O':	/* File name for maregraph data */
					hcum  = &argv[i][2];
					break;
				case 'R':
					error += decode_R (argv[i], &dfXmin, &dfXmax, &dfYmin, &dfYmax);
					got_R = TRUE;
					break;
				case 's':	/* Output velocity grids */ 
					out_velocity = TRUE;
					break;
				case 'S':	/* Output momentum grids */ 
					out_momentum = TRUE;
					break;
				case 'T':	/* File with maregraph positions */
					maregs  = &argv[i][2];
					cumpt = TRUE;
					maregs_in_input = FALSE;
					break;
				default:
					error = TRUE;
					break;
			}
		}
	}

	if (argc == 0 || error) {
		mexPrintf ("SWAN - Um gerador de tsunamis\n\n");
		mexPrintf( "usage: swan(bat,hdr_bat,deform,hdr_deform,params, [maregs], [-G|A|n<name>], [-B<bathy>] [-F<fonte>] [-M] [-N<n_cycles>] [-Rw/e/s/n] [-S] [-s] [-T<mareg>] [-D]\n");
		mexPrintf ("\t-m outputs a 3D grid used to do a movie\n");
		mexPrintf ("\t-A name of ANUGA file\n");
		mexPrintf ("\t-n basename for MOST triplet files (no extension)\n");
		mexPrintf ("\t-B name of bathymetry file. In case it was not transmited in input.\n");
		mexPrintf ("\t-D write grids with the total water depth\n");
		mexPrintf ("\t-F name of source file (default fonte.grd)\n");
		mexPrintf ("\t-G<stem> write grids at the grn intervals. Append file prefix. Files will be called <stem>#.grd\n");
		mexPrintf ("\t-M write grids of max water level [Default wave surface level]\n");
		mexPrintf ("\t-N number of cycles [Default 1010].\n");
		mexPrintf ("\t-R output grids only in the sub-region enclosed by <west/east/south/north>\n");
		mexPrintf ("\t-s write grids with the velocity. Grid names are appended with _U and _V sufixes.\n");
		mexPrintf ("\t-S write grids with the momentum. i.e velocity times water depth.\n");
		mexPrintf ("\t-T name of maregraph file (default mareg.xy)\n");
		mexPrintf("\t-e To be used from the Mirone stand-alone version.\n");
		return;
	}

	if (error) return;

	if (!params_in_input) 		/* If params where not given, read them from file */
		read_params(params);

	if (!bat_in_input) {			/* If bathymetry & source where not given as arguments, load them */
		r_bin_b = read_grd_info_ascii (bathy, &hdr_b);	/* Para saber como alocar a memoria */
		r_bin_f = read_grd_info_ascii (fonte, &hdr_f);	/* e verificar se as grelhas sao compativeis */
		if (r_bin_b < 0 || r_bin_f < 0) {
			mexPrintf ("Invalid grid. Possibly it is in the Surfer 7 format\n"); 
			error++;
		}

		if (hdr_f.nx != hdr_b.nx || hdr_f.ny != hdr_b.ny) {
			mexPrintf ("Bathymetry and source grids have different rows/columns\n"); 
			mexPrintf ("%d %d %d %d\n", hdr_f.nx, hdr_b.nx, hdr_f.ny, hdr_b.ny); 
			error++;
		}
		if (hdr_f.x_min != hdr_b.x_min || hdr_f.x_max != hdr_b.x_max ||
			hdr_f.y_min != hdr_b.y_min || hdr_f.y_max != hdr_b.y_max) {
			mexPrintf ("Bathymetry and source grids do not cover the same region\n"); 
			mexPrintf ("%lf %lf %lf %lf\n", hdr_f.x_min, hdr_b.x_min, hdr_f.x_max, hdr_b.x_max); 
			mexPrintf ("%lf %lf %lf %lf\n", hdr_f.y_min, hdr_b.y_min, hdr_f.y_max, hdr_b.y_max); 
			error++;
		}
		if (error) return;
	}

	if ( !(write_grids || out_sww || out_most || movie || cumpt) ) {
		mexPrintf("Nothing selected for output (grids, movie or mregraphs), exiting\n");
		return;
	}
	if ( water_depth && (out_sww || out_most) ) {
		water_depth = FALSE;
		mexPrintf("SWAN WARNING: Total water option is not compatible with ANUGA|MOST outputs. Ignoring\n");
	}
	if (out_momentum && (out_sww || out_most)) out_momentum = FALSE;
	if (out_velocity && (out_sww || out_most)) out_velocity = FALSE;

	/* dminx and dminy must be in minutes, but dx and dy in meters */
	if (polar == 0) m_per_deg = 1.;
	dminx = (hdr_b.x_max - hdr_b.x_min) / (hdr_b.nx - 1) * 60.;
	dminy = (hdr_b.y_max - hdr_b.y_min) / (hdr_b.ny - 1) * 60.;
	dx = (hdr_b.x_max - hdr_b.x_min) / (hdr_b.nx - 1) * m_per_deg;
	dy = (hdr_b.y_max - hdr_b.y_min) / (hdr_b.ny - 1) * m_per_deg;
	x_min = hdr_b.x_min;		y_min = hdr_b.y_min;
	ip2 = hdr_b.nx;	jp2 = hdr_b.ny;
	ip = ip2 - 2;	jp = jp2 - 2;	/* in original fortran version ip2 = ip + 2 */
	ip1 = ip + 1;	jp1 = jp + 1;
	ncl = ip2 * jp2;
	if (cumpt && !maregs_in_input)
		n_mareg = count_n_maregs(maregs);	/* Count maragraphs number */

	if (cumpt) {
		n_ptmar = n_of_cycles / cumint + 1;
		if ((fp = fopen (hcum, "w")) == NULL) {
			mexPrintf ("%s: Unable to create file %s - exiting\n", "swan", hcum);
			return;
		}
	}

	/* Allocate memory	*/
	if (!bat_in_input) {	/* Bathymetry & Source will be read from files. So need to alloc */
		if ((inicial = (double *) mxMalloc ((size_t)(ncl) * sizeof(double)) ) == NULL) 
			{no_sys_mem("swan --> (inicial)", ncl);	return;}
	}
	if ((dep = (double *) mxCalloc ((size_t)(ncl),	sizeof(double)) ) == NULL) 
		{no_sys_mem("swan --> (dep)", ncl);	return;}
	if ((work = (float *) mxCalloc ((size_t)(ncl),	sizeof(float)) ) == NULL) 
		{no_sys_mem("swan --> (work)", ncl);	return;}
	if ((r = (double *) mxCalloc ((size_t)(ncl),	sizeof(double)) ) == NULL) 
		{no_sys_mem("swan --> (r)", ncl);	return;}
	if ((rn = (double *) mxCalloc ((size_t)(ncl),	sizeof(double)) ) == NULL) 
		{no_sys_mem("swan --> (rn)", ncl);	return;}
	if ((u = (double *) mxCalloc ((size_t)(ncl),	sizeof(double)) ) == NULL) 
		{no_sys_mem("swan --> (u)", ncl);	return;}
	if ((un = (double *) mxCalloc ((size_t)(ncl),	sizeof(double)) ) == NULL) 
		{no_sys_mem("swan --> (un)", ncl);	return;}
	if ((v = (double *) mxCalloc ((size_t)(ncl),	sizeof(double)) ) == NULL) 
		{no_sys_mem("swan --> (v)", ncl);	return;}
	if ((vn = (double *) mxCalloc ((size_t)(ncl),	sizeof(double)) ) == NULL) 
		{no_sys_mem("swan --> (vn)", ncl);	return;}
	if ((h = (double *) mxCalloc ((size_t)(ncl),	sizeof(double)) ) == NULL) 
		{no_sys_mem("swan --> (h)", ncl);	return;}
	if ((hn = (double *) mxCalloc ((size_t)(ncl),	sizeof(double)) ) == NULL) 
		{no_sys_mem("swan --> (hn)", ncl);	return;}
	if ((dxp = (double *) mxCalloc ((size_t)(ncl+1),	sizeof(double)) ) == NULL) 
		{no_sys_mem("swan --> (dxp)", ncl);	return;}
	if ((anglt = (double *) mxCalloc ((size_t)(jp2),	sizeof(double)) ) == NULL) 
		{no_sys_mem("swan --> (anglt)", jp2);	return;}
	/*if (polar != 0) {
		if ((xpp = (double *) mxCalloc ((size_t)(ncl+1), sizeof(double)) ) == NULL) 
			{no_sys_mem("swan --> (xpp)", ncl);	return;}
	} */
	if (max_level) {
		if ((h_bak = (double *) mxCalloc ((size_t)(ncl), sizeof(double)) ) == NULL) 
			{no_sys_mem("swan --> (h_bak)", ncl);	return;}
		if ((zm = (double *) mxCalloc ((size_t)(ncl), sizeof(double)) ) == NULL) 
			{no_sys_mem("swan --> (zm)", ncl);	return;}
	}
	if (cumpt && !maregs_in_input) {
		if ((lcum_p = (mwSize *) mxCalloc ((size_t)(n_mareg), sizeof(mwSize)) ) == NULL) 
			{no_sys_mem("swan --> (lcum_p)", n_mareg);	return;}
		if ((cum_p = (double *) mxCalloc ((size_t)(n_mareg*n_ptmar), sizeof(double)) ) == NULL) 
			{no_sys_mem("swan --> (cum_p)", n_mareg*n_ptmar);	return;}
		if ((time_p = (float *) mxCalloc ((size_t)(n_ptmar), sizeof(float)) ) == NULL) 
			{no_sys_mem("swan --> (time_p)", n_ptmar);	return;}
	}
	if (cumpt && maregs_in_input) {		/* lcum_p is already allocated and filled */
		if ((cum_p = (double *) mxCalloc ((size_t)(n_mareg*n_ptmar), sizeof(double)) ) == NULL) 
			{no_sys_mem("swan --> (cum_p)", n_mareg*n_ptmar);	return;}
		if ((time_p = (float *) mxCalloc ((size_t)(n_ptmar), sizeof(float)) ) == NULL) 
			{no_sys_mem("swan --> (time_p)", n_ptmar);	return;}
	}

	if (cumpt && !maregs_in_input) read_maregs(maregs);	/* Read maregraph locations */
	if (!bat_in_input) {			/* If bathymetry & source where not given as arguments, load them */
		if (!r_bin_b)					/* Read bathymetry */
			read_grd_ascii (bathy, &hdr_b, dep);
		else
			read_grd_bin (bathy, &hdr_b, dep);
		if (r_bin_b)					/* Read source */
			read_grd_bin (fonte, &hdr_f, inicial);
		else
			read_grd_ascii (fonte, &hdr_f, inicial);
	}

	nx = hdr_b.nx;			/* nx has been redeclared before (for counting params) */
	if (bat_in_input) {		/* If bathymetry & source where given as arguments */
		/* Transpose from Matlab orientation to scanline orientation */
		for (i = 0; i < ny; i++)
			for (j = 0; j < nx; j++) dep[i*nx+j] = -dep1[j*ny+i];
	}
	/* Transpose the inicial array from Matlab orientation to scanline orientation */
	for (i = 0; i < ny; i++) {
		for (j = 0; j < nx; j++) {
			dumb = inicial[j*ny+i];		ij = i*nx+j;
			h[ij] = dumb;
			hn[ij] = dumb;
		}
	}

	/* Test that the maregraphs coords -> indeces conversion is correct */
	/*mexPrintf("Adjusted maregraph positions\n");
	mexPrintf("old x\told y\tnew x\tnew y\tdepth\n");
	x_inc = (hdr_b.x_max - hdr_b.x_min) / (ip2 - 1);
	y_inc = (hdr_b.y_max - hdr_b.y_min) / (jp2 - 1);
	for (i = 0; i < n_mareg; i++) {
		x = tmp[i];		y = tmp[i+n_mareg];
		ix = irint((x - hdr_b.x_min) / x_inc);
		jy = irint((y - hdr_b.y_min) / y_inc); 
		x_tmp = hdr_b.x_min + x_inc * ix;	// Adjusted x maregraph pos
		y_tmp = hdr_b.y_min + y_inc * jy;	// Adjusted y maregraph pos
		mexPrintf("%.4f\t%.4f\t%.4f\t%.4f\t%.1f\n",x,y,x_tmp,y_tmp,-dep[lcum_p[i]]);
	}*/

	dangx = dminx / 60.;	dangy = dminy / 60.;
	/*     Polar option (If Polar is 2 then MERCATOR)*/
	if (polar != 0) {
		anglt[0] = y_min;
		for (i = 0; i < jp1; i++) {
			if (polar > 0) anglt[i+1] = anglt[i] + dangy;
			if (polar < 0) anglt[i+1] = anglt[i] - dangy;
			if (polar < 0 && anglt[i+1] < 0.) {
				polar = 1;
				anglt[i+1] = anglt[i] + dangy;
			}
			else {
				angltt = anglt[i] * D2R;
				cang = cos(angltt);
				if (polar == 2) anglt[i+1] = anglt[i] + dangy * cang;
			}
		}
		/*    make dxp array */
		for (j = i = 0; j < jp2; j++) {
			cang = cos(anglt[j] * D2R);
			for (k = 0; k < ip2; k++) {
				dxp[i+1] = dangx * m_per_deg * cang;
				/*xpp[i+1] = xpp[i] + dxp[i+1];*/
				i++;
			}
			/*xpp[i] = 0.;*/
		}
		dxp[0] = dxp[1];	 /* dxp(0) needs to be set to non-zero */
		/*mxFree ((void *) xpp);*/
	}
	else {
		for (i = 0; i < ncl; i++)
			dxp[i] = dx;		/* Do not set if polar */
	}

	if ( cc != 0.) {
		if ((cca = (double *) mxCalloc ((size_t)(ncl),	sizeof(double)) ) == NULL) 
			{no_sys_mem("swan --> (cca)", ncl);	return;}
		for (i = 0; i < ncl; i++)
			cca[i] = cc;
	}

	dtx = (polar == 0) ? dx : dangx;	/* If polar == 0 dx and dy are already in meters */
	dty = (polar == 0) ? dy : dangy;	/* like they must be. Otherwise, they are in degrees */
	lcum = 0;	cycle = 1;	time_h = 0.;

	if (!IamCompiled) {
		rhs[0] = mxCreateDoubleScalar(0.0);
		ptr_wb = mxGetPr(rhs[0]);
		rhs[1] = mxCreateString("title");
		rhs[2] = mxCreateString("Sit and Wait ...");
		mexCallMATLAB(0,NULL,3,rhs,"aguentabar");
	}
	else
		mexEvalString("aguentabar(0,\"title\",\"Sit and Wait...\")");

	/* ---------------- Declarations for the (if) movie option ------------------ */
	if (movie && movie_char) {
		mov_8 = mxCalloc(ncl, sizeof(char));
		mov_8_tmp = mxCalloc(ncl, sizeof(char));
		dims[0] = ny;	dims[1] = nx;	dims[2] = (mwSize)(n_of_cycles / grn + 2);
		plhs[0] = mxCreateNumericArray(3, dims, mxUINT8_CLASS, mxREAL);
		ptr_mov_8 = (unsigned char *)mxGetData(plhs[0]);
	}
	else if (movie && movie_float) {
		mov_32 = (float *)mxCalloc(ncl, sizeof(float));
		dims[0] = ny;	dims[1] = nx;	dims[2] = (mwSize)(n_of_cycles / grn + 2);
		plhs[0] = mxCreateNumericArray(3, dims, mxSINGLE_CLASS, mxREAL);
		ptr_mov_32 = (float *)mxGetData(plhs[0]);
	}

	/* ----------------- Compute vars to use if write grids --------------------- */
	if (!got_R && (write_grids || out_velocity || out_momentum || out_sww || out_most) ) {	
		/* Write grids over the whole region */
		i_start = 0;		i_end = ip2;
		j_start = 0;		j_end = jp2;
		xMinOut = x_min;	yMinOut = y_min;
	}
	else if (got_R && (write_grids || out_velocity || out_momentum || out_sww || out_most) ) {	
		/* Write grids in sub-region */
		i_start = irint((dfXmin - hdr_b.x_min) / dtx);
		j_start = irint((dfYmin - hdr_b.y_min) / dty); 
		i_end   = irint((dfXmax - hdr_b.x_min) / dtx) + 1;
		j_end   = irint((dfYmax - hdr_b.y_min) / dty) + 1;
		xMinOut = hdr_b.x_min + dtx * i_start;	/* Adjustes xMin|yMin to lay on the closest grid node */
		yMinOut = hdr_b.y_min + dty * j_start;
	}
	/* ---------------------------------------------------------------------- */

	if (out_sww) {
		/* ----------------- Open a ANUGA netCDF file for writing --------------- */
		ncid = open_anuga_sww (fname_sww, ids, i_start, j_start, i_end, j_end, ip2,
		dtx, dty, dep, xMinOut, yMinOut, (float)hdr_b.z_min, (float)hdr_b.z_max);
		if (ncid == -1)
			mexErrMsgTxt ("SWAN: failure to create ANUGA SWW file.\n");

		/* To be used when writing the data slices */
		count1_A[0] = 1;	count1_A[1] = (i_end - i_start) * (j_end - j_start);

		stage_range[0] = xmom_range[0] = ymom_range[0] = FLT_MAX;
		stage_range[1] = xmom_range[1] = ymom_range[1] = -FLT_MIN;

		tmp_slice = (float *) mxMalloc (sizeof(float) * (nx * ny));	/* To use inside slice writing */
	}

	if (out_most) {
		/* ----------------- Open a 3 MOST netCDF files for writing ------------- */
		mwSize nx, ny;
		nx = i_end - i_start;		ny = j_end - j_start;
		ncid_most[0] = open_most_nc (basename_most, "HA", ids_ha, nx, ny, dtx, dty, xMinOut, yMinOut);
		ncid_most[1] = open_most_nc (basename_most, "UA", ids_ua, nx, ny, dtx, dty, xMinOut, yMinOut);
		ncid_most[2] = open_most_nc (basename_most, "VA", ids_va, nx, ny, dtx, dty, xMinOut, yMinOut);

		if (ncid_most[0] == -1 || ncid_most[1] == -1 || ncid_most[2] == -1)
			mexErrMsgTxt ("SWAN: failure to create one or more of the MOST files\n");

		ids_most[0] = ids_ha[5];	/* IDs of the Amp, Xmom & Ymom vriables */
		ids_most[1] = ids_ua[5];
		ids_most[2] = ids_va[5];

		/* To be used when writing the data slices */
		count1_M[0] = 1;	count1_M[1] = ny;	count1_M[2] = nx;

		if (!out_sww)
			tmp_slice = (float *)mxMalloc (sizeof(float) * (nx * ny));	/* To use inside slice writing */ 
	}

	/* ---------------------------------------------------------------------- */
	if (time_jump == 0) time_jump = -1;	/* Trick to allow writing zero time grids when jump was not demanded */

	for (k = 0; k < n_of_cycles; k++) {
		if (cycle % 10 == 0) {
			prc = (double)cycle/(double)n_of_cycles;
			if (!IamCompiled) {
				ptr_wb[0] = prc;
				mexCallMATLAB(0,NULL,1,rhs,"aguentabar");
			}
			else {
				sprintf(cmd, "aguentabar(%f)",prc);
				mexEvalString(cmd);
			}
		}

		uvh_(dep, r, rn, u, un, v, vn, h, hn, dxp, cca, m_per_deg);

		if (max_level) {
			change(h, h_bak);	max_z(zm, h_bak);
		}
		if (cumpt) {			/* Want time series at maregraph positions */
			if (cycle % cumint == 0) {	/* Save heights at cumint intervals */
				for (i = 0; i < n_mareg; i++) {
					cum_p[ijc(lcum,i)] = h[lcum_p[i]-1];
					time_p[lcum] = (float)time_h;
				}
				lcum++;
			}
		}
		if ( time_h > time_jump && ( (k % grn) == 0 || k == n_of_cycles - 1) ) {
			if (surf_level) {
				for (i = 0; i < ncl; i++)
					work[i] = (float) h[i];
			}
			else if (max_level) {	/* Max surface level */
				for (i = 0; i < ncl; i++)
					work[i] = (float) zm[i];
			}
			else if (water_depth) {
				for (i = 0; i < ncl; i++)
					if ((work[i] = (float) (h[i] + dep[i])) < 0.) work[i] = 0.;
			}
			if (movie && movie_char) {
				/* Output the "movie" in a 3D uint8 matrix */
				for (i = 0; i < ncl; i++) {
					work_max = MAX(work[i], work_max);
					work_min = MIN(work[i], work_min);
				}
				dz = work_max - work_min;
				for (i = 0; i < ncl; i++)
					mov_8_tmp[i] = (unsigned char)( ((work[i] - work_max) / dz) * 254 + 1);

				/* Transpose to matlab ordering */
				for (i = 0, i2 = ny - 1; i < ny; i++, i2--)
				for (j = 0; j < nx; j++) mov_8[j*ny+i2] = mov_8_tmp[i*nx+j];

				memcpy(ptr_mov_8, mov_8, ncl);
				ptr_mov_8 += ncl;
				n_frames++;
			}
			else if (movie && movie_float) {
				/* Output the "movie" in a 3D float matrix */
				/* Transpose to matlab ordering */
				for (i = 0, i2 = ny - 1; i < ny; i++, i2--)
				for (j = 0; j < nx; j++) mov_32[j*ny+i2] = work[i*nx+j];

				memcpy(ptr_mov_32, mov_32, ncl*4);
				ptr_mov_32 += ncl;
				n_frames++;
			}
			if (write_grids) {
				if (stem[0] == 0)
					sprintf (prenome,"%.5d.grd", irint(time_h) );
				else
					sprintf (prenome, "%s%.5d.grd", stem, irint(time_h) );
				write_grd_bin( prenome, xMinOut, yMinOut, dtx, dty, i_start, j_start, i_end, j_end, ip2, work);
			}
			if (out_velocity) {
				if (stem[0] == 0)
					sprintf (prenome,"%.5d", irint(time_h) );
				else
					sprintf (prenome, "%s%.5d", stem, irint(time_h) );
				for (i = 0; i < ncl; i++) work[i] = (float) u[i];
				write_grd_bin(strcat(prenome,"_U.grd"), xMinOut, yMinOut, dtx, dty, 
						i_start, j_start, i_end, j_end, ip2, work);
				for (i = 0; i < ncl; i++) work[i] = (float) v[i];
				prenome[strlen(prenome) - 6] = '\0';	/* Remove the _U.grd' so that we can add '_V.grd' */
				write_grd_bin( strcat(prenome,"_V.grd"), xMinOut, yMinOut, dtx, dty, 
						i_start, j_start, i_end, j_end, ip2, work);
			}
			if (out_momentum) {
				if (stem[0] == 0)
					sprintf (prenome,"%.5d", irint(time_h) );
				else
					sprintf (prenome, "%s%.5d", stem, irint(time_h) );

				if (water_depth)	/* "work" is already the water depth */ 
					for (i = 0; i < ncl; i++) work[i] = (float) (u[i] * work[i]);
				else
					for (i = 0; i < ncl; i++) {
						if (( work[i] = (float)(h[i] + dep[i]) ) < 0.) work[i] = 0.;
						work[i] *= (float)u[i];
					}

				write_grd_bin( strcat(prenome,"_Uh.grd"), xMinOut, yMinOut, dtx, dty, 
						i_start, j_start, i_end, j_end, ip2, work);

				if (water_depth)
					for (i = 0; i < ncl; i++) work[i] = (float) (v[i] * work[i]);
				else
					for (i = 0; i < ncl; i++) {
						if (( work[i] = (float)(h[i] + dep[i]) ) < 0.) work[i] = 0.;
						work[i] *= (float)v[i];
					}

				prenome[strlen(prenome) - 7] = '\0';	/* Remove the _Uh.grd' so that we can add '_Vh.grd' */
				write_grd_bin( strcat(prenome,"_Vh.grd"), xMinOut, yMinOut, dtx, dty, 
						i_start, j_start, i_end, j_end, ip2, work);
			}

			if (out_sww) {
				if (first_anuga_time) {
					time0 = time_h;
					first_anuga_time = FALSE;
				}
				time_for_anuga = time_h - time0;	/* I think ANUGA wants time starting at zero */
				err_trap (nc_put_vara_double (ncid, ids[6], &start0, &count0, &time_for_anuga));

				write_anuga_slice(ncid, ids[7], i_start, j_start, i_end, j_end, ip2, work,
						h, dep, u, v, tmp_slice, start1_A, count1_A, stage_range, 1, with_land);
				write_anuga_slice(ncid, ids[9], i_start, j_start, i_end, j_end, ip2, work,
						h, dep, u, v, tmp_slice, start1_A, count1_A, xmom_range, 2, with_land);
				write_anuga_slice(ncid, ids[11],i_start, j_start, i_end, j_end, ip2, work,
						h, dep, u, v, tmp_slice, start1_A, count1_A, ymom_range, 3, with_land);

				start1_A[0]++;		/* Increment for the next slice */
			}

			if (out_most) {
				/* Here we'll use the start0 computed above */
				err_trap (nc_put_vara_double (ncid_most[0], ids_ha[4], &start0, &count0, &time_h));
				err_trap (nc_put_vara_double (ncid_most[1], ids_ua[4], &start0, &count0, &time_h));
				err_trap (nc_put_vara_double (ncid_most[2], ids_va[4], &start0, &count0, &time_h));

				write_most_slice(ncid_most, ids_most, i_start, j_start, i_end, j_end, ip2, 
					work, h, dep, u, v, tmp_slice, start1_M, count1_M);
				start1_M[0]++;		/* Increment for the next slice */
			}

			start0++;			/* Only used with netCDF formats */
		}
		time_h += dt;
		cycle++;
	}

	if (out_sww) {		/* Uppdate range values and close SWW file */
		err_trap (nc_put_var_float (ncid, ids[8], stage_range));
		err_trap (nc_put_var_float (ncid, ids[10], xmom_range));
		err_trap (nc_put_var_float (ncid, ids[12], ymom_range));
		err_trap(nc_close (ncid)); 
	}

	if (out_most) {		/* Close MOST files */
		err_trap(nc_close (ncid_most[0])); 
		err_trap(nc_close (ncid_most[1])); 
		err_trap(nc_close (ncid_most[2])); 
	}
	if (out_sww || out_most) mxFree ((void *)tmp_slice);

	if (!IamCompiled) {
		*ptr_wb = 1.0;
		mexCallMATLAB(0,NULL,1,rhs,"aguentabar");
	}
	else
		mexEvalString("aguentabar(1.0)");

	/*     WRITE THE OUTPUT ON FILE CUMHT */
	if (cumpt) {
		for (j = 0; j < lcum; j++) {
			fprintf (fp, "%.3f", time_p[j]);
			for (i = 0; i < n_mareg; i++)
				fprintf (fp, "\t%.4f", cum_p[ijc(j,i)]);
			fprintf (fp, "\n");
		}
	}
	/* Clean up allocated memory. */
	mxDestroyArray(rhs[0]);		mxDestroyArray(rhs[1]);		mxDestroyArray(rhs[2]);
	if (movie && movie_char) {
		mxFree((void *)mov_8);	mxFree((void *)mov_8_tmp);
	}

	if (movie && movie_float) mxFree((void *)mov_32);

	if (cumpt) fclose (fp);
	if (!source_in_input)
		mxFree ((void *) inicial);

	mxFree ((void *) work);	mxFree ((void *) dep);	mxFree ((void *) r);
	mxFree ((void *) h);	mxFree ((void *) rn);	mxFree ((void *) u);
	mxFree ((void *) un);	mxFree ((void *) vn);	mxFree ((void *) hn);
	mxFree ((void *) dxp);	mxFree ((void *) anglt);
	if (max_level) {
		mxFree ((void *) zm);	mxFree ((void *) h_bak); 
	}
	if (cumpt) {
		mxFree ((void *) lcum_p);	mxFree((void *) cum_p);	mxFree ((void *) time_p);	 
	}
	if (cc != 0.) mxFree ((void *) cca);
}

/* --------------------------------------------------------------------------- */
void err_trap(mwSize status) {
	if (status != NC_NOERR)	
		mexPrintf ("swan: error and errorcode = %d\n", status);
}

/* --------------------------------------------------------------------------- */
int bndy_(double *dep, double *r, double *rn, double *u, double *un, double *v, double *vn,
	  double *h, double *hn, double *dxp) {

	double tmp, mul;
	mwSize i, j, ij_0j, ij_1j, ij_2j, ij_i0, ij_i1, ij_i2;
	mwSize ij_ip21_j, ij_ip11_j, ij_i_jp21, ij_i_jp11;

	/*     left boundary */
	for (j = 1; j < jp2; j++) {
		ij_0j = ij(0,j);	ij_1j = ij(1,j);	ij_2j = ij(2,j);
		tmp = d_sqrt (9.8 * dep[ij_0j]);
		if (indl == 1) { 	/*      piston */
			u[ij_0j] = pistal * sin(pistbl * time_h);
			v[ij_0j] = v[ij_1j];
			if (dep[ij_0j] > 0.)
				h[ij_0j] = u[ij_0j] * tmp / 9.8;
			continue;
		}
		else if (indl == 3) { 	/*     reflective */
			u[ij_0j] = 0.;
			v[ij_0j] = v[ij_1j];
			h[ij_0j] = h[ij_1j];
			continue;
		}

		/*     continuative */
		if (polar != 0) dx = dxp[ij_0j];
		if (polar == 2) dy = dxp[ij_0j];	/* MERCATOR */
		/*if (dep[ij(0,j)] < 0.) goto L313;	/* This is bad */
		if (fabs(dep[ij_0j]) > 0.005) dep[ij_0j] = 0.;	/* THIS SEAMS TO PREVENT BORDER INSTABILITIES */
		mul = tmp * dt / dx;
		v[ij_0j] += (v[ij_1j] - v[ij_0j]) * mul;
		if (v[ij_1j] == v[ij_2j])
			v[ij_0j] = v[ij_1j];
		u[ij_0j] += (u[ij_1j] - u[ij_0j]) * mul;
		if (u[ij_1j] == u[ij_2j])
			u[ij_0j] = u[ij_1j];
		h[ij_0j] += (h[ij_1j] - h[ij_0j]) * mul;
		if (h[ij_1j] == h[ij_2j])		/* More useless comparison between floats */
			h[ij_0j] = h[ij_1j];
		if (iopt == 0) continue;
/*L313:*/
		v[ij_0j] = v[ij(1,j+0)];	/* Those 3 would blow te program (v[ij(1,jp2)] doesn't exist) */
		u[ij_0j] = u[ij(1,j+0)];
		h[ij_0j] = h[ij(1,j+0)];
	}
	/*     bottom boundary */
	for (i = 0; i < ip2; i++) {
		ij_i0 = ij(i,0);	ij_i1 = ij(i,1);	ij_i2 = ij(i,2);
		tmp = d_sqrt (9.8 * dep[ij_i0]);
		if (indb == 1) { 	/*     piston */
			v[ij_i0] = pistab * sin(pistbb * time_h);
			u[ij_i0] = u[ij_i1];
			if (dep[ij_i0] > 0.)
				h[ij_i0] = v[ij_i0] * tmp / 9.8;
			continue;
		}
		else if (indb == 3) { 	/*     reflective */
			v[ij_i0] = 0.;
			u[ij_i0] = u[ij_i1];
			h[ij_i0] = h[ij_i1] * 2. - h[ij_i2];
			continue;
		}

		/*     continuative */
		if (polar != 0) dx = dxp[ij_i0];	/* THIS WAS NOT ON THE ORIGINAL CODE. FORGOTTEN? */
		if (polar == 2) dy = dxp[ij_i0];	/* MERCATOR */
		/*if (dep[ij(i,0)] < 0.) goto L323;	/* This is bad */
		if (fabs(dep[ij_i0]) > 0.005) dep[ij_i0] = 0.;	/* THIS SEAMS TO PREVENT BORDER INSTABILITIES */
		mul = tmp * dt / dy;
		v[ij_i0] += (v[ij_i1] - v[ij_i0]) * mul;
		if (v[ij_i1] == v[ij_i2]) 
			v[ij_i0] = v[ij_i1];
		u[ij_i0] += (u[ij_i1] - u[ij_i0]) * mul;
		if (u[ij_i1] == u[ij_i2]) 
			u[ij_i0] = u[ij_i1];
		h[ij_i0] += (h[ij_i1] - h[ij_i0]) * mul;
		if (h[ij_i1] == h[ij_i2]) 
			h[ij_i0] = h[ij_i1];
		if (iopt == 0) continue;
/*L323:*/
		v[ij_i0] = v[ij_i1];
		u[ij_i0] = u[ij_i1];
		h[ij_i0] = h[ij_i1];
	}
	/*     right boundary */
	for (j = 1; j < jp2; j++) {
		ij_ip21_j = ij(ip2-1,j);	ij_ip11_j = ij(ip1-1,j);
		tmp = d_sqrt (9.8 * dep[ij_ip21_j]);
		if (indr == 1) { 	/*     piston */
			u[ij_ip21_j] = -(pistar * sin(pistbr * time_h));
			v[ij_ip21_j] = v[ij_ip11_j];
			/*     Changes of 5/91 */
			u[ij_ip11_j] = u[ij_ip21_j];
			if (dep[ij_ip21_j] > 0.) {
				h[ij_ip21_j] = -(u[ij_ip21_j] * tmp) / 9.8;
				h[ij_ip21_j] = h[ij_ip21_j];
			}
			continue;
		}
		else if (indr == 3) { 	/*     reflective */
			u[ij_ip21_j] = 0.;
			v[ij_ip21_j] = v[ij_ip11_j];
			h[ij_ip21_j] = h[ij_ip11_j];
			continue;
		}

		/*     continuative */
		if (polar != 0) dx = dxp[ij_ip11_j];
		if (polar == 2) dy = dxp[ij_ip11_j];	/* MERCATOR */
		/*if (dep[ij(ip2-1,j)] < 0.) goto L333;	/* This is bad */
		if (fabs(dep[ij_ip21_j]) > 0.005) dep[ij_ip21_j] = 0.;	/* THIS SEAMS TO PREVENT BORDER INSTABILITIES */
		mul = tmp * dt / dx;
		v[ij_ip21_j] += (v[ij_ip11_j] - v[ij_ip21_j]) * mul;
		if (v[ij_ip11_j] == v[ij(ip-1,j)])
			v[ij_ip21_j] = v[ij_ip11_j];
		u[ij_ip21_j] += (u[ij_ip11_j] - u[ij_ip21_j]) * mul;
		if (u[ij_ip11_j] == u[ij(ip-1,j)])
			u[ij_ip21_j] = u[ij_ip11_j];
		h[ij_ip21_j] += (h[ij_ip11_j] - h[ij_ip21_j]) * mul;
		if (h[ij_ip11_j] == h[ij(ip-1,j)])
			h[ij_ip21_j] = h[ij_ip11_j];
		if (iopt == 0) continue;
/*L333:*/
		v[ij_ip21_j] = v[ij_ip11_j];
		u[ij_ip21_j] = u[ij_ip11_j];
		h[ij_ip21_j] = h[ij_ip11_j];
	}
	/*     top boundary */
	/*      do top only to ip instead of ip2    ** */
	/*       Changed from 2 to 1 on 6/90 */
	for (i = 0; i < ip1; i++) {
		ij_i_jp21 = ij(i,jp2-1);	ij_i_jp11 = ij(i,jp1-1);
		tmp = d_sqrt (9.8 * dep[ij_i_jp21]);
		if (indt == 1) { 	/*     piston */
			v[ij_i_jp21] = -(pistat * sin(pistbt * time_h));
			u[ij_i_jp21] = u[ij_i_jp11];
			/*      Changed 5/91 */
			v[ij_i_jp11] = v[ij_i_jp21];
			/*       THE HEIGHT sign  used to be +     4/23/91 */
			if (dep[ij_i_jp11] > 0.)
				h[ij_i_jp21] = -(v[ij_i_jp21] * tmp) / 9.8;
			continue;
		}
		else if (indt == 3) { 	/*     reflective */
			v[ij_i_jp21] = 0.;
			u[ij_i_jp21] = u[ij_i_jp11];
			h[ij_i_jp21] = h[ij_i_jp11];
			continue;
		}

		/*     continuative */
		if (polar != 0) dx = dxp[ij_i_jp21];	/* THIS WAS NOT ON THE ORIGINAL CODE. FORGOTTEN? */
		if (polar == 2) dy = dxp[ij_i_jp21];	/* MERCATOR */
		/*if (dep[ij(i,jp2-1)] < 0.) goto L343;	/* This is bad */
		if (fabs(dep[ij_i_jp21]) > 0.005) dep[ij_i_jp21] = 0.;	/* THIS SEAMS TO PREVENT BORDER INSTABILITIES */
		mul = tmp * dt / dy;
		v[ij_i_jp21] += (v[ij_i_jp11] - v[ij_i_jp21]) * mul;
		if (v[ij_i_jp11] == v[ij(i,jp-1)]) 
			v[ij_i_jp21] = v[ij_i_jp11];
		u[ij_i_jp11] += (u[ij_i_jp11] - u[ij_i_jp21]) * mul;
		if (u[ij_i_jp11] == u[ij_i_jp11]) 
			u[ij_i_jp21] = u[ij_i_jp11];
		h[ij_i_jp11] += (h[ij_i_jp11] - h[ij_i_jp21]) * mul;
		if (h[ij_i_jp11] == h[ij_i_jp11]) 
			h[ij_i_jp21] = h[ij_i_jp11];
		if (iopt == 0) continue;
/*L343:*/
		v[ij_i_jp21] = v[ij_i_jp11];		/* jp doesn't make sense. Perhaps jp1? */
		u[ij_i_jp21] = u[ij_i_jp11];
		h[ij_i_jp21] = h[ij_i_jp11];
	}
	return 0;
}

/* --------------------------------------------------------------------------- */

int uvh_(double *dep, double *r, double *rn, double *u, double *un, double *v, double *vn,
	 double *h, double *hn, double *dxp, double *cca, double m_per_deg) {

	double cang, ck, sa, sb, sang, cang1, tmp, thu, thv;
	double td, tu, tv, angltt, td1, tu1, tv1, tu2, tv2, dph;
	mwSize	i, j, ij_ij, i1_j, i_j1, im1_j, i_jm1;

	sb = 0.;	sa = 0.;
	for (j = 1; j < jp1; j++) {
		for (i = 1; i < ip1; i++) {
			ij_ij = ij(i,j);	i1_j = ij(i+1,j);	i_j1 = ij(i,j+1);
			im1_j = ij(i-1,j);	i_jm1 = ij(i,j-1);
			/* *****5/91*POLAR Each cell has its own dx and cos angle (cang) */
			if (polar != 0) dx = dxp[ij_ij];
			if (polar == 2) dy = dxp[ij_ij];	/* MERCATOR */
			cang1 = 1.;		cang = 1.;
			if (polar != 0) cang = dxp[ij_ij] / (dangx * m_per_deg);
			if (polar != 0) cang1 = dxp[i_j1] / (dangx * m_per_deg);
			/*     donor cell difference */
			/*     will get diffusion if time step is too small */
			td1 = dep[i1_j] + h[i1_j] - r[i1_j];
			td = dep[ij_ij] + h[ij_ij] - r[ij_ij];
			tv1 = dep[i_j1] + h[i_j1] - r[i_j1];
			tv = td;
			if (u[i1_j] > 0.)
				td1 = dep[ij_ij] + h[ij_ij] - r[ij_ij];
			if (u[ij_ij] > 0.)
				td  = dep[im1_j] + h[im1_j] - r[im1_j];
			if (v[i_j1] > 0.)
				tv1 = dep[ij_ij] + h[ij_ij] - r[ij_ij];
			if (v[ij_ij] > 0.)
				tv  = dep[i_jm1] + h[i_jm1] - r[i_jm1];
			/*      Special for Flooding */
			if (td1 < 0.) td1 = 0.;	if (td < 0.) td = 0.;
			if (tv1 < 0.) tv1 = 0.;	if (tv < 0.) tv = 0.;
			hn[ij_ij] = h[ij_ij] - dt * ((u[i1_j] * td1 - u[ij_ij] * td) / 
				dx + (v[i_j1] * cang1 * tv1 - v[ij_ij] * cang * tv) /
				(cang * dy)) + (rn[ij_ij] - r[ij_ij]);
			/*    ROUGH is factor for surface roughness = actual height/ideal height */
			if (rough != 1. && dep[ij_ij] < 0.) 
				hn[ij_ij] *= rough;
		}
	}

	for (j = 1; j < jp1; j++)
		for (i = 1; i < ip1; i++) {
			ij_ij = ij(i,j); 
			h[ij_ij] = hn[ij_ij];
		}

	bndy_(dep, r, rn, u, un, v, vn, h, hn, dxp);

	for (j = 1; j < jp1; j++) {
		angltt = anglt[j] * D2R;
		sang = sin(angltt);
		if (polar != 0) cf = sang * 1.454e-4;
		for (i = 1; i < ip1; i++) {
			ij_ij = ij(i,j);	i1_j = ij(i+1,j);	i_j1 = ij(i,j+1);
			im1_j = ij(i-1,j);	i_jm1 = ij(i,j-1);
			/* **********POLAR Each cell has its own dx  5/91 */
			if (polar != 0) dx = dxp[ij_ij];
			if (polar == 2) dy = dxp[ij_ij];	/* MERCATOR */
			if (cc == 0.) goto L10;
			/*       Special for Flooding */
			dph = dep[ij_ij] + h[ij_ij];
			if (dph <= .1) goto L10;
			ck = cca[ij_ij];
			tmp = d_sqrt(u[ij_ij]*u[ij_ij] + v[ij_ij]*v[ij_ij]) / 
			      (ck * ck * (dep[ij_ij] + h[ij_ij]));
			sb = 9.8 * u[ij_ij] * tmp;
			/*      CORRECTED ll/20/87 */
			sa = 9.8 * v[ij_ij] * tmp;
L10:
			tu1 = u[i1_j] - u[ij_ij];
			tu2 = u[i_j1] - u[ij_ij];
			tv = (v[ij_ij] + v[i_j1] + v[ij(i-1,j+1)] + v[im1_j]) * .25;
			if (u[ij_ij] > 0.) tu1 = u[ij_ij] - u[im1_j];
			if (tv > 0.) tu2 = u[ij_ij] - u[i_jm1];
			thu = h[ij_ij] - h[im1_j];
			tv1 = v[i1_j] - v[ij_ij];
			tv2 = v[i_j1] - v[ij_ij];
			tu = (u[ij_ij] + u[i1_j] + u[i_jm1] + u[ij(i+1,j-1)]) * .25;
			if (tu > 0.) tv1 = v[ij_ij] - v[im1_j];
			if (v[ij_ij] > 0.) tv2 = v[ij_ij] - v[i_jm1];
			/*      CORRECTED ll/20/87 */
			thv = h[ij_ij] - h[ij(i,j-1)];
			un[ij_ij] = u[ij_ij] - dt * (u[ij_ij] * tu1 / dx + tv * tu2 / dy) - 9.8 * 
				dt * thu / dx - dt * (sb - cf * v[ij_ij] - sfx);
			vn[ij_ij] = v[ij_ij] - dt * (tu * tv1 / dx + v[ij_ij] * tv2 / dy) - 9.8 * 
				dt * thv / dy - dt * (sa + cf * u[ij_ij] - sfy);
			/*     FLOODING--Assumes above sea level is a negative depth-- */
		}
	}

	for (j = 1; j < jp1; j++) {
		for (i = 1; i < ip1; i++) {
			ij_ij = ij(i,j);
			r[ij_ij] = rn[ij_ij];
			u[ij_ij] = un[ij_ij];
			v[ij_ij] = vn[ij_ij];
			/*     artificial limits */
			/* Another nonsense. Original code compared if ... < 1e-8, while using float*/
			if (fabs(u[ij_ij]) < 1e-6) u[ij_ij] = 0.;
			if (fabs(v[ij_ij]) < 1e-6) v[ij_ij] = 0.;
			if (fabs(h[ij_ij]) < 5e-3) h[ij_ij] = 0.;
		}
	}
	return 0;
}


/* --------------------------------------------------------------------------- */
int write_grd_bin(char *name, double x_min, double y_min, double x_inc, double y_inc,
		   mwSize i_start, mwSize j_start, mwSize i_end, mwSize j_end, mwSize nX, float *work) {

	/* Writes a grid in the Surfer binary format */
	mwSize i, j;
	double x_max, y_max;
	float work_min = FLT_MAX, work_max = -FLT_MAX, tmp;
	struct srf_header h;
	FILE *fp;

	if ((fp = fopen (name, "wb")) == NULL) {
		mexPrintf("Fatal Error: Could not create file %s!\n", name);
		return (-1);
	}

	x_max = x_min + (i_end - i_start - 1) * x_inc;
	y_max = y_min + (j_end - j_start - 1) * y_inc;

	/* Find zmin/zmax */
	for (j = j_start; j < j_end; j++) {
		for (i = i_start; i < i_end; i++) {
			tmp = work[ijs(i,j,nX)];
			work_max = MAX(tmp, work_max);
			work_min = MIN(tmp, work_min);
		}
	}

	/* store header information and array */
	strcpy (h.id,"DSBB");
	h.nx = (i_end - i_start);	h.ny = (j_end - j_start);
	h.x_min = x_min;		h.x_max = x_max;
	h.y_min = y_min;		h.y_max = y_max;
	h.z_min = (double)work_min;	h.z_max = (double)work_max;

	if (fwrite ((void *)&h, sizeof (struct srf_header), (size_t)1, fp) != 1) {
		fprintf (stderr, "Fatal Error: Error writing file %s!\n", name);
		return (-1);
	}

	for (j = j_start; j < j_end; j++)
		for (i = i_start; i < i_end; i++)
			fwrite ((void *)&work[ijs(i,j,nX)], sizeof(float), (size_t)1, fp);

	fclose(fp);
	return (0);
}

/* ------------------------------------------------------------------------------ */
int read_grd_info_ascii (char *file, struct srf_header *hdr) {

	char line[512], id[5];
	FILE *fp;

	if ((fp = fopen (file, "r")) == NULL) {
		mexPrintf ("%s: Unable to read file %s - exiting\n", "swan", file);
		return (-1);
	}

	fgets (line, 512, fp);
	sscanf (line, "%s", hdr->id);
	fgets (line, 512, fp);
	sscanf (line, "%d %d", &hdr->nx, &hdr->ny);
	fgets (line, 512, fp);
	sscanf (line, "%lf %lf", &hdr->x_min, &hdr->x_max);
	fgets (line, 512, fp);
	sscanf (line, "%lf %lf", &hdr->y_min, &hdr->y_max);
	fgets (line, 512, fp);
	sscanf (line, "%lf %lf", &hdr->z_min, &hdr->z_max);
	fclose(fp);
	sprintf (id, "%.4s", hdr->id);
	if (strcmp (id, "DSAA") == 0)
		return (0);
	else if (strcmp (id, "DSBB") == 0) {
		fp = fopen (file, "rb");
		read_header_bin (fp, hdr);
		fclose(fp);
		return (1);
	}
	else
		return (-1);
}

/* ------------------------------------------------------------------------------ */
int read_grd_ascii (char *file, struct srf_header *hdr, double *work) {

	/* Reads a grid in the Surfer ascii format */
	mwSize i = 0, j, n_field;
	char *p, buffer[512], line[512];
	FILE *fp;

	if ((fp = fopen (file, "r")) == NULL) {
		mexPrintf ("%s: Unable to read file %s - exiting\n", "swan", file);
		return (-1);
	}

	fgets (line, 512, fp);
	sscanf (line, "%s", hdr->id);
	fgets (line, 512, fp);
	sscanf (line, "%d %d", &hdr->nx, &hdr->ny);
	fgets (line, 512, fp);
	sscanf (line, "%lf %lf", &hdr->x_min, &hdr->x_max);
	fgets (line, 512, fp);
	sscanf (line, "%lf %lf", &hdr->y_min, &hdr->y_max);
	fgets (line, 512, fp);
	sscanf (line, "%lf %lf", &hdr->z_min, &hdr->z_max);

	while (fgets (line, 512, fp) != NULL) {
		strcpy (buffer, line);
		n_field = count_col (buffer);	/* Count # of fields in line */
		if (n_field == 0) continue;
		p = (char *)strtok (line, " \t\n\015\032");
		j = 0;
		while (p && j < n_field) {
			sscanf (p, "%lf", &work[i]);
			j++;	i++;
			p = (char *)strtok ((char *)NULL, " \t\n\015\032");
		}
	}
	fclose(fp);
	return (0);
} 

/* -------------------------------------------------------------------- */
int count_col (char *line) {
	/* Count # of fields contained in line */
	mwSize n_col = 0;
	char *p;

	p = (char *)strtok (line, " \t\n\015\032");
	while (p) {	/* Count # of fields */
		n_col++;
		p = (char *)strtok ((char *)NULL, " \t\n\015\032");
	}
	return (n_col);
}

/* -------------------------------------------------------------------- */
int read_header_bin (FILE *fp, struct srf_header *hdr) {
	/* Reads the header of a binary Surfer gridfile */
	fread ((void *)hdr, sizeof (struct srf_header), (size_t)1, fp); 
	return (0);
}

/* -------------------------------------------------------------------- */
int read_grd_bin (char *file, struct srf_header *hdr, double *work) {
	mwSize	i, j, ij, kk;
	float	*tmp;			/* Array pointer for reading in rows of data */
	FILE	*fp;

	if ((fp = fopen (file, "rb")) == NULL) {
		mexPrintf ("%s: Unable to read file %s - exiting\n", "swan", file);
		return (-1);
	}
	fread ((void *)hdr, sizeof (struct srf_header), (size_t)1, fp); 

	/* Allocate memory for one row of data (for reading purposes) */
	tmp = (float *) mxMalloc ((size_t)hdr->nx * sizeof (float));
	for (j = 0; j < (hdr->ny - 1); j++) {
		fread (tmp, sizeof (float), (size_t)hdr->nx, fp);	/* Get one row */
		ij = j * hdr->nx;
		for (i = 0; i < hdr->nx; i++) {
			kk = ij + i;
			work[kk] = tmp[i];
		}
	}
	fclose(fp);
	mxFree ((void *)tmp);
	return (0);
}

/* -------------------------------------------------------------------- */
int read_params(char *file) {
	/* Read parameters that controls SWAN running */
	char line[128];
	FILE *fp;

	if ((fp = fopen (file, "r")) == NULL) {
		mexPrintf ("%s: Unable to open file %s - exiting\n", "swan", file);
		return (-1);
	}

	fgets (line, 128, fp);		/* Comment */
	fgets (line, 128, fp);
	sscanf (line, "%lf", &dt);
	fgets (line, 128, fp);
	sscanf (line, "%d", &grn);
	fgets (line, 128, fp);
	sscanf (line, "%lf", &cf);
	fgets (line, 128, fp);
	sscanf (line, "%lf", &cc);
	fgets (line, 128, fp);
	sscanf (line, "%lf", &sfx);
	fgets (line, 128, fp);
	sscanf (line, "%lf", &sfy);
	fgets (line, 128, fp);
	sscanf (line, "%d", &polar);
	fgets (line, 128, fp);
	sscanf (line, "%lf", &rough);
	fgets (line, 128, fp);
	sscanf (line, "%d", &cumint);
	fgets (line, 128, fp);
	sscanf (line, "%lf", &pistal);
	fgets (line, 128, fp);
	sscanf (line, "%lf", &pistbl);
	fgets (line, 128, fp);
	sscanf (line, "%lf", &pistab);
	fgets (line, 128, fp);
	sscanf (line, "%lf", &pistbb);
	fgets (line, 128, fp);
	sscanf (line, "%lf", &pistar);
	fgets (line, 128, fp);
	sscanf (line, "%lf", &pistbr);
	fgets (line, 128, fp);
	sscanf (line, "%lf", &pistat);
	fgets (line, 128, fp);
	sscanf (line, "%lf", &pistbt);
	fgets (line, 128, fp);		/* Comment */
	fgets (line, 128, fp);
	sscanf (line, "%d", &indl);
	fgets (line, 128, fp);
	sscanf (line, "%d", &indb);
	fgets (line, 128, fp);
	sscanf (line, "%d", &indr);
	fgets (line, 128, fp);
	sscanf (line, "%d", &indt);
	fgets (line, 128, fp);
	sscanf (line, "%d", &iopt);
	fclose (fp);
	return (0);
}

/* -------------------------------------------------------------------- */
int count_n_maregs(char *file) {
	mwSize	i = 0;
	char	line[512];
	FILE	*fp;

	if ((fp = fopen (file, "r")) == NULL) {
		mexPrintf ("%s: Unable to open file %s - exiting\n", "swan", file);
		return (-1);
	}
	while (fgets (line, 512, fp) != NULL) {
		if (line[0] == '#') continue;	/* Jump comment lines */
		i++;
	}
	return(i);
	fclose (fp);
	return (0);
}

/* -------------------------------------------------------------------- */
int read_maregs(char *file) {
	/* Read maregraph positions and convert them to vector indices */
	mwSize	i = 0, ix, jy;
	double	x, y, dx, dy;
	char	line[512];
	FILE	*fp;

	dx = (hdr_b.x_max - hdr_b.x_min) / (hdr_b.nx - 1);
	dy = (hdr_b.y_max - hdr_b.y_min) / (hdr_b.ny - 1);
	if ((fp = fopen (file, "r")) == NULL) {
		mexPrintf ("%s: Unable to open file %s - exiting\n", "swan", file);
		return (-1);
	}

	while (fgets (line, 512, fp) != NULL) {
		if (line[0] == '#') continue;	/* Jump comment lines */
		sscanf (line, "%lf %lf", &x, &y);
		ix = irint((x - hdr_b.x_min) / dx);
		jy = irint((y - hdr_b.y_min) / dy); 
		lcum_p[i] = jy * ip2 + ix; 
		i++;
	}
	fclose (fp);
	return (0);
}

/* -------------------------------------------------------------------- */
void	no_sys_mem (char *where, mwSize n) {	
		mexPrintf ("Fatal Error: %s could not allocate memory, n = %d\n", where, n);
}

/*     ********* CHECK OF MAXIMUM VALUE *********** */
void max_z (double *zm, double *h_bak) {
	mwSize i, ij;

	ij = ip2 * jp2;
	for (i = 0; i < ij; i++)
		if(zm[i] < h_bak[i]) zm[i] = h_bak[i];
}

/* -------------------------------------------------------------------- */
void change (double *h, double *h_bak) {
	mwSize i, ij;

	ij = ip2 * jp2;
	for (i = 0; i < ij; i++)
		h_bak[i] = h[i];
}

/* -------------------------------------------------------------------- */
int decode_R (char *item, double *w, double *e, double *s, double *n) {
	char *text, string[BUFSIZ];
	
	/* Minimalist code to decode option -R extracted from GMT_get_common_args */
	
	mwSize i, error = 0;
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
	mwSize i, colons = 0, suffix;
	double degree, minute, degfrac, second;

	for (i = 0; text[i]; i++) if (text[i] == ':') colons++;
	suffix = (mwSize)text[i-1];	/* Last character in string */
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
int open_most_nc (char *base, char *name_var, mwSize *ids, mwSize nx, mwSize ny, 
		double dtx, double dty, double xMinOut, double yMinOut) {
	/* Open and initialize a MOST netCDF file for writing ---------------- */
	mwSize ncid, m, n, status, dim0[3], dim3[3];
	float dummy = -1e34f;
	double *x, *y;
	char	*long_name = NULL, *units = NULL, *basename = NULL;

	basename = (char *) mxMalloc (strlen(base) * sizeof (char));
	strcpy(basename, base);
	if (!strcmp(name_var,"HA")) {
		strcat(basename,"_ha.nc");
		long_name = "Wave Amplitude";
		units = "CENTIMETERS";
	}
	else if (!strcmp(name_var,"VA")) {
		strcat(basename,"_va.nc");
		long_name = "Velocity Component along Latitude";
		units = "CENTIMETERS/SECOND";
	}
	else if (!strcmp(name_var,"UA")) {
		strcat(basename,"_ua.nc");
		long_name = "Velocity Component along Longitude";
		units = "CENTIMETERS/SECOND";
	}

	if ( (status = nc_create (basename, NC_CLOBBER, &ncid)) != NC_NOERR) {
		mexPrintf ("swan: Unable to create file %s - exiting\n", basename);
		return(-1);
	}
	/* ---- Define dimensions ------------ */
	err_trap (nc_def_dim (ncid, "LON", (size_t) nx, &dim0[0]));
	err_trap (nc_def_dim (ncid, "LAT", (size_t) ny, &dim0[1]));
	err_trap (nc_def_dim (ncid, "TIME", NC_UNLIMITED, &dim0[2]));

	/* ---- Define variables ------------- */
	dim3[0] = dim0[2];	dim3[1] = dim0[1];	dim3[2] = dim0[0];
	err_trap (nc_def_var (ncid, "LON",	NC_DOUBLE,1, &dim0[0], &ids[0]));
	err_trap (nc_def_var (ncid, "LAT",	NC_DOUBLE,1, &dim0[1], &ids[1]));
	err_trap (nc_def_var (ncid, "SLON",	NC_FLOAT,0,  &dim0[0], &ids[2]));
	err_trap (nc_def_var (ncid, "SLAT",	NC_FLOAT,0,  &dim0[1], &ids[3]));
	err_trap (nc_def_var (ncid, "TIME",	NC_DOUBLE,1, &dim0[2], &ids[4]));
	err_trap (nc_def_var (ncid, name_var,	NC_FLOAT,3,  dim3,     &ids[5]));

	/* ---- Variables Attributes --------- */
	err_trap (nc_put_att_text (ncid, ids[0], "units", 12, "degrees_east"));
	err_trap (nc_put_att_text (ncid, ids[0], "point_spacing", 4, "even"));
	err_trap (nc_put_att_text (ncid, ids[1], "units", 13, "degrees_north"));
	err_trap (nc_put_att_text (ncid, ids[1], "point_spacing", 4, "even"));
	err_trap (nc_put_att_text (ncid, ids[2], "units", 12, "degrees_east"));
	err_trap (nc_put_att_text (ncid, ids[2], "long_name", 16, "Source Longitude"));
	err_trap (nc_put_att_text (ncid, ids[3], "units", 13, "degrees_north"));
	err_trap (nc_put_att_text (ncid, ids[3], "long_name", 16, "Source Latitude"));
	err_trap (nc_put_att_text (ncid, ids[4], "units", 7, "SECONDS"));
	err_trap (nc_put_att_text (ncid, ids[5], "long_name", strlen(long_name), long_name));
	err_trap (nc_put_att_text (ncid, ids[5], "units", strlen(units), units));
	err_trap (nc_put_att_float (ncid,ids[5], "missing_value", NC_FLOAT, 1, &dummy));
	err_trap (nc_put_att_float (ncid,ids[5], "_FillValue", NC_FLOAT, 1, &dummy));
	err_trap (nc_put_att_text (ncid, ids[5], "history", 6, "Nikles"));

	/* ---- Global Attributes ------------ */
	err_trap (nc_put_att_text (ncid, NC_GLOBAL, "history", 10, "Mirone Tec"));
	err_trap (nc_put_att_text (ncid, NC_GLOBAL, "title", 22, "Created by Mirone-Swan"));

	err_trap (nc_enddef (ncid));

	/* ---- Write the vector coords ------ */
	x = (double *) mxMalloc (sizeof (double) * nx);
	y = (double *) mxMalloc (sizeof (double) * ny);

	for (n = 0; n < nx; n++) x[n] = xMinOut + n * dtx;
	for (m = 0; m < ny; m++) y[m] = yMinOut + m * dty;
	err_trap (nc_put_var_double (ncid, ids[0], x));
	err_trap (nc_put_var_double (ncid, ids[1], y));
	mxFree ((void *)x); 
	mxFree ((void *)y); 
	mxFree ((void *)basename); 

	return (ncid);
}

/* --------------------------------------------------------------------------- */
void write_most_slice(mwSize *ncid_most, mwSize *ids_most, mwSize i_start, mwSize j_start, mwSize i_end, mwSize j_end,
		mwSize nX, float *work, double *h, double *dep, double *u, double *v, float *tmp,
		size_t *start, size_t *count) {
	/* Write a slice of _ha.nc, _va.nc & _ua.nc MOST netCDF files */
	mwSize i, j, n, ij, k;

	for (n = 0; n < 3; n++) {	/* Loop over, Amplitude, Xmomentum & Ymomentum */
		if (n == 0) {		/* Amplitude */
			for (j = j_start, k = 0; j < j_end; j++)
				for (i = i_start; i < i_end; i++)
					tmp[k++] = work[ijs(i,j,nX)] * 100;
			err_trap (nc_put_vara_float (ncid_most[0], ids_most[0], start, count, tmp));
		}
		else if (n == 1) {		/* X velocity */ 
			for (j = j_start, k = 0; j < j_end; j++)
				for (i = i_start; i < i_end; i++) {
					ij = ijs(i,j,nX);
					if (( tmp[k] = (float)(h[ij] + dep[ij]) ) < 0.) tmp[k] = 0.;
					tmp[k++] *= (float)u[ij] * 100;
				}
			err_trap (nc_put_vara_float (ncid_most[1], ids_most[1], start, count, tmp));
		}
		else {				/* Y velocity */ 
			for (j = j_start, k = 0; j < j_end; j++)
				for (i = i_start; i < i_end; i++) {
					ij = ijs(i,j,nX);
					if (( tmp[k] = (float)(h[ij] + dep[ij]) ) < 0.) tmp[k] = 0.;
					tmp[k++] *= (float)v[ij] * 100;
				}
			err_trap (nc_put_vara_float (ncid_most[2], ids_most[2], start, count, tmp));
		}
	}
}

/* --------------------------------------------------------------------------- */
void write_anuga_slice(mwSize ncid, mwSize z_id, mwSize i_start, mwSize j_start, mwSize i_end, mwSize j_end, mwSize nX, 
		float *work, double *h, double *dep, double *u, double *v, float *tmp,
		size_t *start, size_t *count, float *slice_range, mwSize idx, mwSize with_land) {
	/* Write a slice of either STAGE, XMOMENTUM or YMOMENTUM of a Anuga's .sww netCDF file */
	mwSize i, j, ij, k, ncl;

	ncl = (i_end - i_start)*(j_end - j_start);
	k = 0;
	for (j = j_start; j < j_end; j++) {
		if (idx == 1) 			/* Anuga calls this -> stage */
			if (!with_land)		/* Land nodes are kept = 0 */
				for (i = i_start; i < i_end; i++)
					tmp[k++] = work[ijs(i,j,nX)];
			else {
				for (i = i_start; i < i_end; i++) {
					ij = ijs(i,j,nX);
					if (work[ij] == 0 && dep[ij] < 0) tmp[k++] = (float)-dep[ij];
					else 	tmp[k++] = work[ij];
				}
			}

		else if (idx == 2) {		/* X momentum */
			for (i = i_start; i < i_end; i++) {
				ij = ijs(i,j,nX);
				if (( tmp[k] = (float)(h[ij] + dep[ij]) ) < 0.) tmp[k] = 0.;
				tmp[k++] *= (float)u[ij];
			}
		}
		else {				/* Y momentum */
			for (i = i_start; i < i_end; i++) {
				ij = ijs(i,j,nX);
				if (( tmp[k] = (float)(h[ij] + dep[ij]) ) < 0.) tmp[k] = 0.;
				tmp[k++] *= (float)v[ij];
			}
		}
	}

	/* ----------- Find the min/max of this slice --------- */
	for (k = 0; k < ncl; k++) {
		slice_range[1] = MAX(tmp[k], slice_range[1]);
		slice_range[0] = MIN(tmp[k], slice_range[0]);
	}

	err_trap (nc_put_vara_float (ncid, z_id, start, count, tmp));
}

/* -------------------------------------------------------------------- */
int open_anuga_sww (char *fname_sww, mwSize *ids, mwSize i_start, mwSize j_start, mwSize i_end, mwSize j_end, 
		mwSize nX, double dtx, double dty, double *dep, double xMinOut, double yMinOut, 
		float z_min, float z_max) {

	/* Open and initialize a ANUGA netCDF file for writing ---------------- */
	mwSize ncid, m, n, nx, ny, status, nVolumes, nPoints, dim0[5], dim2[2], dim3[2];
	mwSize i, j, k, m_nx, m1_nx, *volumes, *vertices, v1, v2, v3, v4;
	float dummy2[2], *x, *y, yr, *tmp;
	double dummy, nan, faultPolyX[11], faultPolyY[11], faultSlip[10], faultStrike[10], 
		faultDip[10], faultRake[10], faultWidth[10], faultDepth[10];

	if ( (status = nc_create (fname_sww, NC_CLOBBER, &ncid)) != NC_NOERR) {
		mexPrintf ("swan: Unable to create file %s - exiting\n", fname_sww);
		return(-1);
	}
	/* ---- Define dimensions ------------ */
	nx = i_end - i_start;		ny = j_end - j_start;
	nVolumes = (nx - 1) * (ny - 1) * 2;
	nPoints = nx * ny;
	err_trap (nc_def_dim (ncid, "number_of_volumes", (size_t) nVolumes, &dim0[0]));
	err_trap (nc_def_dim (ncid, "number_of_vertices", (size_t) 3, &dim0[1]));
	err_trap (nc_def_dim (ncid, "numbers_in_range", (size_t) 2, &dim0[2]));
	err_trap (nc_def_dim (ncid, "number_of_points", (size_t) nPoints, &dim0[3]));
	err_trap (nc_def_dim (ncid, "number_of_timesteps", NC_UNLIMITED, &dim0[4]));

	/* ---- Define variables ------------- */
	dim2[0] = dim0[4];		dim2[1] = dim0[3];
	dim3[0] = dim0[0];		dim3[1] = dim0[1];
	err_trap (nc_def_var (ncid, "x",		NC_FLOAT,1, &dim0[3], &ids[0]));
	err_trap (nc_def_var (ncid, "y",		NC_FLOAT,1, &dim0[3], &ids[1]));
	err_trap (nc_def_var (ncid, "z",		NC_FLOAT,1, &dim0[3], &ids[2]));
	err_trap (nc_def_var (ncid, "elevation",	NC_FLOAT,1, &dim0[3], &ids[3]));
	err_trap (nc_def_var (ncid, "elevation_range", 	NC_FLOAT,1, &dim0[2], &ids[4]));
	err_trap (nc_def_var (ncid, "volumes",		NC_INT,  2, dim3, &ids[5]));
	err_trap (nc_def_var (ncid, "time",		NC_DOUBLE,1,&dim0[4], &ids[6]));
	err_trap (nc_def_var (ncid, "stage",		NC_FLOAT,2, dim2, &ids[7]));
	err_trap (nc_def_var (ncid, "stage_range",	NC_FLOAT,1, &dim0[2], &ids[8]));
	err_trap (nc_def_var (ncid, "xmomentum",	NC_FLOAT,2, dim2, &ids[9]));
	err_trap (nc_def_var (ncid, "xmomentum_range", 	NC_FLOAT,1, &dim0[2], &ids[10]));
	err_trap (nc_def_var (ncid, "ymomentum",	NC_FLOAT,2, dim2, &ids[11]));
	err_trap (nc_def_var (ncid, "ymomentum_range", 	NC_FLOAT,1, &dim0[2], &ids[12]));

	/* ---- Global Attributes ------------ */
	err_trap (nc_put_att_text (ncid, NC_GLOBAL, "institution", 10, "Mirone Tec"));
	err_trap (nc_put_att_text (ncid, NC_GLOBAL, "description", 22, "Created by Mirone-Swan"));
	err_trap (nc_put_att_double (ncid, NC_GLOBAL, "xllcorner", NC_DOUBLE, 1, &xMinOut));
	err_trap (nc_put_att_double (ncid, NC_GLOBAL, "yllcorner", NC_DOUBLE, 1, &yMinOut));
	dummy = 29;	err_trap (nc_put_att_double (ncid, NC_GLOBAL, "zone", NC_DOUBLE, 1, &dummy));
	dummy = 0;	err_trap (nc_put_att_double (ncid, NC_GLOBAL, "starttime", NC_DOUBLE, 1, &dummy));
	dummy = 500000;	err_trap (nc_put_att_double (ncid, NC_GLOBAL, "false_easting", NC_DOUBLE, 1, &dummy));
	dummy = 0;	err_trap (nc_put_att_double (ncid, NC_GLOBAL, "false_northing", NC_DOUBLE, 1, &dummy));
	err_trap (nc_put_att_text (ncid, NC_GLOBAL, "datum", 5, "wgs84"));
	err_trap (nc_put_att_text (ncid, NC_GLOBAL, "projection", 3, "UTM"));
	err_trap (nc_put_att_text (ncid, NC_GLOBAL, "units", 1, "m"));
	/* Initialize the following attribs with NaNs. A posterior call will eventualy fill them with the right values */
	nan = mxGetNaN();
	for (i = 0; i < 10; i++) {
		faultPolyX[i] = faultPolyY[i] = faultSlip[i] = faultDip[i] = faultStrike[i] = 
				faultRake[i] = faultWidth[i] = faultDepth[i] = nan;
	}
	faultPolyX[10] = faultPolyY[10] = nan;		/* Those have an extra element */
	err_trap (nc_put_att_double (ncid, NC_GLOBAL, "faultPolyX", NC_DOUBLE, 11, faultPolyX));
	err_trap (nc_put_att_double (ncid, NC_GLOBAL, "faultPolyY", NC_DOUBLE, 11, faultPolyY));
	err_trap (nc_put_att_double (ncid, NC_GLOBAL, "faultStrike", NC_DOUBLE, 10, faultStrike));
	err_trap (nc_put_att_double (ncid, NC_GLOBAL, "faultSlip", NC_DOUBLE, 10, faultSlip));
	err_trap (nc_put_att_double (ncid, NC_GLOBAL, "faultDip", NC_DOUBLE, 10, faultDip));
	err_trap (nc_put_att_double (ncid, NC_GLOBAL, "faultRake", NC_DOUBLE, 10, faultRake));
	err_trap (nc_put_att_double (ncid, NC_GLOBAL, "faultWidth", NC_DOUBLE, 10, faultWidth));
	err_trap (nc_put_att_double (ncid, NC_GLOBAL, "faultDepth", NC_DOUBLE, 10, faultDepth));

	/* ---- Write the vector coords ------ */
	x = (float *) mxMalloc (sizeof (float) * (nx * ny));
	y = (float *) mxMalloc (sizeof (float) * (nx * ny));
	vertices = (mwSize *) mxMalloc (sizeof (mwSize) * (nx * ny));
	volumes = (mwSize *) mxMalloc (sizeof (mwSize) * (nVolumes * 3));

	/* Construct 2 triangles per 'rectangular' element */
	i = 0;
	for (m = 0; m < ny; m++) {		/* By rows    - Y */
		yr = (float)(m * dty);
		for (n = 0; n < nx; n++) {	/* By columns - X */
			x[i] = (float)(n * dtx);
			y[i] = yr;
			vertices[i] = i;
			i++;
		}
	}

	for (n = i = 0; n < nx - 1; n++) {			/* X */
		for (m = 0; m < ny - 1; m++) {			/* Y */
			m_nx = m * nx;		m1_nx = (m + 1) * nx;
			v1 = vertices[n + m_nx];	v2 = vertices[n + 1 + m_nx];
			v3 = vertices[n + 1 + m1_nx];	v4 = vertices[n + m1_nx];
			volumes[i] = v1;	volumes[i+1] = v2;	volumes[i+2] = v3;
			i += 3;
			volumes[i] = v1;	volumes[i+1] = v3;	volumes[i+2] = v4;
			i += 3;
		}
	}

	mxFree ((void *)vertices);

	err_trap (nc_enddef (ncid));
	err_trap (nc_put_var_float (ncid, ids[0], x));
	err_trap (nc_put_var_float (ncid, ids[1], y));
	mxFree ((void *)x); 
	mxFree ((void *)y); 

	tmp = (float *) mxMalloc (sizeof (float) * (nx * ny));
	for (j = j_start, k = 0; j < j_end; j++) {
		for (i = i_start; i < i_end; i++)
			tmp[k++] = (float)-dep[ijs(i,j,nX)];
	}

	err_trap (nc_put_var_float (ncid, ids[2], tmp));	/* z */

	err_trap (nc_put_var_float (ncid, ids[3], tmp));	/* elevation */
	dummy2[0] = z_min;		dummy2[1] = z_max;
	err_trap (nc_put_var_float (ncid, ids[4], dummy2));	/* elevation_range */
	err_trap (nc_put_var_int (ncid, ids[5], volumes));

	mxFree ((void *)tmp);
	mxFree ((void *)volumes);

	return (ncid);
}
