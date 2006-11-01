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
 *	Translated to C & mexified By
 *	Joaquim Luis - 2005
 *
 *	version WITH waitbar
 */

#include "mex.h"
#include <float.h>
#include <math.h>

#define	FALSE	0
#define	TRUE	1
#define M_PI	3.14159265358979323846
#define D2R		M_PI / 180.
#define LF		0x0A
#define CR		0x0D

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
#define ijs(i,j) ((i) + (j)*i_end)
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

void no_sys_mem (char *where, int n);
int count_col (char *line);
int read_grd_info_ascii (char *file, struct srf_header *hdr);
int read_header_bin (FILE *fp, struct srf_header *hdr);
int write_grd_ascii(int kk, double x_min, double y_min, double dtx, double dty, int i_end, int j_end, float *work);
int write_grd_bin (int kk, double x_min, double y_min, double dtx, double dty, int i_end, int j_end, float *work);
int read_grd_ascii (char *file, struct srf_header *hdr, double *work);
int read_grd_bin (char *file, struct srf_header *hdr, double *work);
int read_params(char *file);
int read_maregs(char *file);
int count_n_maregs(char *file);
int uvh_(double *dep, double *r, double *rn, double *u, double *un, double *v, double *vn, double *h, double *hn, double *dxp, double *cca);
int bndy_(double *dep, double *r, double *rn, double *u, double *un, double *v, double *vn, double *h, double *hn, double *dxp);
void max_z (double *zm, double *h_bak);
void change (double *h, double *h_bak);

typedef int BOOLEAN;              /* BOOLEAN used for logical variables */
int	ip, jp, polar, dumb, indl, indb, indr, indt, iopt;
int	ip1, jp1, ip2, jp2, grn, cumint, *lcum_p;
float	grav = (float)9.8, m_per_deg = (float)111317.1;
double	dx, dy, dt, cf, cc, sfx, sfy, time_h, rough;
double	pistal, pistbl, pistab, pistbb, pistar, pistbr, pistat, pistbt;
double	dangx, dangy, *anglt;
char	stem[80], *fonte;			/* Name pointer for tsunami source file */
BOOLEAN first_in_uvh = TRUE;
struct	srf_header hdr_b;
struct	srf_header hdr_f;

/* --------------------------------------------------------------------------- */
/* Matlab Gateway routine */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	float	*work, dz, *ptr_mov_32, *mov_32, *time_p;
	float	work_min = FLT_MAX, work_max = -FLT_MAX;
	double	*inicial, *r, *rn, *u, *un, *v, *vn, *h, *h_bak, *hn, *zm;
	double	*dxp, *xpp, *cca, x, y, small = 1e-6;
	double	x_min, y_min, dminx, dminy, dtx, dty, *head, *tmp;
	double	*dep, *dep1, *cum_p, cang, angltt;
	double	*ptr, *ptr1, *h_bar, tmp_ptr[1];		/* Pointers to be used in the waitbar */
	int	i, i2, j, k, ncl, lcum, n_mareg, n_ptmar, cycle;
	int	w_ascii = FALSE, w_bin = TRUE, cumpt = FALSE, error = FALSE;
	int	r_bin_b, r_bin_f, surf_level = TRUE, max_level = FALSE, water_depth = FALSE;
	int	argc, n_arg_no_char = 0, nx, ny, dims[3], n_frames = 0;
	int	n_of_cycles = 1010;	/* Numero de ciclos a calcular */
	char	*bathy = NULL;		/* Name pointer for bathymetry file */
	char 	*params = NULL;		/* Name pointer for parameters file */
	char 	*hcum = NULL;		/* Name pointer for cumulative hight file */
	char 	*maregs = NULL;		/* Name pointer for maregraph positions file */
	char	**argv, w_bar_title[] = "Aguenta ai";
	unsigned char	*ptr_mov_8, *mov_8, *mov_8_tmp;
	BOOLEAN	params_in_input = FALSE, bat_in_input = FALSE, source_in_input = FALSE;
	BOOLEAN	write_grids = FALSE, movie = FALSE, movie_char = FALSE, movie_float = FALSE;
	BOOLEAN	maregs_in_input = FALSE;
	mxArray *rhs[2], *rhs1[2], *lhs[1];
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
		dt = tmp[0];		grn = (int)tmp[1];	cf = tmp[2];
		cc = tmp[3];		sfx = tmp[4];		sfy = tmp[5];
		polar = (int)tmp[6];	rough = tmp[7];		cumint = (int)tmp[8];
		pistal = tmp[9];	pistbl = tmp[10];	pistab = tmp[11];
		pistbb = tmp[12];	pistar = tmp[13];	pistbr = tmp[14];
		pistat = tmp[15];	pistbt = tmp[16];	indl = (int)tmp[17];
		indb = (int)tmp[18];	indr = (int)tmp[19];	indt = (int)tmp[20];
		iopt = (int)tmp[21];
		params_in_input = TRUE;
	}
	if(n_arg_no_char == 6) {		/* A maregraph vector was given as the sixth argument*/
		tmp = mxGetPr(prhs[5]);
		n_mareg = mxGetM(prhs[5]);
		dx = head[7];		dy = head[8];
		lcum_p = (int *) calloc ((size_t)(n_mareg), sizeof(int));
		for (i = 0; i < n_mareg; i++) {
			x = tmp[i];		y = tmp[i+n_mareg];	/* Matlab vectors are stored by columns */
			lcum_p[i] = (irint((y - hdr_b.y_min) / dy) ) * hdr_b.nx + irint((x - hdr_b.x_min) / dx);
		}
		maregs_in_input = TRUE;
		cumpt = TRUE;
	}

	/* get the length of the input string */
	argv=(char **)mxCalloc(argc, sizeof(char *));
	for (i = 0; i < argc; i++) {
		argv[i] = (char *)mxArrayToString(prhs[i+n_arg_no_char]);
	}

	for (i = 0; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {
				case 'a':	/* Write ascii grids */
					w_ascii = TRUE;
					w_bin = FALSE;
					break;
				case 'f':	/* Movie */
					movie = TRUE;
					movie_char = FALSE;
					movie_float = TRUE;
					break;
				case 'm':	/* Movie */
					movie = TRUE;
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
				case 'M':
					max_level = TRUE;
					surf_level = FALSE;
					break;
				case 'N':	/* Numero de ciclos a calcular */
					n_of_cycles = atoi(&argv[i][2]);
					break;
				case 'O':	/* File name for maregraph data */
					hcum  = &argv[i][2];
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
		mexPrintf( "usage: swan [-B<bathy>] [-F<fonte>] [-M] [-N<n_cycles>] [-T<mareg>] [-D] [-a]\n");
		mexPrintf ("\t-a write ascii Surfer grids [Default is binary]\n");
		mexPrintf ("\t-m outputs a 3D grid used to do a movie\n");
		mexPrintf ("\t-B name of bathymetry file (default lap.grd)\n");
		mexPrintf ("\t-D INACABADO write grids with the total water depth\n");
		mexPrintf ("\t-F name of source file (default fonte.grd)\n");
		mexPrintf ("\t-G<stem> write grids at the grn intervals. Append file prefix. Files will be called <stem>#.grd\n");
		mexPrintf ("\t-M write grids of max water level [Default wave surface level]\n");
		mexPrintf ("\t-N number of cycles [Default 1010].\n");
		mexPrintf ("\t-T name of maregraph file (default mareg.xy)\n");
	}

	if (error) return;

	if (!params_in_input) {			/* If params where not given, read them from file */
		read_params(params);
	}
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
	if (!write_grids && !movie && !cumpt) {
		mexPrintf("Nothing selected for output (grids, movie or mregraphs), exiting\n");
		return;
	}

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
	if (cumpt && !maregs_in_input) {
		n_mareg = count_n_maregs(maregs);	/* Count maragraphs number */
	}
	if (cumpt) {
		n_ptmar = n_of_cycles / cumint + 1;
		if ((fp = fopen (hcum, "w")) == NULL) {
			mexPrintf ("%s: Unable to create file %s - exiting\n", "swan", hcum);
			return;
		}
	}

	/* Allocate memory	*/
	if (!bat_in_input) {	/* Bathymetry & Source will be read from files. So need to alloc */
		if ((inicial = (double *) calloc ((size_t)(ncl), sizeof(double)) ) == NULL) 
			{no_sys_mem("swan --> (inicial)", ncl);	return;}
	}
	if ((dep = (double *) calloc ((size_t)(ncl), sizeof(double)) ) == NULL) 
		{no_sys_mem("swan --> (dep)", ncl);	return;}
	if ((work = (float *) calloc ((size_t)(ncl), sizeof(float)) ) == NULL) 
		{no_sys_mem("swan --> (work)", ncl);	return;}
	if ((r = (double *) calloc ((size_t)(ncl), sizeof(double)) ) == NULL) 
		{no_sys_mem("swan --> (r)", ncl);	return;}
	if ((rn = (double *) calloc ((size_t)(ncl), sizeof(double)) ) == NULL) 
		{no_sys_mem("swan --> (rn)", ncl);	return;}
	if ((u = (double *) calloc ((size_t)(ncl), sizeof(double)) ) == NULL) 
		{no_sys_mem("swan --> (u)", ncl);	return;}
	if ((un = (double *) calloc ((size_t)(ncl), sizeof(double)) ) == NULL) 
		{no_sys_mem("swan --> (un)", ncl);	return;}
	if ((v = (double *) calloc ((size_t)(ncl), sizeof(double)) ) == NULL) 
		{no_sys_mem("swan --> (v)", ncl);	return;}
	if ((vn = (double *) calloc ((size_t)(ncl), sizeof(double)) ) == NULL) 
		{no_sys_mem("swan --> (vn)", ncl);	return;}
	if ((h = (double *) calloc ((size_t)(ncl), sizeof(double)) ) == NULL) 
		{no_sys_mem("swan --> (h)", ncl);	return;}
	if ((hn = (double *) calloc ((size_t)(ncl), sizeof(double)) ) == NULL) 
		{no_sys_mem("swan --> (hn)", ncl);	return;}
	if ((dxp = (double *) calloc ((size_t)(ncl+1), sizeof(double)) ) == NULL) 
		{no_sys_mem("swan --> (dxp)", ncl);	return;}
	if ((cca = (double *) calloc ((size_t)(ncl), sizeof(double)) ) == NULL) 
		{no_sys_mem("swan --> (cca)", ncl);	return;}
	if ((anglt = (double *) calloc ((size_t)(jp2), sizeof(double)) ) == NULL) 
		{no_sys_mem("swan --> (anglt)", jp2);	return;}
	/*if (polar != 0) {
		if ((xpp = (double *) calloc ((size_t)(ncl+1), sizeof(double)) ) == NULL) 
			{no_sys_mem("swan --> (xpp)", ncl);	return;}
	} */
	if (max_level) {
		/* Not shure that h_bak is necessary */
		if ((h_bak = (double *) calloc ((size_t)(ncl), sizeof(double)) ) == NULL) 
			{no_sys_mem("swan --> (h_bak)", ncl);	return;}
		if ((zm = (double *) calloc ((size_t)(ncl), sizeof(double)) ) == NULL) 
			{no_sys_mem("swan --> (zm)", ncl);	return;}
	}
	if (cumpt && !maregs_in_input) {
		if ((lcum_p = (int *) calloc ((size_t)(n_mareg), sizeof(int)) ) == NULL) 
			{no_sys_mem("swan --> (lcum_p)", n_mareg);	return;}
		if ((cum_p = (double *) calloc ((size_t)(n_mareg*n_ptmar), sizeof(double)) ) == NULL) 
			{no_sys_mem("swan --> (cum_p)", n_mareg*n_ptmar);	return;}
		if ((time_p = (float *) calloc ((size_t)(n_ptmar), sizeof(float)) ) == NULL) 
			{no_sys_mem("swan --> (time_p)", n_ptmar);	return;}
	}
	if (cumpt && maregs_in_input) {		/* lcum_p is already allocated and filled */
		if ((cum_p = (double *) calloc ((size_t)(n_mareg*n_ptmar), sizeof(double)) ) == NULL) 
			{no_sys_mem("swan --> (cum_p)", n_mareg*n_ptmar);	return;}
		if ((time_p = (float *) calloc ((size_t)(n_ptmar), sizeof(float)) ) == NULL) 
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
			h[i*nx+j] = inicial[j*ny+i];
			hn[i*nx+j] = inicial[j*ny+i];
		}
	}

	/* Test that the maregraphs coords -> indeces conversion is correct */
	/*mexPrintf("Adjusted maregraph positions\n");
	mexPrintf("old x\told y\tnew x\tnew y\tdepth\n");
	for (i = 0; i < n_mareg; i++)
		mexPrintf("%.4f\t%.4f\t%.1f\n",x,y,-dep[lcum_p[i]]);*/

	dangx = dminx / 60.;	dangy = dminy / 60.;
	/*     Polar option (If Polar is 2 then MERCATOR)*/
	if (polar != 0) {
		anglt[0] = y_min;
		for (i = 0; i < jp1; i++) {
			if (polar > 0) anglt[i+1] = anglt[i] + dangy;
			if (polar < 0) anglt[i+1] = anglt[i] - dangy;
			if (polar < 0 && anglt[i+1] < 0.) goto L54;
			angltt = anglt[i] * D2R;
			cang = cos(angltt);
			if (polar == 2) anglt[i+1] = anglt[i] + dangy * cang;
			goto L53;
L54:
			polar = 1;
			anglt[i+1] = anglt[i] + dangy;
L53:;
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
		/*free ((void *) xpp);*/
	}
	for (i = 0; i < ncl; i++) {
		cca[i] = cc;
		if (polar == 0) dxp[i] = dx;	/* Do not set if polar */
	}

	dtx = (polar == 0) ? dx : dangx;	/* If polar == 0 dx and dy are already in meters */
	dty = (polar == 0) ? dy : dangy;	/* like they must be. Otherwise, they are in degrees */
	lcum = 0;	cycle = 1;	time_h = 0.;

	lhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
	rhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
	rhs1[0] = mxCreateDoubleMatrix(1, 1, mxREAL);	/* For updating the bar in waitbar */
	ptr = mxGetPr(rhs[1]);
	ptr1 = mxGetPr(rhs1[0]);
	tmp_ptr[0] = 0;					/* Start the waitbar with zero length */
	memcpy(ptr, tmp_ptr, 8);
	rhs[0] = mxCreateString(w_bar_title);	/* multiwaitbar message */
	rhs1[1] = mxCreateString(w_bar_title);	/* Waitbar message */
	mexCallMATLAB(1,lhs,2,rhs1,"waitbar");
	h_bar = mxGetPr(lhs[0]);				/* Save the waitbar handle */

	/* Declarations for the (if) movie option */
	if (movie && movie_char) {
		mov_8 = mxCalloc(ncl, sizeof(char));
		mov_8_tmp = mxCalloc(ncl, sizeof(char));
		dims[0] = ny;	dims[1] = nx;	dims[2] = (int)(n_of_cycles / grn + 2);
		plhs[0] = mxCreateNumericArray(3, dims, mxUINT8_CLASS, mxREAL);
		ptr_mov_8 = (unsigned char *)mxGetData(plhs[0]);
	}
	else if (movie && movie_float) {
		mov_32 = (float *)mxCalloc(ncl, sizeof(float));
		dims[0] = ny;	dims[1] = nx;	dims[2] = (int)(n_of_cycles / grn + 2);
		plhs[0] = mxCreateNumericArray(3, dims, mxSINGLE_CLASS, mxREAL);
		ptr_mov_32 = (float *)mxGetData(plhs[0]);
	}

	for (k = 0; k < n_of_cycles; k++) {
		if (cycle % 10 == 0) {
			tmp_ptr[0] = (double)cycle/(double)n_of_cycles;
			memcpy(ptr1, tmp_ptr, 8);
			mexCallMATLAB(0,NULL,1,rhs1,"waitbar");
			h_bar = mxGetPr(lhs[0]);			/* Save the waitbar handle */
			/*mexPrintf("Ciclo %d \tde %d\t %lg\r", cycle, n_of_cycles, *h_bar); */
		}
		time_h += dt;
		uvh_(dep, r, rn, u, un, v, vn, h, hn, dxp, cca);
		if (max_level) {
			change(h, h_bak);	max_z(zm, h_bak);
		}
		if (cumpt) {			/* Want time series at maregraph positions */
			if (cycle % cumint == 0) {	/* Save heights at cumint intervals */
				for (i = 0; i < n_mareg; i++) {
					cum_p[ijc(lcum,i)] = h[lcum_p[i]-1];
					time_p[lcum] = time_h;
				}
				lcum++;
			}
		}
		if ((k % grn) == 0 || k == n_of_cycles - 1) {
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
				for (i = 0; i < ncl; i++) {
					mov_8_tmp[i] = (unsigned char)( ((work[i] - work_max) / dz) * 254 + 1);
				}
				/* Transpose to matlab ordering */
				for (i = 0, i2 = ny - 1; i < ny; i++, i2--)
				for (j = 0; j < nx; j++) mov_8[j*ny+i2] = mov_8_tmp[i*nx+j];

				/*mexPrintf("n_frames = %d\tplanos = %d\n",n_frames,dims[2]); */
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
			if (write_grids && w_ascii) {
				write_grd_ascii ((int)time_h, x_min, y_min, dtx, dty, ip2, jp2, work);
				mexPrintf ("\t\t\t\tWrote grelha %d\n", (int)time_h);
			}
			else if (write_grids && !w_ascii) {
				write_grd_bin ((int)time_h, x_min, y_min, dtx, dty, ip2, jp2, work);
				mexPrintf ("\t\t\t\tWrote grelha %d\n", (int)time_h);
			}
		}
		cycle++;
	}

	mxSetPr(rhs[0],h_bar);	mxSetPr(rhs[1],NULL);	mxSetPr(lhs[0],NULL);
	mexCallMATLAB(0,lhs,1,rhs,"close");

	/*     WRITE THE OUTPUT ON FILE CUMHT */
	if (cumpt) {
		for (j = 0; j < lcum; j++) {
			fprintf (fp, "%lg", time_p[j]);
			for (i = 0; i < n_mareg; i++)
				fprintf (fp, "\t%.4f", cum_p[ijc(j,i)]);
			fprintf (fp, "\n");
		}
	}
	/* Clean up allocated memory. */
	mxDestroyArray(rhs[0]);
	mxDestroyArray(rhs[1]);
	mxDestroyArray(rhs1[0]);
	mxDestroyArray(rhs1[1]);
	mxDestroyArray(lhs[0]);
	if (movie && movie_char) {
		mxFree(mov_8);	mxFree(mov_8_tmp);
	}
	if (movie && movie_float) {
		mxFree(mov_32);
	}

	if (cumpt) fclose (fp);
	if (!source_in_input) {
		free ((void *) inicial);
	}
	free ((void *) work);	free ((void *) dep);
	free ((void *) r);	free ((void *) h);
	free ((void *) rn);	free ((void *) u);	free ((void *) un);
	free ((void *) vn);	free ((void *) hn);	free ((void *) dxp);
	free ((void *) cca);	free ((void *) anglt);
	if (max_level) {
		free ((void *) zm);	free ((void *) h_bak); 
	}
	if (cumpt) {
		free ((void *) lcum_p);	free((void *) cum_p);	free ((void *) time_p);	 
	}
}

/* --------------------------------------------------------------------------- */

int bndy_(double *dep, double *r, double *rn, double *u, double *un, double *v, double *vn, double *h, double *hn, double *dxp) {
	double tmp, mul;
	int i, j;

	/*     left boundary */
	for (j = 1; j < jp2; j++) {
		tmp = d_sqrt (grav * dep[ij(0,j)]);
		if (indl == 1) goto L111;
		if (indl == 2) goto L112;
		/*     reflective */
		u[ij(0,j)] = 0.;
		v[ij(0,j)] = v[ij(1,j)];
		h[ij(0,j)] = h[ij(1,j)];
		goto L110;
		/*     continuative */
L112:
		if (polar != 0) dx = dxp[ij(0,j)];
		if (polar == 2) dy = dxp[ij(0,j)];	/* MERCATOR */
		/*if (dep[ij(0,j)] < 0.) goto L313;	/* This is bad */
		if (fabs(dep[ij(0,j)]) > 0.005) dep[ij(0,j)] = 0.;	/* THIS SEAMS TO PREVENT BORDER INSTABILITIES */
		mul = tmp * dt / dx;
		v[ij(0,j)] += (v[ij(1,j)] - v[ij(0,j)]) * mul;
		if (v[ij(1,j)] == v[ij(2,j)])
			v[ij(0,j)] = v[ij(1,j)];
		u[ij(0,j)] += (u[ij(1,j)] - u[ij(0,j)]) * mul;
		if (u[ij(1,j)] == u[ij(2,j)])
			u[ij(0,j)] = u[ij(1,j)];
		h[ij(0,j)] += (h[ij(1,j)] - h[ij(0,j)]) * mul;
		if (h[ij(1,j)] == h[ij(2,j)])		/* More useless comparison between floats */
			h[ij(0,j)] = h[ij(1,j)];
		if (iopt == 0) goto L210;
L313:
		v[ij(0,j)] = v[ij(1,j+0)];	/* Those 3 would blow te program (v[ij(1,jp2)] doesn't exist) */
		u[ij(0,j)] = u[ij(1,j+0)];
		h[ij(0,j)] = h[ij(1,j+0)];
L210:
		goto L110;
		/*      piston */
L111:
		u[ij(0,j)] = pistal * sin(pistbl * time_h);
		v[ij(0,j)] = v[ij(1,j)];
		if (dep[ij(0,j)] > 0.)
			h[ij(0,j)] = u[ij(0,j)] * tmp / grav;
L110:;
	}
	/*     bottom boundary */
	for (i = 0; i < ip2; i++) {
		tmp = d_sqrt (grav * dep[ij(i,0)]);
		if (indb == 1) goto L121;
		if (indb == 2) goto L122;
		/*     reflective */
		v[ij(i,0)] = 0.;
		u[ij(i,0)] = u[ij(i,1)];
		h[ij(i,0)] = h[ij(i,1)] * 2. - h[ij(i,2)];
		goto L120;
		/*     continuative */
L122:
		if (polar != 0) dx = dxp[ij(i,0)];	/* THIS WAS NOT ON THE ORIGINAL CODE. FORGOTTEN? */
		if (polar == 2) dy = dxp[ij(i,0)];	/* MERCATOR */
		/*if (dep[ij(i,0)] < 0.) goto L323;	/* This is bad */
		if (fabs(dep[ij(i,0)]) > 0.005) dep[ij(i,0)] = 0.;	/* THIS SEAMS TO PREVENT BORDER INSTABILITIES */
		mul = tmp * dt / dy;
		v[ij(i,0)] += (v[ij(i,1)] - v[ij(i,0)]) * mul;
		if (v[ij(i,1)] == v[ij(i,2)]) 
			v[ij(i,0)] = v[ij(i,1)];
		u[ij(i,0)] += (u[ij(i,1)] - u[ij(i,0)]) * mul;
		if (u[ij(i,1)] == u[ij(i,2)]) 
			u[ij(i,0)] = u[ij(i,1)];
		h[ij(i,0)] += (h[ij(i,1)] - h[ij(i,0)]) * mul;
		if (h[ij(i,1)] == h[ij(i,2)]) 
			h[ij(i,0)] = h[ij(i,1)];
		if (iopt == 0) goto L211;
L323:
		v[ij(i,0)] = v[ij(i,1)];
		u[ij(i,0)] = u[ij(i,1)];
		h[ij(i,0)] = h[ij(i,1)];
L211:
		goto L120;
		/*     piston */
L121:
		v[ij(i,0)] = pistab * sin(pistbb * time_h);
		u[ij(i,0)] = u[ij(i,1)];
		if (dep[ij(i,0)] > 0.)
			h[ij(i,0)] = v[ij(i,0)] * tmp / grav;
L120:;

	}
	/*     right boundary */
	for (j = 1; j < jp2; j++) {
		tmp = d_sqrt (grav * dep[ij(ip2-1,j)]);
		if (indr == 1) goto L131;
		if (indr == 2) goto L132;
		/*     reflective */
		u[ij(ip2-1,j)] = 0.;
		v[ij(ip2-1,j)] = v[ij(ip1-1,j)];
		h[ij(ip2-1,j)] = h[ij(ip1-1,j)];
		goto L130;
		/*     continuative */
L132:
		if (polar != 0) dx = dxp[ij(ip1-1,j)];
		if (polar == 2) dy = dxp[ij(ip1-1,j)];	/* MERCATOR */
		/*if (dep[ij(ip2-1,j)] < 0.) goto L333;	/* This is bad */
		if (fabs(dep[ij(ip2-1,j)]) > 0.005) dep[ij(ip2-1,j)] = 0.;	/* THIS SEAMS TO PREVENT BORDER INSTABILITIES */
		mul = tmp * dt / dx;
		v[ij(ip2-1,j)] += (v[ij(ip1-1,j)] - v[ij(ip2-1,j)]) * mul;
		if (v[ij(ip1-1,j)] == v[ij(ip-1,j)])
			v[ij(ip2-1,j)] = v[ij(ip1-1,j)];
		u[ij(ip2-1,j)] += (u[ij(ip1-1,j)] - u[ij(ip2-1,j)]) * mul;
		if (u[ij(ip1-1,j)] == u[ij(ip-1,j)])
			u[ij(ip2-1,j)] = u[ij(ip1-1,j)];
		h[ij(ip2-1,j)] += (h[ij(ip1-1,j)] - h[ij(ip2-1,j)]) * mul;
		if (h[ij(ip1-1,j)] == h[ij(ip-1,j)])
			h[ij(ip2-1,j)] = h[ij(ip1-1,j)];
		if (iopt == 0) goto L212;
L333:
		v[ij(ip2-1,j)] = v[ij(ip1-1,j)];
		u[ij(ip2-1,j)] = u[ij(ip1-1,j)];
		h[ij(ip2-1,j)] = h[ij(ip1-1,j)];
L212:
		goto L130;
		/*     piston */
L131:
		u[ij(ip2-1,j)] = -(pistar * sin(pistbr * time_h));
		v[ij(ip2-1,j)] = v[ij(ip1-1,j)];
		/*     Changes of 5/91 */
		u[ij(ip1-1,j)] = u[ij(ip2-1,j)];
		if (dep[ij(ip2-1,j)] > 0.) {
			h[ij(ip2-1,j)] = -(u[ij(ip2-1,j)] * tmp) / grav;
			h[ij(ip1-1,j)] = h[ij(ip2-1,j)];
		}
L130:;
	}
	/*     top boundary */
	/*      do top only to ip instead of ip2    ** */
	/*       Changed from 2 to 1 on 6/90 */
	for (i = 0; i < ip1; i++) {
		tmp = d_sqrt (grav * dep[ij(i,jp2-1)]);
		if (indt == 1) goto L141;
		if (indt == 2) goto L142;
		/*     reflective */
		v[ij(i,jp2-1)] = 0.;
		u[ij(i,jp2-1)] = u[ij(i,jp1-1)];
		h[ij(i,jp2-1)] = h[ij(i,jp1-1)];
		goto L140;
		/*     continuative */
L142:
		if (polar != 0) dx = dxp[ij(i,jp2-1)];	/* THIS WAS NOT ON THE ORIGINAL CODE. FORGOTTEN? */
		if (polar == 2) dy = dxp[ij(i,jp2-1)];	/* MERCATOR */
		/*if (dep[ij(i,jp2-1)] < 0.) goto L343;	/* This is bad */
		if (fabs(dep[ij(i,jp2-1)]) > 0.005) dep[ij(i,jp2-1)] = 0.;	/* THIS SEAMS TO PREVENT BORDER INSTABILITIES */
		mul = tmp * dt / dy;
		v[ij(i,jp2-1)] += (v[ij(i,jp1-1)] - v[ij(i,jp2-1)]) * mul;
		if (v[ij(i,jp1-1)] == v[ij(i,jp-1)]) 
			v[ij(i,jp2-1)] = v[ij(i,jp1-1)];
		u[ij(i,jp2-1)] += (u[ij(i,jp1-1)] - u[ij(i,jp2-1)]) * mul;
		if (u[ij(i,jp1-1)] == u[ij(i,jp-1)]) 
			u[ij(i,jp2-1)] = u[ij(i,jp1-1)];
		h[ij(i,jp2-1)] += (h[ij(i,jp1-1)] - h[ij(i,jp2-1)]) * mul;
		if (h[ij(i,jp1-1)] == h[ij(i,jp-1)]) 
			h[ij(i,jp2-1)] = h[ij(i,jp1-1)];
		if (iopt == 0) goto L213;
L343:
		v[ij(i,jp2-1)] = v[ij(i,jp-1)];		/* jp doesn't make sense. Perhaps jp1? */
		u[ij(i,jp2-1)] = u[ij(i,jp-1)];
		h[ij(i,jp2-1)] = h[ij(i,jp-1)];
L213:
		goto L140;
		/*     piston */
L141:
		v[ij(i,jp2-1)] = -(pistat * sin(pistbt * time_h));
		u[ij(i,jp2-1)] = u[ij(i,jp1-1)];
		/*      Changed 5/91 */
		v[ij(i,jp1-1)] = v[ij(i,jp2-1)];
		/*       THE HEIGHT sign  used to be +     4/23/91 */
		if (dep[ij(i,jp2-1)] > 0.)
			h[ij(i,jp2-1)] = -(v[ij(i,jp2-1)] * tmp) / grav;
L140:;
	}
	return 0;
}

/* --------------------------------------------------------------------------- */

int uvh_(double *dep, double *r, double *rn, double *u, double *un, double *v, double *vn, double *h, double *hn, double *dxp, double *cca) {

	double cang, ck, sa, sb, sang, cang1, tmp, thu, thv;
	double td, tu, tv, angltt, td1, tu1, tv1, tu2, tv2, dph;
	int	i, j;

	sb = 0.;	sa = 0.;
	for (j = 1; j < jp1; j++) {
		for (i = 1; i < ip1; i++) {
			/* *****5/91*POLAR Each cell has its own dx and cos angle (cang) */
			if (polar != 0) dx = dxp[ij(i,j)];
			if (polar == 2) dy = dxp[ij(i,j)];	/* MERCATOR */
			cang1 = 1.;		cang = 1.;
			if (polar != 0) cang = dxp[ij(i,j)] / (dangx * m_per_deg);
			if (polar != 0) cang1 = dxp[ij(i,j+1)] / (dangx * m_per_deg);
			/*     donor cell difference */
			/*     will get diffusion if time step is too small */
			td1 = dep[ij(i+1,j)] + h[ij(i+1,j)] - r[ij(i+1,j)];
			td = dep[ij(i,j)] + h[ij(i,j)] - r[ij(i,j)];
			tv1 = dep[ij(i,j+1)] + h[ij(i,j+1)] - r[ij(i,j+1)];
			tv = td;
			if (u[ij(i+1,j)] > 0.)
				td1 = dep[ij(i,j)] + h[ij(i,j)] - r[ij(i,j)];
			if (u[ij(i,j)] > 0.)
				td = dep[ij(i-1,j)] + h[ij(i-1,j)] - r[ij(i-1,j)];
			if (v[ij(i,j+1)] > 0.)
				tv1 = dep[ij(i,j)] + h[ij(i,j)] - r[ij(i,j)];
			if (v[ij(i,j)] > 0.)
				tv = dep[ij(i,j-1)] + h[ij(i,j-1)] - r[ij(i,j-1)];
			/*      Special for Flooding */
			if (td1 < 0.) td1 = 0.;	if (td < 0.) td = 0.;
			if (tv1 < 0.) tv1 = 0.;	if (tv < 0.) tv = 0.;
			hn[ij(i,j)] = h[ij(i,j)] - dt * ((u[ij(i+1,j)] * td1 - u[ij(i,j)] * td) / 
				dx + (v[ij(i,j+1)] * cang1 * tv1 - v[ij(i,j)] * cang * tv) /
				(cang * dy)) + (rn[ij(i,j)] - r[ij(i,j)]);
			/*    ROUGH is factor for surface roughness = actual height/ideal height */
			if (dep[ij(i,j)] < 0.) 
				hn[ij(i,j)] *= rough;
		}
	}
	/*      SPECIAL */
	/* This is completly idiot. Comparing dep to 0. ??? I'll comment this section. JL 12-11-04 */ 
	/*for (j = 1; j < jp1; j++) {
		for (i = 1; i < ip1; i++) { */
			/*     zero depth boundary treated as reflective  10/1/89 method */
			/*if (dep[ij(i-1,j)] == 0. && dep[ij(i,j)] != 0.)
				hn[ij(i-1,j)] = hn[ij(i,j)];
			if (dep[ij(i+1,j)] == 0. && dep[ij(i,j)] != 0.)
				hn[ij(i+1,j)] = hn[ij(i,j)];
			if (dep[ij(i,j-1)] == 0. && dep[ij(i,j)] != 0.)
				hn[ij(i,j-1)] = hn[ij(i,j)];
			if (dep[ij(i,j+1)] == 0. && dep[ij(i,j)] != 0.)
				hn[ij(i,j+1)] = hn[ij(i,j)];
			if (dep[ij(i-1,j)] == 0. && dep[ij(i,j)] == 0. && dep[ij(i+1,j)] == 0.
		    	&& dep[ij(i,j-1)] == 0. && dep[ij(i,j+1)] == 0.)
				hn[ij(i,j)] = 0.;
		}
	} */
	for (j = 1; j < jp1; j++) {
		for (i = 1; i < ip1; i++) {
		    	h[ij(i,j)] = hn[ij(i,j)];
		}
	}
	bndy_(dep, r, rn, u, un, v, vn, h, hn, dxp);
	for (j = 1; j < jp1; j++) {
		angltt = anglt[j] * D2R;
		sang = sin(angltt);
		if (polar != 0) cf = sang * 1.454e-4;
		for (i = 1; i < ip1; i++) {
			/* **********POLAR Each cell has its own dx  5/91 */
			if (polar != 0) dx = dxp[ij(i,j)];
			if (polar == 2) dy = dxp[ij(i,j)];	/* MERCATOR */
			if (cc == 0.) goto L10;
			/*       Special for Flooding */
			dph = dep[ij(i,j)] + h[ij(i,j)];
			if (dph <= .1) goto L10;
			ck = cca[ij(i,j)];
			tmp = d_sqrt(u[ij(i,j)]*u[ij(i,j)] + v[ij(i,j)]*v[ij(i,j)]) / 
			      (ck * ck * (dep[ij(i,j)] + h[ij(i,j)]));
			sb = grav * u[ij(i,j)] * tmp;
			/*      CORRECTED ll/20/87 */
			sa = grav * v[ij(i,j)] * tmp;
L10:
			tu1 = u[ij(i+1,j)] - u[ij(i,j)];
			tu2 = u[ij(i,j+1)] - u[ij(i,j)];
			tv = (v[ij(i,j)] + v[ij(i,j+1)] + v[ij(i-1,j+1)] + v[ij(i-1,j)]) / 4.;
			if (u[ij(i,j)] > 0.) tu1 = u[ij(i,j)] - u[ij(i-1,j)];
			if (tv > 0.) tu2 = u[ij(i,j)] - u[ij(i,j-1)];
			thu = h[ij(i,j)] - h[ij(i-1,j)];
			tv1 = v[ij(i+1,j)] - v[ij(i,j)];
			tv2 = v[ij(i,j+1)] - v[ij(i,j)];
			tu = (u[ij(i,j)] + u[ij(i+1,j)] + u[ij(i,j-1)] + u[ij(i+1,j-1)]) / 4.;
			if (tu > 0.) tv1 = v[ij(i,j)] - v[ij(i-1,j)];
			if (v[ij(i,j)] > 0.) tv2 = v[ij(i,j)] - v[ij(i,j-1)];
			/*      CORRECTED ll/20/87 */
			thv = h[ij(i,j)] - h[ij(i,j-1)];
			un[ij(i,j)] = u[ij(i,j)] - dt * (u[ij(i,j)] * tu1 / dx + tv * tu2 / dy) - grav * 
				dt * thu / dx - dt * (sb - cf * v[ij(i,j)] - sfx);
			vn[ij(i,j)] = v[ij(i,j)] - dt * (tu * tv1 / dx + v[ij(i,j)] * tv2 / dy) - grav * 
				dt * thv / dy - dt * (sa + cf * u[ij(i,j)] - sfy);
			/*     FLOODING--Assumes above sea level is a negative depth-- */
		}
	}
	/*     update */
	/* This is completly idiot. Comparing dep to 0. ??? I'll comment this section. JL 12-11-04 */ 
	/*for (j = 1; j < jp1; j++) {
		for (i = 0; i < ip1; i++) { */
		/*     zero depth boundary treated as reflective  10/1/89 method */
		/*if (dep[ij(i-1,j)] == 0. && dep[ij(i,j)] != 0.)
			un[ij(i-1,j)] = 0.;
		if (dep[ij(i+1,j)] == 0. && dep[ij(i,j)] != 0.)
			un[ij(i+1,j)] = 0.;
		if (dep[ij(i,j-1)] == 0. && dep[ij(i,j)] != 0.)
			vn[ij(i,j-1)] = 0.;
		if (dep[ij(i,j+1)] == 0. && dep[ij(i,j)] != 0.)
			vn[ij(i,j+1)] = 0.;
		if (dep[ij(i-1,j)] == 0. && dep[ij(i,j)] != 0.)
			vn[ij(i-1,j)] = v[ij(i,j)];
		if (dep[ij(i+1,j)] == 0. && dep[ij(i,j)] != 0.)
			vn[ij(i+1,j)] = v[ij(i,j)];
		if (dep[ij(i,j-1)] == 0. && dep[ij(i,j)] != 0.)
			un[ij(i,j-1)] = un[ij(i,j)];
		if (dep[ij(i,j+1)] == 0. && dep[ij(i,j)] != 0.)
			un[ij(i,j+1)] = un[ij(i,j)];
		}
	} */
	for (j = 1; j < jp1; j++) {
		for (i = 1; i < ip1; i++) {
			r[ij(i,j)] = rn[ij(i,j)];
			u[ij(i,j)] = un[ij(i,j)];
			v[ij(i,j)] = vn[ij(i,j)];
			/*     artificial limits */
			/* Another idiotic comparison. Original code compared if ... < 1e-8, while using float*/
			if (fabs(u[ij(i,j)]) < 1e-6) u[ij(i,j)] = 0.;
			if (fabs(v[ij(i,j)]) < 1e-6) v[ij(i,j)] = 0.;
			if (fabs(h[ij(i,j)]) < 1e-6) h[ij(i,j)] = 0.;
		}
	}
	return 0;
}

/* --------------------------------------------------------------------------- */

int write_grd_ascii (int kk, double x_min, double y_min, double dtx, double dty, int i_end, int j_end, float *work) {

	/* Writes a grid in the Surfer ascii format */
	int i, j;
	double x_max, y_max;
	float work_min = FLT_MAX, work_max = -FLT_MAX;
	char name[80];
	FILE *fp;

	if (stem[0] == 0)
		sprintf (name,"%d%s\0", kk,".grd");
	else
		sprintf (name, "%s%d%s", stem, kk,".grd");

	/*sprintf (name, "%d\0", kk); */
	/*strcat (name, ".grd"); */

	if ((fp = fopen (name, "w")) == NULL) {
		fprintf (stderr, "%s: Unable to create file %s - exiting\n", "swan", name);
		return (-1);
	}

	x_max = x_min + (i_end - 1) * dtx;
	y_max = y_min + (j_end - 1) * dty;

	/* Find zmin/zmax */
	for(i = 0; i < i_end; i++) {
		for(j = 0; j < j_end; j++) {
			work_max = MAX(work[ijs(i,j)], work_max);
			work_min = MIN(work[ijs(i,j)], work_min);
		}
	}

	fprintf (fp, "DSAA\n");
	fprintf (fp, "%d %d\n", i_end, j_end);
	fprintf (fp, "%lg %lg\n", x_min, x_max);
	fprintf (fp, "%lg %lg\n", y_min, y_max);
	fprintf (fp, "%lg %lg\n", work_min, work_max);
	for(j = 0; j < j_end; j++) {
		for(i = 0; i < i_end; i++) {
			fprintf (fp, "%lg\n", work[ijs(i,j)]);
		}
	}
	fclose(fp);
	return (0);
} 

int write_grd_bin (int kk, double x_min, double y_min, double dtx, double dty, int i_end, int j_end, float *work) {

	/* Writes a grid in the Surfer binary format */
	int i, j;
	double x_max, y_max;
	float work_min = FLT_MAX, work_max = -FLT_MAX;
	char name[80];
	struct srf_header h;
	FILE *fp;

	if (stem[0] == 0)
		sprintf (name,"%d%s\0", kk,".grd");
	else
		sprintf (name, "%s%d%s", stem, kk,".grd");

	/*sprintf (name, "%d\0", kk); */
	/*strcat (name, ".grd"); */

	if ((fp = fopen (name, "wb")) == NULL) {
		mexPrintf("Fatal Error: Could not create file %s!\n", name);
		return (-1);
	}

	x_max = x_min + (i_end - 1) * dtx;
	y_max = y_min + (j_end - 1) * dty;

	/* Find zmin/zmax */
	for (i = 0; i < i_end*j_end; i++) {
		work_max = MAX(work[i], work_max);
		work_min = MIN(work[i], work_min);
	}

	/* store header information and array */
	strcpy (h.id,"DSBB");
	h.nx = i_end;	 	h.ny = j_end;
	h.x_min = x_min;	h.x_max = x_max;
	h.y_min = y_min;	h.y_max = y_max;
	h.z_min = (double)work_min;	h.z_max = (double)work_max;

	if (fwrite ((void *)&h, sizeof (struct srf_header), (size_t)1, fp) != 1) {
		fprintf (stderr, "Fatal Error: Error writing file %s!\n", name);
		return (-1);
	}
	for(j = 0; j < j_end; j++) {
		for(i = 0; i < i_end; i++) {
			fwrite ((void *)&work[ijs(i,j)], sizeof(float), (size_t)1, fp);
		}
	}
	fclose(fp);
	return (0);
}

int read_grd_info_ascii (char *file, struct srf_header *hdr) {

	char line[512], id[5];
	FILE *fp;

	if ((fp = fopen (file, "r")) == NULL) {
		mexPrintf ("%s: Unable to read file %s - exiting\n", "swan", file);
		return (-1);
	}

	fgets (line, 512, fp);
	sscanf (line, "%s", &hdr->id);
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

int read_grd_ascii (char *file, struct srf_header *hdr, double *work) {

	/* Reads a grid in the Surfer ascii format */
	int i = 0, j, n_field;
	char *p, buffer[512], line[512];
	FILE *fp;

	if ((fp = fopen (file, "r")) == NULL) {
		mexPrintf ("%s: Unable to read file %s - exiting\n", "swan", file);
		return (-1);
	}

	fgets (line, 512, fp);
	sscanf (line, "%s", &hdr->id);
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

int count_col (char *line) {
	/* Count # of fields contained in line */
	int n_col = 0;
	char *p;

	p = (char *)strtok (line, " \t\n\015\032");
	while (p) {	/* Count # of fields */
		n_col++;
		p = (char *)strtok ((char *)NULL, " \t\n\015\032");
	}
	return (n_col);
}

int read_header_bin (FILE *fp, struct srf_header *hdr) {
	/* Reads the header of a binary Surfer gridfile */
	fread ((void *)hdr, sizeof (struct srf_header), (size_t)1, fp); 
	return (0);
}

int read_grd_bin (char *file, struct srf_header *hdr, double *work) {
	int	i, j, ij, kk;
	float	*tmp;			/* Array pointer for reading in rows of data */
	FILE	*fp;

	if ((fp = fopen (file, "rb")) == NULL) {
		mexPrintf ("%s: Unable to read file %s - exiting\n", "swan", file);
		return (-1);
	}
	fread ((void *)hdr, sizeof (struct srf_header), (size_t)1, fp); 

	/* Allocate memory for one row of data (for reading purposes) */
	tmp = (float *) calloc ((size_t)hdr->nx, sizeof (float));
	for (j = 0; j < (hdr->ny - 1); j++) {
		fread (tmp, sizeof (float), (size_t)hdr->nx, fp);	/* Get one row */
		ij = j * hdr->nx;
		for (i = 0; i < hdr->nx; i++) {
			kk = ij + i;
			work[kk] = tmp[i];
		}
	}
	fclose(fp);
	free ((void *)tmp);
	return (0);
}

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

int count_n_maregs(char *file) {
	int	i = 0;
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

int read_maregs(char *file) {
	/* Read maregraph positions and convert them to vector indices */
	int	i = 0, ix, jy;
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
		sscanf (line, "%f %f", &x, &y);
		ix = irint((x - hdr_b.x_min) / dx);
		jy = irint((y - hdr_b.y_min) / dy); 
		lcum_p[i] = jy * ip2 + ix; 
		i++;
	}
	fclose (fp);
	return (0);
}

void	no_sys_mem (char *where, int n) {	
		mexPrintf ("Fatal Error: %s could not allocate memory, n = %d\n", where, n);
}

/*     ********* CHECK OF MAXIMUM VALUE *********** */
void max_z (double *zm, double *h_bak) {
	int i, ij;

	ij = ip2 * jp2;
	for (i = 0; i < ij; i++)
		if(zm[i] < h_bak[i]) zm[i] = h_bak[i];
}

void change (double *h, double *h_bak) {
	int i, ij;

	ij = ip2 * jp2;
	for (i = 0; i < ij; i++)
		h_bak[i] = h[i];
}
