/*
This program extracts data from level-3 HDF bin files

Regions of interest can be specified by longitude and latitude
boundaries or by a central coordinate and a radius in kilometers.

Norman Kuring		14-Feb-2003	Original development
Norman Kuring		11-Dec-2007	Fix memory-overrun bug and add a
					couple of calls to free().

 *	Mexified (and simplified to output only sum / weights) By
 *	Joaquim Luis - Dec 2009
*/

#include "mex.h"
#include <math.h>
#include "hdf.h"

#define PI	3.1415926535897932384626433832795029L
#define D2R	PI / 180.

#define NUMROWS		2160
#define EARTH_RADIUS	6371.229

static int	basebin[NUMROWS];
static short int	numbin[NUMROWS];
static float	latbin[NUMROWS];

#define BLIST_FIELDS "bin_num,weights"
#define BLIST_SIZE	8
static unsigned char	blistrec[BLIST_SIZE];
static int	bin_num;
static float	weights;
static	VOIDP	bufptrs[] = {&bin_num,&weights};

#define PREC_SIZE	8

int	initbin(void);
int	latlon2bin(double lat, double lon);
void	bin2latlon(int bin, double *clat, double *clon);
short int	lat2row(double lat);
int	rowlon2bin(short int row, double lon);
int	binsearch(int bin, int vdata_id, int numrecs);
double	constrain_lat(double lat);
double	constrain_lon(double lon);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	int	i, n, b, mode = 0, numbins = 0, error = 0, *binnums = NULL, n_alloc = 10000;
	int	file_id,vdata_ref,vdata_id,numrecs,pvd_id, nTrue, *recno, totbins;
	char	*infile, *param, *param_fields;
	static unsigned char	paramrec[PREC_SIZE];
	float	*pdata_s;
	static float	sum;
	double	radius, clon, clat, *pdata_d;
	double	north, south, west, east;
	static	VOIDP	paramptrs[] = {&sum};

	if(!mxIsChar(prhs[0])) {
		mexPrintf("SWREADL3_M ERROR: First input must contain the HDF filename.\n");
		error++;
	}
	if(!mxIsChar(prhs[1])) {
		mexPrintf("SWREADL3_M ERROR: Second input must contain the name of the variable to read.\n");
		error++;
	}
	if (nrhs == 3 && !(mxGetN(prhs[2]) == 3 || mxGetN(prhs[2]) == 4)) { 
		mexPrintf("SWREADL3_M ERROR: Thirth argument must be of the form:\n");
		mexPrintf("\t\t[lon lat radius]\tor\n");
		mexPrintf("\t\t[west east shouth north]\n");
		error++;
	}

	if (error) {
		mexPrintf("Usage:\n\tout = swreadl3_m(main_file, parameter, [lon lat radius])\n or\n");
		mexPrintf("\tout = swreadl3_m(main_file, parameter, [west east south north])\n");
		mexErrMsgTxt("\n");
	}

	infile = (char *)mxArrayToString(prhs[0]);
	param  = (char *)mxArrayToString(prhs[1]);

	if (nrhs == 3) {
		pdata_d = mxGetPr(prhs[2]);
		if (mxGetN(prhs[2]) == 3) {	/* Input arguments are: main_file parameter lon lat radius */
			clon    = constrain_lon(pdata_d[0]);
			clat    = constrain_lat(pdata_d[1]);
			radius  = pdata_d[3];
			mode    = 1;
		}
		else {			/* Input arguments are: main_file parameter west east south north */
			west  = constrain_lon(pdata_d[0]);
			east  = constrain_lon(pdata_d[1]);
			south = constrain_lat(pdata_d[2]);
			north = constrain_lat(pdata_d[3]);
		}
	}
	else {
		west = -180;	east = 180;	south = -90;	north = 90;
	}

	totbins = initbin();

	if (mode == 1) {	/* Input arguments are: main_file_path parameter lon lat radius. */
		double	radius_degrees;

		radius_degrees = (radius/EARTH_RADIUS) / D2R;
		if (radius_degrees > 180) {		/* The entire globe has been selected. */
			binnums = (int *)realloc(binnums,totbins*sizeof(int));
			if (binnums == NULL) {
				mexPrintf("-E- %s line %d: Memory allocation failed.\n", __FILE__,__LINE__);
				mexErrMsgTxt("\n");
			}
			for (b = 1; b <= totbins; b++) binnums[numbins++] = b;
		}
		else {
			double	sin_c_over_2;
			short int	n_row, s_row, row;

			sin_c_over_2 = sin((radius/EARTH_RADIUS)/2);

			north  = clat + radius_degrees;
			south  = clat - radius_degrees;
			binnums = (int *)mxMalloc(n_alloc * sizeof(int));

			if (north > 90) {
				/* Chosen radius extends north to pole and over the other side to some
				latitude south of the north pole.  Set the new north to this new
				latitude.  All bins north of this latitude get selected. */

				north = 180 - north;
				n_row = lat2row(north);
				n = totbins - basebin[n_row] + 1;
				if ((numbins + n) > n_alloc) {
					n_alloc = numbins + n; 
					if ((binnums = (int *)mxRealloc(binnums, n_alloc * sizeof(int))) == NULL)
						mexErrMsgTxt("SWREADL3B_M: Could not realocate more memory\n");
				}
				for (b = basebin[n_row]; b <= totbins; b++) binnums[numbins++] = b;
			}
			else
				n_row = lat2row(north) + 1;

			if (south < -90) {
				/* Chosen radius extends south to pole and over the other side to some
				latitude north of the south pole.  Set the new south to this new
				latitude.  All bins south of this latitude get selected.  */
				int	b;

				south = -180 - south;
				s_row = lat2row(south);
				n = basebin[s_row] + numbin[s_row];
				if ((numbins + n) > n_alloc) {
					n_alloc = numbins + n; 
					if ((binnums = (int *)mxRealloc(binnums, n_alloc * sizeof(int))) == NULL)
						mexErrMsgTxt("SWREADL3B_M: Could not realocate more memory\n");
				}
				for(b = 1; b <= n; b++) binnums[numbins++] = b;
			}
			else
				s_row = lat2row(south) - 1;

			for (row = s_row + 1; row < n_row; row++) {
				double	deltalon, sin_deltalat_over_2;

				sin_deltalat_over_2 = sin((latbin[row] - clat)*D2R / 2);

				/* The following equation is a rearranged version of
				Equation 5-3a on page 30 of
				Map Projections -- A Working Manual
				by John P. Snyder
				U.S. Geological Survey Professional Paper 1395 1987
				(Fourth printing 1997) */
				deltalon = 2 * asin( sqrt( fabs( (sin_c_over_2*sin_c_over_2 - sin_deltalat_over_2*sin_deltalat_over_2)
				/( cos(latbin[row]*D2R) * cos(clat*D2R) ) )));

				deltalon /= D2R;

				west = constrain_lon(clon - deltalon);
				east = constrain_lon(clon + deltalon);

				if (east < west) {	/* User's region of interest spans the 180-degree meridian. */
					int   b,b1,b2,b3,b4,n1,n2;
					b1 = rowlon2bin(row,west);
					b2 = rowlon2bin(row, 180);
					b3 = rowlon2bin(row,-180);
					b4 = rowlon2bin(row,east);
					n1 = b2 - b1 + 1;
					n2 = b4 - b3 + 1;
					if ((numbins + n1 + n2) > n_alloc) {
						n_alloc += (numbins + n1 + n2); 
						if ((binnums = (int *)mxRealloc(binnums, n_alloc * sizeof(int))) == NULL)
							mexErrMsgTxt("SWREADL3B_M: Could not realocate more memory\n");
					}
					for(b = b1; b <= b2; b++) binnums[numbins++] = b;
					for(b = b3; b <= b4; b++) binnums[numbins++] = b;
				}
				else {	/* User's region of interest does not span the 180-degree meridian. */
					int   b,b1,bn;
 
					b1 = rowlon2bin(row,west);
					bn = rowlon2bin(row,east);
					n = bn - b1 + 1;
					if ((numbins + n) > n_alloc) {
						n_alloc += (numbins + n); 
						if ((binnums = (int *)mxRealloc(binnums, n_alloc * sizeof(int))) == NULL)
							mexErrMsgTxt("SWREADL3B_M: Could not realocate more memory\n");
					}
					for (b = b1; b <= bn; b++) binnums[numbins++] = b;
				}
			}
		}
	}

	else {		/* Input arguments are: main_file_path parameter west east south north. */
		short int	n_row, s_row, row;

		if (south > north) {
			double	tmp;
			mexPrintf( "Specified south latitude is greater than specified north latitude.\n");
			mexPrintf("I will swap the two.\n");
			tmp = north;	north = south;	south = tmp;
		}

		n_row = lat2row(north);
		s_row = lat2row(south);
		binnums = (int *)mxMalloc(n_alloc * sizeof(int));

		for (row = s_row; row <= n_row; row++) {

			if (east < west) {	/* User's region of interest spans the 180-degree meridian. */
				int	b1,b2,b3,b4,n1,n2;
        			b1 = rowlon2bin(row,west);
				b2 = rowlon2bin(row, 180);
				b3 = rowlon2bin(row,-180);
				b4 = rowlon2bin(row,east);
				n1 = b2 - b1 + 1;
				n2 = b4 - b3 + 1;
				if ((numbins + n1 + n2) > n_alloc) {
					n_alloc += (numbins + n1 + n2); 
					if ((binnums = (int *)mxRealloc(binnums, n_alloc * sizeof(int))) == NULL)
						mexErrMsgTxt("SWREADL3B_M: Could not realocate more memory\n");
				}
				for (b = b1; b <= b2; b++) binnums[numbins++] = b;
				for (b = b3; b <= b4; b++) binnums[numbins++] = b;
			}
			else {	/* User's region of interest does not span the 180-degree meridian. */
				int	b1,bn;

				b1 = rowlon2bin(row,west);
				bn = rowlon2bin(row,east);
				n = bn - b1 + 1;
				if ((numbins + n) > n_alloc) {
					n_alloc += (numbins + n); 
					if ((binnums = (int *)mxRealloc(binnums, n_alloc * sizeof(int))) == NULL)
						mexErrMsgTxt("SWREADL3B_M: Could not realocate more memory\n");
				}
				for (b = b1; b <= bn; b++) binnums[numbins++] = b;
			}
		}
	}

	/* Now that I have a list of desired bins, I can extract the corresponding data from the level-3 bin files */

	/* Open the HDF file. */
	file_id = Hopen(infile, DFACC_READ, 0);
	if(file_id == FAIL){
		mexPrintf("-E- %s line %d: Hopen(%s,DFACC_READ,0) failed.\n", __FILE__,__LINE__,infile);
		mexErrMsgTxt("\n");
	}
	/* Initialize the Vdata interface. */
	if(Vstart(file_id) == FAIL){
		mexPrintf("-E- %s line %d: Vstart(%d) failed.\n", __FILE__,__LINE__,file_id);
		mexErrMsgTxt("\n");
	}

	/* Open the "BinList" Vdata. */
	vdata_ref = VSfind(file_id,"BinList");
	if(vdata_ref == FAIL){
		mexPrintf("-E- %s line %d: VSfind(%d,\"BinList\") failed.\n", __FILE__,__LINE__,file_id);
		mexErrMsgTxt("\n");
	}
	vdata_id = VSattach(file_id, vdata_ref, "r");
	if(vdata_id == FAIL){
		mexPrintf("-E- %s line %d: VSattach(%d,%d,\"r\") failed.\n", __FILE__,__LINE__,file_id,vdata_ref);
		mexErrMsgTxt("\n");
	}
	/* Find out how many bins are stored in this file. */
	numrecs = VSelts(vdata_id);
	if(numrecs == FAIL){
		mexPrintf("-E- %s line %d: VSelts(%d) failed.\n", __FILE__,__LINE__,vdata_id);
		mexErrMsgTxt("\n");
	}
	/* Set up to read the fields in the BinList Vdata records. */
	if(VSsetfields(vdata_id,BLIST_FIELDS) == FAIL){
		mexPrintf("-E- %s line %d: VSsetfields(%d,%s) failed.\n", __FILE__,__LINE__,vdata_id,BLIST_FIELDS);
		mexErrMsgTxt("\n");
	}

	/* Open the parameter-specific Vdata. */
	vdata_ref = VSfind(file_id,param);
	if(vdata_ref == 0){
		mexPrintf("-E- %s line %d: VSfind(%d,\"%s\") failed.\n", __FILE__,__LINE__,file_id,param);
		mexErrMsgTxt("\n");
	}
	pvd_id = VSattach(file_id, vdata_ref, "r");
	if(pvd_id == FAIL){
		mexPrintf("-E- %s line %d: VSattach(%d,%d,\"r\") failed.\n", __FILE__,__LINE__,file_id,vdata_ref);
		mexErrMsgTxt("\n");
	}
	/* Set up to read the fields in the parameter-specific Vdata records. */
	{
		int	len;
		len = strlen(param) + strlen("_sum,") + 1;
		param_fields = (char *)mxMalloc(len);
		strcpy(param_fields,param);
		strcat(param_fields,"_sum");
	}
	if (VSsetfields(pvd_id,param_fields) == FAIL) {
		mexPrintf("-E- %s line %d: VSsetfields(%d,%s) failed.\n", __FILE__,__LINE__,pvd_id,param_fields);
		mexErrMsgTxt("\n");
	}

	recno = (int *)mxCalloc (numbins, sizeof (int));
	for (i = 0, nTrue = 0; i < numbins; i++) {
		recno[i] = binsearch(binnums[i],vdata_id,numrecs);
		if (recno[i] >= 0) nTrue++;
	}
	plhs[0] = mxCreateNumericMatrix (nTrue,3,mxSINGLE_CLASS,mxREAL);
	pdata_s = mxGetData(plhs[0]);

	for (i = 0, n = 0; i < numbins; i++) {

		if (recno[i] >= 0) {
			/* Read the sum for the the specified parameter for this bin. */

			if (VSseek(pvd_id,recno[i]) == FAIL){
				mexPrintf("-E- %s line %d: VSseek(%d,%d) failed.\n", __FILE__,__LINE__,pvd_id,recno[i]);
				mexErrMsgTxt("\n");
			}
			if (VSread(pvd_id,paramrec,1,FULL_INTERLACE) != 1){
				mexPrintf("-E- %s line %d: ",__FILE__,__LINE__);
				mexPrintf("VSread(%d,paramrec,1,FULL_INTERLACE) failed.\n", pvd_id);
				mexErrMsgTxt("\n");
			}

			/* VSfpack() sets the global sum variable via the paramptrs pointer array. */

			if (VSfpack( pvd_id,_HDF_VSUNPACK,param_fields,paramrec,PREC_SIZE,1,NULL,paramptrs ) == FAIL){
				mexPrintf("-E- %s line %d: VSfpack(%d, ...) failed.\n",__FILE__,__LINE__, pvd_id);
				mexErrMsgTxt("\n");
			}

			/* Get the geographical coordinates associated with this bin. */
			bin2latlon(binnums[i],&clat,&clon);

			/* Output the results. */
			pdata_s[n]         = (float)clon;
			pdata_s[n+nTrue]   = (float)clat;
			pdata_s[n+2*nTrue] = (float)(sum / weights);
			n++;
		}
	}

	if(VSdetach(pvd_id) == FAIL){
		mexPrintf("-E- %s line %d: VSdetach(%d) failed.\n", __FILE__,__LINE__,pvd_id);
		mexErrMsgTxt("\n");
	}
	if(VSdetach(vdata_id) == FAIL){
		mexPrintf("-E- %s line %d: VSdetach(%d) failed.\n", __FILE__,__LINE__,vdata_id);
		mexErrMsgTxt("\n");
	}
	if(Vend(file_id) == FAIL){
		mexPrintf("-E- %s line %d: Vend(%d) failed.\n", __FILE__,__LINE__,file_id);
		mexErrMsgTxt("\n");
	}
	if(Hclose(file_id) == FAIL){
		mexPrintf("-E- %s line %d: Hclose(%d) failed.\n", __FILE__,__LINE__,file_id);
		mexErrMsgTxt("\n");
	}

	mxFree((void *)param_fields);
	mxFree((void *)binnums);
}

int binsearch(int bin, int vdata_id, int numrecs){
	int	lo, hi, mid;

	lo = 0;
	hi = numrecs - 1;
	while (lo <= hi) {
		mid = (lo + hi)/2;
		if (VSseek(vdata_id,mid) == FAIL) {
			mexPrintf("-E- %s line %d: VSseek(%d,%d) failed.\n", __FILE__,__LINE__,vdata_id,mid);
			mexErrMsgTxt("\n");
		}
		if (VSread(vdata_id,blistrec,1,FULL_INTERLACE) != 1) {
			mexPrintf("-E- %s line %d: ",__FILE__,__LINE__);
			mexPrintf("VSread(%d,blistrec,1,FULL_INTERLACE) failed.\n", vdata_id);
			mexErrMsgTxt("\n");
		}

		/* VSfpack() sets the global bin_num variable (and others) via the bufptrs pointer array. */

		if ( VSfpack(vdata_id,_HDF_VSUNPACK,BLIST_FIELDS,blistrec,BLIST_SIZE,1,NULL,bufptrs) == FAIL) {
			mexPrintf("-E- %s line %d: VSfpack(%d, ...) failed.\n",__FILE__,__LINE__, vdata_id);
			mexErrMsgTxt("\n");
		}
		if     (bin < bin_num) hi = mid - 1;
		else if(bin > bin_num) lo = mid + 1;
		else                   return(mid);
	}
	return(-1);
}

/* The following functions are based on the pseudocode found in Appendix A of:

	Campbell, J.W., J.M. Blaisdell, and M. Darzi, 1995:
	Level-3 SeaWiFS Data Products: Spatial and Temporal Binning Algorithms.
	NASA Tech. Memo. 104566, Vol. 32,
	S.B. Hooker, E.R. Firestone, and J.G. Acker, Eds.,
	NASA Goddard Space Flight Center, Greenbelt, Maryland
*/

int initbin(void) {
	int row;

	basebin[0] = 1;
	for (row = 0; row < NUMROWS; row++) {
		latbin[row] = ((row + 0.5)*180.0/NUMROWS) - 90.0;
		numbin[row] = (short int)(2*NUMROWS*cos(latbin[row]*D2R) + 0.5);
		if (row > 0)
			basebin[row] = basebin[row - 1] + numbin[row - 1];
	}
	return(basebin[NUMROWS - 1] + numbin[NUMROWS - 1] - 1);
}

short int lat2row(double lat) {
	short int row;

	row = (short int)((90 + lat)*NUMROWS/180.0);
	if (row >= NUMROWS) row = NUMROWS - 1;
	return(row);
}

int rowlon2bin(short int row, double lon){
	short int col;
	int	bin;

	lon = constrain_lon(lon);
	col = (short int)((lon + 180.0)*numbin[row]/360.0);
	if(col >= numbin[row]) col = numbin[row] - 1;
	bin = basebin[row] + col;
	return(bin);
}

int latlon2bin(double lat, double lon){
	short int row, col;
	int	bin;

	/* Constrain latitudes to [-90,90] and longitudes to [-180,180]. */
	lat = constrain_lat(lat);
	lon = constrain_lon(lon);

	row = lat2row(lat);
	col = (short int)((lon + 180.0)*numbin[row]/360.0);
	if (col >= numbin[row]) col = numbin[row] - 1;
	bin = basebin[row] + col;
	return(bin);
}

void bin2latlon(int bin, double *clat, double *clon){
	short int row;

	row = NUMROWS - 1;
	if (bin < 1) bin = 1;
	while (bin < basebin[row]) row--;
	*clat = latbin[row];
	*clon = 360.0*(bin - basebin[row] + 0.5)/numbin[row] - 180.0;
}

double constrain_lat(double lat){
	if (lat >  90) lat =  90;
	if (lat < -90) lat = -90;
	return(lat);
}

double constrain_lon(double lon){
	while (lon < -180) lon += 360;
	while (lon >  180) lon -= 360;
	return(lon);
}
