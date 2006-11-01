/*
 * Output a uint8 or uint16 matrix scaled to the [0-255] or [0-65535] range. The input matrix can have
 * any of the following types: double, single, Int32, Int16 and UInt16 
 *
 * Alternatively give two outputs to recover only z_min & z_max
 *
 * Note: In fact the scaling is done in [0-254] (or [0-65534]) and add 1 to the result. In this way the
 * the first color entry is reserved for the background color (NaNs). Apparently IVS uses the same
 * technique when scaling to uint16.
 *
 * Revision 1.0	18-Sep-2004 Joaquim Luis
 * Changes 	26-Oct-2004 - Added a test for uint8 with NaNs
 *		17-Apr-2006 - Added a noDataValue for ints
 *
 */

#ifndef MIN
#define MIN(x, y) (((x) < (y)) ? (x) : (y))	/* min and max value macros */
#endif
#ifndef MAX
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#endif

#include "mex.h"
#include "float.h"

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	double  range8, *z_min, *z_max, *z_8, min8 = DBL_MAX, max8 = -DBL_MAX;
	double	*pNodata, *which_scale;
	float	range, min = FLT_MAX, max = -FLT_MAX, *z_4, new_range, nodata;
	int     nx, ny, i, is_double = 0, is_single = 0, is_int32 = 0, is_int16 = 0;
	int     is_uint16 = 0, scale_range = 1, *i_4, scale8, scale16;
	int	got_nodata = 0;
	short int *i_2;
	unsigned short int *ui_2, *out16;
	unsigned char *out8;
  
  	/*  check for proper number of arguments */
	if((nrhs < 1 || nrhs > 3) || (nlhs < 1 || nlhs > 2)) {
		mexPrintf ("usage: img8  = scaleto8(Z);\n");
		mexPrintf (" 	   img8  = scaleto8(Z,8,noDataValue);\n");
		mexPrintf (" 	   img16 = scaleto8(Z,16);\n");
		mexPrintf (" 	   img16 = scaleto8(Z,16,noDataValue);\n");
		mexPrintf (" 	   [z_min,z_max] = scaleto8(Z);\n");
		return;
	}
  
	if(mxGetN(prhs[0])*mxGetM(prhs[0]) == 1)
		mexErrMsgTxt("SCALETO8 ERROR: Input must be a matrix.");
 
	if(nrhs == 1) {
		new_range = 254;		/* scale to uint8 */
		scale8 = 1;	scale16 = 0;
	}
	else {
		nx = mxGetN (prhs[1]);
		if (nx > 1)
			mexErrMsgTxt("SCALETO8 ERROR: Second argument must be a scalar == 8 or 16.");
		which_scale = mxGetPr(prhs[1]);
		if (*which_scale == 8) {
			new_range = 254;		/* scale to uint8 */
			scale8 = 1;	scale16 = 0;
		}
		else if (*which_scale == 16) { 
			new_range = 65534;		/* scale to uint16 */
			scale16 = 1;	scale8 = 0;
		}
		else
			mexErrMsgTxt("SCALETO8 ERROR: Second argument must be a scalar == 8 or 16.");
		if (nrhs == 3) {
			nx = mxGetN (prhs[2]);
			if (nx > 1)
				mexErrMsgTxt("SCALETO8 ERROR: Thirth argument must be a scalar.");
			pNodata = mxGetPr(prhs[2]);
			nodata = (float)pNodata[0];
			got_nodata = 1;
		}
	}

	if(nlhs == 2) 
		scale_range = 0;	/* Output min/max */

	/* Find out in which data type was given the input array */
	if (mxIsDouble(prhs[0])) {
		z_8 = mxGetPr(prhs[0]);
		is_double = 1;
	}
	else if (mxIsSingle(prhs[0])) {
		z_4 = mxGetData(prhs[0]);
		is_single = 1;
	}
	else if (mxIsInt32(prhs[0])) {
		i_4 = mxGetData(prhs[0]);
		is_int32 = 1;
	}
	else if (mxIsInt16(prhs[0])) {
		i_2 = mxGetData(prhs[0]);
		is_int16 = 1;
	}
	else if (mxIsUint16(prhs[0])) {
		ui_2 = mxGetData(prhs[0]);
		is_uint16 = 1;
	}
	else {
		mexPrintf("SCALETO8 ERROR: Unknown input data type.\n");
		mexErrMsgTxt("Valid types are:double, single, Int32, Int16 and UInt16.\n");
	}

	/*  get the dimensions of the matrix input */
	ny = mxGetM(prhs[0]);
	nx = mxGetN(prhs[0]);

	/*  set the output pointer to the output matrix */
	if (scale_range) {	/* If scale the output into new_range */
		/*  create a C pointer to a copy of the output matrix */
		if (scale8) {
			plhs[0] = mxCreateNumericMatrix (ny,nx,mxUINT8_CLASS ,mxREAL);
			out8 = mxGetData(plhs[0]);
		}
		else {
			plhs[0] = mxCreateNumericMatrix (ny,nx,mxUINT16_CLASS ,mxREAL);
			out16 = mxGetData(plhs[0]);
		}
	}
	else {
		plhs[0] = mxCreateDoubleMatrix (1,1,mxREAL);
		plhs[1] = mxCreateDoubleMatrix (1,1,mxREAL);
		/*  create a C pointer to a copy of the output matrix */
		z_min = mxGetPr(plhs[0]);
		z_max = mxGetPr(plhs[1]);
	}
  

	/* Find min and max value */
	if (is_double) {
		for (i = 0; i < nx*ny; i++) {
			if (mxIsNaN(z_8[i])) continue;
			min8 = MIN(min8,z_8[i]);
			max8 = MAX(max8,z_8[i]);
		}
		if (scale_range) {	/* Scale data into the new_range ([0 255] or [0 65535]) */ 
			range8 = max8 - min8;	/* This range8 means input has 8 bytes */
			if (scale8) {	/* Scale to uint8 */
				for (i = 0; i < nx*ny; i++) {
					if (mxIsNaN(z_8[i])) continue;
					out8[i] = (char)(((z_8[i] - min8) / range8) * new_range + 1);
				} 
			} 
			else if (scale16) {	/* Scale to uint16 */
				for (i = 0; i < nx*ny; i++) {
					if (mxIsNaN(z_8[i])) continue;
					out16[i] = (unsigned short int)(((z_8[i] - min8) / range8) * new_range + 1);
				} 
			} 
		}
		else {			/* min/max required */
			*z_min = min8;	*z_max = max8;
		}
	}
	else if (is_single) {
		for (i = 0; i < nx*ny; i++) {
			if (mxIsNaN((double)z_4[i])) continue;
			min = MIN(min,z_4[i]);
			max = MAX(max,z_4[i]);
		}
		if (scale_range) {	/* Scale data into the new_range ([0 255] or [0 65535]) */ 
			range = max - min;
			if (scale8) {	/* Scale to uint8 */
				for (i = 0; i < nx*ny; i++) {
					if ( mxIsNaN( (double)z_4[i] ) ) continue;
					out8[i] = (char)(((z_4[i] - min) / range) * new_range + 1);
				} 
			} 
			else if (scale16) {	/* Scale to uint16 */
				for (i = 0; i < nx*ny; i++) {
					if ( mxIsNaN( (double)z_4[i] ) ) continue;
					out16[i] = (unsigned short int)(((z_4[i] - min) / range) * new_range + 1);
				} 
			} 
		}
		else {			/* min/max required */
			*z_min = min;	*z_max = max;
		}
	}
	else if (is_int32) {
		if (!got_nodata) {
			for (i = 0; i < nx*ny; i++) {
				min = MIN(min,i_4[i]);
				max = MAX(max,i_4[i]);
			}
		}
		else {
			for (i = 0; i < nx*ny; i++) {
				if ((double)i_4[i] == nodata) continue;
				min = MIN(min,i_4[i]);
				max = MAX(max,i_4[i]);
			}
		}
		if (scale_range) {	/* Scale data into the new_range ([0 255] or [0 65535]) */ 
			range = max - min;
			if (scale8) {	/* Scale to uint8 */
				for (i = 0; i < nx*ny; i++)
					out8[i] = (char)((((float)i_4[i] - min) / range) * new_range + 1);
			}
			else if (scale16) {	/* Scale to uint16 */
				for (i = 0; i < nx*ny; i++)
					out16[i] = (unsigned short int)((((float)i_4[i] - min) / range) * new_range + 1);
			}
		}
		else {			/* min/max required */
			*z_min = min;	*z_max = max;
		}
	}
	else if (is_int16) {
		if (!got_nodata) {
			for (i = 0; i < nx*ny; i++) {
				min = MIN(min,i_2[i]);
				max = MAX(max,i_2[i]);
			}
		}
		else {
			for (i = 0; i < nx*ny; i++) {
				if ((double)i_2[i] == nodata) continue;
				min = MIN(min,i_2[i]);
				max = MAX(max,i_2[i]);
			}
		}
		if (scale_range) {	/* Scale data into the new_range ([0 255] or [0 65535]) */ 
			range = max - min;
			if (scale8) {	/* Scale to uint8 */
				for (i = 0; i < nx*ny; i++)
					out8[i] = (char)((((float)i_2[i] - min) / range) * new_range + 1);
			}
			else if (scale16) {	/* Scale to uint16 */
				for (i = 0; i < nx*ny; i++)
					out16[i] = (unsigned short int)((((float)i_2[i] - min) / range) * new_range + 1);
			}
		}
		else {			/* min/max required */
			*z_min = min;	*z_max = max;
		}
	}
	else if (is_uint16) {
		if (!got_nodata) {
			for (i = 0; i < nx*ny; i++) {
				min = MIN(min,ui_2[i]);
				max = MAX(max,ui_2[i]);
			}
		}
		else {
			for (i = 0; i < nx*ny; i++) {
				if ((double)ui_2[i] == nodata) continue;
				min = MIN(min,ui_2[i]);
				max = MAX(max,ui_2[i]);
			}
		}
		if (scale_range) {	/* Scale data into the new_range ([0 255] or [0 65535]) */ 
			range = max - min;
			if (scale8) {	/* Scale to uint8 */
				for (i = 0; i < nx*ny; i++)
					out8[i] = (char)((((float)ui_2[i] - min) / range) * new_range + 1);
			}
			else if (scale16) {	/* Scale to uint16 */
				for (i = 0; i < nx*ny; i++)
					out16[i] = (unsigned short int)((((float)ui_2[i] - min) / range) * new_range + 1);
			}
		}
		else {			/* min/max required */
			*z_min = min;	*z_max = max;
		}
	}
}
