/* Copyright 1993-2007 The MathWorks, Inc. */

/* $Revision: 5.19.4.3 $  $Date: 2007/06/04 21:09:41 $ */

/*
 *	imhistc.c
 *
 *	Y = IMHISTC( A, N, ISSCALED, TOP ) makes an N bin histogram of 
 *      the image matrix A.  This function is written as an auxiliary to 
 *      the M-file IMHIST.M.
 *
 *	This is a MEX file for MATLAB.
 *
 */

#include <math.h>
#include "mex.h"
#include "mwsize.h"
#include "iptutil.h"

void ScaledHistDouble(double *a, double top, mwSize n, mwSize length, double *y);
void HistDouble(double *a, mwSize n, mwSize length, double *y);
void ScaledHistSingle(float *a, double top, mwSize n, mwSize length, double *y);
void HistSingle(float *a, mwSize n, mwSize length, double *y);
void HistUint8(uint8_T *a, mwSize n, mwSize length, double *y);
void ScaledHistUint8(uint8_T *a, double top, mwSize n, mwSize length, double *y);
void HistUint16(uint16_T *a, mwSize n, mwSize length, double *y);
void ScaledHistUint16(uint16_T *a, double top, mwSize n, mwSize length, double *y);
void ValidateInputs( const mxArray *prhs[], int nrhs, const mxArray **A, mwSize *n, bool *isScaled, double *top);

/******************************************************************/
/*   The ScaledHist___ functions are used for an intensity image. */
/*   The Hist___ functions are used if a colormap is provided.    */
/******************************************************************/

void ScaledHistDouble(double *a, double top, mwSize n, mwSize length, double *y) {
    mwSize i;
    double scale;
    double z;
    int warnNanFlag = 0;

    scale = (double) (n-1) / top;


    for (i = 0; i < length; i++) {
	if(mxIsNaN(a[i])) {
	    warnNanFlag = 1;
	    z = 0;
	}
	else
	    z = floor(a[i] * scale + .5);

        if (z < 0.0) 
            y[0]++;
        else if (z > (n-1)) 
            y[n-1]++;
        else 
            y[(mwSize) z]++;
    }

    if (warnNanFlag == 1)
        mexWarnMsgIdAndTxt("Images:imhistc:inputHasNaNs",
			   "Converting NaN inputs to 0.");

}

/***************************************************************/
void ScaledHistSingle(float *a, double top, mwSize n, mwSize length, double *y) {
    mwSize i;
    float scale;
    float z;
    int warnNanFlag = 0;

    scale = (float) (n-1) / (float) top;

    for (i = 0; i < length; i++) {

	if(mxIsNaN(a[i])) {
	    warnNanFlag = 1;
	    z = 0;
	}
	else
	    z = (float) floor(a[i] * scale + .5);

        if (z < 0.0)
            y[0]++;
 
        else if (z > (n-1)) 
            y[n-1]++;
     
        else 
            y[(mwSize) z]++;
    }

    if (warnNanFlag == 1) {
        mexWarnMsgIdAndTxt("Images:imhistc:inputHasNaNs",
			   "Converting NaN inputs to 0.");
    }
}

/***************************************************************/
void HistDouble(double *a, mwSize n, mwSize length, double *y) {
    mwSize i;
    mwSize idx;
    int warnRngFlag = 0;
    int warnNanFlag = 0;

    for (i = 0; i < length; i++) {

	if(mxIsNaN(a[i])) {
	    warnNanFlag = 1;
	    idx = 0; /* Treat NaNs like zeros */
	}
	else {
            double a_minus_1 = a[i] - 1.0;
            if (a_minus_1 < 0) {
                warnRngFlag = 1;
                idx = 0;
            }
            else if (a_minus_1 >= n) {
                warnRngFlag = 1;
                idx = n - 1;
            }
            else
                idx = (mwSize) a_minus_1;
	}

        y[idx]++;
    }

    if (warnRngFlag == 1) {
        mexWarnMsgIdAndTxt("Images:imhistc:outOfRange",
			   "Input has out-of-range values");
    }
    else if (warnNanFlag == 1) {
        mexWarnMsgIdAndTxt("Images:imhistc:inputHasNaNs",
			   "Converting NaN inputs to 0.");
    }
    
}

/***************************************************************/
void HistSingle(float *a, mwSize n, mwSize length, double *y) {
    mwSize i;
    mwSize idx;
    int warnRngFlag = 0;
    int warnNanFlag = 0;

    for (i = 0; i < length; i++) {

	if(mxIsNaN(a[i])) {
	    warnNanFlag = 1;
	    idx = 0; /* Treat NaNs like zeros */
	}
	else {
            double a_minus_1 = a[i] - 1.0;
            if (a_minus_1 < 0.0) {
                warnRngFlag = 1;
                idx = 0;
            }
            else if (a_minus_1 >= n) {
                warnRngFlag = 1;
                idx = n - 1;
            }
            else
                idx = (mwSize) a_minus_1;
        }

        y[idx]++;
    }

    if (warnRngFlag == 1) {
        mexWarnMsgIdAndTxt("Images:imhistc:outOfRange",
			   "Input has out-of-range values");
    }
    else if (warnNanFlag == 1) {
        mexWarnMsgIdAndTxt("Images:imhistc:inputHasNaNs",
			   "Converting NaN inputs to 0.");
    }
    
}

/***************************************************************/
void HistUint8(uint8_T *a, mwSize n, mwSize length, double *y) {
    mwSize i;
    mwSize idx;
    int warnFlag = 0;

    if(n==256) {
        for (i = 0; i < length; i++) {
            y[a[i]]++;
        }
    }
    else {
        for (i = 0; i < length; i++) {
            idx = a[i];
            if (idx > (n-1)) {
                warnFlag = 1;
                y[n-1]++;
            } else
                y[idx]++;
        }
    }
    
    if (warnFlag == 1) {
        mexWarnMsgIdAndTxt("Images:imhistc:outOfRange",
			   "Input has out-of-range values");
    }
    
}

/***************************************************************/
void ScaledHistUint8(uint8_T *a, double top, mwSize n, mwSize length, double *y) {
    mwSize i;
    double scale;
    mwSize z;

    scale = (double) (n-1) / top;

    for (i = 0; i < length; i++) {
        z = (mwSize)(a[i] * scale + 0.5);
        
	if (z > (n-1))
            y[n-1]++;
        else
            y[z]++;
    }

}

/***************************************************************/
void HistUint16(uint16_T *a, mwSize n, mwSize length, double *y) {
    mwSize i;
    mwSize idx;
    int warnFlag = 0;

    if(n==65536) {
        for (i = 0; i < length; i++) {
            y[a[i]]++;
        }
    }
    else {
	for (i = 0; i < length; i++) {
	    idx = a[i];
	    if (idx > (n-1)) {
		warnFlag = 1;
		y[n-1]++;
	    } else {
		y[idx]++;
	    }
	}
    }

    if (warnFlag == 1) {
        mexWarnMsgIdAndTxt("Images:imhistc:outOfRange", 
			   "Input has out-of-range values");
    }
    
}

/***************************************************************/
void ScaledHistUint16(uint16_T *a, double top, mwSize n, mwSize length, double *y) {
    mwSize i;
    double scale;
    mwSize z;

    scale = (double) (n-1) / top;

    for (i = 0; i < length; i++) {
        z = (mwSize)(a[i] * scale + 0.5);

        if (z > (n-1))
            y[n-1]++;
        else
            y[z]++;
    }
}

/***************************************************************/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    mwSize length;
    const mxArray *A;
    mwSize n;
    bool isScaled;
    double top;
    double *a_real;
    uint8_T *a_int8;
    uint16_T *a_int16;
    double *y;
    mxArray *Y;
    float *a_single;

    (void) nlhs; /* unused parameter */

    ValidateInputs(prhs, nrhs, &A, &n, &isScaled, &top);
    length = mxGetM(A) * mxGetN(A);

    Y = mxCreateDoubleMatrix(n, 1, mxREAL);
    y = (double *) mxGetData(Y);

    if (mxIsDouble(A)) {
        a_real = (double *) mxGetData(A);
        if (isScaled)
            ScaledHistDouble(a_real, top, n, length, y);
        else 
            HistDouble(a_real, n, length, y);

    } 
    else if (mxIsSingle(A)) {
        a_single = (float *) mxGetData(A);
        if (isScaled) 
            ScaledHistSingle(a_single, top, n, length, y);
        else 
            HistSingle(a_single, n, length, y);

    } 
    else if (mxIsUint8(A)) {
        a_int8 = (uint8_T *) mxGetData(A);
        if (isScaled) {
            if ((n == 256) && (top == 255.0)) 
                HistUint8(a_int8, n, length, y);
            else 
                ScaledHistUint8(a_int8, top, n, length, y);
        } 
        else 
            HistUint8(a_int8, n, length, y);
    } 
    else if (mxIsUint16(A)) {
        a_int16 = (uint16_T *) mxGetData(A);
        if (isScaled) {
            if ((n == 65536) && (top == 65535.0)) {
                HistUint16(a_int16, n, length, y);
            } else {
                ScaledHistUint16(a_int16, top, n, length, y);
            }
        } else {
            HistUint16(a_int16, n, length, y);
        }
    }
        
    /* Done! Give the answer back */
    plhs[0] = Y;
}


void ValidateInputs( const mxArray *prhs[], int nrhs, const mxArray **A, mwSize *n, bool *isScaled, double *top) {
    int i;
    mwSize length;
    double n_real;
    double isScaled_real;

    if(nrhs != 4){
        mexErrMsgIdAndTxt("Images:imhistc:invalidNumOfInputs",
			  "Four inputs are required");
    }

    for (i = 0; i < nrhs; i++) {
        if (mxIsComplex(prhs[i])) {
            mexWarnMsgIdAndTxt("Images:imhistc:ignoreImaginaryPart",
			       "Ignoring imaginary part of complex inputs");
        }
        if (!mxIsNumeric(prhs[i])) {
            mexErrMsgIdAndTxt("Images:imhistc:mustBeNumeric",
			      "Inputs to IMHISTC must be numeric");
        }
    }

    *A = prhs[0];

    for (i = 1; i < nrhs; i++) {
        if (!mxIsDouble(prhs[i])) {
            mexErrMsgIdAndTxt("Images:imhistc:mustBeDouble",
			      "IMHISTC inputs 2, 3, and 4 must be DOUBLE");
        }
    }

    length = mxGetM(prhs[1]) * mxGetN(prhs[1]);
    if (length == 0) {
        n_real = 0;
        mexErrMsgIdAndTxt("Images:imhistc:invalidSecondInput",
			  "Second input to IMHISTC must not be empty");
    } else if (length == 1) {
        n_real = *((double *) mxGetData(prhs[1]));
    } else {
        mexWarnMsgIdAndTxt("Images:imhistc:invalidSecondInput",
			   "Second input to IMHISTC should be a scalar");
        n_real = *((double *) mxGetData(prhs[1]));
    }
    *n = (mwSize) n_real;
    if (((double) *n) != n_real) {
        mexWarnMsgIdAndTxt("Images:imhistc:invalidSecondInput",
			   "Second input to IMHISTC should be an integer");
    }
    if (*n <= 0) {
        mexErrMsgIdAndTxt("Images:imhistc:invalidSecondInput",
			  "Second input to IMHISTC should be positive");
    }

    length = mxGetM(prhs[2]) * mxGetN(prhs[2]);
    if (length == 0) {
        isScaled_real = 0;
        mexErrMsgIdAndTxt("Images:imhistc:invalidThirdInput",
			  "Third input to IMHISTC must not be empty");
    } else if (length == 1) {
        isScaled_real = *((double *) mxGetData(prhs[2]));
    } else {
        mexWarnMsgIdAndTxt("Images:imhistc:invalidThirdInput",
			   "Third input to IMHISTC should be a scalar");
        isScaled_real = *((double *) mxGetData(prhs[2]));
    }
    if (isScaled_real != 0.0) {
        *isScaled = true;
    } else {
        *isScaled = false;
    }

    length = mxGetM(prhs[3]) * mxGetN(prhs[3]);
    if (length == 0) {
        mexErrMsgIdAndTxt("Images:imhistc:invalidFourthInput",
			  "Fourth input to IMHISTC must not be empty");
    } else if (length == 1)
        *top = *((double *) mxGetData(prhs[3]));
    else {
        mexWarnMsgIdAndTxt("Images:imhistc:invalidFourthInput",
			   "Fourth input to IMHISTC should be a scalar");
        *top = *((double *) mxGetData(prhs[3]));
    }
    if (*top <= 0.0) {
        mexErrMsgIdAndTxt("Images:imhistc:invalidFourthInput",
			  "Fourth input to IMHISTC must be positive");
    }
}    
