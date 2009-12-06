/*
 * GRAYTO8 MEX-file
 *
 * B = GRAYTO8(A) converts the double array A to uint8 by scaling A by 255
 * and then rounding.  NaN's in A are converted to 0.  Values in A greater 
 * than 1.0 are converted to 255; values less than 0.0 are converted to 0.
 *
 * B = GRAYTO8(A) converts the uint16 array A by scaling the elements of A
 * 255/65535, rounding, and then casting to uint8.
 *
 * Copyright 1993-2007 The MathWorks, Inc.
 *
 */

#include "mex.h"
#include "math.h"
#include "mwsize.h"

void ConvertFromDouble(double *pr, uint8_T *qr, mwSize numElements) {
    mwSize k;
    double val;

    for (k = 0; k < numElements; k++) {
        val = *pr++;
        if (mxIsNaN(val)) {
            *qr++ = 0;
        }
        else {
            val = val * 255.0 + 0.5;
            if (val > 255.0) val = 255.0;
            if (val < 0.0)   val = 0.0;
            *qr++ = (uint8_T) val;
        }
    }
}

void ConvertFromSingle(float *pr, uint8_T *qr, mwSize numElements) {
    mwSize k;
    float val;

    for (k = 0; k < numElements; k++) {
        val = *pr++;
        if (mxIsNaN(val))
            *qr++ = 0;
        else {
            val = val * 255.0f + 0.5f;
            if (val > 255.0) val = 255.0;
            if (val < 0.0)   val = 0.0;
            *qr++ = (uint8_T) val;
        }
    }
}

void ConvertFromUint16(uint16_T *pr, uint8_T *qr, mwSize numElements) {
    mwSize k;
    double factor = 1.0 / 257.0;

    for (k = 0; k < numElements; k++)
        *qr++ = (uint8_T) ( (double) (*pr++) * factor + 0.5 );
}

void ValidateInputs(int nrhs, const mxArray *prhs[]) {
    if (nrhs < 1) {
        mexErrMsgIdAndTxt("Images:grayto8:tooFewInputs",
                          "%s","Too few input arguments.");
    }
    if (nrhs > 1) {
        mexErrMsgIdAndTxt("Images:grayto8:tooManyInputs",
                          "%s","Too many input arguments.");
    }
    if (!mxIsDouble(prhs[0]) && !mxIsUint16(prhs[0]) && !mxIsSingle(prhs[0])) {
        mexErrMsgIdAndTxt("Images:grayto8:invalidType",
                          "%s","Input must be double, single, or uint16.");
    }
    if (mxIsComplex(prhs[0])) {
        mexWarnMsgIdAndTxt("Images:grayto8:ignoringImaginaryPartOfInput",
                           "%s","Ignoring imaginary part of input.");
    }
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    uint8_T *qr;

    (void) nlhs;  /* unused parameter */

    ValidateInputs(nrhs, prhs);

    plhs[0] = mxCreateNumericArray(mxGetNumberOfDimensions(prhs[0]),
                                   mxGetDimensions(prhs[0]), 
                                   mxUINT8_CLASS, mxREAL);
    qr = (uint8_T *) mxGetData(plhs[0]);

    if (mxIsDouble(prhs[0]))
        ConvertFromDouble((double *) mxGetData(prhs[0]), qr, mxGetNumberOfElements(prhs[0]));

    else if (mxIsUint16(prhs[0]))
        ConvertFromUint16((uint16_T *) mxGetData(prhs[0]), qr, mxGetNumberOfElements(prhs[0]));

    else
        ConvertFromSingle((float *) mxGetData(prhs[0]), qr, mxGetNumberOfElements(prhs[0]));

}
