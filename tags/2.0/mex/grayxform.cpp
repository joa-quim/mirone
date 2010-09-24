//Copyright 1993-2007 The MathWorks, Inc.


//    grayxform.cpp .MEX file

//                out = grayxform(in, T)

//                Apply a graylevel transform to an image.
               
//                in should be of class uint8, uint16, or double.
//                T should be double, in the range [0,1].
               
//                Chris Griffin
//                February 1998



#include "mex.h"
#include "mwsize.h"
#include "grayxform.h"

//////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////
void checkInputs(int nlhs, int nrhs, const mxArray *prhs[]) {
    if(nrhs != 2) 
        mexErrMsgIdAndTxt("Images:grayxform:invalidInputArgs",
                          "%s","Invalid input arguments.");

    if(nlhs > 1) 
        mexErrMsgIdAndTxt("Images:grayxform:tooManyOutputArgs",
                          "%s","Too many output arguments.");

    // check image
    if(mxIsComplex(prhs[0]))
        mexErrMsgIdAndTxt("Images:grayxform:inputImageMustBeReal",
                          "%s","Input image is complex, only real"
                          " images are supported.");

    // check transform array.
    if (!mxIsDouble(prhs[1])) 
        mexErrMsgIdAndTxt("Images:grayxform:trasformArrayMustBeDouble",
                          "%s","Transformation array should be of"
                          " class double.");

    if(mxIsEmpty(prhs[1]))
        mexErrMsgIdAndTxt("Images:grayxform:emptyTransformArray",
                          "%s","Empty transformation array.");

    if(mxIsComplex(prhs[1])) 
        mexWarnMsgIdAndTxt("Images:grayxform:"
                           "imaginaryPartOfTransformArrayIgnored",
                           "%s","Complex transformation array, imaginary part"
                           " ignored.");

    double *t = mxGetPr(prhs[1]);

    for(mwSize i=0; i<mxGetNumberOfElements(prhs[1]); i++) {
        if(t[i]<0 || t[i]>1)
            mexErrMsgIdAndTxt("Images:grayxform:"
                              "outOfRangeElementsOfTransformArray",
                              "%s","Elements of transformation array"
                              " outside the range [0,1].");
    }
}

//////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////
extern "C"
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    checkInputs(nlhs, nrhs, prhs);

    // get information about input image.
    const mxArray *in       = prhs[0];
    void          *inData   = mxGetData(in);
    mwSize        nelements = mxGetNumberOfElements(in);
    mwSize        ndims     = mxGetNumberOfDimensions(in);
    const mwSize  *dims     = mxGetDimensions(in);
    mxClassID     classID   = mxGetClassID(in);

    // get information about transformation array.
    const mxArray *T        = prhs[1];
    double        *tData    = (double *)mxGetData(T);
    mwSize        maxTidx   = mxGetNumberOfElements(T)-1;

    // initialize output array
    mxArray *out = mxCreateNumericArray(ndims, dims, classID, mxREAL);
    void *outData = mxGetData(out);

    // perform transformation
    switch (classID) {
    case mxUINT8_CLASS:
        transformUintArray((uint8_T *)inData, tData, maxTidx, nelements, (uint8_T *)outData, 255.0);
        break;
    case mxUINT16_CLASS:
        transformUintArray((uint16_T *)inData, tData, maxTidx, nelements, (uint16_T *)outData, 65535.0);
        break;
    case mxDOUBLE_CLASS:
        transformDoubleSingleArray((double *)inData, tData, maxTidx, nelements, (double *)outData);
        break;
    case mxSINGLE_CLASS:
        transformDoubleSingleArray((float *)inData, tData, maxTidx, nelements, (float *)outData);
        break;

    default:
        mexErrMsgIdAndTxt("Images:grayxform:unsupportedInputClass",
                          "%s","Unsupported input data class.");
        break;
    }
    plhs[0] = out;
}
