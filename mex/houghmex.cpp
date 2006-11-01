//////////////////////////////////////////////////////////////////////////////
// HOUGHMEX.<mex>
//
// H = houghmex(BW,RHO,THETA)
//
// This mex function implements Standard Hough Transform (SHT). It requires
// a binary input image.  It also requires two vectors of RHO and THETA
// values over which the SHT is evaluated. The resulting Hough transform is 
// returned in a double Matlab array.
//
// BW
//     logical image
// RHO
//     rho vector
// THETA
//     theta vector in radians
// H
//     output hough matrix
//
// $Revision $  $Date: 2004/08/10 01:45:11 $
// Copyright 1993-2003 The MathWorks, Inc.
//
// $Revision $  $Date: 2006/08/29 J. LUIS
//		Changed output to SINGLE (& it probably can be further changed to uint16)
//////////////////////////////////////////////////////////////////////////////

#include <math.h>
#include "mex.h"
#include "typeconv.h"

static char rcsid[] = "$Id: houghmex.cpp,v 1.1.8.1 2004/08/10 01:45:11 batserve Exp $";

//////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////
#define FCN_NAME "houghmex"
void checkInputs(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs != 3) {
        mexErrMsgIdAndTxt("Images:houghmex:invalidNumInputs",
                          "HOUGHMEX requires three input arguments.");
    }

    // Perform some rudimentary checks
    const mxArray *BW       = prhs[0]; 
    const mxArray *RHO      = prhs[1];
    const mxArray *THETA    = prhs[2];

    // Check class of the inputs
    if ( mxGetClassID(BW) != mxLOGICAL_CLASS ) {
        mexErrMsgIdAndTxt("Images:houghmex:invalidType",
                          "%s %s %s\n%s",
                          "Function", FCN_NAME, "expected its"
                          " first input argument, BW,",
                          "to be logical.");
    }

    if ( mxGetClassID(RHO) != mxDOUBLE_CLASS ) {
        mexErrMsgIdAndTxt("Images:houghmex:invalidType",
                          "%s %s %s\n%s",
                          "Function", FCN_NAME, "expected its"
                          " second input argument, RHO,",
                          "to be double.");
    }

    if ( mxGetClassID(THETA) != mxDOUBLE_CLASS ) {
        mexErrMsgIdAndTxt("Images:houghmex:invalidType",
                          "%s %s %s\n%s",
                          "Function", FCN_NAME, "expected its"
                          " third input argument, THETA,",
                          "to be double.");
    }
}

//////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////
void hough(uint8_T *inImage, int inRows, int inCols, double *rho,
               int rhoLen, double *theta, int thetaLen, float *outImage) {
    int i, n_inRows;
    int rhoIdx;
    double firstRho = rho[0];

    // Allocate space for the cos/sin lookup tables
    double *cost = (double *)mxCalloc(thetaLen,sizeof(double));
    double *sint = (double *)mxCalloc(thetaLen,sizeof(double));

    // Precompute the sin and cos tables
    for(i=0; i<thetaLen; i++)
    {
        cost[i] = cos(theta[i]);
        sint[i] = sin(theta[i]);
    }

    // Compute the factor for converting back to the rho matrix index
    double slope = (rhoLen - 1)/(rho[rhoLen-1] - firstRho);

    // Compute the hough transform
    for(int n=0; n < inCols; n++) {
	n_inRows = n * inRows;
        for(int m=0; m < inRows; m++) {
            if(inImage[n_inRows+m]) { // if pixel is on
                for(int thetaIdx=0; thetaIdx<thetaLen; thetaIdx++) {
                    // x*cos(theta)+y*sin(theta)=rho
                    double myrho = n*cost[thetaIdx]+m*sint[thetaIdx];
                    // convert to bin index
                    roundAndCast(&rhoIdx, slope*(myrho - firstRho));
                    outImage[thetaIdx*rhoLen+rhoIdx]++; //accumulate
                }
            }
        }
    }   
    
    // Clean up the memory!
    if(cost) mxFree(cost);
    if(sint) mxFree(sint);
}

//////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////
extern "C"
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    checkInputs(nlhs, plhs, nrhs, prhs);

    // Input arrays
    const mxArray *BW       = prhs[0]; 
    uint8_T       *inImage  = (uint8_T *)mxGetData(BW);
    const mxArray *RHO      = prhs[1];
    double        *rho      = (double *)mxGetData(RHO);
    const mxArray *THETA    = prhs[2];
    double        *theta    = (double *)mxGetData(THETA);
    
    // Input parameters
    //const int *dims = mxGetDimensions(BW);
    int inRows      = mxGetM(BW);
    int inCols      = mxGetN(BW);
    int rhoLen      = mxGetNumberOfElements(RHO);
    int thetaLen    = mxGetNumberOfElements(THETA);
    int dims[3];

    // Create output array
    //mxArray       *OUT      = mxCreateDoubleMatrix(rhoLen,thetaLen,mxREAL);
    dims[0] = rhoLen;
    dims[1] = thetaLen;
    dims[2] = 1;
    mxArray       *OUT      = mxCreateNumericArray(2,dims,mxSINGLE_CLASS,mxREAL);
    float         *outImage = (float *)mxGetData(OUT);
    
    // Compute and return the results
    hough(inImage, inRows, inCols, rho, rhoLen, theta, thetaLen, outImage);
    plhs[0] = OUT;
}
