// Copyright 1993-2007 The MathWorks, Inc.
// $Revision.1 $  $Date: 2007/06/04 21:09:46 $

// IMRECONSTRUCT MEX-file
//
// K = IMRECONSTRUCT(J,I,CONN) performs grayscale reconstruction with
// J as the marker image and I as the mask image.  CONN specifies
// connectivity.
//
// Input and output specs
// ----------------------
// J:    N-D real matrix, uint8, uint16, or double
//       empty allowed
//       Infs allowed
//       NaNs not allowed
//
// I:    N-D real matrix, same size and class as J
//       Infs allowed
//       NaNs not allowed
//       Elements of I must be >= corresponding elements of J
// 
// CONN: See connectivity spec.
//
// K:    N-D real matrix, same size and class as I and J.
//       logical if and only if both inputs are logical.

#include "neighborhood.h"
#include "reconstruct.h"
#include "mex.h"
#include "mwsize.h"

//////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////
void CheckInputs(int nrhs, const mxArray *prhs[]) {
    mwSize ndims;
    const mwSize *size_0;
    const mwSize *size_1;
    mwSize k;

    if (nrhs < 2)
        mexErrMsgIdAndTxt("Images:imreconstructmex:tooFewInputs",
                          "IMRECONSTRUCTMEX needs at least two input "
                          "arguments.");

    else if (nrhs > 3)
        mexErrMsgIdAndTxt("Images:imreconstructmex:tooManyInputs",
                          "IMRECONSTRUCTMEX takes at most two input arguments.");

    if (mxGetClassID(prhs[0]) != mxGetClassID(prhs[1]))
        mexErrMsgIdAndTxt("Images:imreconstruct:notSameClass","%s",
                          "Function imreconstruct expected MARKER and MASK "
                          "to have the same class.");

    ndims = mxGetNumberOfDimensions(prhs[0]);
    if (ndims != mxGetNumberOfDimensions(prhs[1]))
        mexErrMsgIdAndTxt("Images:imreconstruct:notSameSize",
                          "%s",
                          "Function imreconstruct expected MARKER and MASK "
                          "to be the same size.");

    size_0 = mxGetDimensions(prhs[0]);
    size_1 = mxGetDimensions(prhs[1]);
    for (k = 0; k < ndims; k++) {
        if (size_0[k] != size_1[k])
            mexErrMsgIdAndTxt("Images:imreconstruct:notSameSize",
                              "%s",
                              "Function imreconstruct expected MARKER and "
                              "MASK to be the same size.");
    }

}
//////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////
extern "C"
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    mxArray *J;
    const mxArray *I;
    mwSize num_dims;
    mwSize num_elements;
    const mwSize *input_size;
    void *prJ;
    void *prI;
    Neighborhood_T nhood;
    NeighborhoodWalker_T trailing_walker;
    NeighborhoodWalker_T leading_walker;
    NeighborhoodWalker_T walker;

    (void) nlhs;  /* unused parameter */

    CheckInputs(nrhs, prhs);

    // The reconstruction algorithm works in-place on a copy of the
    //input marker image.  At the end this copy will hold the result.
    J = mxDuplicateArray(prhs[0]);
    plhs[0] = J;
    if (mxIsEmpty(J))
        return;

    I = prhs[1];
    num_elements = mxGetNumberOfElements(prhs[0]);
    num_dims = mxGetNumberOfDimensions(prhs[0]);
    input_size = mxGetDimensions(prhs[0]);

    if (nrhs < 3) 
        nhood = nhMakeDefaultConnectivityNeighborhood(num_dims);
    else
        nhood = nhMakeNeighborhood(prhs[2],NH_CENTER_MIDDLE_ROUNDDOWN);
    
    trailing_walker = nhMakeNeighborhoodWalker(nhood, input_size, num_dims, 
                                      NH_SKIP_CENTER | NH_SKIP_LEADING);
    leading_walker = nhMakeNeighborhoodWalker(nhood, input_size, num_dims, 
                                      NH_SKIP_CENTER | NH_SKIP_TRAILING);
    walker = nhMakeNeighborhoodWalker(nhood, input_size, num_dims,
                                      NH_SKIP_CENTER);
    nhDestroyNeighborhood(nhood);

    prJ = mxGetData(J);
    prI = mxGetData(I);
    
    switch (mxGetClassID(J)) {
    case mxLOGICAL_CLASS:
        compute_reconstruction((mxLogical *)prJ, (mxLogical *)prI, num_elements, 
                               walker, trailing_walker, leading_walker);
        break;

    case mxUINT8_CLASS:
        compute_reconstruction((uint8_T *)prJ, (uint8_T *)prI, num_elements, 
                               walker, trailing_walker, leading_walker);
        break;
        
    case mxUINT16_CLASS:
        compute_reconstruction((uint16_T *)prJ, (uint16_T *)prI, num_elements, 
                               walker, trailing_walker, leading_walker);
        break;
        
    case mxUINT32_CLASS:
        compute_reconstruction((uint32_T *)prJ, (uint32_T *)prI, num_elements, 
                               walker, trailing_walker, leading_walker);
        break;
        
    case mxINT8_CLASS:
        compute_reconstruction((int8_T *)prJ, (int8_T *)prI, num_elements, 
                               walker, trailing_walker, leading_walker);
        break;
        
    case mxINT16_CLASS:
        compute_reconstruction((int16_T *)prJ, (int16_T *)prI, num_elements, 
                               walker, trailing_walker, leading_walker);
        break;
        
    case mxINT32_CLASS:
        compute_reconstruction((int32_T *)prJ, (int32_T *)prI, num_elements, 
                               walker, trailing_walker, leading_walker);
        break;
        
    case mxSINGLE_CLASS:
        do_nan_check((float *)prJ,(float *)prI,num_elements);
        compute_reconstruction((float *)prJ, (float *)prI, num_elements, 
                               walker, trailing_walker, leading_walker);
        break;
        
    case mxDOUBLE_CLASS:
        do_nan_check((double *)prJ,(double *)prI,num_elements);
        compute_reconstruction((double *)prJ, (double *)prI, num_elements, 
                               walker, trailing_walker, leading_walker);
        break;
        
    default:
        mexErrMsgIdAndTxt("Images:imreconstruct:badClass",
                          "%s",
                          "Unsupported input class.");
        break;
    }

    nhDestroyNeighborhoodWalker(trailing_walker);
    nhDestroyNeighborhoodWalker(leading_walker);
    nhDestroyNeighborhoodWalker(walker);
}
