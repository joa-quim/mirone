// Copyright 1993-2007 The MathWorks, Inc.
  
//////////////////////////////////////////////////////////////////////////////
//  Helper MEX-file for IMFILTER.  
//  
//  Inputs:
//  prhs[0] -         - Padded image
//  prhs[1] - int32_T - Size of unpadded image
//  prhs[2] - double  - Filter coefficients (all)
//  prhs[3] - double  - Filter coefficients (non-zero)
//  prhs[4] - double  - Connectivity matrix used for neighborhood creation.
//  prhs[5] - int32_T - Starting point in padded image
//  prhs[6] - int     - OR'd flags for user defined options
//////////////////////////////////////////////////////////////////////////////

#include <math.h>
#include "mex.h"
#include "mwsize.h"
#include "imfilter_mex.h"
#include "neighborhood.h"
#include "iptutil.h"

//////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////
mxArray *Filter::evaluate(void) {
    mxAssert((fFlags != 0), ERR_STRING("Filter::fFlags","evaluate()"));
    mxAssert((fImageSize != NULL), ERR_STRING("Filter::fImageSize","evaluate()"));
    mxAssert((fKernel != NULL), ERR_STRING("Filter::fKernel","evaluate()"));
    mxAssert((fPadImage != NULL), 
             ERR_STRING("Filter::fPadImage","evaluate()"));
    mxAssert((fConn != NULL), ERR_STRING("Filter::fConn","evaluate()"));
    mxAssert((fStart != NULL), ERR_STRING("Filter::fStart","evaluate()"));
    mxAssert((fNonZeroKernel != NULL), 
             ERR_STRING("Filter::fNonZeroKernel","evaluate()"));


    void *In[2], *Out[2];
    mxClassID image_class = mxGetClassID(fPadImage);
    mwSize   *image_size  = fImageSize;
    mwSize    numInDims   = fImageNumDims;
    bool      cmplx       = mxIsComplex(fPadImage);

    // The filter will always be real (prhs[2]) but the image may be complex
    // (prhs[0]). If the filter is complex, the M-code will call imfilter_mex 
    // once each for the real and imaginary parts of the filter.
    mxArray *mxOutput;
    In[R] = mxGetData(fPadImage);

    if(cmplx)
    {
	In[I] = mxGetImagData(fPadImage);
	mxOutput =  mxCreateNumericArray(numInDims, image_size, image_class, mxCOMPLEX);
        Out[I] = mxGetImagData(mxOutput);
    }
    else
	mxOutput =  mxCreateNumericArray(numInDims, image_size, image_class, mxREAL);
    
    Out[R] = mxGetData(mxOutput);



    //Filter real and imaginary (if available) parts of the image.
    for (int i=0; i < (cmplx ? 2 : 1); i++) {    
        switch (image_class) {
          case mxLOGICAL_CLASS:
            ndim_filter((mxLogical *)In[i], (mxLogical *)Out[i]);
#if defined(__APPLE__) && defined(__ppc__) //Because mxLogical = uint8_T on the MAC
            {
                mwSize numInPixels = 1;   
                
                for(mwSize k=0; k < numInDims; k++) {
                    //image_size holds the dimensions of the input image
                    //before padding
                    numInPixels *= image_size[k];
                }
                
                for(mwSize k=0; k < numInPixels; k++) {
                    //Anything greater than 0 becomes true thus converting
                    //uint8 back to logical
                    ((mxLogical *)Out[i])[k] = (((mxLogical *)Out[i])[k] > 0);
                }
            }
#endif
            break;
            
          case mxUINT8_CLASS:

            if(useIPPL()) {
                ipplFilterFlpKern((uint8_T *)In[i], (uint8_T *)Out[i],
                                  ippiFilter32f_8u_C1R);
                break;
            }

            ndim_filter((uint8_T *)In[i], (uint8_T *)Out[i]);
            break;
            
          case mxINT8_CLASS:
            ndim_filter((int8_T *)In[i], (int8_T *)Out[i]);
            break;
            
          case mxUINT16_CLASS:
            ndim_filter((uint16_T *)In[i], (uint16_T *)Out[i]);
            break;
            
          case mxINT16_CLASS:

            if(useIPPL()) {
                ipplFilterFlpKern((int16_T *)In[i], (int16_T *)Out[i], ippiFilter32f_16s_C1R);
                break;
            }

            ndim_filter((int16_T *)In[i], (int16_T *)Out[i]);
            break;
            
          case mxUINT32_CLASS:
            ndim_filter((uint32_T *)In[i], (uint32_T *)Out[i]);
            break;
            
          case mxINT32_CLASS:
            ndim_filter((int32_T *)In[i], (int32_T *)Out[i]);
            break;
                        
          case mxSINGLE_CLASS:
            if(useIPPL()) {
                ipplFilterFlpKern((float *)In[i], (float *)Out[i], ippiFilter_32f_C1R);
                break;
            }

            ndim_filter((float *)In[i], (float *)Out[i]);
            break;
            
          case mxDOUBLE_CLASS:
            ndim_filter((double *)In[i], (double *)Out[i]);
            break;

          default:
            mexErrMsgIdAndTxt("Images:imfilter_mex:unsupportedInputClass",
                              "%s","Unsupported input class.");
            break;
        }
    }

    return(mxOutput);
}

//////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////
mwSize Filter::update_linear_index(mwSize *current_pos, mwSize *image_size,
                                mwSize *padded_cumprod, mwSize num_dims, mwSize *start)
{
    current_pos[0]++;
    for(mwSize i=0;i<num_dims-1;i++) {
        if( static_cast<mwSignedIndex>(current_pos[i]) - 
            static_cast<mwSignedIndex>(start[i]) > 
            (static_cast<mwSignedIndex>(image_size[i])-1) ) {
            current_pos[i] = start[i];
            current_pos[i+1]++;
        }
        else
            break;
    }
    
    return(sub_to_ind(current_pos,padded_cumprod,num_dims));
}

//////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////
bool Filter::useIPPL(void) {
    const mwSize maxNumelForIPP = 0xFFFFFFFFUL;

    return fUseIPPL && (mxGetNumberOfDimensions(fKernel) == 2) &&
        (mxGetNumberOfElements(fPadImage) <= maxNumelForIPP) &&
        (mxGetNumberOfElements(fKernel)   <= maxNumelForIPP);
}

//////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////
Filter::Filter() {
    //Initialize member variables
    fFlags        = 0;
    fImageNumDims = 0;
    fImageSize    = NULL;
    fKernel       = NULL;
    fPadImage     = NULL;
    fConn         = NULL;
    fStart        = NULL;
    fNonZeroKernel = NULL;

#ifdef __i386__  //Don't use this code on non-Intel platforms

    fMwIppl  = new MwIppl();
    fUseIPPL = true;
    
    ippiFilter32f_8u_C1R = (ippiFilter32f_8u_C1R_T)
        fMwIppl->getMethodPointer("ippiFilter32f_8u_C1R", ippI);
    if(!ippiFilter32f_8u_C1R) fUseIPPL = false;

    ippiFilter32f_16s_C1R = (ippiFilter32f_16s_C1R_T)
        fMwIppl->getMethodPointer("ippiFilter32f_16s_C1R", ippI);
    if(!ippiFilter32f_16s_C1R) fUseIPPL = false;

    ippiFilter_32f_C1R = (ippiFilter_32f_C1R_T)
        fMwIppl->getMethodPointer("ippiFilter_32f_C1R", ippI);
    if(!ippiFilter_32f_C1R) fUseIPPL = false;

#else
    fUseIPPL = false;

#endif

}

//////////////////////////////////////////////////////////////////////////////
// Note:  The destructor will be called on matlab exit or
//        when using matlab's "clear" function
//////////////////////////////////////////////////////////////////////////////
Filter::~Filter() {
    mxFree(fImageSize);
    
#ifdef __i386__  //Don't use this code on non-Intel platforms
    if(fMwIppl  != NULL) delete fMwIppl;
#endif

}

//////////////////////////////////////////////////////////////////////////////
// Size of unpadded image (int32_T)
//////////////////////////////////////////////////////////////////////////////
void Filter::setInSize(const mxArray *inSize) {
    mxFree(fImageSize);

    fImageNumDims = mxGetNumberOfElements(inSize);
    fImageSize = (mwSize *) mxMalloc(fImageNumDims * sizeof(*fImageSize));
    mexMakeMemoryPersistent(fImageSize);

    double *pr = (double *) mxGetData(inSize);
    for (mwSize k = 0; k < fImageNumDims; k++)
        fImageSize[k] = static_cast<mwSize>(pr[k]);
}

//////////////////////////////////////////////////////////////////////////////
// OR'd flags for user defined options (int)
//////////////////////////////////////////////////////////////////////////////
void Filter::setFlags(const mxArray *flags) {
    fFlags = (int)mxGetScalar(flags);
}
//////////////////////////////////////////////////////////////////////////////
// Full kernel (double)
//////////////////////////////////////////////////////////////////////////////
void Filter::setKernel(const mxArray *kernel) {
    fKernel = kernel;
}
//////////////////////////////////////////////////////////////////////////////
// Non-zero filter coefficients (double)
//////////////////////////////////////////////////////////////////////////////
void Filter::setNonZeroKernel(const mxArray *nonZeroKernel) {
    fNonZeroKernel = nonZeroKernel;
}
//////////////////////////////////////////////////////////////////////////////
// padded image
//////////////////////////////////////////////////////////////////////////////
void Filter::setPadImage(const mxArray *padImage) {
    fPadImage = padImage;
}
//////////////////////////////////////////////////////////////////////////////
// Connectivity matrix used for neighborhood creation (double)
//////////////////////////////////////////////////////////////////////////////
void Filter::setConn(const mxArray *conn) {
    // Use the neighborhood function nhCheckDomain to check
    // the connectivity matrix.
    nhCheckDomain(conn);
    fConn = conn;
}
//////////////////////////////////////////////////////////////////////////////
// Starting point in padded image (int32_T)
//////////////////////////////////////////////////////////////////////////////
void Filter::setStart(const mxArray *start) {
    mxFree(fStart);

    mwSize numel = mxGetNumberOfElements(start);
    fStart = (mwSize *) mxMalloc(numel * sizeof(*fStart));
    mexMakeMemoryPersistent(fStart);

    double *pr = (double *) mxGetData(start);
    for (mwSize k = 0; k < numel; k++) {
        fStart[k] = static_cast<mwSize>(pr[k]);
    }
}

//////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////
Filter filter;

extern "C"
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    (void) nlhs;  // unused parameter
 
    if (nrhs != 7) {
        mexErrMsgIdAndTxt("Images:imfilter_mex:invalidNumInputs", "%s",
                          "IMFILTER_MEX needs 7 input arguments.");
    }

    filter.setPadImage(prhs[0]);
    filter.setInSize(prhs[1]);
    filter.setKernel(prhs[2]);
    filter.setNonZeroKernel(prhs[3]);
    filter.setConn(prhs[4]);
    filter.setStart(prhs[5]);
    filter.setFlags(prhs[6]);

    //Filter the image
    plhs[0] = filter.evaluate();
}
