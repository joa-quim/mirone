// $Revision: 1.2.4.4 $
// Copyright 1993-2007 The MathWorks, Inc.

#ifndef _IMFILTER_MEX_H
#define _IMFILTER_MEX_H

#include <string.h>
#include "mex.h"
#include "mwsize.h"
#include "neighborhood.h"
#include "typeconv.h"
#include "mwippl.h"
#include "iptutil_cpp.h"

#define CONV 2
#define SAME 4
#define FULL 8

#define R 0
#define I 1

class Filter {

  private:

    //Member variables
    //////////////////

    MwIppl        *fMwIppl;
    bool           fUseIPPL;
    int            fFlags;
    mwSize        *fImageSize;
    mwSize        fImageNumDims;
    mwSize        *fStart;
    const mxArray *fKernel;
    const mxArray *fPadImage;
    const mxArray *fNonZeroKernel;
    const mxArray *fConn;

    //Pointers to IPPL methods

    ippiFilter32f_8u_C1R_T  ippiFilter32f_8u_C1R;   //float kernel
    ippiFilter32f_16s_C1R_T ippiFilter32f_16s_C1R;
    ippiFilter_32f_C1R_T    ippiFilter_32f_C1R;

    //Other IPPL methods to consider:
    //ippiConvFull_8u_C1R_T   ippiConvFull_8u_C1R;    //full convolution
    //ippiConvFull_16s_C1R_T  ippiConvFull_16s_C1R;
    //ippiConvFull_32f_C1R_T  ippiConvFull_32f_C1R;
    //ippiConvValid_8u_C1R_T  ippiConvValid_8u_C1R;   //valid convolution
    //ippiConvValid_16s_C1R_T ippiConvValid_16s_C1R;
    //ippiConvValid_32f_C1R_T ippiConvValid_32f_C1R;
    //ippiFilter_8u_C1R_T     ippiFilter_8u_C1R;      //32 bit int kernel
    //ippiFilter_16s_C1R_T    ippiFilter_16s_C1R;

    //Methods
    /////////

    mwSize update_linear_index(mwSize  *current_pos,
                               mwSize  *image_size,
                               mwSize  *padded_cumprod, 
                               mwSize   num_dims,
                               mwSize  *start);

    bool useIPPL(void);

    //Templetized methods
    /////////////////////

    //////////////////////////////////////////////////////////////////////////
    //  ndim_filter performs N-dimensional filtering on input buffer inBuf
    /////////////////////////////////////////////////////////////////////////
    template<typename _T1>
        void ndim_filter(_T1 *inBuf, _T1 *outBuf) {
            // The neighborhood and neighborhood walker is used to reflect 
            // the neighborhood for convolution and to calculate the filter
            // offset values.            
            Neighborhood_T nhood;
            NeighborhoodWalker_T walker;
            
            const mwSize *padSize = mxGetDimensions(fPadImage); 
            mwSize    *inSize     = fImageSize;
            mwSize     numInDims  = fImageNumDims;
            mwSize     numPadDims = mxGetNumberOfDimensions(fPadImage);
            
            // Do convolution only if specified by the user, otherwise, do
            // correlation. If doing convolution, the center needs to 
            // be placed "opposite" of where we want it after the 
            // reflection.
            int centerLoc;
            if( fFlags & FULL )
                centerLoc = (fFlags & CONV) ? NH_CENTER_UL : NH_CENTER_LR;
            else 
                centerLoc = (fFlags & CONV) ? NH_CENTER_MIDDLE_ROUNDUP : NH_CENTER_MIDDLE_ROUNDDOWN;

            nhood = nhMakeNeighborhood(fConn, centerLoc);
            if(fFlags & CONV) nhReflectNeighborhood(nhood);
    
            walker = nhMakeNeighborhoodWalker(nhood,padSize, numPadDims,0);

            mwSignedIndex *filter_offsets = nhGetWalkerNeighborOffsets(walker);
            
            mwSize *padded_cumprod = (mwSize *)mxMalloc(sizeof(*padded_cumprod) * (numPadDims+1));
            
            //Assuming at least one pixel in image. This works for 
            //empty inputs too.
            mwSize num_image_pixels = 1;
            padded_cumprod[0] = 1;
            
            for(mwSize i=0; i < numPadDims; i++)
                padded_cumprod[i+1] = padded_cumprod[i]*padSize[i];
            
            for(mwSize i=0; i < numInDims; i++)
                num_image_pixels *= inSize[i];
            
            /* Initialize current position */
            mwSize *current_pos = (mwSize *) mxMalloc(sizeof(*current_pos) *
                                                numPadDims);
            for (mwSize k = 0; k < numPadDims; k++)
                current_pos[k] = fStart[k];
            mwSize p = sub_to_ind(current_pos,padded_cumprod,
                                  numPadDims);
            
            // Convolve
            mwSize     numKernElem = mxGetNumberOfElements(fNonZeroKernel);
            double *nonZeroKernel  = mxGetPr(fNonZeroKernel);
            for(mwSize i=0; i<num_image_pixels; i++) {
                double sum = 0;
                for(mwSize j=0; j < numKernElem; j++)
                    sum += inBuf[ p+filter_offsets[j] ]*nonZeroKernel[j];
                saturateRoundAndCast(&outBuf[i],sum);
                
                p = update_linear_index(current_pos,inSize,padded_cumprod, numPadDims,fStart);
            }
    
            // Clean up
            nhDestroyNeighborhood(nhood);
            nhDestroyNeighborhoodWalker(walker);
            mxFree(padded_cumprod);
            mxFree(current_pos); 
        }

    //////////////////////////////////////////////////////////////////////////
    //
    //////////////////////////////////////////////////////////////////////////
    template<typename _T1, typename _T2>
        void ipplFilterFlpKern(_T1 *inBuf, _T1 *outBuf, _T2 funcPtr)
        {
#ifdef __i386__
            double *dblKernel = (double *)mxGetData(fKernel);
            mwSize kRows = mxGetM(fKernel);
            mwSize kCols = mxGetN(fKernel);
            
            mwSize numKernelElements = kRows*kCols;                  
            float *fltKernel = (float *)mxCalloc(numKernelElements, 
                                                 sizeof(float));
            
            //Repackage the kernel into a float array
            if(fFlags & CONV) {
                for (mwSize j=0; j<numKernelElements; j++)
                    fltKernel[j] = (float)dblKernel[j];
            }
            else
            {   //rotate the kernel by 180 degrees to get correlation
                for (mwSize j=0; j<numKernelElements; j++)
                    fltKernel[j] = (float)dblKernel[numKernelElements-j-1];
            }
            
            //IN DIMS
            IppiSize kernelSize = {static_cast<int>(kRows), 
                                   static_cast<int>(kCols)};
            const mwSize *inDims = mxGetDimensions(fPadImage);
            int inRows = static_cast<int>(inDims[0]);
            int inCols = static_cast<int>(inDims[1]);
            mwSize inPlanes = (mxGetNumberOfElements(fPadImage) == 0) ? 
                0 : mxGetNumberOfElements(fPadImage) / inRows / inCols;
            
            //OUT DIMS
            int outRows = static_cast<int>(fImageSize[0]);
            int outCols = static_cast<int>(fImageSize[1]);
            
            //Intel's anchor point with respect to UL corner of the kernel.
            //Notice that for 'same' conv or corr, the anchor in some of 
            //the cases may be outside of the kernel dimensions; although 
            //I don't believe that Intel software developers intended
            //their code being used with such inputs, the resulting
            //matrices match the desired Matlab output.  For full correlation
            //and convolution, it will always work out to be the lower 
            //right corner of the kernel.
            IppiPoint anchor = {inRows - outRows, inCols-outCols};
            
            IppiSize dstRoiSize = {outRows, outCols};                  
            int srcStep = inRows*sizeof(_T1);
            int dstStep = dstRoiSize.width*sizeof(_T1);

            for (unsigned int k = 0; k < inPlanes; k++) {
                //Invoke Intel's code
                IppStatus statusCode = (*funcPtr)
                    (inBuf + k*inRows*inCols,srcStep,
                     outBuf + k*outRows*outCols,dstStep,
                     dstRoiSize,fltKernel,kernelSize,
                     anchor);
                
                fMwIppl->ipplCheckStatus(statusCode);
            }

            if (fltKernel) mxFree(fltKernel);
#else
            // unused parameters
            (void) inBuf;
            (void) outBuf;
            (void) funcPtr;

            mexErrMsgIdAndTxt("Images:imfilter:badIppPlatform",
                              "Internal problem: reached ipplFilterFlpKern() unexpectedly.");

#endif
        }

  public:

    Filter();
    virtual ~Filter();

    void setInSize(const mxArray *inSize);
    void setFlags(const mxArray *kernel);
    void setKernel(const mxArray *kernel);
    void setPadImage(const mxArray *padImage);
    void setNonZeroKernel(const mxArray *nonZeroKernel);
    void setConn(const mxArray *conn);
    void setStart(const mxArray *start);

    mxArray *evaluate(void);
};

#endif
