// Copyright 1993-2007 The MathWorks, Inc.
// $Revision.1 $  $Date: 2007/06/04 21:09:33 $

#ifndef _GRAYXFORM_H
#define _GRAYXFORM_H

#include "mex.h"
#include "mwsize.h"

//////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////
template<typename _T>
void transformDoubleSingleArray(_T *a, double *t, mwSize maxTidx, mwSize nelements, _T *out) {
    _T val;
    mwSize index;
    
    for(mwSize i=0; i < nelements; i++) {
        // Clip a[i] to the range [0,1].
        val = a[i];
        if ((val >= 0.0) && (val <= 1.0)) {
            index = (mwSize) (val * maxTidx + 0.5);
        }
        else if (val > 1.0) {
            index = maxTidx;
        }
        else { // if ((val < 0.0) || mxIsNaN((double)val)) 
            index = 0;
        }
        out[i] = (_T)t[index];
    }
}

//////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////
template<typename _T>
void transformUintArray(_T *a, double *t, mwSize maxTidx, mwSize nelements, _T *out, double scaleFactor) {
    if(maxTidx == (mwSize)scaleFactor) {
        // Perfect fit, we don't need to scale the index.
        for(mwSize i=0; i < nelements; i++) 
            out[i] = (_T) (scaleFactor * t[a[i]] + 0.5);
    }
    else {
        mwSize index;

        // Scale the index by maxTidx.
        double scale = maxTidx / scaleFactor;

        for(mwSize i=0; i < nelements; i++) {
            index = (mwSize) (a[i] * scale + 0.5);
            out[i] = (_T) (scaleFactor * t[index] + 0.5);
        }
    }
}

#endif //_GRAYXFORM_H
