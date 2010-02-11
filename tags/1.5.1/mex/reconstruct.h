// Copyright 1993-2007 The MathWorks, Inc.
// $Revision: 1.1.6.5 $  $Date: 2007/06/04 21:10:17 $

#ifndef _RECONSTRUCT_H
#define _RECONSTRUCT_H

#include "mex.h"
#include "mwsize.h"
#include "neighborhood.h"
#include "queue.h"

//////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////
template<typename _T>
inline void do_nan_check(_T *J,_T *I,mwSize num_elements) {
    for (mwSize p = 0; p < num_elements; p++) {
        if (mxIsNaN(J[p]) || mxIsNaN(I[p]))
            mexErrMsgIdAndTxt("Images:imreconstruct:containsNaN",
                              "%s",
                              "MARKER and MASK may not contain NaNs.");
    }
}     

//////////////////////////////////////////////////////////////////////////////
//
// This file contains the function body for the image reconstruction
// algorithm.  
//
// Algorithm reference: Luc Vincent, "Morphological Grayscale Reconstruction
// in Image Analysis: Applications and Efficient Algorithms," IEEE
// Transactions on Image Processing, vol. 2, no. 2, April 1993, pp. 176-201.
//
// The algorithm used here is "fast hybrid grayscale reconstruction,"
// described as follows on pp. 198-199:
//
// I: mask image (binary or grayscale)
// J: marker image, defined on domain D_I, J <= I.
//   Reconstruction is determined directly in J.
//
// Scan D_I in raster order:
// Let p be the current pixel;
// J(p) <- (max{J(q),q member_of N_G_plus(p) union {p}}) ^ I(p)
//      [Note that ^ here refers to "pointwise minimum.]
//
// Scan D_I in antiraster order:
//  Let p be the current pixel;
//  J(p) <- (max{J(q),q member_of N_G_minus(p) union {p}}) ^ I(p)
//      [Note that ^ here refers to "pointwise minimum.]
//  If there exists q member_of N_G_minus(p) such that J(q) < J(p) and
//      J(q) < I(q), then fifo_add(p)
//
// Propagation step:
//  While fifo_empty() is false
//  p <- fifo_first()
//  For every pixel q member_of N_G(p):
//    If J(q) < J(p) and I(q) ~= J(q), then
//      J(q) <- min{J(p),I(q)}
//      fifo_add(q)
//
//////////////////////////////////////////////////////////////////////////////
template<typename _T>
void compute_reconstruction(_T *J, _T *I, mwSize numElements, 
                         NeighborhoodWalker_T walker,
                         NeighborhoodWalker_T trailingWalker,
                         NeighborhoodWalker_T leadingWalker) {
    _T Jp, Jq, Iq;
    _T maxPixel;
    mwSize p, q, k;
    mwSize pp;

    // Enforce the requirement that J <= I.  We need to check this here,
    // because if it isn't true the algorithm might not terminate.

    for (k = 0; k < numElements; k++) {
        if (J[k] > I[k])
            mexErrMsgIdAndTxt("Images:imreconstruct:markerGreaterThanMask",
                              "%s", "MARKER pixels must be <= MASK pixels.");
    }
    
    Queue<mwSize> Queue;
    Queue.initialize(0);

    // First pass,  Scan D_I in raster order (upper-left to lower-right,
    // along the columns):

    for (p = 0; p < numElements; p++) {
        // "Let p be the current pixel"
        // "J(p) <- (max{J(q),q member_of N_G_plus(p) union {p}}) ^ I(p)"

        // Find the maximum value of the (y,x) pixel
        // plus all the pixels in the "plus" neighborhood 
        // of (y,x).

        maxPixel = J[p];
        nhSetWalkerLocation(trailingWalker, p);
        while (nhGetNextInboundsNeighbor(trailingWalker, &q, NULL)) {
            if (J[q] > maxPixel)
                maxPixel = J[q];
        }

        // Now set the (y,x) pixel of image J to the minimum
        // of maxPixel and the (y,x) pixel of image I.

        J[p] = (maxPixel < I[p]) ? maxPixel : I[p]; 
    }

    // Second pass,  Scan D_I in antiraster order (lower-right to upper-left,
    // along the columns):

    for (pp = 0; pp < numElements; pp++) {
        p = numElements - 1 - pp;

        // "Let p be the current pixel"
        // "J(p) <- (max{J(q),q member_of N_G_minus(p) union {p}}) ^ I(p)"

        // Find the maximum value of the (y,x) pixel
        // plus all the pixels in the "minus" neighborhood 
        // of (y,x).

        maxPixel = J[p];
        nhSetWalkerLocation(leadingWalker, p);
        while (nhGetNextInboundsNeighbor(leadingWalker, &q, NULL)) {
            if (J[q] > maxPixel)
                maxPixel = J[q];
        }
    

        // Now set the (y,x) pixel of image J to the minimum
        // of maxPixel and the (y,x) pixel of image I.

        J[p] = (maxPixel < I[p]) ? maxPixel : I[p];
        
        // "If there exists q member_of N_G_minus(p) 
        // such that J(q) < J(p) and J(q) < I(q), then fifo_add(p)"

        nhSetWalkerLocation(leadingWalker, p);
        while (nhGetNextInboundsNeighbor(leadingWalker, &q, NULL)) {
            if ((J[q] < J[p]) && (J[q] < I[q])) {
                Queue.put(p);
                break;
            }
        }
    }

    // "Propagation step:
    // While fifo_empty() is false"

    while (Queue.getSequenceLength() > 0) {
        // "p <= fifo_first()"

        p = Queue.get();
        Jp = J[p];

        // "For every pixel q member_of N_g(p):"

        nhSetWalkerLocation(walker, p);
        while (nhGetNextInboundsNeighbor(walker, &q, NULL)) {
            
            // "If J(q) < J(p) and I(q) ~= J(q), then
            //  J(q) <- min{J(p),I(q)}
            //  fifo_add(q)"

            Jq = J[q];
            Iq = I[q];
            if ((Jq < Jp) && (Iq != Jq)) {
                J[q] = (Jp < Iq) ? Jp : Iq;
                Queue.put(q);
            }
        }
    }

    Queue.freeSequence();
}

#endif
