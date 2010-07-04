// 
// Copyright 1993-2007 The MathWorks, Inc.
// $Revision: 1.1.6.3 $
//

#ifndef RESAMPSEP_TYPES_H
#define RESAMPSEP_TYPES_H

//////////////////////////////////////////////////////////////////////////////
// Typedefs
//////////////////////////////////////////////////////////////////////////////

/*============================ Iterator Struct =============================*/

typedef struct Iterator {
    mwSize   ndims;         // number of dimensions
    mwSignedIndex   offset; // current offset value
    mwSize   length;        // number subscript values
    mwSignedIndex*  size;   // size of each dimension
    mwSignedIndex*  subs;   // current subscript vector
}
Iterator;

/*========================== Config Struct =================================*/

/*
 * The transform dimensions in tdims may be out of sequence. That's
 * fine, but the values in tsize and cpTrans must
 * be consistent with the order used in tdims. cpTrans and cpOther
 * each contain elements of the cumulative product vector of size,
 * with a one inserted at the front. The elements corresponding
 * to tdims are ordered the same way as the positions listed in
 * tdims -- so if the elements of tdims are out of sequence, then
 * the values in cpTrans will not be monotonically increasing.
 * No problem.
 */

typedef struct Config {
    mwSize  ndims;    // total number of dimensions
    mwSize  nTrans;   // number of transform dimensions
    mwSize  nOther;   // number of non-transform dimensions
    mwSize  tlength;  // total number of transform dimension elements
    mwSize* size;     // size of each dimension, in order
    mwSize* tdims;    // position of each transform dimension
    mwSignedIndex* tsize;  // size of each transform dimension (sometimes set to -1)
    mwSignedIndex* osize;    // size of each non-transform dimension
    mwSize* cpTrans;  // cumulative product at the position of each transform 
                      // dimension
    mwSize* cpOther;  // cumulative product at the position of each 
                      // non-transform dimension
}
Config;

/*========================== Kernel Struct =================================*/

typedef struct Kernel {
    double      halfwidth;
    mwSize      nSamplesInPositiveHalf;
    double*     positiveHalf;
    mwSize      stride;      // The maximum number of integer values that 
                             // the kernel can span, given any possible shift.
    double      indexFactor; // Precomputing this value will save us time in
                             // EvaluateKernel function
}
Kernel;

/*=================== Pad Method and Convolver Struct =====================*/

typedef enum { Fill, Bound, Replicate, Circular, Symmetric } PadMethod;

typedef struct Convolver {
    PadMethod padmethod;          // Method for defining values for points that map
                                  // outside A
    mwSize           ndims;       // Number of input transform dimensions
    mwSize           nPoints;     // Total number of interpolating points
    mwSignedIndex*   size;        // Length of the weight array for each transform 
                                  // dimension
    mwSize*          cumsum;      // Cumulative sum of sizes, with a zero inserted 
                                  // at the front
    mwSize*          cumprod;     // Cumulative product of sizes, with a one inserted
                                  // at the front
    double**         weights;     // Array of pointers to the weight arrays for each
                                  // transform dimension
    mwSignedIndex*   tsize;       // Input tsize

    mwSignedIndex**  tsub;        // Array pointers to the subscript arrays for each
                                  // transform dimension
    mwSignedIndex**  subs;        // A list of transform subscripts for each local 
                                  // point
    bool*            useFill;     // Boolean array, one value for each local point

    double*          weight_data; // Storage for the arrays pointed to by weights
    mwSignedIndex*   tsub_data;   // Storage for the arrays pointed to by tsub
    mwSignedIndex*   subs_data;   // Storage for the arrays pointed to by subs

    double*          lo;          // Lower bounds for range checking
    double*          hi;          // Upper bounds for range checking

    Kernel**         kernelset;   // An interpolating kernel for each input
                                  // transform dimension
    Iterator*        localIt;     // Iterator for local grid points in transform
                                  // space
 }
Convolver;

#endif
