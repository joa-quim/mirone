/*
 * Copyright 1993-2007 The MathWorks, Inc.
 * $Revision: 1.5.4.3 $  $Date: 2007/06/04 21:10:05 $
 */

#ifndef MORPHMEX_H
#define MORPHMEX_H

#include <math.h>
#include "neighborhood.h"
#include "typeconv.h"

#define BITS_PER_WORD 32
#define LEFT_SHIFT(x,shift) (shift == 0 ? x : (shift == BITS_PER_WORD ? 0 : x << shift))
#define RIGHT_SHIFT(x,shift) (shift == 0 ? x : (shift == BITS_PER_WORD ? 0 : x >> shift))

void dilate_logical(mxLogical *In, mxLogical *Out, mwSize num_elements,
                    NeighborhoodWalker_T walker);

void erode_logical(mxLogical *In, mxLogical *Out, mwSize num_elements,
                   NeighborhoodWalker_T walker);

//////////////////////////////////////////////////////////////////////////////
// Perform flat grayscale dilation on a uint8 array.
//
// Inputs
// ======
// In             - pointer to first element of input array
// num_elements   - number of elements in input and output arrays
// walker         - neighborhood walker corresponding to reflected structuring
//                  element
//
// Output
// ======
// Out            - pointer to first element of output array
//
// Note that implementing dilation properly with this function requires 
// passing in a neighborhood walker constructed from a reflected neighborhood.
//////////////////////////////////////////////////////////////////////////////
template<typename _T>
void dilateGrayFlat(_T *In, _T *Out, mwSize num_elements, 
                    NeighborhoodWalker_T walker)
{
    for (mwSize p = 0; p < num_elements; p++)
    {
        _T val;
        _T new_val;
        mwSize q;

        setMin(&val);
        nhSetWalkerLocation(walker, p);
        while (nhGetNextInboundsNeighbor(walker, &q, NULL))
        {
            new_val = In[q];
            if (new_val > val)
            {
                val = new_val;
            }
        }
        Out[p] = val;
    }
}

//////////////////////////////////////////////////////////////////////////////
// Perform flat grayscale dilation on input array.
//
// Inputs
// ======
// In             - pointer to first element of input array
// num_elements   - number of elements in input and output arrays
// walker         - neighborhood walker corresponding to reflected structuring
//                  element
//
// Output
// ======
// Out            - pointer to first element of output array
//////////////////////////////////////////////////////////////////////////////
template<typename _T>
void erodeGrayFlat(_T *In, _T *Out, mwSize num_elements, 
                   NeighborhoodWalker_T walker)
{
    for (mwSize p = 0; p < num_elements; p++)
    {
        mwSize q;
        _T val;
        _T new_val;

        setMax(&val);
        nhSetWalkerLocation(walker, p);
        while (nhGetNextInboundsNeighbor(walker, &q, NULL))
        {
            new_val = In[q];
            if (new_val < val)
            {
                val = new_val;
            }
        }
        Out[p] = val;
    }
}

void dilate_gray_nonflat_uint8(uint8_T *In, uint8_T *Out, mwSize num_elements,
                               NeighborhoodWalker_T walker, double *heights);

void dilate_gray_nonflat_uint16(uint16_T *In, uint16_T *Out, mwSize num_elements,
                                NeighborhoodWalker_T walker, double *heights);

void dilate_gray_nonflat_uint32(uint32_T *In, uint32_T *Out, mwSize num_elements,
                                NeighborhoodWalker_T walker, double *heights);

void dilate_gray_nonflat_int8(int8_T *In, int8_T *Out, mwSize num_elements,
                              NeighborhoodWalker_T walker, double *heights);

void dilate_gray_nonflat_int16(int16_T *In, int16_T *Out, mwSize num_elements,
                               NeighborhoodWalker_T walker, double *heights);

void dilate_gray_nonflat_int32(int32_T *In, int32_T *Out, mwSize num_elements,
                               NeighborhoodWalker_T walker, double *heights);

void dilate_gray_nonflat_single(float *In, float *Out, mwSize num_elements,
                                NeighborhoodWalker_T walker, double *heights);

void dilate_gray_nonflat_double(double *In, double *Out, mwSize num_elements,
                                NeighborhoodWalker_T walker, double *heights);

void erode_gray_nonflat_uint8(uint8_T *In, uint8_T *Out, mwSize num_elements,
                              NeighborhoodWalker_T walker, double *heights);

void erode_gray_nonflat_uint16(uint16_T *In, uint16_T *Out, mwSize num_elements,
                               NeighborhoodWalker_T walker, double *heights);

void erode_gray_nonflat_uint32(uint32_T *In, uint32_T *Out, mwSize num_elements,
                               NeighborhoodWalker_T walker, double *heights);

void erode_gray_nonflat_int8(int8_T *In, int8_T *Out, mwSize num_elements,
                             NeighborhoodWalker_T walker, double *heights);

void erode_gray_nonflat_int16(int16_T *In, int16_T *Out, mwSize num_elements,
                              NeighborhoodWalker_T walker, double *heights);

void erode_gray_nonflat_int32(int32_T *In, int32_T *Out, mwSize num_elements,
                              NeighborhoodWalker_T walker, double *heights);

void erode_gray_nonflat_single(float *In, float *Out, mwSize num_elements,
                               NeighborhoodWalker_T walker, double *heights);

void erode_gray_nonflat_double(double *In, double *Out, mwSize num_elements,
                               NeighborhoodWalker_T walker, double *heights);

void dilate_packed_uint32(uint32_T *In, uint32_T *Out, mwSize M, mwSize N,
                          mwSignedIndex *rc_offsets, mwSize num_neighbors);

void erode_packed_uint32(uint32_T *In, uint32_T *Out, mwSize M, mwSize N,
                         mwSignedIndex *rc_offsets, mwSize num_neighbors,
                         mwSize unpacked_M);


#endif /* MORPHMEX_H */
