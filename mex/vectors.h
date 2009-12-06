/* Copyright 2003-2007 The MathWorks, Inc. */

/*
 * $Revision: 1.1.8.3 $
 */

#include "mex.h"
#include "mwsize.h"

/*
 * dimReshape computes a "pseudo-size" (M1-by-M2-by-M3) for an array.
 *
 * Input parameter numdims is the return value of mxGetNumberOfDimensions().
 *
 * Input parameter size is the return value of mxGetDimensions().
 *
 * Input parameter dim is the desired action dimension, zero-based.
 * dim may be any nonnegative integer.  There is no requirement that dim be
 * less than numdims.
 *
 * Output parameter M1 is the product of the sizes of the original array 
 * along all dimensions less than dim.
 *
 * Output parameter M2 is the size of the original array along the specified 
 * dimension, dim.
 *
 * Output parameter M3 is the product of the sizes of the original array 
 * along all dimensions greater than dim.
 *
 * Example: numdims = 3, size[] = {5, 3, 2}
 *
 *     dim    M1    M2    M3
 *     ---    --    --    --
 *      0     1     5     6
 *      1     5     3     2
 *      2     15    2     1
 *      3     30    1     1
 *      4     30    1     1
 *      ...   ...   ...   ...
 *      
 */
void dimReshape(mwSize numdims, const mwSize *size, mwSize dim, mwSize *M1, mwSize *M2, mwSize *M3);

/*
 * copyVector copies a vector of data from the input buffer to the output
 * buffer.
 *
 * Input parameter in_buffer points to the beginning of the input buffer.
 *
 * Input parameter out_buffer points to the beginning of the output buffer.
 * The contents of the output buffer are modified by this operation.
 *
 * Input parameter element_size is the size, in bytes, of each element
 * of the vector.
 *
 * Input parameter in_vector_stride is the number of bytes from one
 * input vector element to the next.
 *
 * Input parameter out_vector_stride is the number of bytes from one
 * output vector element to the next.
 *
 * Input parameter vector_length is the number of vector elements.
 *
 * Implementation note: if in_vector_stride and out_vector_stride are
 * both 1, the operation is performed with a single call to memcpy.
 */
void copyVector(uint8_T *in_buffer,
                uint8_T *out_buffer,
                size_t   element_size,
                mwSize   in_vector_stride,
                mwSize   out_vector_stride,
                mwSize   vector_length);
