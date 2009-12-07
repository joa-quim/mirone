/*
 * Copyright 1993-2007 The MathWorks, Inc.
 * $Revision: 1.1.6.2 $  $Date: 2007/06/04 21:09:21 $
 *
 * Packed 2-D binary dilation and erosion.
 */

#include "morphmex.h"

#define BITS_PER_WORD 32
#define LEFT_SHIFT(x,shift) (shift == 0 ? x : (shift == BITS_PER_WORD ? 0 : x << shift))
#define RIGHT_SHIFT(x,shift) (shift == 0 ? x : (shift == BITS_PER_WORD ? 0 : x >> shift))

/*
 * dilate_packed_uint32
 * Packed binary dilation
 *
 * Inputs
 * ======
 * In            - pointer to first element of input array
 * M             - number of rows of packed input array
 * N             - number of columns of packed input array
 * rc_offsets    - Row-column offset locations corresponding to
 *                 each element of the structuring element; storage
 *                 order is that for a MATLAB array, num_neighbors-by-2
 * num_neighbors - number of neighbors in the structuring element
 *
 * Output
 * ======
 * Out           - pointer to first element of output array
 */
void dilate_packed_uint32(uint32_T *In, uint32_T *Out, mwSize MM, mwSize NN,
                          mwSignedIndex *rc_offsets, mwSize num_neighbors) {
    mwSignedIndex M = static_cast<mwSignedIndex>(MM);
    mwSignedIndex N = static_cast<mwSignedIndex>(NN);
    mwSignedIndex *column_offset;
    mwSignedIndex *row_offset1;
    mwSignedIndex *row_offset2;
    mwSignedIndex *bit_shift1;
    mwSignedIndex *bit_shift2;
    uint32_T val;
    uint32_T shifted_val;
    
    column_offset = (mwSignedIndex *) mxCalloc(num_neighbors, sizeof(*column_offset));
    row_offset1 = (mwSignedIndex *) mxCalloc(num_neighbors, sizeof(*row_offset1));
    row_offset2 = (mwSignedIndex *) mxCalloc(num_neighbors, sizeof(*row_offset2));
    bit_shift1 = (mwSignedIndex *) mxCalloc(num_neighbors, sizeof(*bit_shift1));
    bit_shift2 = (mwSignedIndex *) mxCalloc(num_neighbors, sizeof(*bit_shift2));
    
    for (mwSize k = 0; k < num_neighbors; k++) {
        column_offset[k] = rc_offsets[k + num_neighbors];
        mwSignedIndex r = rc_offsets[k];
        row_offset1[k] = (mwSignedIndex) floor((double) r / BITS_PER_WORD);
        row_offset2[k] = row_offset1[k] + 1;
        bit_shift1[k] = r - BITS_PER_WORD*row_offset1[k];
        bit_shift2[k] = BITS_PER_WORD - bit_shift1[k];
    }

    for (mwSignedIndex c = 0; c < N; c++) {
        for (mwSignedIndex r = 0; r < M; r++) {
            val = *In++;
            if (val != 0) {
                for (mwSize k = 0; k < num_neighbors; k++) {
                    mwSignedIndex cc = c + column_offset[k];
                    if ((cc >= 0) && (cc < N)) {
                        mwSignedIndex rr = r + row_offset1[k];
                        if ((rr >= 0) && (rr < M)) {
                            shifted_val = LEFT_SHIFT(val,bit_shift1[k]);
                            Out[cc*M + rr] |= shifted_val;
                        }
                        rr = r + row_offset2[k];
                        if ((rr >= 0) && (rr < M)) {
                            shifted_val = RIGHT_SHIFT(val,bit_shift2[k]);
                            Out[cc*M + rr] |= shifted_val;
                        }
                    }
                }
            }
        }
    }

    mxFree(column_offset);
    mxFree(row_offset1);
    mxFree(row_offset2);
    mxFree(bit_shift1);
    mxFree(bit_shift2);
}


/*
 * erode_packed_uint32
 * Packed binary erosion
 *
 * Inputs
 * ======
 * In            - pointer to first element of input array
 * M             - number of rows of packed input array
 * N             - number of columns of packed input array
 * rc_offsets    - Row-column offset locations corresponding to
 *                 each element of the structuring element; storage
 *                 order is that for a MATLAB array, num_neighbors-by-2
 * num_neighbors - number of neighbors in the structuring element
 * unpacked_M    - number of rows of unpacked input array
 *
 * Output
 * ======
 * Out           - pointer to first element of output array
 */
void erode_packed_uint32(uint32_T *In, uint32_T *Out, mwSize MM, mwSize NN,
                         mwSignedIndex *rc_offsets, mwSize num_neighbors, mwSize unpacked_M) {
    mwSignedIndex M = static_cast<mwSignedIndex>(MM);
    mwSignedIndex N = static_cast<mwSignedIndex>(NN);

    mwSignedIndex *column_offset;
    mwSignedIndex *row_offset1;
    mwSignedIndex *row_offset2;
    mwSignedIndex *bit_shift1;
    mwSignedIndex *bit_shift2;
    uint32_T val;
    uint32_T shifted_val;
    uint32_T last_row_mask;
    mwSize num_real_bits_in_last_row;
    
    column_offset = (mwSignedIndex *) mxCalloc(num_neighbors, sizeof(*column_offset));
    row_offset1 = (mwSignedIndex *) mxCalloc(num_neighbors, sizeof(*row_offset1));
    row_offset2 = (mwSignedIndex *) mxCalloc(num_neighbors, sizeof(*row_offset2));
    bit_shift1 = (mwSignedIndex *) mxCalloc(num_neighbors, sizeof(*bit_shift1));
    bit_shift2 = (mwSignedIndex *) mxCalloc(num_neighbors, sizeof(*bit_shift2));

    for (mwSize k = 0; k < 2*num_neighbors; k++)
        rc_offsets[k] = -rc_offsets[k];
    
    for (mwSize k = 0; k < num_neighbors; k++) {
        column_offset[k] = (mwSignedIndex) rc_offsets[k + num_neighbors];
        mwSignedIndex r = (mwSignedIndex) rc_offsets[k];
        row_offset1[k] = (mwSignedIndex) floor((double) r / BITS_PER_WORD);
        row_offset2[k] = row_offset1[k] + 1;
        bit_shift1[k] = r - BITS_PER_WORD*row_offset1[k];
        bit_shift2[k] = BITS_PER_WORD - bit_shift1[k];
    }

    num_real_bits_in_last_row = unpacked_M % BITS_PER_WORD;
    if (num_real_bits_in_last_row == 0)
        num_real_bits_in_last_row = BITS_PER_WORD;

    last_row_mask = 0;
    for (mwSize k = 0; k < num_real_bits_in_last_row; k++)
        last_row_mask |= 1 << k;

    for (mwSignedIndex c = 0; c < N; c++) {
        for (mwSignedIndex r = 0; r < M; r++) {
            val = ~(*In++);
            if (r == (M - 1))
                val = val & last_row_mask;
            if (val != 0) {
                for (mwSize k = 0; k < num_neighbors; k++) {
                    mwSignedIndex cc = c + column_offset[k];
                    if ((cc >= 0) && (cc < N)) {
                        mwSignedIndex rr = r + row_offset1[k];
                        if ((rr >= 0) && (rr < M)) {
                            shifted_val = LEFT_SHIFT(val,bit_shift1[k]);
                            Out[cc*M + rr] |= shifted_val;
                        }
                        rr = r + row_offset2[k];
                        if ((rr >= 0) && (rr < M)) {
                            shifted_val = RIGHT_SHIFT(val,bit_shift2[k]);
                            Out[cc*M + rr] |= shifted_val;
                        }
                    }
                }
            }
        }
    }

    for (mwSignedIndex k = 0; k < M*N; k++)
        Out[k] = ~Out[k];

    /*
     * Mask out extraneous bits in the last row.
     */
    for (mwSignedIndex c = 0; c < N; c++)
        Out[c*M + M - 1] &= last_row_mask;

    mxFree(column_offset);
    mxFree(row_offset1);
    mxFree(row_offset2);
    mxFree(bit_shift1);
    mxFree(bit_shift2);
}
