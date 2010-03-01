/*
 * Copyright 1993-2003 The MathWorks, Inc.
 * $Revision: 1.1.8.2 $  $Date: 2007/06/04 21:10:08 $
 */

#ifndef NEIGHBORHOOD_H
#define NEIGHBORHOOD_H

#include "mex.h"
#include "mwsize.h"

#define NH_USE_ALL       0
#define NH_SKIP_TRAILING 1
#define NH_SKIP_LEADING  2
#define NH_SKIP_CENTER   4

#define NH_CENTER_MIDDLE_ROUNDUP   0
#define NH_CENTER_MIDDLE_ROUNDDOWN 2
#define NH_CENTER_UL               4
#define NH_CENTER_LR               8

typedef struct Neighborhood_tag {
    /* 
     * Conceptually, array_coords is a num_neighbors-by-num_dims array 
     * containing relative offsets.  num_dims is the array dimension
     * that varies quickest.
     */
    mwSignedIndex *array_coords;  
    mwSize num_neighbors;
    mwSize num_dims;
} *Neighborhood_T;

typedef struct NeighborhoodWalker_tag {
    /* 
     * Conceptually, array_coords is a num_neighbors-by-num_dims array 
     * containing relative offsets.
     */
    mwSignedIndex *array_coords;

    /*
     * neighbor_offsets is an an array containing linear neighbor
     * offsets, computed with respect to a given image size.
     */
    mwSignedIndex *neighbor_offsets;
    mwSize *image_size;

    /*
     * Contains the array coordinates for the image pixel whose
     * neighborhood we are about to walk.
     */
    mwSize *center_coords;

    /*
     * Contains the cumulative product of the image_size array;
     * used in the sub_to_ind and ind_to_sub calculations.
     */
    mwSize *cumprod;

    /*
     * Linear index of image pixel whose neighborhood we are about to
     * walk.
     */
    mwSize pixel_offset;

    /*
     * Used to filter out certain neighbors in a neighborhood walk.
     */
    bool    *use;

    /*
     * Index of the next neighbor in the walk.
     */
    mwSignedIndex next_neighbor_idx;
    mwSize num_neighbors;
    mwSize num_dims;

    /*
     * Flag to make sure walker isn't used until nhSetWalkerLocation is called.
     */
    bool ready_for_use;

} *NeighborhoodWalker_T;

Neighborhood_T nhMakeNeighborhood(const mxArray *D,int center_location);
NeighborhoodWalker_T nhMakeNeighborhoodWalker(Neighborhood_T nhood,
                                              const mwSize *input_size,
                                              mwSize input_dims,
                                              unsigned int flags);
Neighborhood_T nhMakeDefaultConnectivityNeighborhood(mwSize num_dims);
void nhReflectNeighborhood(Neighborhood_T nhood);
void nhDestroyNeighborhoodWalker(NeighborhoodWalker_T walker);
void nhDestroyNeighborhood(Neighborhood_T nhood);
void nhSetWalkerLocation(NeighborhoodWalker_T walker, mwSize p);
bool nhGetNextInboundsNeighbor(NeighborhoodWalker_T walker, mwSize *p, mwSize *idx);
void nhCheckDomain(const mxArray *D);

extern void nhCheckConnectivityDomain(const mxArray *D,
                                      const char *function_name,
                                      const char *variable_name,
                                      int argument_position);

mwSize sub_to_ind(mwSize *coords, mwSize *cumprod, mwSize num_dims);
mwSignedIndex sub_to_ind_signed(mwSignedIndex *coords, mwSize *cumprod, mwSize num_dims);
void ind_to_sub(mwSize p, mwSize num_dims, mwSize *cumprod, mwSize *coords);
void ind_to_sub(mwSize p, mwSize num_dims, mwSize *cumprod, mwSignedIndex *coords);

mwSignedIndex *nhGetWalkerNeighborOffsets(NeighborhoodWalker_T walker);
mwSize num_nonzeros(const mxArray *D);
Neighborhood_T allocate_neighborhood(mwSize num_neighbors, mwSize num_dims);
mwSize nhGetNumNeighbors(NeighborhoodWalker_T walker);

/////////////////////////////////////////////////////////////////////////
//
/////////////////////////////////////////////////////////////////////////
template<typename _T1>
Neighborhood_T create_neighborhood_general_template(_T1 *pr, const mxArray *D, int center_location) {
    mwSize num_neighbors;
    mwSize num_dims;
    mwSize num_elements;
    mwSize *cumprod;
    Neighborhood_T result;
    mwSize count = 0;

    mxAssert(D != NULL, "");

    num_neighbors = num_nonzeros(D);
    num_dims = mxGetNumberOfDimensions(D);
    num_elements = mxGetNumberOfElements(D);
    const mwSize *size_D = mxGetDimensions(D);

    result = allocate_neighborhood(num_neighbors, num_dims);
    cumprod = (mwSize *)mxMalloc(num_dims * sizeof(*cumprod));
    cumprod[0] = 1;
    for (mwSize p = 1; p < num_dims; p++)
        cumprod[p] = cumprod[p-1] * size_D[p-1];
    
    for (mwSize p = 0; p < num_elements; p++) {
        if (*pr) {
            mwSignedIndex *coords = result->array_coords + count*num_dims;
            ind_to_sub(p, num_dims, cumprod, coords);
            if(center_location == NH_CENTER_MIDDLE_ROUNDDOWN) {
                /*
                 * Subtract the location of the center pixel from
                 * the neighbor coordinates.
                 */
                for (mwSize q = 0; q < num_dims; q++)
                    coords[q] -= (size_D[q] - 1) / 2;
            }
            else if(center_location == NH_CENTER_UL) {
                /*
                 * No change required for center in Upper Left
                 * assuming that ind_to_sub returns a zero based
                 * subscript with the top in the upper left
                 */
            }
            else if(center_location == NH_CENTER_LR) {
                for (mwSize q = 0; q < num_dims; q++)
                    coords[q] -= (size_D[q] - 1);
            }
            if(center_location == NH_CENTER_MIDDLE_ROUNDUP) {
                /*
                 * Subtract the location of the center pixel from
                 * the neighbor coordinates.
                 */
                for (mwSize q = 0; q < num_dims; q++)
                    coords[q] -= (size_D[q] - 1) / 2 + ((size_D[q] - 1) % 2 ? 1:0);
            }
            count++;
        }
        pr++;
    }
    
    mxFree(cumprod);
    return result;
}

#endif
