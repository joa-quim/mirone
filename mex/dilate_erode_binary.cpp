/* 
 * Copyright 1993-2007 The MathWorks, Inc.
 * $Revision: 1.1.6.2 $  $Date: 2007/06/04 21:09:17 $
 * 
 * Core algorithms for binary dilation and erosion.
 */

#include "morphmex.h"

/*
 * dilate_logical
 * Perform binary dilation on a logical array.
 *
 * Inputs
 * ======
 * In            - pointer to first element of input array
 * num_elements  - number of elements in input and output arrays
 * walker        - neighborhood walker representing the structuring element
 *
 * Output
 * ======
 * Out           - pointer to first element of output array
 */
void dilate_logical(mxLogical *In, mxLogical *Out, mwSize num_elements, NeighborhoodWalker_T walker) {
    for (mwSize p = 0; p < num_elements; p++) {
        if (In[p]) {
            mwSize q;

            nhSetWalkerLocation(walker, p);
            while (nhGetNextInboundsNeighbor(walker, &q, NULL))
                Out[q] = 1;
        }
    }
}


/*
 * erode_logical
 * Perform binary erosion on a logical array.
 *
 * Inputs
 * ======
 * In            - pointer to first element of input array
 * num_elements  - number of elements in input and output arrays
 * walker        - neighborhood walker representing the structuring element
 *
 * Output
 * ======
 * Out           - pointer to first element of output array
 */
void erode_logical(mxLogical *In, mxLogical *Out, mwSize num_elements, NeighborhoodWalker_T walker) {
    for (mwSize p = 0; p < num_elements; p++) {
        mwSize q;

        Out[p] = 1;
        nhSetWalkerLocation(walker, p);
        while (nhGetNextInboundsNeighbor(walker, &q, NULL)) {
            if (In[q] == 0) {
                Out[p] = 0;
                break;
            }
        }
    }
}
