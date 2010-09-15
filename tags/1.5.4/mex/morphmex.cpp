/*
 * Copyright 1993-2007 The MathWorks, Inc.
 * $Revision: 1.1.6.6 $  $Date: 2007/06/04 21:10:04 $
 *
 * MORPHMEX(MEX_METHOD, B, NHOOD, HEIGHT, UNPACKED_M)
 *
 * MEX_METHOD:
 *            'dilate_binary'
 *            'erode_binary'
 *            'dilate_binary_packed'
 *            'erode_binary_packed'
 *            'dilate_gray_flat'
 *            'erode_gray_flat'
 *            'dilate_gray_nonflat'
 *            'erode_gray_nonflat'
 *
 * B: input image array
 *
 * NHOOD: neighborhood of structuring element; N-D double array of 0s and 1s
 *
 * HEIGHT: double array of heights; same size as NHOOD
 *
 * UNPACKED_M: row size of original input image before it was packed;
 *             only used for the erode_binary_packed method; otherwise
 *             it is ignored.
 */

#include <string.h>
#include <math.h>
#include "morphmex.h"
#include "erode_linear.h"
#include "dilate_linear.h"
#include "vectors.h"

/*
 * When performing flat gray scale erosion or dilation with a line segment,
 * use a specialized algorithm if the line segment is at least this long.
 */
#define MIN_LINE_SEGMENT_LENGTH 4

void dilate_binary(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    Neighborhood_T nhood;
    NeighborhoodWalker_T walker;
    const mwSize *input_size;
    mwSize input_dims;
    const mxArray *input_image = prhs[1];
    const mxArray *input_nhood = prhs[2];
    mxArray *output_image;

    // unused parameters
    (void) nlhs;
    (void) nrhs;

    if (!mxIsLogical(input_image))
        mexErrMsgIdAndTxt("Images:morphmex:"
                          "inputImageMustBeLogicalForDilateMethod",
                          "%s","Input image must be logical for"
                          " dilate_binary method.");
    
    input_dims = mxGetNumberOfDimensions(input_image);
    input_size = mxGetDimensions(input_image);

    nhood = nhMakeNeighborhood(input_nhood,NH_CENTER_MIDDLE_ROUNDDOWN);
    walker = nhMakeNeighborhoodWalker(nhood, input_size, input_dims, 0U);

    output_image = mxCreateLogicalArray(input_dims, input_size);

    dilate_logical((mxLogical *) mxGetData(input_image),
                   (mxLogical *) mxGetData(output_image),
                   mxGetNumberOfElements(input_image),
                   walker);
    
    nhDestroyNeighborhood(nhood);
    nhDestroyNeighborhoodWalker(walker);

    plhs[0] = output_image;
}

void erode_binary(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    Neighborhood_T nhood;
    NeighborhoodWalker_T walker;
    const mwSize *input_size;
    mwSize input_dims;
    const mxArray *input_image = prhs[1];
    const mxArray *input_nhood = prhs[2];
    mxArray *output_image;

    // unused parameters
    (void) nlhs;
    (void) nrhs;

    if (!mxIsLogical(input_image))
        mexErrMsgIdAndTxt("Images:morphmex:"
                          "inputImageMustBeLogicalForErodeMethod",
                          "%s","Input image must be logical for "
                          "erode_binary method.");

    input_dims = mxGetNumberOfDimensions(input_image);
    input_size = mxGetDimensions(input_image);

    nhood = nhMakeNeighborhood(input_nhood,NH_CENTER_MIDDLE_ROUNDDOWN);
    walker = nhMakeNeighborhoodWalker(nhood, input_size, input_dims, 0U);

    output_image = mxCreateLogicalArray(input_dims, input_size);
    
    erode_logical((mxLogical *) mxGetData(input_image),
                  (mxLogical *) mxGetData(output_image),
                  mxGetNumberOfElements(input_image),
                  walker);
    
    nhDestroyNeighborhood(nhood);
    nhDestroyNeighborhoodWalker(walker);

    plhs[0] = output_image;
}

/*
 * is_linear_strel() determines if input neighborhood is a contiguous line
 * of pixels oriented along a single dimension of the array.
 *
 * Input parameter input_nhood is a logical array specifying the strel neighborhood.
 *
 * Return value is true if input neighborhood is a linear strel oriented along one of the
 *     array dimensions.  Otherwise returns false.  Returns false for an empty input.
 *
 * Output parameter dim is the dimension along which the linear strel is oriented. It is
 *     zero-based.
 *
 * Output parameter lambda is the length (number of nonzero pixels) of the linear strel.
 *
 * Output parameter origin is the origin offset of the linear strel, with respect to
 *     first nonzero pixel location in the input neighborhood.
 *
 * Output parameters are only set if the return value is true.
 */
bool is_linear_strel(const mxArray *input_nhood, mwSize *dim, mwSize *lambda, mwSignedIndex *origin) {
    bool      result                      = false;
    mwSize    numdims                     = mxGetNumberOfDimensions(input_nhood);
    const mwSize *size                    = mxGetDimensions(input_nhood);
    mwSize    numel                       = mxGetNumberOfElements(input_nhood);
    double    *pr                         = (double *) mxGetData(input_nhood);
    mwSize    num_nonsingleton_dimensions = 0;
    mwSize    nonsingleton_dimension      = 0;
    mwSize    first_nonzero_value         = 0;
    mwSize    last_nonzero_value          = 0;
    mwSize    k;

    mxAssert(mxIsDouble(input_nhood), "Input neighborhood expected to be double.");

    /*
     * Check to see if there is only one nonsingleton dimension.
     * Remember which dimension is the nonsingleton dimension.
     */
    if (! mxIsEmpty(input_nhood)) {
        for (k = 0; k < numdims; k++) {
            if (size[k] != 1) {
                num_nonsingleton_dimensions++;
                nonsingleton_dimension = k;
            }
        }

        if (num_nonsingleton_dimensions == 1)
            /* Not a linear strel. */
            result = true;
    }

    /*
     * Find the linear offset of the first nonzero value
     * in the neighborhood.  If there are no nonzero values
     * then it's not a linear strel.
     */
    if (result) {
        bool found = false;

        for (k = 0; k < numel; k++) {
            if (pr[k]) {
                first_nonzero_value = k;
                found = true;
                break;
            }
        }

        if (! found)
            result = false;
    }

    if (result) {
        /*
         * Find the linear offset of the last nonzero value
         * in the neighborhood.
         */
        for (mwSize k_up = 0; k_up < numel; k_up++) {
            k = numel - 1 - k_up;
            if (pr[k]) {
                last_nonzero_value = k;
                break;
            }
        }

        /*
         * If there are any zero-valued pixels between the first 
         * and the last nonzero-valued pixels, then it's not a linear
         * strel.
         */
        for (k = first_nonzero_value + 1; k < last_nonzero_value; k++) {
            if (! pr[k]) {
                result = false;
                break;
            }
        }
    }

    if (result) {
        *dim = nonsingleton_dimension;
        *lambda = last_nonzero_value - first_nonzero_value  + 1;
        *origin = ((size[nonsingleton_dimension] - 1) / 2) - first_nonzero_value;
    }

    return(result);
                
}

/*
 * dilate_gray_linear() performs grayscale flat dilation with a line segment.
 *
 * Input parameter A is the input image.
 *
 * Output parameter B is the output image.
 *
 * Input parameter dim is the dimension (zero-based) along which to perform
 *     the dilation.
 *
 * Input parameter lambda is the length (in pixels) of the line segment.
 *
 * Input parameter origin is the offset of the structuring element origin,
 *     with respect to the first pixel in the line segment.
 *
 */
void dilate_gray_linear(const mxArray *A, mxArray *B, mwSize dim, mwSize lambda, mwSignedIndex origin) {
    mwSize        numdims      = mxGetNumberOfDimensions(A);
    const mwSize  *size        = mxGetDimensions(A);
    size_t        elem_size    = mxGetElementSize(A);
    uint8_T       *in          = (uint8_T *) mxGetData(A);
    uint8_T       *out         = (uint8_T *) mxGetData(B);

    mwSize working_length;

    uint8_T *h;
    uint8_T *g;
    uint8_T *f;
    uint8_T *r;

    mwSize M1;
    mwSize M2;
    mwSize M3;

    /*
     * Compute a pseudo-size, M1-by-M2-by-M3, for the input array.  M2 is the
     * size of the array along the specified dimension.  M1 is the product
     * of the sizes of the array along all dimensions less than the specified
     * dimension.  M3 is the product of the sizes of the array along all
     * dimensions greater than the specified dimension.
     */
    dimReshape(numdims, size, dim, &M1, &M2, &M3);

    /* Working vectors must be a multiple of the structuring element length. */
    working_length = ((M2%lambda) == 0) ? M2 : (lambda * ((M2 / lambda) + 1));

    h = (uint8_T *) mxMalloc(working_length * elem_size);
    g = (uint8_T *) mxMalloc(working_length * elem_size);
    f = (uint8_T *) mxMalloc(working_length * elem_size);
    r = (uint8_T *) mxMalloc(working_length * elem_size);

    for (mwSize k3 = 0; k3 < M3; k3++) {
        for (mwSize k1 = 0; k1 < M1; k1++) {
            /*
             * In this call to copyVector, M1 is the input vector stride,
             * 1 is the output vector stride, and M2 is the number of 
             * elements to copy.
             */
            copyVector(in, f, elem_size, M1, 1, M2);
            
            switch (mxGetClassID(A)) {
              case mxDOUBLE_CLASS:
                dilateWithLine((double *) f, (double *) g, (double *) h, (double *) r, 
                              -((double) mxGetInf()), lambda, origin, M2, working_length);
                break;

              case mxSINGLE_CLASS:
                dilateWithLine((float *) f, (float *) g, (float *) h, (float *) r,
                              -((float) mxGetInf()), lambda, origin, M2, working_length);
                break;

              case mxUINT8_CLASS:
                dilateWithLine((uint8_T *) f, (uint8_T *) g, (uint8_T *) h, (uint8_T *) r,
                              MIN_uint8_T, lambda, origin, M2, working_length);
                break;

              case mxUINT16_CLASS:
                dilateWithLine((uint16_T *) f, (uint16_T *) g, (uint16_T *) h, (uint16_T *) r,
                              MIN_uint16_T, lambda, origin, M2, working_length);
                break;

              case mxUINT32_CLASS:
                dilateWithLine((uint32_T *) f, (uint32_T *) g, (uint32_T *) h, (uint32_T *) r,
                              MIN_uint32_T, lambda, origin, M2, working_length);
                break;

              case mxUINT64_CLASS:
                dilateWithLine((uint64_T *) f, (uint64_T *) g, (uint64_T *) h, (uint64_T *) r,
                              MIN_uint64_T, lambda, origin, M2, working_length);
                break;

              case mxINT8_CLASS:
                dilateWithLine((int8_T *) f, (int8_T *) g, (int8_T *) h, (int8_T *) r,
                              MIN_int8_T, lambda, origin, M2, working_length);
                break;

              case mxINT16_CLASS:
                dilateWithLine((int16_T *) f, (int16_T *) g, (int16_T *) h, (int16_T *) r,
                              MIN_int16_T, lambda, origin, M2, working_length);
                break;

              case mxINT32_CLASS:
                dilateWithLine((int32_T *) f, (int32_T *) g, (int32_T *) h, (int32_T *) r,
                              MIN_int32_T, lambda, origin, M2, working_length);
                break;

              case mxINT64_CLASS:
                dilateWithLine((int64_T *) f, (int64_T *) g, (int64_T *) h, (int64_T *) r,
                              MIN_int64_T, lambda, origin, M2, working_length);
                break;

              default:
                mexErrMsgIdAndTxt("Images:morphmex:unexpectedCase", "Reached unexpected switch case.");
                break;

            }

            /*
             * In this call to copyVector, 1 is the input vector stride,
             * M1 is the output vector stride, and M2 is the number of 
             * elements to copy.
             */
            copyVector(r, out, elem_size, 1, M1, M2);

            /*
             * When moving along the first pseudo-dimension (M1), vectors are separated
             * by one element.
             */
            in  += elem_size;
            out += elem_size;
        }

        /*
         * When moving along the third pseudo-dimension (M3), vectors are separated
         * by M1*M2 elements.  Since we just moved once through the first dimension,
         * we move ahead by (M1*M2)-M1 elements.
         */
        in  += (M1*M2 - M1) * elem_size;
        out += (M1*M2 - M1) * elem_size;
    }

    mxFree(f);
    mxFree(g);
    mxFree(h);
    mxFree(r);
}

void dilate_gray_flat(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    Neighborhood_T nhood;
    NeighborhoodWalker_T walker;
    const mwSize *input_size;
    mwSize input_dims;
    int input_class;
    const mxArray *input_image = prhs[1];
    const mxArray *input_nhood = prhs[2];
    mxArray *output_image;
    void *In;
    void *Out;
    mwSize num_elements;
    mwSize dim;           /* Orientation dimension for linear strel.                          */
    mwSize lambda;        /* Length of linear strel.                                          */
    mwSignedIndex origin; /* Origin of linear strel, with respect to left-most nonzero pixel. */

    // unused parameters
    (void) nlhs;
    (void) nrhs;

    input_size = mxGetDimensions(input_image);
    input_dims = mxGetNumberOfDimensions(input_image);
    input_class = mxGetClassID(input_image);
    output_image = mxCreateNumericArray(input_dims, 
                                        input_size,
                                        mxGetClassID(input_image),
                                        mxREAL);
    plhs[0] = output_image;

    if (is_linear_strel(input_nhood, &dim, &lambda, &origin) && (lambda >= MIN_LINE_SEGMENT_LENGTH))
    {
        dilate_gray_linear(input_image, output_image, dim, lambda, origin);
        return;
    }

    nhood = nhMakeNeighborhood(input_nhood,NH_CENTER_MIDDLE_ROUNDDOWN);
    nhReflectNeighborhood(nhood);
    walker = nhMakeNeighborhoodWalker(nhood, input_size, input_dims, 0U);

    num_elements = mxGetNumberOfElements(input_image);
    In = mxGetData(input_image);
    Out = mxGetData(output_image);

    switch (input_class) {
    case mxUINT8_CLASS:
        dilateGrayFlat((uint8_T *)In, (uint8_T *)Out, num_elements, walker);
        break;
        
    case mxUINT16_CLASS:
        dilateGrayFlat((uint16_T *)In, (uint16_T *)Out, num_elements, walker);
        break;
        
    case mxUINT32_CLASS:
        dilateGrayFlat((uint32_T *)In, (uint32_T *)Out, num_elements, walker);
        break;
        
    case mxINT8_CLASS:
        dilateGrayFlat((int8_T *)In, (int8_T *)Out, num_elements, walker);
        break;
        
    case mxINT16_CLASS:
        dilateGrayFlat((int16_T *)In, (int16_T *)Out, num_elements, walker);
        break;
        
    case mxINT32_CLASS:
        dilateGrayFlat((int32_T *)In, (int32_T *)Out, num_elements, walker);
        break;
        
    case mxSINGLE_CLASS:
        dilateGrayFlat((float *)In, (float *)Out, num_elements, walker);
        break;
        
    case mxDOUBLE_CLASS:
        dilateGrayFlat((double *)In, (double *)Out, num_elements, walker);
        break;
        
    default:
        mexErrMsgIdAndTxt("Images:morphmex:internalInvalidClass",
                          "%s","Internal problem: invalid input image class.");
    }

    nhDestroyNeighborhood(nhood);
    nhDestroyNeighborhoodWalker(walker);

}

void erode_gray_linear(const mxArray *A, mxArray *B, mwSize dim, mwSize lambda, mwSignedIndex origin) {
    mwSize        numdims      = mxGetNumberOfDimensions(A);
    const mwSize  *size        = mxGetDimensions(A);
    size_t        elem_size    = mxGetElementSize(A);
    uint8_T       *in          = (uint8_T *) mxGetData(A);
    uint8_T       *out         = (uint8_T *) mxGetData(B);

    mwSize working_length;

    uint8_T *h;
    uint8_T *g;
    uint8_T *f;
    uint8_T *r;

    mwSize M1;
    mwSize M2;
    mwSize M3;

    dimReshape(numdims, size, dim, &M1, &M2, &M3);

    /* Working vectors must be a multiple of the structuring element length. */
    working_length = ((M2%lambda) == 0) ? M2 : (lambda * ((M2 / lambda) + 1));

    h = (uint8_T *) mxMalloc(working_length * elem_size);
    g = (uint8_T *) mxMalloc(working_length * elem_size);
    f = (uint8_T *) mxMalloc(working_length * elem_size);
    r = (uint8_T *) mxMalloc(working_length * elem_size);

    for (mwSize k3 = 0; k3 < M3; k3++) {
        for (mwSize k1 = 0; k1 < M1; k1++) {
            copyVector(in, f, elem_size, M1, 1, M2);
            
            switch (mxGetClassID(A)) {
              case mxDOUBLE_CLASS:
                erodeWithLine((double *) f, (double *) g, (double *) h, (double *) r, 
                              (double) mxGetInf(), lambda, origin, M2, working_length);
                break;

              case mxSINGLE_CLASS:
                erodeWithLine((float *) f, (float *) g, (float *) h, (float *) r,
                              (float) mxGetInf(), lambda, origin, M2, working_length);
                break;

              case mxUINT8_CLASS:
                erodeWithLine((uint8_T *) f, (uint8_T *) g, (uint8_T *) h, (uint8_T *) r,
                              MAX_uint8_T, lambda, origin, M2, working_length);
                break;

              case mxUINT16_CLASS:
                erodeWithLine((uint16_T *) f, (uint16_T *) g, (uint16_T *) h, (uint16_T *) r,
                              MAX_uint16_T, lambda, origin, M2, working_length);
                break;

              case mxUINT32_CLASS:
                erodeWithLine((uint32_T *) f, (uint32_T *) g, (uint32_T *) h, (uint32_T *) r,
                              MAX_uint32_T, lambda, origin, M2, working_length);
                break;

              case mxUINT64_CLASS:
                erodeWithLine((uint64_T *) f, (uint64_T *) g, (uint64_T *) h, (uint64_T *) r,
                              MAX_uint64_T, lambda, origin, M2, working_length);
                break;

              case mxINT8_CLASS:
                erodeWithLine((int8_T *) f, (int8_T *) g, (int8_T *) h, (int8_T *) r,
                              MAX_int8_T, lambda, origin, M2, working_length);
                break;

              case mxINT16_CLASS:
                erodeWithLine((int16_T *) f, (int16_T *) g, (int16_T *) h, (int16_T *) r,
                              MAX_int16_T, lambda, origin, M2, working_length);
                break;

              case mxINT32_CLASS:
                erodeWithLine((int32_T *) f, (int32_T *) g, (int32_T *) h, (int32_T *) r,
                              MAX_int32_T, lambda, origin, M2, working_length);
                break;

              case mxINT64_CLASS:
                erodeWithLine((int64_T *) f, (int64_T *) g, (int64_T *) h, (int64_T *) r,
                              MAX_int64_T, lambda, origin, M2, working_length);
                break;

              default:
                mexErrMsgIdAndTxt("Images:morphmex:unexpectedCase", "Reached unexpected switch case.");
                break;

            }

            copyVector(r, out, elem_size, 1, M1, M2);

            /*
             * When moving along the first pseudo-dimension (M1), vectors are separated
             * by one element.
             */
            in  += elem_size;
            out += elem_size;
        }

        /*
         * When moving along the third pseudo-dimension (M3), vectors are separated
         * by M1*M2 elements.  Since we just moved once through the first dimension,
         * we move ahead by (M1*M2)-M1 elements.
         */
        in  += (M1*M2 - M1) * elem_size;
        out += (M1*M2 - M1) * elem_size;
    }

    mxFree(f);
    mxFree(g);
    mxFree(h);
    mxFree(r);
}

void erode_gray_flat(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    Neighborhood_T nhood;
    NeighborhoodWalker_T walker;
    const mwSize *input_size;
    mwSize input_dims;
    int input_class;
    const mxArray *input_image = prhs[1];
    const mxArray *input_nhood = prhs[2];
    mxArray *output_image;
    void *In;
    void *Out;
    mwSize num_elements;

    mwSize dim;           /* Orientation dimension for linear strel.                          */
    mwSize lambda;        /* Length of linear strel.                                          */
    mwSignedIndex origin; /* Origin of linear strel, with respect to left-most nonzero pixel. */

    // unused parameters
    (void) nlhs;
    (void) nrhs;

    input_size = mxGetDimensions(input_image);
    input_dims = mxGetNumberOfDimensions(input_image);
    input_class = mxGetClassID(input_image);
    output_image = mxCreateNumericArray(input_dims, 
                                        input_size,
                                        mxGetClassID(input_image),
                                        mxREAL);
    plhs[0] = output_image;
    
    if (is_linear_strel(input_nhood, &dim, &lambda, &origin) && (lambda >= MIN_LINE_SEGMENT_LENGTH))
    {
        erode_gray_linear(input_image, output_image, dim, lambda, origin);
        return;
    }

    nhood = nhMakeNeighborhood(input_nhood,NH_CENTER_MIDDLE_ROUNDDOWN);
    walker = nhMakeNeighborhoodWalker(nhood,
                                      input_size,
                                      input_dims,
                                      0U);

    num_elements = mxGetNumberOfElements(input_image);
    In = mxGetData(input_image);
    Out = mxGetData(output_image);

    switch (input_class) {
      case mxUINT8_CLASS:
        erodeGrayFlat((uint8_T *)In, (uint8_T *)Out, num_elements, walker);
        break;
        
      case mxUINT16_CLASS:
        erodeGrayFlat((uint16_T *)In, (uint16_T *)Out, num_elements, walker);
        break;
        
      case mxUINT32_CLASS:
        erodeGrayFlat((uint32_T *)In, (uint32_T *)Out, num_elements, walker);
        break;
        
      case mxINT8_CLASS:
        erodeGrayFlat((int8_T *)In, (int8_T *)Out, num_elements, walker);
        break;
        
      case mxINT16_CLASS:
        erodeGrayFlat((int16_T *)In, (int16_T *)Out, num_elements, walker);
        break;
        
      case mxINT32_CLASS:
        erodeGrayFlat((int32_T *)In, (int32_T *)Out, num_elements, walker);
        break;
        
      case mxSINGLE_CLASS:
        erodeGrayFlat((float *)In, (float *)Out, num_elements, walker);
        break;
        
      case mxDOUBLE_CLASS:
        erodeGrayFlat((double *)In, (double *)Out, num_elements, walker);
        break;
        
      default:
        mexErrMsgIdAndTxt("Images:morphmex:internalInvalidClass",
                          "%s","Internal problem: invalid input image class.");
    }

    nhDestroyNeighborhood(nhood);
    nhDestroyNeighborhoodWalker(walker);

}

double *get_heights(const mxArray *input_nhood, const mxArray *input_height) {
    mwSize num_neighbors;
    mwSize num_elements;
    double *heights;
    double *pr;
    double *input_nhood_pr;
    double *input_height_pr;
    
    num_elements = mxGetNumberOfElements(input_nhood);
    input_nhood_pr = (double *) mxGetData(input_nhood);
    num_neighbors = 0;
    for (mwSize k = 0; k < num_elements; k++) {
        if (input_nhood_pr[k] != 0)
            num_neighbors++;
    }
    
    heights = (double *) mxCalloc(num_neighbors, sizeof(double));
    input_height_pr = (double *) mxGetData(input_height);
    pr = heights;
    for (mwSize k = 0; k < num_elements; k++) {
        if (input_nhood_pr[k] != 0.0) {
            *pr = input_height_pr[k];
            pr++;
        }
    }

    return heights;
}

void dilate_gray_nonflat(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    Neighborhood_T nhood;
    NeighborhoodWalker_T walker;
    const mwSize *input_size;
    mwSize input_dims;
    int input_class;
    const mxArray *input_image = prhs[1];
    const mxArray *input_nhood = prhs[2];
    const mxArray *input_height = prhs[3];
    mxArray *output_image;
    void *In;
    void *Out;
    mwSize num_elements;
    double *heights;

    // unused parameters
    (void) nlhs;
    (void) nrhs;

    input_size = mxGetDimensions(input_image);
    input_dims = mxGetNumberOfDimensions(input_image);
    input_class = mxGetClassID(input_image);
    output_image = mxCreateNumericArray(input_dims, input_size, mxGetClassID(input_image), mxREAL);
    
    nhood = nhMakeNeighborhood(input_nhood,NH_CENTER_MIDDLE_ROUNDDOWN);
    nhReflectNeighborhood(nhood);
    walker = nhMakeNeighborhoodWalker(nhood, input_size, input_dims, 0U);

    num_elements = mxGetNumberOfElements(input_image);
    In = mxGetData(input_image);
    Out = mxGetData(output_image);

    heights = get_heights(input_nhood, input_height);

    switch (input_class) {
    case mxUINT8_CLASS:
        dilate_gray_nonflat_uint8((uint8_T *)In, (uint8_T *)Out, num_elements, walker, heights);
        break;
        
    case mxUINT16_CLASS:
        dilate_gray_nonflat_uint16((uint16_T *)In, (uint16_T *)Out, num_elements, walker, heights);
        break;
        
    case mxUINT32_CLASS:
        dilate_gray_nonflat_uint32((uint32_T *)In, (uint32_T *)Out, num_elements, walker, heights);
        break;
        
    case mxINT8_CLASS:
        dilate_gray_nonflat_int8((int8_T *)In, (int8_T *)Out, num_elements, walker, heights);
        break;
        
    case mxINT16_CLASS:
        dilate_gray_nonflat_int16((int16_T *)In, (int16_T *)Out, num_elements, walker, heights);
        break;
        
    case mxINT32_CLASS:
        dilate_gray_nonflat_int32((int32_T *)In, (int32_T *)Out, num_elements, walker, heights);
        break;
        
    case mxSINGLE_CLASS:
        dilate_gray_nonflat_single((float *)In, (float *)Out, num_elements, walker, heights);
        break;
        
    case mxDOUBLE_CLASS:
        dilate_gray_nonflat_double((double *)In, (double *)Out, num_elements, walker, heights);
        break;
        
    default:
        mexErrMsgIdAndTxt("Images:morphmex:internalInvalidClass",
                          "%s","Internal problem: invalid input image class.");
    }

    nhDestroyNeighborhood(nhood);
    nhDestroyNeighborhoodWalker(walker);
    mxFree(heights);

    plhs[0] = output_image;
}

void erode_gray_nonflat(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    Neighborhood_T nhood;
    NeighborhoodWalker_T walker;
    const mwSize *input_size;
    mwSize input_dims;
    int input_class;
    const mxArray *input_image = prhs[1];
    const mxArray *input_nhood = prhs[2];
    const mxArray *input_height = prhs[3];
    mxArray *output_image;
    void *In;
    void *Out;
    mwSize num_elements;
    double *heights;

    // unused parameters
    (void) nlhs;
    (void) nrhs;

    input_size = mxGetDimensions(input_image);
    input_dims = mxGetNumberOfDimensions(input_image);
    input_class = mxGetClassID(input_image);
    output_image = mxCreateNumericArray(input_dims, input_size, mxGetClassID(input_image), mxREAL);
    
    nhood = nhMakeNeighborhood(input_nhood,NH_CENTER_MIDDLE_ROUNDDOWN);
    walker = nhMakeNeighborhoodWalker(nhood, input_size, input_dims, 0U);

    num_elements = mxGetNumberOfElements(input_image);
    In = mxGetData(input_image);
    Out = mxGetData(output_image);
    
    heights = get_heights(input_nhood, input_height);

    switch (input_class)
    {
    case mxUINT8_CLASS:
        erode_gray_nonflat_uint8((uint8_T *)In, (uint8_T *)Out, num_elements, walker, heights);
        break;
        
    case mxUINT16_CLASS:
        erode_gray_nonflat_uint16((uint16_T *)In, (uint16_T *)Out, num_elements, walker, heights);
        break;
        
    case mxUINT32_CLASS:
        erode_gray_nonflat_uint32((uint32_T *)In, (uint32_T *)Out, num_elements, walker, heights);
        break;
        
    case mxINT8_CLASS:
        erode_gray_nonflat_int8((int8_T *)In, (int8_T *)Out, num_elements, walker, heights);
        break;
        
    case mxINT16_CLASS:
        erode_gray_nonflat_int16((int16_T *)In, (int16_T *)Out, num_elements, walker, heights);
        break;
        
    case mxINT32_CLASS:
        erode_gray_nonflat_int32((int32_T *)In, (int32_T *)Out, num_elements, walker, heights);
        break;
        
    case mxSINGLE_CLASS:
        erode_gray_nonflat_single((float *)In, (float *)Out, num_elements, walker, heights);
        break;
        
    case mxDOUBLE_CLASS:
        erode_gray_nonflat_double((double *)In, (double *)Out, num_elements, walker, heights);
        break;
        
    default:
        mexErrMsgIdAndTxt("Images:morphmex:internalInvalidClass",
                          "%s","Internal problem: invalid input image class.");
    }

    nhDestroyNeighborhood(nhood);
    nhDestroyNeighborhoodWalker(walker);
    mxFree(heights);

    plhs[0] = output_image;
}

void get_rc_offsets(const mxArray *nhood, mwSize *num_neighbors, mwSignedIndex **rc_offsets) {
    mwSize M;
    mwSize N;
    mwSize center_row;
    mwSize center_col;
    double *pr;
    mwSize num_elements;
    mwSize counter;
    mwSize c;
    mwSize r;
    
    if (mxGetNumberOfDimensions(nhood) != 2) {
        mexErrMsgIdAndTxt("Images:morphmex:nhoodMustBe2DForPackedMethods",
                          "%s","Neighborhood must be 2-D for packed methods.");
    }

    M = mxGetM(nhood);
    N = mxGetN(nhood);
    num_elements = M*N;

    center_row = (M-1)/2;
    center_col = (N-1)/2;
    *num_neighbors = 0;
    pr = (double *) mxGetData(nhood);
    
    for (mwSize k = 0; k < num_elements; k++) {
        if (pr[k] != 0.0)
            (*num_neighbors)++;
    }
    
    *rc_offsets = (mwSignedIndex *) mxCalloc(*num_neighbors * 2, sizeof(**rc_offsets));
    
    counter = 0;
    for (c = 0; c < N; c++) {
        for (r = 0; r < M; r++) {
            if (*pr != 0.0) {
                (*rc_offsets)[counter] = static_cast<mwSignedIndex>(r) - center_row;
                (*rc_offsets)[counter + *num_neighbors] = 
                    static_cast<mwSignedIndex>(c) - center_col;
                counter++;
            }
            pr++;
        }
    }
}

void dilate_packed(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    const mxArray *input_image = prhs[1];
    const mxArray *input_nhood = prhs[2];
    mxArray *output_image;
    uint32_T *In;
    uint32_T *Out;
    mwSize M;
    mwSize N;
    mwSize num_neighbors;
    mwSignedIndex *rc_offsets;

    // unused parameters
    (void) nlhs;
    (void) nrhs;

    if (!mxIsUint32(input_image)) {
        mexErrMsgIdAndTxt("Images:morphmex:"
                          "inputImageMustBeUint32ForPackedMethods",
                          "%s","Input image must be uint32 for packed methods.");
    }

    if (mxGetNumberOfDimensions(input_image) != 2) {
        mexErrMsgIdAndTxt("Images:morphmex:inputImageMustBe2DForPackedMethods",
                          "%s","Input image must be 2-D for packed methods.");
    }

    M = mxGetM(input_image);
    N = mxGetN(input_image);
    output_image = mxCreateNumericMatrix(M,N,mxUINT32_CLASS,mxREAL);

    get_rc_offsets(input_nhood, &num_neighbors, &rc_offsets);
    
    In = (uint32_T *) mxGetData(input_image);
    Out = (uint32_T *) mxGetData(output_image);

    dilate_packed_uint32(In, Out, M, N, rc_offsets, num_neighbors);
    
    mxFree(rc_offsets);

    plhs[0] = output_image;
}

void erode_packed(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    const mxArray *input_image = prhs[1];
    const mxArray *input_nhood = prhs[2];
    mxArray *output_image;
    uint32_T *In;
    uint32_T *Out;
    mwSize M;
    mwSize N;
    mwSize num_neighbors;
    mwSize unpacked_M;
    mwSignedIndex *rc_offsets;

    (void) nlhs;  // unused parameter

    if (!mxIsUint32(input_image)) {
        mexErrMsgIdAndTxt("Images:morphmex:"
                          "inputImageMustBeUint32ForPackedMethods","%s",
                          "Input image must be uint32 for packed methods.");
    }

    if (mxGetNumberOfDimensions(input_image) != 2) {
        mexErrMsgIdAndTxt("Images:morphmex:inputImageMustBe2DForPackedMethods",
                          "%s","Input image must be 2-D for packed methods.");
    }

    if (nrhs < 5) {
        mexErrMsgIdAndTxt("Images:morphmex:missingMForPackedErosion",
                          "%s","M must be provided for packed erosion.");
    }

    if (mxGetScalar(prhs[4]) < 0.0) {
        mexErrMsgIdAndTxt("Images:morphmex:inputMMustBeNonnegative",
                          "%s","M must be nonnegative.");
    }
    unpacked_M = (mwSize) mxGetScalar(prhs[4]);
    
    M = mxGetM(input_image);
    N = mxGetN(input_image);
    output_image = mxCreateNumericMatrix(M,N,mxUINT32_CLASS,mxREAL);

    get_rc_offsets(input_nhood, &num_neighbors, &rc_offsets);
    
    In = (uint32_T *) mxGetData(input_image);
    Out = (uint32_T *) mxGetData(output_image);

    erode_packed_uint32(In, Out, M, N, rc_offsets, num_neighbors, 
                        unpacked_M);
    
    mxFree(rc_offsets);

    plhs[0] = output_image;
}

void check_for_nans(const mxArray *input_image) {
    double *double_ptr;
    float *single_ptr;
    mwSize k;
    mwSize num_elements;
    
    num_elements = mxGetNumberOfElements(input_image);
    
    if (mxIsDouble(input_image)) {
        double_ptr = (double *) mxGetData(input_image);
        for (k = 0; k < num_elements; k++) {
            if (mxIsNaN(double_ptr[k])) {
                mexErrMsgIdAndTxt("Images:morphmex:expectedNonnan",
                                  "%s",
                                  "Input image may not contain NaNs.");
            }
        }
    }
    else if (mxIsSingle(input_image))
    {
        single_ptr = (float *) mxGetData(input_image);
        for (k = 0; k < num_elements; k++) {
            if (mxIsNaN(single_ptr[k])) {
                mexErrMsgIdAndTxt("Images:morphmex:expectedNonnan",
                                  "%s","Input image may not contain NaNs.");
            }
        }
    }
    else {
        mexErrMsgIdAndTxt("Images:morphmex:internalBadInputClass",
                          "%s","Internal problem: unexpected input"
                          " class in check_for_nans.");
    }
}

typedef void (matlab_fcn)(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

void check_inputs_generic(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    const mxArray *method;
    const mxArray *input_image;
    const mxArray *nhood;
    const mxArray *height;
    const mxArray *unpacked_M;
    mwSize num_nhood_dims;
    const mwSize *nhood_size;
    const mwSize *height_size;
    mwSize i;
    double scalar;

    // unused parameters
    (void) nlhs;
    (void) plhs;
    
    if (nrhs < 4) {
        mexErrMsgIdAndTxt("Images:morphmex:tooFewInputs", "%s","Not enough input arguments.");
    }
    
    method = prhs[0];
    input_image = prhs[1];
    nhood = prhs[2];
    height = prhs[3];
    
    if (!mxIsChar(method)) {
        mexErrMsgIdAndTxt("Images:morphmex:firstInputMustBeMethodString","%s",
                          "First input argument must be a method string.");
    }
    
    if (!mxIsNumeric(input_image) && !mxIsLogical(input_image)) {
        mexErrMsgIdAndTxt("Images:morphmex:inputImageMustBeNumericOrLogical",
                          "%s","Input image must be numeric or logical.");
    }
    if (mxIsSparse(input_image)) {
        mexErrMsgIdAndTxt("Images:morphmex:inputImageMustBeNonsparse",
                          "%s","Input image must not be sparse.");
    }
    if (mxIsComplex(input_image)) {
        mexErrMsgIdAndTxt("Images:morphmex:inputImageMustBeReal",
                          "%s","Input image must be real.");
    }
    
    if (mxIsDouble(input_image) || mxIsSingle(input_image)) {
        check_for_nans(input_image);
    }

    nhCheckDomain(nhood);
    
    if (mxIsSparse(height)) {
        mexErrMsgIdAndTxt("Images:morphmex:heightMustBeNonsparse",
                          "%s","Height must not be sparse.");
    }
    if (!mxIsDouble(height)) {
        mexErrMsgIdAndTxt("Images:morphmex:heightMustBeDoubleArray",
                          "%s","Height must be a double array.");
    }
    if (mxIsComplex(height)) {
        mexErrMsgIdAndTxt("Images:morphmex:heightMustBeReal",
                          "%s","Height must be real.");
    }

    num_nhood_dims = mxGetNumberOfDimensions(nhood);
    if (num_nhood_dims != mxGetNumberOfDimensions(height)) {
        mexErrMsgIdAndTxt("Images:morphmex:nhoodAndHeightMustHaveSameSize",
                          "%s",
                          "Neighborhood and height must have the same size.");
    }
    nhood_size = mxGetDimensions(nhood);
    height_size = mxGetDimensions(height);
    for (i = 0; i < num_nhood_dims; i++) {
        if (nhood_size[i] != height_size[i]) {
            mexErrMsgIdAndTxt("Images:morphmex:nhoodAndHeightMustHaveSameSize",
                              "%s","Neighborhood and height must have"
                              " the same size.");
        }
    }

    if (nrhs == 5) {
        unpacked_M = prhs[4];
        if (! mxIsDouble(unpacked_M) || (mxGetNumberOfElements(unpacked_M) != 1)) {
            mexErrMsgIdAndTxt("Images:morphmex:inputMMustBeDoubleScalar",
                              "%s","M must be a double scalar.");
        }
        scalar = mxGetScalar(unpacked_M);
        if (floor(scalar) != scalar) {
            mexErrMsgIdAndTxt("Images:morphmex:inputMMustBeInteger",
                              "%s","M must be an integer.");
        }
    }
}

extern "C"
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    struct 
    {
        char *name;
        matlab_fcn *func;
    }
    morph_methods[] = {
        {"dilate_binary",        dilate_binary},
        {"erode_binary",         erode_binary},
        {"dilate_binary_packed", dilate_packed},
        {"erode_binary_packed",  erode_packed},
        {"dilate_gray_flat",     dilate_gray_flat},
        {"erode_gray_flat",      erode_gray_flat},
        {"dilate_gray_nonflat",  dilate_gray_nonflat},
        {"erode_gray_nonflat",   erode_gray_nonflat},
        {"",                     NULL}
    };

    int i = 0;
    char *method_string;
    int method_idx = 0;

    check_inputs_generic(nlhs, plhs, nrhs, prhs);
    
    method_string = mxArrayToString(prhs[0]);
        
    while (morph_methods[i].func != NULL) {
        if (strcmp(method_string, morph_methods[i].name) == 0) {
            method_idx = i;
            break;
        }
        i++;
    }
    
    mxFree(method_string);
    
    if (method_idx >= 0) {
        (*(morph_methods[method_idx].func))(nlhs, plhs, nrhs, prhs);
    }
    else {
        mexErrMsgIdAndTxt("Images:morphmex:unknownMethodString",
                          "%s","Unknown method string.");
    }
}

