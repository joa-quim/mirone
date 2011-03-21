/*
 * Copyright 1993-2007 The MathWorks, Inc.
 * $Revision: 1.1.6.4 $  $Date: 2007/06/04 21:10:22 $
 */

/*
 * Defines the MEX function for separable resamplers.  MATLAB function
 * syntax:
 *
 *   B = resampsep( A, M, tdims_A, tdims_B, fsize_A, fsize_B, F, R );
 *
 * resampsep() is an example of a resample_fcn as decribed in
 * makeresampler.m.
 */

/* 
 * Implementation remarks:
 * -----------------------
 * 
 * Iterator example: Path through a 4 x 3 x 2 array
 * 
 *    0   4   8           12  16  20
 *    1   5   9           13  17  21
 *    2   6  10           14  18  22
 *    3   7  11           15  19  23
 * 
 * 
 * Computing array offset from permuted and restructured cumulative
 * product array
 * 
 *   The basic idea is to perform all permutations, and as many 
 *   multiplications as possible, up front, before starting the processing 
 *   loops.  So we translate array subscripts to offsets in an unusual order.
 *   The essence is to construct the cumulative product of the array 
 *   dimensions, shifted forward with a one inserted at the front, then take
 *   its dot product with the subscript array.  However, we partition this
 *   dot product into contributions from transform dimensions and 
 *   contributions from other dimensions, permute the former, and compute
 *   the dot products at different times.
 *   
 *   Output array example:
 *     
 *     size_B = [300 400   3   7  10]
 *   
 *     tdims_B = [2 1 5]
 *   
 *   Cumulative product of sizes, with 1 inserted at the front:
 *   
 *     cp = [1   300   300*400  300*400*3   300*400*3*7]
 *   
 *   Transform dimensions: Take the 2nd, 1st, and 5th elements of cp:
 *   
 *     oCon->cpTrans = [300  1  300*400*3*7]
 *     
 *   Other dimensions: Take what's left, in order:
 *   
 *     oCon->cpOther = [300*400   300*400*3]
 *     
 *   Total output offset:
 *   
 *     (sum over k)( oCon->cpTrans[k] * oTransIt->subs[k] )
 *     
 *     + (sum over j)( oCon->cpOther[j] * otherIt->subs[j] )
 *   
 *   The sums are computed in Subs2Offset(), which is very simple and efficient
 *   because of all the work done up front.
 *   
 *   In ResampSep(), the outer loop goes over the transform dimensions and
 *   the inner loop goes over the other dimensions, so we compute and save
 *   the transform dimension contributions before starting the inner loop.
 *   
 *   The input array offsets are handled using the same general idea. But the
 *   contribution of the other dimensions to the offset is re-used from the
 *   output array computation.  And the transform dimension contributions are
 *   computed in the (innermost) convolution loop.
 *
 */

#include "mex.h"
#include "mwsize.h"
#include "matrix.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "iptutil.h"
#include "resampsep.h"

/*===================== Array Count Utility ====================*/

mwSize  ResampSep::GetCount( const mxArray* A ) {
    return mxGetNumberOfElements(A);
}

/*================ Cumulative Product Utility ===================*/
/*
 * Return the cumulative products of the sizes of an mxArray,
 * with a one inserted at the front.
 */
mwSize*  ResampSep::GetCumProd( const mxArray* A ) {
    mwSize        n = mxGetNumberOfDimensions(A);
    const mwSize *d = mxGetDimensions(A);

    mwSize*        cumprod = (mwSize *)mxCalloc(n + 1, sizeof(mwSize));

    cumprod[0] = 1;
    for( mwSize k = 0; k < n; k++ )
        cumprod[k+1] = d[k] * cumprod[k];;

    return cumprod;
}

/*==================== Iterator utilities ====================== */

Iterator*  ResampSep::NewIterator( const mwSize ndims, const mwSignedIndex* size ) {
    Iterator* i = (Iterator*)mxCalloc(1, sizeof(*i));
    i->length = 1;
    i->size = (mwSignedIndex*)mxCalloc(ndims,sizeof(mwSignedIndex));
    i->subs = (mwSignedIndex*)mxCalloc(ndims,sizeof(mwSignedIndex));

    i->ndims  = ndims;
    i->offset = 0;
    for( mwSize k = 0; k < ndims; k++ ) {
        i->subs[k] = 0;
        i->size[k] = size[k];
        i->length *= size[k];
    }
    return i;
}
/*
 *
 *
 */
void  ResampSep::DestroyIterator(Iterator* i) {
    if( i != NULL ) {
        mxFree(i->size);
        mxFree(i->subs);
        mxFree(i);
    }
}
/*
 *
 *
 */
inline void ResampSep::IncrementIterator(Iterator* i) {
    bool done = false;
    mwSize k = 0;
    while( !done && k < i->ndims ) {
        (i->subs[k])++;
        if( i->subs[k] < i->size[k] )
            done = true;
        else
            i->subs[k++] = 0;
    }
    (i->offset)++;
}
/*
 *
 *
 */
inline bool ResampSep::DoneIterating(const Iterator* i) {
    return i->offset >= static_cast<mwSignedIndex>(i->length);
}
/*
 *
 *
 */
void ResampSep::ResetIterator(Iterator* i) {
    i->offset = 0;
    for( mwSize k = 0; k < i->ndims; k++ )
        i->subs[k] = 0;
}

/*====================== Config related methods ====================*/

/*
 *
 *
 */
Config* ResampSep::NewConfig( const mxArray* fsize, const mxArray* tdims )
{
    mwSize* cumprod;
    bool* isTrans;

    Config* c = (Config*)mxCalloc(1,sizeof(*c));

    c->ndims  = GetCount(fsize);
    c->nTrans = GetCount(tdims);
    c->nOther = (c->ndims) - (c->nTrans);

    c->size    = (mwSize*)mxCalloc(c->ndims, sizeof(mwSize));
    c->tdims   = (mwSize*)mxCalloc(c->nTrans,sizeof(mwSize));
    c->tsize   = (mwSignedIndex*)mxCalloc(c->nTrans,sizeof(mwSignedIndex));
    c->cpTrans = (mwSize*)mxCalloc(c->nTrans,sizeof(mwSize));
    c->osize   = (mwSignedIndex*)mxCalloc(c->nOther,sizeof(mwSignedIndex));
    c->cpOther = (mwSize*)mxCalloc(c->nOther,sizeof(mwSize));

    isTrans = (bool*)mxCalloc(c->ndims,sizeof(bool));
    cumprod = (mwSize*)mxCalloc(1 + c->ndims,sizeof(mwSize));
    cumprod[0] = 1;
    for( mwSize k = 0; k < c->ndims; k++ )
    {
        c->size[k] = static_cast<mwSize>(((mxGetPr(fsize))[k]));
        cumprod[k+1] = cumprod[k] * c->size[k];
    }

    c->tlength = 1;
    for( mwSize k = 0; k < c->nTrans; k++ )
    {
        /* Subtract one to change to zero-based indexing */
        c->tdims[k]   = -1 + static_cast<mwSize>((mxGetPr(tdims))[k]);
        c->tsize[k]   = c->size[c->tdims[k]];
        c->cpTrans[k] = cumprod[c->tdims[k]];
        isTrans[c->tdims[k]] = true;
        (c->tlength) *= (c->tsize[k]);
    }

    /*
     * c->cpTrans contains the cumulative product components corresponding
     * to the transform dimensions, listed in the same order as c->tdims.
     * isTrans is true for all transform dimensions and false for the others.
     * Now copy the remaining sizes to c->osize and the remaining
     * cumulative product components to c->cpOther.
     */

    mwSize j = 0;
    for(mwSize k = 0; k < c->ndims; k++ )
    {
        if( !isTrans[k] )
        {
            c->osize[j]   = c->size[k];
            c->cpOther[j] = cumprod[k];
            j++;
        }
    }

    mxFree(cumprod);
    mxFree(isTrans);

    return c;
}


void ResampSep::DestroyConfig(Config* c)
{
    if( c != NULL )
    {
        mxFree(c->size);
        mxFree(c->tdims);
        mxFree(c->tsize);
        mxFree(c->cpTrans);
        mxFree(c->osize);
        mxFree(c->cpOther);
        mxFree(c);
    }
}

/*====================== Kernel related methods ====================*/


int ResampSep::IsKernelDefArray( const mxArray* C )
{
    return !( (mxGetClassID(C) != mxCELL_CLASS)
              || (GetCount(C) != 2)
              || (mxGetClassID(mxGetCell(C,0)) != mxDOUBLE_CLASS)
              || (mxGetClassID(mxGetCell(C,1)) != mxDOUBLE_CLASS)
              || (GetCount(mxGetCell(C,0)) != 1)
              || (mxGetScalar(mxGetCell(C,0)) <= 0.0)
              || (GetCount(mxGetCell(C,1)) < 2) );
}


Kernel* ResampSep::NewKernel( const mxArray* kernelDef )
{
    Kernel* k = NULL;



    if( kernelDef != NULL && IsKernelDefArray(kernelDef) )
    {
        const mxArray* halfwidth    = mxGetCell(kernelDef,0);
        const mxArray* positiveHalf = mxGetCell(kernelDef,1);
        k = (Kernel *)mxCalloc(1,sizeof(Kernel));
        k->halfwidth = mxGetScalar(halfwidth);
        k->stride = (mwSize)(ceil(2.0 * (k->halfwidth)));
        k->nSamplesInPositiveHalf = GetCount(positiveHalf);
        k->positiveHalf = 
            (double *)mxCalloc(k->nSamplesInPositiveHalf,sizeof(double));
        k->indexFactor = (k->nSamplesInPositiveHalf - 1) / k->halfwidth;
        for( mwSize i = 0; i < k->nSamplesInPositiveHalf; i++ )
        {
            k->positiveHalf[i] = (mxGetPr(positiveHalf))[i];
        }
    }
    else if( kernelDef == NULL || !mxIsEmpty(kernelDef) )
    {
        mxAssert(0, "Kernel definition must be either empty or a cell"
                 " array with form {halfWidth, positiveHalf}.");
    }

    return k;
}

void ResampSep::DestroyKernel( Kernel* k )
{
    if( k != NULL )
    {
       mxFree(k->positiveHalf);
       mxFree(k);
    }
}


/*
 * Apply bilinear interpolation to the resampling kernel if in range,
 * return zero otherwise.
 *
 * This function performs a lookup operation.  It computes the distance
 * from the new pixel to neighboring pixels and looks up a weight for that
 * distance in the interpolation kernel.  Interpolation kernel is
 * essentialy a lookup table.
 */

double ResampSep::EvaluateKernel( const Kernel* k, const double t )
{
    double result = 0.0;

    // To do: Decide which side should have equality. 
    // (We're doing a convolution.)
    if( -(k->halfwidth) < t && t <= (k->halfwidth) )
    {
        double x = k->indexFactor * fabs(t);
        mwSize index = (mwSize)x;  /* This is equivalent to (mwSize)floor(x) if x>0 */
        if( index >= k->nSamplesInPositiveHalf - 1 )
        {
            result = k->positiveHalf[k->nSamplesInPositiveHalf - 1];
        }
        else
        {
            /* WJ Surprisingly, removing this operation, by replacing it with 
               result = k->positiveHalf[index]; did not provide much of a 
               speedup */
            double w1 = x - index;
            double w0 = 1.0 - w1;
            result = w0 * (k->positiveHalf[index]) + 
                     w1 * (k->positiveHalf[index + 1]);
        }
    }

    return result;
}

Kernel** ResampSep::CreateKernelset( const mwSize nTrans, const mxArray* K )
{
    const bool singleton = (GetCount(K) == 1);
    const mxArray* sharedKernelDef = NULL;
    Kernel** kernelset = (Kernel**)mxCalloc(nTrans,sizeof(Kernel*));
    if( !mxIsEmpty(K) && mxGetClassID(K) != mxCELL_CLASS )
    {
        mexErrMsgIdAndTxt("Images:resampsep:inputKMustBeEmptyOrCellArray",
                          "%s","K (interpolating kernel array) should be"
                          " either empty or a cell array.");
    }

    if( singleton || mxIsEmpty(K) || IsKernelDefArray(K) )
    {
        /* Non-null K0 => All kernels are the same. */
        sharedKernelDef = (singleton ? mxGetCell(K,0) : K);
    }

    for( mwSize j = 0; j < nTrans; j++ )
    {
        kernelset[j] =
            (sharedKernelDef != NULL ? NewKernel(sharedKernelDef) :
             NewKernel(mxGetCell(K,j)));
    }
    return kernelset;
}

void ResampSep::DestroyKernels( const mwSize nTrans, Kernel* kernelset[] )
{
    for( mwSize j = 0; j < nTrans; j++ )
    {
        DestroyKernel(kernelset[j]);
    }

    mxFree(kernelset);
}

/*=========================== Pad Method helper function ===================*/

PadMethod ResampSep::PadMethodFromString( const mxArray* padmethod )
{
    const char* usage =
        "padmethod must be 'fill', 'bound', 'replicate',"
        " 'circular', or 'symmetric'.";
    
    char buf[32];
    if( mxGetString(padmethod, buf, 32) )
    {
        mexErrMsgIdAndTxt("Images:resampsep:badPadmethod",
                          "%s",usage);
    }

    if( strcmp(buf,"fill")      == 0 ) return Fill;
    if( strcmp(buf,"bound")     == 0 ) return Bound;
    if( strcmp(buf,"replicate") == 0 ) return Replicate;
    if( strcmp(buf,"circular")  == 0 ) return Circular;
    if( strcmp(buf,"symmetric") == 0 ) return Symmetric;

    mexErrMsgIdAndTxt("Images:resampsep:badPadMethod",
                      "%s",usage);
    return Replicate;
}

/*========================== Convolver related methods =====================*/

Convolver* ResampSep::NewConvolver( const Config* iCon, const mxArray* K,
                                    const PadMethod padmethod  )
{
    Convolver* c = (Convolver*)mxCalloc(1,sizeof(*c));

    c->padmethod  = padmethod;
    c->ndims      = iCon->nTrans;
    c->size       = (mwSignedIndex*)mxCalloc(c->ndims, sizeof(mwSignedIndex));
    c->tsize      = (mwSignedIndex*)mxCalloc(c->ndims, sizeof(mwSignedIndex));
    c->cumsum     = (mwSize*)mxCalloc(c->ndims + 1, sizeof(mwSize));
    c->cumprod    = (mwSize*)mxCalloc(c->ndims + 1, sizeof(mwSize));
    c->weights    = (double**)mxCalloc(c->ndims, sizeof(double*));
    c->tsub       = (mwSignedIndex**)mxCalloc(c->ndims, sizeof(mwSignedIndex*));
    c->kernelset  = CreateKernelset(c->ndims, K);

    c->cumprod[0] = 1;
    c->cumsum[0]  = 0;
    for( mwSize k = 0; k < c->ndims; k++ )
    {
        c->size[k] = (c->kernelset[k] == NULL ? 1 : c->kernelset[k]->stride);
        c->tsize[k] = iCon->tsize[k];
        c->cumsum[k+1]  = c->size[k] + c->cumsum[k];
        c->cumprod[k+1] = c->size[k] * c->cumprod[k];
    }

    c->weight_data = (double*)mxCalloc(c->cumsum[c->ndims],sizeof(double));
    c->tsub_data   = (mwSignedIndex*)mxCalloc(c->cumsum[c->ndims],sizeof(mwSignedIndex));
    for( mwSize k = 0; k < c->ndims; k++ )
    {
        c->weights[k] = &(c->weight_data[c->cumsum[k]]);
        c->tsub[k]    = &(c->tsub_data[c->cumsum[k]]);
    }

    c->nPoints     = c->cumprod[c->ndims];
    c->useFill     = (bool*)mxCalloc(c->nPoints,sizeof(bool));
    c->subs        = (mwSignedIndex**)mxCalloc(c->nPoints, sizeof(mwSignedIndex*));
    c->subs_data   = (mwSignedIndex*)mxCalloc((c->nPoints) * (c->ndims),sizeof(mwSignedIndex));
    for( mwSize j = 0; j < c->nPoints; j++ )
    {
        c->subs[j] = &(c->subs_data[j * (c->ndims)]);
    }

    c->lo = (double*)mxCalloc(c->ndims,sizeof(double));
    c->hi = (double*)mxCalloc(c->ndims,sizeof(double));
    for( mwSize k = 0; k < c->ndims; k++ )
    {
        double h = ( (c->kernelset[k] != NULL && padmethod == Fill) ?
                      c->kernelset[k]->halfwidth : 0.5 );
        if( c->tsize[k] >= 1)
        {
            c->lo[k] = -h;
            c->hi[k] = (c->tsize[k] - 1) + h;
        }
        else /* Never in bounds if tsize is zero. */
        {
            c->lo[k] =  1;
            c->hi[k] = -1;
        }
    }

    c->localIt   = NewIterator(c->ndims, c->size);

    return c;
}


void ResampSep::DestroyConvolver( Convolver* c)
{
    if( c != NULL )
    {
        mxFree(c->size);
        mxFree(c->tsize);
        mxFree(c->cumprod);
        mxFree(c->cumsum);
        mxFree(c->weights);
        mxFree(c->tsub);
        mxFree(c->subs);
        mxFree(c->weight_data);
        mxFree(c->tsub_data);
        mxFree(c->useFill);
        mxFree(c->subs_data);
        mxFree(c->lo);
        mxFree(c->hi);
        DestroyKernels(c->ndims,c->kernelset);
        DestroyIterator(c->localIt);
        mxFree(c);
    }
}

inline
bool ResampSep::UseConvolution( const Convolver* c, const double* p )
{
    if( c->padmethod == Fill || c->padmethod == Bound )
    {
        for( mwSize k = 0; k < c->ndims; k++ )
        {
            if( !(c->lo[k] <= p[k] && p[k] < c->hi[k]) ) return false;
        }
    }

    return true;
}

/*
 * Adjust subscript according to the pad method. If a sub falls
 * outside the interval [0 tsize-1], follow a rule to adjust it to
 * fall within the interval or (if padMethod == Fill) set it to -1.
 * Also, adjust the subscript to -1 in the case of an empty input
 * array.
 */


mwSignedIndex ResampSep::AdjustSubscript( const mwSignedIndex subscript, const mwSignedIndex tsize, 
                                          const PadMethod padMethod )
{
    mwSignedIndex sub = subscript;
    mwSignedIndex tsize2;

    if( tsize <= 0 ) return -1; // Avoid potential zero-divide with
                                // empty input array

    switch( padMethod )
    {
     case Fill:
        sub = (sub < 0 ? -1 : (sub >= tsize ? -1 : sub));
        break;

     case Bound:
     case Replicate:
        sub = (sub < 0 ? 0 : (sub >= tsize ? tsize - 1 : sub));
        break;

     case Circular:
        sub = (sub >= 0 ? sub % tsize : tsize - 1 - ((-sub - 1) % tsize) );
        break;

     case Symmetric:
        tsize2 = 2 * tsize;
        sub = (sub >= 0 ? sub % (tsize2) : 
               tsize2 - 1 - ((-sub - 1) % tsize2) );
        sub = (sub >= tsize ? (tsize2) - sub - 1  : sub);
        break;
    }

    return sub;
}

/*
 * Given a specific point, p, construct subscript and weight arrays
 * needed to compute the value of the convolution at a single point.
 * (Use Convolve() to perform the actual computation.)
 * Set up nearest-neighbor for each dimension with a null kernel.
 * Adjust the subscripts to fall within the interval [0 tsize[k]]
 * or set to -1. A null kernel signifies nearest-neighbor interpolation.
 */

void ResampSep::InitConvolver( Convolver* c, const double* p )
{
    for( mwSize k = 0; k < c->ndims; k++ )
    {
        mwSignedIndex s0;
        mwSignedIndex s;

        // Either use the kernel or simply round p[k] for nearest
        // neighbor resampling.
        if(c->kernelset[k] != NULL)
        {
            /*
             * Translate the kernel so that its center is at center, then
             * return the lowest integer for which the kernel is defined.
             */
            s0 = (mwSignedIndex)(ceil(p[k] - c->kernelset[k]->halfwidth));

            for( mwSignedIndex j = 0; j < c->size[k]; j++ )
            {
                s = s0 + j;

                /* use the kernel */
                c->weights[k][j] = EvaluateKernel(c->kernelset[k], p[k] - s);
                c->tsub[k][j] = AdjustSubscript( s, c->tsize[k],
                                                 c->padmethod );
            }
        }
        else
        {
            /* use MATLAB-Compatible Rounding Function */
	    /* this used to be (int)floor(p[k]+0.5); */
            s0 = (p[k] < 0.0 ? (mwSignedIndex)(p[k]-0.5) : (mwSignedIndex)(p[k] + 0.5));

            for( mwSignedIndex j = 0; j < c->size[k]; j++ )
            {
                s = s0 + j;

                // Set the single weight to unity for nearest 
                // neighbor resampling.
                c->weights[k][j] = 1.0;
                c->tsub[k][j] = AdjustSubscript( s, c->tsize[k], 
                                                 c->padmethod );
            }
        }
    }

    /* Save the outer product set of c->tsub in t->localsubs */
    for( ResetIterator(c->localIt); !DoneIterating(c->localIt); 
         IncrementIterator(c->localIt) )
    {
        c->useFill[c->localIt->offset] = false;
        for( mwSize k = 0; k < c->ndims; k++ )
        {
            mwSignedIndex s = c->tsub[k][c->localIt->subs[k]];
            c->subs[c->localIt->offset][k] = s;

            /* Turn on useFill if the adjusted subscript is less than
               one in any of the dimensions. */
            if( s == -1 ) c->useFill[c->localIt->offset] = true;
        }
    }
}

/*
 * Convolve
 *
 * Logically, v is a multidimensional array of input samples with
 * size given by c->size. Loop over each dimension of v from highest
 * to lowest and apply the weights, collapsing v successively until 
 * the weighted sum is fully accumluated in v[0].
 * (This function takes advantage of the separability of the
 * resampling kernels.)
 */

double ResampSep::Convolve( const Convolver* c, double* v )
{
    for (mwSize k_up = 0; k_up < c->ndims; k_up++)
    {
        mwSize k = c->ndims - 1 - k_up;

        mwSize      n = c->size[k];
        double* w = c->weights[k];

        // Use the cumulative size product to iterate over the first
        // k-1 dimensions of v -- all the uncollapsed dimensions except 
        // for the highest. (For each value of q we visit a different 
        // point in this subspace.)
        mwSize s = c->cumprod[k];
        for( mwSize q = 0; q < s; q++ )
        {
            // Take the inner product of the k-th weight array and a (1-D) 
            // profile across the highest (i.e., k-th) uncollapsed dimension
            // of v.  Re-use memory by storing the result back in the first 
            // element of the profile.
            double t = 0.0;
            for( mwSize r = 0; r < n; r++ )
            {
                t += w[r] * v[q + r * s];
            }
            v[q] = t;
        }
    }

    /* After repeated iterations, everything is collapsed into the */
    /* first element of v. It is now zero-dimensional (a scalar).  */
    return v[0];
}


/*======= Special 'Cumulative Product' for the Fill Value Array ============*/


mwSize* ResampSep::CreateCPFill( const mwSize nOther, const mxArray* F )
{
    /*
     * Prepare for efficient computation of offsets into the fill
     * array F. Create a special cumulative product array for F.
     * Its length is nOther (which may be greater than NDIMS(F)) and
     * elements corresponding to a value of 1 in SIZE(F) are set to zero.
     * The result can be multiplied (dot product via Subs2Offset) with
     * a set of non-transform subscripts to produce the correct offset into F.
     */

    const mwSize   ndims_F = mxGetNumberOfDimensions(F);
    const mwSize*  size_F  = mxGetDimensions(F);
    mwSize*        cpFill  = (mwSize*)mxCalloc( nOther, sizeof(mwSize) );

    mwSize partialprod = 1;
    for( mwSize k = 0; k < ndims_F && k < nOther; k++ )
    {
        if( size_F[k] == 1 )
            cpFill[k] = 0;
        else
            cpFill[k] = partialprod;

        partialprod *= size_F[k];
    }

    for( mwSize k = ndims_F; k < nOther; k++ )
    {
        cpFill[k] = 0;
    }

    return cpFill;
}

/*========================== Subs2Offset =========================*/

/*
 * The cumulative product values in cumprod must be pre-ordered
 * to be consistent with the order of the subscripts in subs.
 */

mwSignedIndex ResampSep::Subs2Offset( const mwSize ndims, const mwSize* cumprod, 
                                      const mwSignedIndex* subs )
{
    mwSignedIndex offset = 0;

    mxAssert(subs    != NULL, "");
    mxAssert(cumprod != NULL, "");

    for (mwSize k = 0; k < ndims; k++)
    {
        offset += subs[k] * cumprod[k];
    }

    return offset;
}

//////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////
mxArray* ResampSep::evaluate()
{
    mxClassID mClassID  = mxGetClassID(f_A);
    bool      bComplexA = mxIsComplex(f_A);


    //Allocate space for the output array and get pointer to its data
    mxArray    *B = mxCreateNumericArray(f_oCon->ndims, 
                                         f_oCon->size,
                                         mClassID,
                                         (bComplexA?mxCOMPLEX:mxREAL) );
    void *ptrBr = mxGetData(B);

    void *ptrBi = NULL;
    if(bComplexA) ptrBi = mxGetImagData(B);


    switch(mClassID)
    {
      case(mxLOGICAL_CLASS):
        resample((uint8_T *)ptrBr, (uint8_T *)ptrBi);
        break;

      case(mxUINT8_CLASS):
        resample((uint8_T *)ptrBr, (uint8_T *)ptrBi);
        break;

      case(mxINT8_CLASS):
        resample((int8_T *)ptrBr, (int8_T *)ptrBi);
        break;

      case(mxUINT16_CLASS):
        resample((uint16_T *)ptrBr, (uint16_T *)ptrBi);
        break;

      case(mxINT16_CLASS):
        resample((int16_T *)ptrBr, (int16_T *)ptrBi);
        break;

      case(mxUINT32_CLASS):
        resample((uint32_T *)ptrBr, (uint32_T *)ptrBi);
        break;

      case(mxINT32_CLASS):
        resample((int32_T *)ptrBr, (int32_T *)ptrBi);        
        break;

      case(mxSINGLE_CLASS):
        resample((float *)ptrBr, (float *)ptrBi);
        break;

      case(mxDOUBLE_CLASS):
        resample((double *)ptrBr, (double *)ptrBi);
        break;

      default:
        //Should never get here
        mexErrMsgIdAndTxt("Images:resampsep:unsuppDataType", "%s",
                          "Unsupported data type.");
        break;
    }

    return(B);
}

//////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////
ResampSep::ResampSep( 
        const mxArray*  A,
        const mxArray*  M,
        const mxArray*  tdims_A,
        const mxArray*  tdims_B,
        const mxArray*  fsize_A,
        const mxArray*  fsize_B,
        const mxArray*  F,
        const mxArray*  padstr,
        const mxArray*  K
        )
{
    /* Working objects and storage:
     *   cumprodM    Cumulative product array for computing offsets into M
     *   iCon:       Configuration of input dimensions
     *   oCon:       Configuration of output dimensions
     *   convolver:  Weights and subscripts needed for convolution
     *   oTransIt:   Iterator for output transform space
     *   otherIt:    Iterator for non-transform space
     *   cpFill:     Cumulative products for computing offsets into the 
     *               fill array
     *   p:          Current point in input transform space
     *   vReal:      Input values to be used in convolution (real parts)
     *   vImag:      Input values to be used in convolution (imaginary parts)
     *
     * Return value:
     *   B:          Output array (allocated in the main routine)
     */

    f_cumprodM  = GetCumProd(M);
    f_iCon      = NewConfig( fsize_A, tdims_A );
    f_oCon      = NewConfig( fsize_B, tdims_B );
    f_convolver = NewConvolver( f_iCon, K, PadMethodFromString(padstr) );
    f_oTransIt  = NewIterator( f_oCon->nTrans, f_oCon->tsize );
    f_otherIt   = NewIterator( f_oCon->nOther, f_oCon->osize );
    f_cpFill    = CreateCPFill( f_oCon->nOther, F );
    f_p         = (double*)mxCalloc( f_iCon->nTrans, sizeof(double) );
    f_vReal     = (double*)mxCalloc( f_convolver->nPoints,
                                     sizeof(double) );
    f_vImag     = (double*)mxCalloc( f_convolver->nPoints, 
                                     sizeof(double) );


    f_M = M;
    f_F = F;
    f_A = A;
}

//////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////
ResampSep::~ResampSep()
{
    mxFree(f_cumprodM);
    DestroyConfig(f_iCon);
    DestroyConfig(f_oCon);
    DestroyConvolver(f_convolver);
    DestroyIterator(f_oTransIt);
    DestroyIterator(f_otherIt);
    mxFree(f_cpFill);
    mxFree(f_p);
    mxFree(f_vReal);
    mxFree(f_vImag);
}

/*========================== mexFunction ===========================*/

extern "C"
void mexFunction( int nlhs, mxArray       *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
    const mxArray* A       = prhs[0];  // Input array
    const mxArray* M       = prhs[1];  // Inverse mapping from output to input
    const mxArray* tdims_A = prhs[2];  // List of input transform dimensions
    const mxArray* tdims_B = prhs[3];  // List of output transform dimensions
    const mxArray* fsize_A = prhs[4];  // Full size of input array
    const mxArray* fsize_B = prhs[5];  // Full size of output block
    const mxArray* F       = prhs[6];  // Fill value array
    const mxArray* R       = prhs[7];  // Resampler

    const mxArray* rdata;   // Resampler's data struct
    const mxArray* padstr;  // Pad method string
    const mxArray* K;       // Interplating kernel cell array

    if( nrhs != 8 )
    {
        mexErrMsgIdAndTxt("Images:resampsep:eightInputsAreRequired",
                          "%s","Eight input arguments are required.");
    }
    else if( nlhs > 1 )
    {
        mexErrMsgIdAndTxt("Images:resampsep:tooManyOutputs",
                          "%s","Too many output arguments.");
    }

    if( !mxIsStruct(R) )
    {
        mexErrMsgIdAndTxt("Images:resampsep:inputRMustBeStruct",
                          "%s","R must be a struct.");
    }

    if( (rdata = mxGetField(R,0,"rdata")) == NULL )
    {
        mexErrMsgIdAndTxt("Images:resampsep:missingRdataFieldInR",
                          "%s","R must have an rdata field.");
    }

    if( !mxIsStruct(rdata) )
    {
        mexErrMsgIdAndTxt("Images:resampsep:rdataMustBeAStruct",
                          "%s","R.rdata must be a struct.");
    }

    if( (padstr = mxGetField(R,0,"padmethod")) == NULL )
    {
        mexErrMsgIdAndTxt("Images:resampsep:missingPadmethodFieldInR",
                          "%s","R must have a padmethod field.");
    }

    if( (K = mxGetField(rdata,0,"K")) == NULL )
    {
         mexErrMsgIdAndTxt("Images:resampsep:rdataMustHaveKField",
                           "%s","R.rdata must have a K (interpolating kernel)"
                           " field.");
    }


    //Instantiate the class
    ResampSep resampSep( A, M, tdims_A, tdims_B, fsize_A, 
                         fsize_B, F, padstr, K );

    plhs[0] = resampSep.evaluate();

}
