// $Revision: 1.1.8.1 $
// Copyright 1993-2007 The MathWorks, Inc.
/*
 * Implementation of Sequence abstract data type.
 *
 * A sequence is a homogeneous dynamic array.  When you initialize the
 * sequence, you specify the desired number of items in the array.
 *
 * You can set and get items at a particular zero-based array offset using
 * setItem() and getItem().
 *
 * You can remove array items from either the beginning of the array
 * (removeLow) or the end (removeHigh).  You can add array items from either
 * the beginning of the array (addLow) or the end (addHigh).
 *
 * Use freeSequence to destroy the array and reset sequence's member values to
 * their initial state when class is instantiated.
 *
 * The array used for storing the sequence dynamically grows itself as
 * necessary when items are added.  It does not (in the current
 * implementation) ever shrink itself when items are removed.
 *
 */

#ifndef _SEQUENCEMEX_H
#define _SEQUENCEMEX_H

#include "mex.h"
#include "iptutil_cpp.h"
#include <string.h>

template< typename _T >
class Sequence
{
 private:

    //Member variables
    _T *fArray;
    mwSize fSequenceLength; 
    mwSize fArrayLength;    
    mwSize fHead;            // item # of start of sequence

    //function template declaration
    void expand();
    inline _T *getItemAddress(mwSize i);
    inline _T get(mwSize i);

 public:

    Sequence()
    {
        fArray = NULL;
        fSequenceLength = 0;
        fArrayLength = 0;
        fHead = 0;

    }
    ~Sequence()
    {
        if (fArray)
        {
            mxFree(fArray);
        }
    }

    mwSize getSequenceLength() 
    { 
        return fSequenceLength; 
    }

    // function template declaration
    inline _T getItem(mwSize i); 
    inline _T removeLow();
    inline _T removeHigh();

    void initialize(mwSize h);
    inline void setItem(mwSize i, _T x);
    inline void addHigh(_T x);
    void addLow(_T x); 
    void freeSequence();
    void copyToBuf(_T *x);
};


// Private methods

//////////////////////////////////////////////////////////////////////////////
// getItemAddress
// Return the memory address of the sequence item whose array index is i.
//////////////////////////////////////////////////////////////////////////////

template< typename _T>
inline _T *Sequence<_T>::getItemAddress(mwSize i)
{
    mxAssert(fArray != NULL,ERR_STRING("Sequence::getItemAddress()",
                                       "fArray is NULL."));

    return (_T *)(fArray + ((fHead + i) % fArrayLength));
}


//////////////////////////////////////////////////////////////////////////////
// Expand
// Grow the sequence storage array by a factor of 2.
//////////////////////////////////////////////////////////////////////////////

template< typename _T>
void Sequence<_T>::expand()
{
    mxAssert(fArray != NULL,ERR_STRING("Sequence::getItemAddress()",
                                       "fArray is NULL."));
    mwSize oldLength;
    void *newArray;

    oldLength = fArrayLength;
    mwSize newLength = 2 * fArrayLength;
    newArray = mxRealloc(fArray, newLength * sizeof(_T));
    
    if (!newArray)
    {
        mexErrMsgIdAndTxt("Images:sequencemex:outOfMemory",
                          "%s", "Out of memory.");
    }
    fArray = (_T *)newArray;
    fArrayLength = newLength;
    if (fHead > 0)
    {
        // Move items that were at the tail end of the array to the
        // tail end of the reallocated array.
        void *old = getItemAddress(0);
        memcpy((_T *)old + oldLength, old, (oldLength - fHead)*sizeof(_T));
        fHead += oldLength;
    }
}

//////////////////////////////////////////////////////////////////////////////
// Get
// Get item stored in address
//////////////////////////////////////////////////////////////////////////////
template<typename _T>
inline _T Sequence<_T>::get(mwSize i)
{
    mxAssert(fArray != NULL,ERR_STRING("Sequence::getItem()",
                                       "fArray is NULL."));
    return(*getItemAddress(i));
}

// Public methods

//////////////////////////////////////////////////////////////////////////////
// Initialize sequence. Input is hint which is in the initally allocated
// length in number of items (not bytes) of the sequence array.
//////////////////////////////////////////////////////////////////////////////

template< typename _T>
void Sequence<_T>::initialize(mwSize hint)
{
     if (hint == 0)
     {
         hint = 16;
     }
     fArray = (_T *) mxMalloc(sizeof(_T)*hint);
     if (!fArray)
     {
         mexErrMsgIdAndTxt("Images:sequencemex:outOfMemory",
                          "Out of memory.");
     }
     fArrayLength = hint;
 }

//////////////////////////////////////////////////////////////////////////////
// Get the i-th sequence item.
//////////////////////////////////////////////////////////////////////////////

template<typename _T>
inline _T Sequence<_T>::getItem(mwSize i)
{

    mxAssert(fArray != NULL,ERR_STRING("Sequence::getItem()",
                                       "fArray is NULL."));
    mxAssert(i < fSequenceLength,
             ERR_STRING("Sequence::getItem()", "i is invalid."));

    return(get(i));
}

//////////////////////////////////////////////////////////////////////////////
// Set the i-th sequence item. 
//////////////////////////////////////////////////////////////////////////////

template<typename _T>
inline void Sequence<_T>::setItem(mwSize i, _T x)
{
    mxAssert(fArray != NULL,ERR_STRING("Sequence::setItem()",
                                       "fArray is NULL."));
    mxAssert(i < fSequenceLength,
             ERR_STRING("Sequence::setItem()", "i is invalid."));

    *getItemAddress(i) = x;
}

//////////////////////////////////////////////////////////////////////////////
// Add item to end of the sequence.  This increments the sequence length by 1.
//////////////////////////////////////////////////////////////////////////////

template<typename _T>
inline void Sequence<_T>::addHigh(_T x)
{
    if (fSequenceLength == fArrayLength)
    {
        expand();
    }
    mwSize len = fSequenceLength++;
    setItem(len, x);
}

//////////////////////////////////////////////////////////////////////////////
// Add item to the beginning of the sequence.  This increments the sequence
// length by 1. Item i becomes item i+1 after this operation.
//////////////////////////////////////////////////////////////////////////////

template<typename _T>
void Sequence<_T>::addLow(_T x)
{
    if (fSequenceLength == fArrayLength)
    {
        expand();
    }
    if (--fHead < 0)
    {
        fHead = fArrayLength - 1;
    }
    fSequenceLength++;
    setItem(0,x);        //start of sequence
}
//////////////////////////////////////////////////////////////////////////////
// Remove the first item in the sequence.  This has the effect of decrementing
// the sequence length.  Item i becomes item i-1 after this operation. fHead
// is also incremented by 1. X contains the value that is removed from the
// sequence at the end of the function.
//////////////////////////////////////////////////////////////////////////////

template<typename _T>
inline _T Sequence<_T>::removeLow()
{
    mxAssert(fArray != NULL,ERR_STRING("Sequence::removeLow()",
                                       "fArray is NULL."));
    
    _T x = get(0);      //start of sequence
    
    fHead = (fHead + 1) % fArrayLength;
    --fSequenceLength;

    return(x);
}

//////////////////////////////////////////////////////////////////////////////
// Remove the last item in the sequence.  This has the effect of decrementing
// the sequence length. X contains the value that is removed from the sequence
// at the end of the function.
//////////////////////////////////////////////////////////////////////////////

template<typename _T>
inline _T Sequence<_T>::removeHigh()
{
    mxAssert(fArray != NULL,ERR_STRING("Sequence::removeHigh()",
                                       "fArray is NULL."));

    return(get(--fSequenceLength));
}

//////////////////////////////////////////////////////////////////////////////
// Remove the last item in the sequence.  This has the effect of decrementing
// the sequence length. X contains the value that is removed from the sequence
// at the end of the function.
//////////////////////////////////////////////////////////////////////////////

template<typename _T>
void Sequence<_T>::freeSequence()
{
    if (fArray)
    {
        mxFree(fArray);
        fArray = NULL;
    }
    fSequenceLength = 0;
    fArrayLength = 0;
    fHead = 0;
}

//////////////////////////////////////////////////////////////////////////////
// Copy sequence's contents to a buffer.
//////////////////////////////////////////////////////////////////////////////

template<typename _T>
void Sequence<_T>::copyToBuf(_T *a)
{
    if (fHead == 0)
    {
        memcpy(a, fArray, fSequenceLength*sizeof(_T));
    }
    else
    {
        for (mwSize k = 0; k < fSequenceLength; k++)
        {
            a[k] = getItem(k);
        }
    }
}

#endif
