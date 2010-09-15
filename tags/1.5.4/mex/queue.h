// Copyright 1993-2007 The MathWorks, Inc.
// $Revision: 1.1.8.1 $  $Date: 2008/12/22 23:43:43 $

#ifndef TOOLBOX_IMAGES_IMAGES_QUEUE_H
#define TOOLBOX_IMAGES_IMAGES_QUEUE_H

#include "sequencemex.h"

template<typename _T>
class Queue : public Sequence<_T>
{
 public:
    inline void put(_T x);
    inline _T get(void);
};

//////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////
template<typename _T>
inline void Queue<_T>::put(_T x)
{
    this->addHigh(x);
}

//////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////
template<typename _T>
inline _T Queue<_T>::get()
{
    return(this->removeLow());
}


#endif // QUEUE_H 
