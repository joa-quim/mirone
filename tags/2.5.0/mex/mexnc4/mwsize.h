/*
 * Minimalist thing to define mwSize and let newer MEXs compiler in R13
 * This should be harmless to ML versions where mwSize is natively defined.
 */

#ifndef MX_COMPAT_32

typedef int mwSize;
typedef int mwIndex;
typedef int mwSignedIndex;

#define MWSIZE_MAX    2147483647UL
#define MWINDEX_MAX   2147483647UL
#define MWSINDEX_MAX  2147483647L
#define MWSINDEX_MIN -2147483647L
#define MWSIZE_MIN    0UL
#define MWINDEX_MIN   0UL

#endif
