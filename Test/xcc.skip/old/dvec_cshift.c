
/*
 * This routine performs a circular shift of elements in an
 * M_REAL floating-point vector. The signs of shift are defined:
 *  (-) shift right (or down)
 *  (+) shift left  (or up  )
 */

#include <stdlib.h>
#include <string.h>
#include "aces.h" /* defines the M_REAL floating-point type */
#define _DTMP_DIM 8 /* the length of the dTmp buffer */

void
#ifdef C_SUFFIX
  dvec_cshift_
#else
  dvec_cshift
#endif /* C_SUFFIX */
( M_REAL * pdVec, long * plLen, long * plShift)
{
    long lShift2; /* changeable arguments */

    /* quick returns */
#ifdef _DEBUG
    if (0 > *plLen)
    {
        printf("@dvec_cshift: The vector cannot have a negative length.\n");
        abort();
    }
#endif /* _DEBUG */
    if (2 > *plLen) return;
    if ( 0 == (lShift2 = *plShift % *plLen) ) return;

    /* reduce calls to memmove by shifting the other way fewer times */
    if ( labs(lShift2<<1) > *plLen )
    {
        long lTmp = (0 < lShift2) ? -(*plLen) : *plLen ;
        lShift2 += lTmp;
    }

    /* Note:
     *    If memmove fails, the code will abort without warning with
     * a bus error. Attempts at checking the value of *memmove on exit
     * are meaningless since memmove will not return after failing.
     */

    /* shift right/down */
    if (0 > lShift2)
    {
        M_REAL *pdVec_save, *pdVec_dest, dTmp[_DTMP_DIM];

        long lAwait = -lShift2;
        int  iNum   = (lAwait % _DTMP_DIM);
        if (0 < iNum)
        {
            long   lSrcLen = *plLen - iNum;
            size_t Bytes   = lSrcLen*sizeof(dTmp[1]);
            pdVec_save = pdVec + lSrcLen;
            pdVec_dest = pdVec + iNum;
            {
                int i;
                for ( i=0; i<iNum; i++) dTmp[i] = *(pdVec_save+i);
                memmove( pdVec_dest, pdVec, Bytes);
                for ( i=0; i<iNum; i++) *(pdVec+i) = dTmp[i];
                lAwait -= iNum;
            }
        }
        iNum = _DTMP_DIM;
        if ( 0 < (lAwait /= _DTMP_DIM) )
        {
            long   lSrcLen = *plLen - iNum;
            size_t Bytes   = lSrcLen*sizeof(dTmp[1]);
            pdVec_save = pdVec + lSrcLen;
            pdVec_dest = pdVec + iNum;
            for ( lAwait=(-lAwait); 0>lAwait; lAwait++)
            {
                int i;
                for ( i=0; i<iNum; i++) dTmp[i] = *(pdVec_save+i);
                memmove( pdVec_dest, pdVec, Bytes);
                for ( i=0; i<iNum; i++) *(pdVec+i) = dTmp[i];
            }
        }
    }
    /* shift left/up */
    else
    {
        M_REAL *pdVec_src, *pdVec_restore, dTmp[_DTMP_DIM];

        long lAwait = lShift2;
        int  iNum   = (lAwait % _DTMP_DIM);
        if (0 < iNum)
        {
            long   lSrcLen = *plLen - iNum;
            size_t Bytes   = lSrcLen*sizeof(dTmp[1]);
            pdVec_src     = pdVec + iNum;
            pdVec_restore = pdVec + lSrcLen;
            {
                int i;
                for ( i=0; i<iNum; i++) dTmp[i] = *(pdVec+i);
                memmove( pdVec, pdVec_src, Bytes);
                for ( i=0; i<iNum; i++) *(pdVec_restore+i) = dTmp[i];
                lAwait -= iNum;
            }
        }
        iNum = _DTMP_DIM;
        if ( 0 < (lAwait /= _DTMP_DIM) )
        {
            long   lSrcLen = *plLen - iNum;
            size_t Bytes   = lSrcLen*sizeof(dTmp[1]);
            pdVec_src     = pdVec + iNum;
            pdVec_restore = pdVec + lSrcLen;
            for ( lAwait=(-lAwait); 0>lAwait; lAwait++)
            {
                int i;
                for ( i=0; i<iNum; i++) dTmp[i] = *(pdVec+i);
                memmove( pdVec, pdVec_src, Bytes);
                for ( i=0; i<iNum; i++) *(pdVec_restore+i) = dTmp[i];
            }
        }
    }

    return;
}

