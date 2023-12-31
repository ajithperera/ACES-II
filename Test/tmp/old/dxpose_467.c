
/*
 * This routine transposes a rectangular matrix in situ.
 * The algorithm was blatantly ripped from netlib/TOMS.
 * (cf. http://www.netlib.org/tomspdf/467.pdf)
 *
 * NOTE:
 *    This algorithm was originally a copy of algorithm 513
 * which is a rework of algorithm 380. Matrices of certain sizes
 * did not completely transpose (e.g., 34x19, 34x13, 34x4, 33x17,
 * 31x11, 31x7, 31x6, 29x17, 29x9, 29x8, 29x5, 26x11, 23x12, 19x7,
 * and 19x4). It is possible an error arose in reordering the
 * logic since both the C and Fortran versions show the error.
 * However, algorithm 467 is faster and [Anthony Yau] didn't want
 * to chase down a logic error.
 *
 * ARGUMENTS:
 *   a(,) - the double precision array to transpose
 *   n1   - THE 'FASTEST-RUNNING' INDEX
 *          [In Fortran, this is the number of rows in a(,), but
 *           in C, this is the number of columns.]
 *   n2   - THE 'SLOWEST-RUNNING' INDEX
 *          [In Fortran, this is the number of columns in a(,), but
 *           in C, this is the number of rows.]
 */

#include <stdio.h>
#include <stdlib.h>

#define bool  short
#define TRUE  1
#define FALSE 0

#define PRIMES 8 /* the maximum number of unique primes for factoring */

void
#ifdef C_SUFFIX
  dxpose_467_
#else
  dxpose_467
#endif /* C_SUFFIX */
(double * a, long * n1, long * n2)
{

    if ((*n1 < 2) || (*n2 < 2)) return;

    if (*n1 == *n2)
    {
        long j;
        for ( j=0; j<(*n1-1); j++)
        {
            long i;
            for ( i=(j+1); i<*n1; i++)
            {
                double b       = *(a+i+(*n1*j));
                *(a+i+(*n1*j)) = *(a+j+(*n1*i));
                *(a+j+(*n1*i)) = b;
            }
        }
    }
    else
    {

        long n = *n1;
        long m = ((*n1)*(*n2))-1;
        long lCount, lDiv, lStart;

        long lFactor[PRIMES], lPower[PRIMES], lExponent[PRIMES], lExp[PRIMES];
        long lNumber; /* to reinforce the point, lNumber must be <= PRIMES */

        unsigned char *work; long workdim;

        long p; /* a common array pointer */

        bool repeat_main_loop;

        {
            size_t size;
            workdim = ((*n1)+(*n2))>>1;
            size = 1 + (workdim<<3);
            work = malloc(size);
            if (work == NULL)
            {
                printf("@dmat_xpose: Not enough free system memory.\n");
                exit(-1);
            }
            else
            {
                for ( p=0; p<size; p++) *(work+p) &= 0;
                workdim = size<<3;
            }
        }

        {
            /* factor modulus m into prime powers */
            /*
             * EXAMPLE:
             *    m = 1960 = (2**3)*(5**1)*(7**2)
             *    lNumber   = 3
             *    lFactor   = 2, 5,  7
             *    lExponent = 3, 1,  2
             *    lPower    = 8, 5, 49
             */

            long p = -1;      /* start the index before the array   */
            long lCurFac = 0; /* the current prime factor           */
            long lEm = m;     /* a temp argument to abuse liberally */
            long lDiv = 2;    /* the first prime divisor to test    */

            bool still_going = TRUE;
            while (still_going)
            {
                long lQuot = lEm/lDiv;
                still_going = FALSE;
                if (lEm == lDiv*lQuot)
                {
                    if (lDiv > lCurFac)
                    {
                        p += 1;
                        lFactor[p]   = lDiv;
                        lPower[p]    = lDiv;
                        lCurFac      = lDiv;
                        lExponent[p] = 1;
                    }
                    else
                    {
                        lPower[p]    *= lDiv;
                        lExponent[p] += 1;
                    }
                    lEm = lQuot;
                    still_going = TRUE;
                }
                else
                {
                    if (lQuot > lDiv)
                    {
                        if (lDiv == 2) lDiv  = 3;
                        else           lDiv += 2;
                        still_going = TRUE;
                    }
                }
            } /* while (still_going) */

            if (lEm > 1)
            {
                if (lEm > lCurFac)
                {
                    p += 1;
                    lFactor[p]   = lEm;
                    lPower[p]    = lEm;
                    lExponent[p] = 1;
                }
                else
                {
                    lPower[p]    *= lEm;
                    lExponent[p] += 1;
                }
            }

            lNumber = 1 + p; /* the number of unique primes that divide m */
        }

        for ( p=0; p<lNumber; p++) lExp[p] = 0;
        lCount = m;
        for ( p=0; p<lNumber; p++)
        {
            lCount /= lFactor[p];
            lCount *= lFactor[p]-1;
        }
        lDiv   = 1;
        lStart = 1;

        repeat_main_loop = TRUE;
        while (repeat_main_loop)
        {

            bool cycle_elements;
            long mmlSt = m - lStart;

            if (lStart == lDiv)
            {
                cycle_elements = TRUE;
            }
            else
            {
                if (    (lStart <= workdim)
                     && ((*(work+(lStart>>3)) >> (lStart % 8)) & 1)
                   )
                {
                    cycle_elements = FALSE;
                }
                else
                {
                    long lSolD = lStart/lDiv;
                    cycle_elements = TRUE;
                    for ( p=0; p<lNumber; p++)
                    {
                        if (    (lExp[p] != lExponent[p])
                             && ((lSolD % lFactor[p]) == 0)
                           ) cycle_elements = FALSE;
                    }
                    if ((lStart > workdim) && cycle_elements)
                    {
                        long lTest = lStart;
                        for (;;)
                        {
                            lTest = (n*lTest) % m;
                            if ((lTest < lStart) || (lTest > mmlSt))
                            {
                                cycle_elements = FALSE; break;
                            }
                            else
                            {
                                if ((lTest <= lStart) || (lTest >= mmlSt))
                                {
                                    cycle_elements = TRUE; break;
                                }
                            }
                        } /* for (;;) */
                    } /* if ((lStart > workdim) && cycle_elements) */
                } /* if ((lStart <= workdim) && (...)) */
            } /* if (lStart == lDiv) */

            if (cycle_elements)
            {
                double atemp = *(a+lStart);
                double btemp = *(a+mmlSt);
                long   a1    = lStart;
                bool repeat_inner_loop = TRUE;
                while (repeat_inner_loop)
                {
                    long a2   = (n*a1) % m;
                    long mma1 = m - a1;
                    long mma2 = m - a2;
                    if (a1   <= workdim) *(work+(a1  >>3)) |= (1 << (a1   % 8));
                    if (mma1 <= workdim) *(work+(mma1>>3)) |= (1 << (mma1 % 8));
                    lCount -= 2;
                    repeat_inner_loop = FALSE;
                    if (lStart == a2)
                    {
                        *(a+a1)   = atemp;
                        *(a+mma1) = btemp;
                    }
                    else
                    {
                        if (lStart == mma2)
                        {
                            *(a+a1)   = btemp;
                            *(a+mma1) = atemp;
                        }
                        else
                        {
                            *(a+a1)   = *(a+a2);
                            *(a+mma1) = *(a+mma2);
                            a1 = a2;
                            repeat_inner_loop = TRUE;
                        }
                    } /* if (lStart == a2) */
                } /* while (repeat_inner_loop) */
            } /* if (cycle_elements) */

            repeat_main_loop = FALSE;
            lStart += lDiv;
            if (lCount > 0)
            {
                repeat_main_loop = TRUE;
            }
            else
            {
                for ( p=0; p<lNumber; p++)
                {
                    if (!repeat_main_loop)
                    {
                        if (lExp[p] == lExponent[p])
                        {
                            lExp[p]  = 0;
                            lDiv    /= lPower[p];
                        }
                        else
                        {
                            lExp[p] += 1;
                            lDiv    *= lFactor[p];
                            if (lDiv >= (m>>1)) return;
                            lCount = m/lDiv;
                            {
                                long i;
                                for ( i=0; i<lNumber; i++)
                                {
                                    if (lExp[i] != lExponent[i])
                                    {
                                        lCount /= lFactor[i];
                                        lCount *= lFactor[i]-1;
                                    }
                                }
                                for ( i=0; i<(workdim>>3); i++) *(work+i) &= 0;
                            }
                            lStart = lDiv;
                            repeat_main_loop = TRUE;
                        }
                    } /* if (!repeat_main_loop) */
                } /* for ( p=0; p<lNumber; p++) */
            } /* if (lCount > 0) */

        } /* while (repeat_main_loop) */

        free(work);
    } /* if (*n1 == *n2) */

    return;
}

