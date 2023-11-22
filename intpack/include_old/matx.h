
#ifndef _MATX_H_
#define _MATX_H_

c A block-diagonal symmetric matrix can be stored in several different forms.
c These include:
c   full:  All elements (including zeroes) are stored.
c   sqr :  Each square block is stored (both upper and lower triangles), but
c          none of the zeros.
c   tri :  Only the upper triangular part of each block is stored.
c
c An individual block can be referenced in one of several ways including:
c   triblk :  One of the triangular blocks.
c   triseg :  A specific block in a "tri" matrix.
c   sqrblk :  One of the square blocks.
c   sqrseg :  A specific block in a "sqr" matrix.
c   fullseg:  A specific block in a "full" matrix (the zeroes are ignored).
c
c To make working with these matrices and blocks easier, the following
c constants are defined:

#define MAT_FULL     0
#define MAT_SQR      1
#define MAT_TRI      2

#define MAT_FULLSEG  10
#define MAT_SQRBLK   20
#define MAT_SQRSEG   21
#define MAT_TRIBLK   30
#define MAT_TRISEG   31

#endif /* _MATX_H_ */

