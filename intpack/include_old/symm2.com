
#ifndef _SYMM2_COM_
#define _SYMM2_COM_

c  irpsz1 - this vector contains the sizes of each of the (aa|aa)
c           and (ab|ab) symmetry combinations.  the element of
c           the vector is indexed by indx(a,b), where "a" and "b"
c           are the irrep numbers.
c
c  irpsz2 - same as irpsz1, except for (aa|bb) symmetry combinations.
c           the indexing is indx(a-1,b).
c
c  irpds1 - vector which contains the distribution size of the
c           (aa|aa) and (ab|ab) symmetry combinations.  indexed
c           as indx(a,b), where "a" and "b" are irrep numbers.
c
c  irpds2 - vector which contains the distribution sizes (both
c           bra and ket) of the (aa|bb) symmetry combinations.
c           irpds2(odd) contain the bra distribution sizes,
c           irpds2(even) contain the ket distribution sizes.  these
c           entries are accessed by 2*indx(a-1,b)-1 and 2*indx(a-1,b),
c           or an equivalent counting.
c
c  inewvc - vector which contains the relative orbital number within
c           an irrep for each of the orbitals (no dropped mos).
c
c  idxvec - lookup vector for determining irreps of orbitals

#include "maxbasfn.par"

      integer irpsz1(36),irpsz2(28),irpds1(36),
     &    irpds2(56),inewvc(maxbasfn),idxvec(maxbasfn)

      common /symm2/ irpsz1,irpsz2,irpds1,irpds2,
     &    inewvc,idxvec

c  ipkoff -

      integer ipkoff(73)

      common /pkoff/ ipkoff

#endif /* _SYMM2_COM_ */

