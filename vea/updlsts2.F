C
C  UPDATES THE ACTUAL LENGTHS OF AN S2 LIST.
C  INPUT IS orb_m, num_s, AND irp_s, WHICH ALLOW DETERMINATION OF
C  MOIODS(irp_bra,iFamily) AND DistType_bra WHICH ALLOWS DETERMINATION OF
C  DISTRIBUTION SIZE
C
      SUBROUTINE UPDLSTS2(iFamily, orb_m, num_s, DistType_bra, irp_s)
      IMPLICIT NONE
      integer iFamily, orb_m(*), num_s, DistType_bra, irp_s
      integer irp_bra, irp_m, nCols(8), i
C
#include "syminf.com" /* for nirrep and dirprd */
#include "sympop.com" /* for irpdpd */
#ifdef _ASSERT
      i = 0
      if ((DistType_bra.lt.1).or.(22.lt.DistType_bra)) then
         print *, '@UPDLSTS2: Assertion failed.'
         print *, '   DistType_bra = ',DistType_bra
         i = 1
      end if
      if ((irp_s.lt.1).or.(nirrep.lt.irp_s)) then
         print *, '@UPDLSTS2: Assertion failed.'
         print *, '   irp_s = ',irp_s
         i = 1
      end if
      if (i.ne.0) call aces_exit(i)
#endif
C
C FIRST DETERMINE ROW-SIZE FROM num_s AND irp_s
C
      do irp_bra = 1, nirrep
         irp_m = dirprd(irp_bra,irp_s)
         nCols(irp_bra) = orb_m(irp_m) * num_s
      end do
C
C  THE RIGHT HAND TYPE [I,S] DOES NOT HAVE A LABEL
C
      call aces_family_danger(nirrep,                iFamily,
     &                        irpdpd(1,DistType_bra),nCols,
     &                        DistType_bra,          1)
      RETURN
      END

