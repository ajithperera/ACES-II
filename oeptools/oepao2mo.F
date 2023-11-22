      subroutine oepao2mo (
     &     nbas, nocc, evec,
     &     aa, ma, pp, ph, hh)
c
c This routine transforms a packed one-electron integral from AO to MO
c representation.
c
c In: nbas, nocc, evec, aa
c Scr: ma
c Out: pp, ph, hh
c 
c IVS, $Id: oepao2mo.FPP,v 1.1.1.1 2008/05/29 19:35:39 taube Exp $
c 
      implicit none
c
c     Arguments
c 
      integer
     &     nbas, nocc
c     
      double precision
     &     evec (nbas, nbas), aa (*), ma (*), pp (*), ph (*), hh (*)
c
c     Local variables
c     
      integer
     &     n1, n2, i1, i2, i3
c
c     Define the one- and two-electron offset operators
c     
#define DEFINE_OPERATORS
#include <oepintind.com>
#define DEFINE_VARIABLES
#define N_BAS nbas
#define N_OCC nocc
#include <oepintind.com>
c
c     Initialize the arrays with zeroth.
c     
      call zero (ma, n_ma)
      call zero (pp, n_pp)
      call zero (ph, n_ph)
      call zero (hh, n_hh)
c
c     AA->MA. Transform the first AO index into a MO index.
c     
      do n1 = 1, n_aa
         i2 = i2_aa (n1)
         i1 = i1_aa (n1, i2)
         do i3 = 1, n_a
            n2 = i_ma (i3, i2)
            ma (n2) = ma (n2) + evec (i1, i3) * aa (n1)
            if (i1.ne.i2) then
               n2 = i_ma (i3, i1)
               ma (n2) = ma (n2) + evec (i2, i3) * aa (n1)
            endif
         enddo
      enddo
c
c     MA->PP,PH,HH. Transform the second AO index into a MO index.
c     
      do n1 = 1, n_pp
         i2 = i2_pp (n1)
         i1 = i1_pp (n1, i2) + n_h
         i2 = i2 + n_h
         do i3 = 1, n_a
            n2 = i_ma (i1, i3)
            pp (n1) = pp (n1) + ma (n2) * evec (i3, i2)
         enddo
      enddo
c     
      do n1 = 1, n_ph
         i2 = i2_ph (n1)
         i1 = i1_ph (n1, i2) + n_h
         do i3 = 1, n_a
            n2 = i_ma (i1, i3)
            ph (n1) = ph (n1) + ma (n2) * evec (i3, i2)
         enddo
      enddo
c     
      do n1 = 1, n_hh
         i2 = i2_hh (n1)
         i1 = i1_hh (n1, i2)
         do i3 = 1, n_a
            n2 = i_ma (i1, i3)
            hh (n1) = hh (n1) + ma (n2) * evec (i3, i2)
         enddo
      enddo
c     
      return
      end


