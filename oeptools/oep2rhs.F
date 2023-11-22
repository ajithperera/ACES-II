      subroutine oep2rhs (
     &     nbas, nocc, naux,
     &     isdff, isdfw, isdww, isfdf, iswdw, isdffdiag, hybrid,
     &     ener, eval, oep1au, 
     &     auxpp, auxph, auxhh, f1pp, f1ph, f1hh,
     &     ppph, pphh, phph, phhh, intpp, intph, inthh, intr12, rhs,
     &     dor12)
c
c This routine calculate the right hand side of the second-order OEP
c equation and the second-order energy.
c
c For the detail on the first- and second-order OEP equation see
c document REPORTS/OEP2.DERIVATION.TEX
c
c To make the routine more universial the contribution of the each
c diagram can be neglected by setting the corresponding flag ISxxx to
c FALSE. E.g., by the exchange-only calculation can be emulating by
c setting all the flags to FALSE.
c
c In: nbas, nocc, naux, eval, auxpp, auxph, auxhh, f1pp, f1ph, f1hh,
c     ppph, pphh, phph, phhh
c Scr: intpp, intph, inthh,
c Out: rhs
c
c Igor Schweigert, Jan 2004
c $Id: oep2rhs.FPP,v 1.4 2008/06/06 18:09:17 taube Exp $
c
      implicit none
c
c     Arguments
c
      logical
     &     isdff, isdfw, isdww, isfdf, iswdw, isdffdiag, dor12
c     
      integer
     &     nbas, nocc, naux
c
      double precision
     &     hybrid, ener, eval (nbas), oep1au (naux), 
     &     auxpp (*), auxph (*), auxhh (*),
     &      f1pp (*),  f1ph (*),  f1hh (*),
     &      ppph (*),  pphh (*),  phph (*), phhh (*),
     &     intpp (*), intph (*), inthh (*), rhs (*)
c
c     Local variables
c
      integer
     &     n, n1, n2, n3, n4, a, b, c, i1, i2, i3, i4
c
      double precision
     &     e1, e2, x
c     R12 variables
c     Need 4 byte integers to connect with Noga's code
      integer*4
     &     nh4,np4
c
      integer
     &     unit_mos 
c     
      double precision           
     &     er12,r12fact,intr12(*)
c
c     Define the one- and two-electron array offset operators
c
#define DEFINE_OPERATORS
#include <oepintind.com>
#define DEFINE_VARIABLES
#define N_BAS nbas
#define N_OCC nocc
#include <oepintind.com>
c
c     Messages to the user
c     
 5905 format (
     &     t3,'@OEP2RHS-W: Reference eigenvalues are degenerate.',/,
     &     t3,'The HOMO-LUMO gap is ',f20.10)
c     
c     Calculate the MOMEs of the one-particle part of the perturbation
c     by subtracting the first-order OEP MOMEs from non-local exchange
c     MOMEs and scale by the 1 minus the fraction of the non-local
c     exchange.
c
      call dgemm (
     &     'n', 'n', n_pp, 1, naux,
     &     -1.d0+hybrid, auxpp, n_pp,
     &     oep1au, naux,
     &     1.d0-hybrid, f1pp, n_pp)
      call dgemm (
     &     'n', 'n', n_ph, 1, naux,
     &     -1.d0+hybrid, auxph, n_ph,
     &     oep1au, naux,
     &     1.d0-hybrid, f1ph, n_ph)
      call dgemm (
     &     'n', 'n', n_hh, 1, naux,
     &     -1.d0+hybrid, auxhh, n_hh,
     &     oep1au, naux,
     &     1.d0-hybrid, f1hh, n_hh)
c     
c     Calculate the second-order contribution to the total energy: <>
c     Calculate the singles contibutions. <> Calculate the doubles
c     contribution.
c     
      e1 = ener
      do n1 = 1, n_ph
         i2 = i2_ph (n1)
         i1 = i1_ph (n1, i2)
         ener = ener +
     &        2.d0 * f1ph (n1) * f1ph (n1) / (eval (i2) - eval (i1+n_h))
      enddo
c$$$      print *,'E2s = ', ener - e1, ' GAP = ', eval (n_h) - eval (n_h+1)
      if (eval (n_h) - eval (n_h+1).gt.-1.d-3)
     &     write (6, 5905) eval (n_h) - eval (n_h+1)
c      e1r12=ener
c      do n1 = 1, n_hh
c         ener = ener + 2.d0 * v1(n1) * c1(n1)
c      end do
c$$$      print *,'E2r12S = ', ener - e1r12
c     
      e2 = ener
      do n1 = 1, n_phph
         i4 = i4_phph (n1)
         i3 = i3_phph (n1, i4)
         i2 = i2_phph (n1, i3, i4)
         i1 = i1_phph (n1, i2, i3, i4)
         e1 = eval (i2) + eval (i4) - eval (i1+n_h) - eval (i3+n_h)
         n2 = i_phph_cr (i1, i4, i3, i2)
         ener = ener +
     &        phph (n1) * (2.d0 * phph (n1) - phph (n2)) / e1
         if (i_ph (i1, i2) .ne. i_ph (i3, i4)) then
            n2 = i_phph_cr (i3, i2, i1, i4)
            ener = ener +
     &           phph (n1) * (2.d0 * phph (n1) - phph (n2)) /e1
         endif
      enddo
c$$$      print *, 'E2d = ', ener - e2
c      e2r12 = ener
c      do n1 = 1, n_hhhh
c      do n1 = 1, n_phph
c         i4 = i4_phph (n1)
c         i3 = i3_phph (n1, i4)
c         i2 = i2_phph (n1, i3, i4)
c         i1 = i1_phph (n1, i2, i3, i4)
c         e1 = eval (i2) + eval (i4) - eval (i1+n_h) - eval (i3+n_h)
c         n2 = i_phph_cr (i1, i4, i3, i2)
c         ener = ener +
c     &        phph (n1) * (2.d0 * phph (n1) - phph (n2)) / e1
c         if (i_ph (i1, i2) .ne. i_ph (i3, i4)) then
c            n2 = i_phph_cr (i3, i2, i1, i4)
c            ener = ener +
c     &           phph (n1) * (2.d0 * phph (n1) - phph (n2)) /e1
c         endif
c      enddo
c$$$      print *, 'E2d = ', ener - e2
c     
c     Initialize the intermediates.
c
      call zero (intpp, n_pp)
      call zero (intph, n_ph)
      call zero (inthh, n_hh)
c
c     DFF - SINGLES PH. Calculate the contribution of the DFF diagram.
c     
      if (.not.isdff .or. hybrid.eq.1d0) goto 190
      do n = 1, n_ph
         i = i2_ph (n)
         b = i1_ph (n, i)
         do a = 1, n_p
            if (.not.isdffdiag.and.a.eq.b) goto 110
            e1 = eval (i) - eval (b+n_h)
            e2 = eval (i) - eval (a+n_h)
            n1 = i_pp_cr (b, a)
            n2 = i_ph    (a, i)
            intph (n) = intph (n) + 4.d0 * 
     &           f1pp (n1) * f1ph (n2) /
     &           (e1 * e2)
 110        continue
         enddo
      enddo
c
      do n = 1, n_ph
         j = i2_ph (n)
         a = i1_ph (n, j)
         do i = 1, n_h
            if (.not.isdffdiag.and.i.eq.j) goto 120
            e1 = eval (j) - eval (a+n_h)
            e2 = eval (i) - eval (a+n_h)
            n1 = i_hh_cr (j, i)
            n2 = i_ph    (a, i)
            intph (n) = intph (n) - 4.d0 *
     &           f1hh (n1) * f1ph (n2) /
     &           (e1 * e2)
 120        continue
         enddo
      enddo     
 190  continue
c
c     DFW - SINGLE and DOUBLES PH. Calculate the contribution of the DFW
c     and DWF diagrams
c     
      if (.not.isdfw .or. hybrid.eq.1d0) goto 290
      do n = 1, n_ph
         j = i2_ph (n)
         b = i1_ph (n, j)
         do a = 1, n_p
            do i = 1, n_h
               e1 = eval (j) - eval (b+n_h)
               e2 = eval (i) - eval (a+n_h)
               n1 = i_ph (a, i)
               n2 = i_pphh_cr (b, a, i, j)
               n3 = i_phph_cr (b, i, a, j)
               intph (n) = intph (n) - 4.d0 *
     &              f1ph (n1) * (pphh (n2) + phph (n3)) /
     &              (e1 * e2)
            enddo
         enddo
      enddo
 290  continue
c
c     DWW - DOUBLES PH. Calculate the contribution of the DWW diagram.
c
      if (.not.isdww) goto 390
      do n = 1, n_ph
         i = i2_ph (n)
         c = i1_ph (n, i)
         do a = 1, n_p
            do b = 1, a
               do j = 1, n_h
                  e1 = eval (i) - eval (c+n_h)
                  e2 = eval (i) + eval (j) - eval (a+n_h) - eval (b+n_h)
                  n1 = i_ppph_cr (a, c, b, j)
                  n2 = i_phph_cr (a, i, b, j)
                  n3 = i_ppph_cr (b, c, a, j)
                  n4 = i_phph_cr (a, j, b, i)
                  if (b.eq.a) e1 = 2.d0 * e1
                  intph (n) = intph (n) + 4.d0 * (
     &                 ppph (n1) * phph (n2) * 2.d0 -
     &                 ppph (n1) * phph (n4) -
     &                 ppph (n3) * phph (n2) +
     &                 ppph (n3) * phph (n4) * 2.d0) /
     &                 (e1 * e2)
               enddo
            enddo
         enddo
      enddo
c
      do n = 1, n_ph
         k = i2_ph (n)
         a = i1_ph (n, k)
         do b = 1, n_p
            do i = 1, n_h
               do j = 1, i
                  e1 = eval (k) - eval (a+n_h)
                  e2 = eval (i) + eval (j) - eval (a+n_h) - eval (b+n_h)
                  n1 = i_phhh_cr (b, j, k, i)
                  n2 = i_phph_cr (a, i, b, j)
                  n3 = i_phhh_cr (b, i, k, j)
                  n4 = i_phph_cr (a, j, b, i)
                  if (j.eq.i) e1 = 2.d0 * e1
                  intph (n) = intph (n) - 4.d0 * (
     &                 phhh (n1) * phph (n2) * 2.d0 -
     &                 phhh (n1) * phph (n4) -
     &                 phhh (n3) * phph (n2) +
     &                 phhh (n3) * phph (n4) * 2.d0) /
     &                 (e1 * e2)
               enddo
            enddo
         enddo
      enddo
c     
 390  continue

c
c     FDF - SINGLES PP, HH. Calculate the contribution of the FDF diagram
c     
      if (.not.isfdf .or. hybrid.eq.1d0) goto 490
      do n = 1, n_pp
         a = i2_pp (n)
         b = i1_pp (n, a)
         do i = 1, n_h
            e1 = eval (i) - eval (b+n_h)
            e2 = eval (i) - eval (a+n_h)
            n1 = i_ph (b, i)
            n2 = i_ph (a, i)
            intpp (n) = intpp (n) + 2.d0 * 
     &           f1ph (n1) * f1ph (n2) /
     &           (e1 * e2)
         enddo
      enddo
c     
      do n = 1, n_hh
         i = i2_hh (n)
         j = i1_hh (n, i)
         do a = 1, n_p
            e1 = eval (j) - eval (a+n_h)
            e2 = eval (i) - eval (a+n_h)
            n1 = i_ph (a, j)
            n2 = i_ph (a, i)
            inthh (n) = inthh (n) - 2.d0 * 
     &           f1ph (n1) * f1ph (n2) /
     &           (e1 * e2)
         enddo
      enddo
 490  continue

c
c     WDW - DOUBLES PP, HH. Calculate the contribution of the WDW diagram.
c
      if (.not.iswdw) goto 590
      do n = 1, n_pp
         a = i2_pp (n)
         c = i1_pp (n, a)
         do b = 1, n_p
            do i = 1, n_h
               do j = 1, i
                  e1 = eval (i) + eval (j) - eval (c+n_h) - eval (b+n_h)
                  e2 = eval (i) + eval (j) - eval (a+n_h) - eval (b+n_h)
                  n1 = i_phph_cr (c, i, b, j)
                  n2 = i_phph_cr (a, i, b, j)
                  n3 = i_phph_cr (c, j, b, i)
                  n4 = i_phph_cr (a, j, b, i)
                  if (i.eq.j) e1 = 2.d0 * e1
                  intpp (n) = intpp (n) + 2.d0 * (
     &                 phph (n1) * phph (n2) * 2.d0 -
     &                 phph (n1) * phph (n4) -
     &                 phph (n3) * phph (n2) +
     &                 phph (n3) * phph (n4) * 2.d0) /
     &                 (e1 * e2)
               enddo
            enddo
         enddo
      enddo
c
      do n = 1, n_hh
         i = i2_hh (n)
         k = i1_hh (n, i)
         do j = 1, n_h
            do a = 1, n_p
               do b = 1, a
                  e1 = eval (i) + eval (j) - eval (a+n_h) - eval (b+n_h)
                  e2 = eval (k) + eval (j) - eval (a+n_h) - eval (b+n_h)
                  n1 = i_phph_cr (a, k, b, j)
                  n2 = i_phph_cr (a, i, b, j)
                  n3 = i_phph_cr (a, j, b, k)
                  n4 = i_phph_cr (a, j, b, i)
                  if (a.eq.b) e1 = 2.d0 * e1
                  inthh (n) = inthh (n) - 2.d0 * (
     &                 phph (n1) * phph (n2) * 2.d0 -
     &                 phph (n1) * phph (n4) -
     &                 phph (n3) * phph (n2) +
     &                 phph (n3) * phph (n4) * 2.d0) /
     &                 (e1 * e2)
               enddo
            enddo
         enddo
      enddo
c$$$      write (6,*) 'INTHH 1: ', inthh (1)
c$$$      inthh (1) = 0.d0
c$$$      do a = 1, n_p
c$$$         do b = 1, n_p
c$$$            n = i_phph_cr (a, 1, b, 1)
c$$$            e1 = 2.d0 * eval (1) - eval (a+1) - eval (b+1)
c$$$            inthh (1) = inthh (1) - 2.d0 * phph (n)**2 / e1**2
c$$$         enddo
c$$$      enddo
c$$$      write (6,*) 'INTHH 1: ', inthh (1)
 590  continue
      
      if (dor12) then
C       R12 contribution
C       Completely separate read of R12 OV_OO File test      
        unit_mos=8
        open(unit=unit_mos,file='OV_OO_BLOCK_F12',ACCESS='SEQUENTIAL',
     &       STATUS='OLD',FORM='UNFORMATTED')
        read(unit=unit_mos) er12,nh4,np4
        n1=nh4
        n2=np4
        read(unit=unit_mos) (intr12(i),i=1,n1*n2)
        close(unit=unit_mos)
C       R12 array is stored hole-particle
C       OEP expects particle-hole
C       Need to check proper factor... don't know exactly how to check this yet
        r12fact=1.0d0      
        do n = 1, n_ph
C n1 is the index of the transposed matrix               
          i = i2_ph (n)
          a = i1_ph (n, i)
          n1 = i + n_h * (a-1)   
          intph(n) = intph(n)+r12fact*intr12(n1)
        enddo
      endif       
       
c
c     Calculate the OEP r.h.s. in auxiliary basis by contracting the RHS
c     intermediate with the auxiliary MOMEs: <> Scale the diagonal
c     elements of auxiliary PP and HH MOMEs by a factor of one-half. <>
c     Contract the intermediates with the auxiliary MOMEs, type by
c     type. <> Scale the diagonal elements of auxiliary PP and HH MOMEs
c     by a factor of two.
c
c     Note that the fact the full matrix is twice as long as PP, PH, HH
c     triangular matrices is taken into account by factor of two in the
c     DGEMMs. Hence, in order to prevent the double-counting of the
c     diagonal elements (PP and HH), they are scaled by a factor of
c     one-half.
c
      do i = 1, naux
         do j = 1, n_p
            n = i_pp (j, j) + n_pp * (i-1)
            auxpp (n) = .5d0 * auxpp (n)
         enddo
         do j = 1, n_h
            n = i_hh (j, j) + n_hh * (i-1)
            auxhh (n) = .5d0 * auxhh (n)
         enddo
      enddo
c
      call dgemm (
     &     'n', 'n', 1, naux, n_pp,
     &     2.d0, intpp, 1,
     &     auxpp, n_pp,
     &     0.d0, rhs, 1)
c
      call dgemm (
     &     'n', 'n', 1, naux, n_ph,
     &     1.d0, intph, 1,
     &     auxph, n_ph,
     &     1.d0, rhs, 1)
c
      call dgemm (
     &     'n', 'n', 1, naux, n_hh,
     &     2.d0, inthh, 1,
     &     auxhh, n_hh,
     &     1.d0, rhs, 1)
c     
      do i = 1, naux
         do j = 1, n_p
            n = i_pp (j, j) + n_pp * (i-1)
            auxpp (n) = 2.d0 * auxpp (n)
         enddo
         do j = 1, n_h
            n = i_hh (j, j) + n_hh * (i-1)
            auxhh (n) = 2.d0 * auxhh (n)
         enddo
      enddo
c     
      return
      end
