      subroutine oep1homo (
     &     nbas, nocc, naux, hybrid, 
     &     scfevc, hfxaa, oep1au, auxaa, homo)
c
c This routine calculates the difference between the expectation values
c of non-local and local exchange operators with respect to the HOMO.
c Since these values must be equal, the difference gives an estimation of
c the error arising from the incompleteness of the AO and auxiliary basis
c sets. HYBRID gives the value for the fraction of the nonlocal exchange. 
c
c Igor Schweigert, Sep 2004
c $Id: oep1homo.FPP,v 1.1.1.1 2008/05/29 19:35:39 taube Exp $
c
      implicit none
c
c     Arguments
c
      integer
     &     nbas, nocc, naux
c
      double precision
     &     scfevc (nbas, nbas), hfxaa (*), oep1au (naux),
     &     auxaa (*), homo, hybrid
c
c     Local variables
c
      integer
     &     n1, i1, i2
c
c     Define the two-electron integral offset operators and the
c     associated variables.
c
#define DEFINE_OPERATORS     
#include "oepintind.com"
#define DEFINE_VARIABLES
#define N_BAS nbas
#define N_OCC nocc
#include "oepintind.com"
c
c     Calculate the difference between the OEP and non-local exchange
c     HOMO MOMEs. <> Subtract the exchange OEP AOMEs from the non-local
c     exchange AOMEs. <> Calculate the difference and scale it by 1 -
c     the fraction of nonlocal exchange. <> Restore the non-local
c     exchange AOMEs.
c
      call dgemm (
     &     'n', 'n', n_aa, 1, naux,
     &     -1.d0, auxaa, n_aa,
     &     oep1au, naux,
     &     1.d0, hfxaa, n_aa)
c     
      homo = 0.d0
      do n1 = 1, n_aa
         i2 = i2_aa (n1)
         i1 = i1_aa (n1, i2)
         if (i1 .eq. i2) then
            homo=homo+scfevc(i1,n_h)*scfevc(i2,n_h)*hfxaa(n1)
         else
            homo=homo+scfevc(i1,n_h)*scfevc(i2,n_h)*hfxaa(n1)*2.d0
         endif
      enddo
      homo = homo * (1.0d0-hybrid)
c
      call dgemm (
     &     'n', 'n', n_aa, 1, naux,
     &     1.d0, auxaa, n_aa,
     &     oep1au, naux,
     &     1.d0, hfxaa, n_aa)
c
      return
      end

         
         
      
      
