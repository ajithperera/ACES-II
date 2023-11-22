      subroutine oeprefham (
     &     nbas, nocc, naux, reference,
     &     scfevl, scfevc, hrtaa, hfxaa, refevl)
c
c This routine calculates the reference one-particle Hamiltonian
c eigenvalues according to the type of the reference chosen.
c
c In:  nbas, nocc, naux, scfevl, scfevc, hrtaa, hfxpp, hfxhh
c Out: refevl
c 
c Igor Schweigert, Mar 2004
c $Id: oeprefham.FPP,v 1.1.1.1 2008/05/29 19:35:39 taube Exp $
c 
      implicit none
c
c     Arguments
c     
      integer
     &     nbas, nocc, naux
c 
      character*(*)
     &     reference
c     
      double precision
     &     scfevl (nbas), scfevc (nbas, nbas),
     &     hrtaa (*), hfxaa (*), refevl (nbas)
c     
c     Local variables
c     
      integer
     &     n
c 
      double precision
     &     x, OEPORBSHIFT
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
     &     t3,'@OEPREFHAM-W: The reference eigenvalues are degenerate.'/
     &     t3,'The HOMO-LUMO gap is ',f20.10,' a.u.',/)
 5910 format (
     &     t3,'@OEPREFHAM-I: The virtual orbital energies will be',/
     &     t3,'shifted by ',f10.5,' a.u.',/)
c
c     Set the orbital shift (later will be set by a keyword)
c     
      OEPORBSHIFT = 2.d-1
c     
c     If the OEP Hamiltonian is used as the reference Hamiltonian, copy
c     the OEP eigenvalues into the reference eigenvalues.
c 
      if (reference (1:4).eq.'OEP ') then
         call dcopy (
     &        nbas,
     &        scfevl, 1,
     &        refevl, 1)
      endif
c
c     If the OEP Hamiltonian with the diagonal Fock elemenets is used as
c     the reference Hamiltonian, calculate the Fock diagonal elements
c     and copy them into the reference eigenvalues: <> Initialize the
c     eigenvalues with zeros. <> Add the contribution from the Hartree
c     Hamiltonian. <> Add the contribution from the non-local exchange
c     potential. <> Return.
c     
      if (reference (1:4).eq.'DIAG'.or.reference (1:4).eq.'SEMI') then
         call zero (refevl, nbas)
         do n = 1, n_aa
            j = i2_aa (n)
            i = i1_aa (n, j)
            x = 2.d0
            if (i.eq.j) x = 1.d0
            do k = 1, n_a
               refevl (k) = refevl (k) + 
     &              x * (hrtaa (n) + hfxaa (n)) *
     &              scfevc (i, k) * scfevc (j, k)
            enddo
         enddo
      endif
c     
      if (refevl (n_h) - refevl (n_h+1).gt.-1.d-3) then
         write (6, 5905) refevl (n_h) - refevl (n_h+1)
         write (6, 5910) OEPORBSHIFT
         do i = 1, n_p
            refevl (n_h+i) = refevl (n_h+i) + OEPORBSHIFT
         enddo
      endif
c 
      return
      end
      
      
