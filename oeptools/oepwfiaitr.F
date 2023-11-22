      subroutine oepwfiaitr (
     &     nbas, nocc, naux,
     &     icore, maxmem, usdmem,
     &     auxao, scfham, scfevl, scfevc, hfxao)
c
c This routine solves the exchange Optimized Effective Potential
c equation APPROXIMATELY finding a local potential whose hole-particle
c MEs weighted with the energy denominators are equal to those of
c non-local exchange potential.
c
c The difference with OEPFIA method is that we scale both l.h.s. and
c r.h.s with the energy denominators, which would not affect the
c analytical solution, but does affect the numerical solution.
c
c This is done by projecting the corresponding equation to an auxiliary
c basis set and solving the resulting set of linear equations. Note that
c matrix of auxiliary MOMEs is, in general, rectangular and
c singular. The approximate solution to these equation can be obtained
c using the Singular Value Decomposition. 
c
c The interface with VSCF is the same as for VSCF/OEPITER, see comments
c therein.
c
c Igor Schweigert, Sep 2003
c $Id: oepwfiaitr.FPP,v 1.1.1.1 2008/05/29 19:35:39 taube Exp $
c
      implicit none
c
c     Common blocks
c
#include <symm2.com>
#include <machsp.com>
c
c     The input parameters                                         
c
      integer
     $     nbas, nocc, naux, 
     $     icore (*), maxmem, usdmem
c
      double precision
     $     auxao (itriln (nirrep+1), naux),
     $     scfham (itriln (nirrep+1)), scfevl (nbas),
     &     scfevc (nbas, nbas), hfxao (itriln (nirrep+1))
c
c     Local variables
c
      integer
     $     x, i, a,
     $     i010, i020, i030, i040, i050, i060, itop, coretop,
     $     auxhp, oepau, hfxhp, nsvd, ierr, ndep
c
      double precision
     $     svdthr, dsvd, ddot
c
      integer ibufln, luint
      parameter (ibufln=600)
      parameter (luint=10)
c
c     Messages to the user
c     
 1005 format (
     $     t3,'@OEPWFIAITR-F: DGESVD failed to converge')
 1010 format (
     $     t3,'@OEPWFIAITR-F: Illegal argument to DGESVD, argument #',
     $     i2)
 1020 format (
     $     t3,'@OEPWFIAITR-F: ',i2,' singularities of the auxiliary',/,
     $     t3,'hole-particle MEs have been removed.')
c
c     Set the SVD threshold. This is temporary, later it will be set
c     by a keyword
c
      svdthr = 1.d-8
c
c     Set the top of the core to the begging of the ICORE array.
c     (It has been already been shifted in VSCF)
c
      coretop = 1
c
c     Transform auxiliary AO MEs to MO MEs:
c       1/ allocate memory for MO MEs and a scratch array;
c       2/ check if there is enough memory available;
c       3/ transform auxilary AO MEs to occ.-virt. MEs. and
c          scale them with the energy denominators.
c
      auxhp = coretop
      coretop = auxhp + iintfp * nocc * (nbas - nocc) * naux
      i010 = coretop
      i020 = i010 + iintfp * nbas * nbas
      itop = i020 + iintfp * nocc * nbas
c
      if (itop .gt. maxmem) call insmem (
     &     'OEPFIAITR-F: Transformation of auxiliary b.s. AO MEs',
     &     itop+usdmem, maxmem+usdmem)
c
      do x=1, naux
         call oepao2mo (
     &        nbas, nocc, 'IA  ', scfevc, auxao (1, x), 
     &        icore (i010), icore (i020),
     &        icore (auxhp+iintfp*nocc*(nbas-nocc)*(x-1)))
         call oepenrscl (
     &        nbas, nocc, 'IA  ', scfevl,
     &        icore (auxhp+iintfp*nocc*(nbas-nocc)*(x-1)))
      enddo
c
c     Transform non-local exchange AO MEs to MO MEs:
c       1/ allocate memory for MO MEs and a scratch array;
c       2/ check if there is enough memory available;
c       3/ transform n.-l. exch. AOMEs to occ.-virt. MEs;
c       4/ scale the n.-l. exch. MOMEs with the energy
c          denominators.
c
      hfxhp = coretop
      coretop = hfxhp + iintfp * nocc * (nbas - nocc) * naux
      i010 = coretop
      i020 = i010 + iintfp * nbas * nbas
      itop = i020 + iintfp * nocc * nbas
c
      if (itop .gt. maxmem) call insmem (
     &     'OEPFIAITR-F: Transformation of n.-l. exch. AO MEs',
     &     itop+usdmem, maxmem+usdmem)
c
      call oepao2mo (
     &     nbas, nocc, 'IA  ', scfevc, hfxao, 
     &     icore (i010), icore (i020),
     &     icore (hfxhp))
c     
      call oepenrscl (
     &     nbas, nocc, 'IA  ', scfevl, icore (hfxhp))
c
c     Build the inverse of the auxiliary MOMEs and find the local pontential
c     whose hole-particle MEs are equal to those of n.-l. exch. To do this:
c       1/ figure out what is the least dimension of the auxiliary MOMEs;
c       2/ allocate memory for the local potential auxiliary expansion, SVD
c          eigenvalues, left and right eigenvectors, and a scratch array;
c       3/ check if there is enough memory available;
c       4/ find SVD of the auxiliary MOMEs;
c       5/ check if SVD converged
c       6/ invert the singular eigenvalues;
c       7/ let user know how many singularities have been removed
c       8/ contract n.-l. exch. hole-particle MEs with the SVD's right
c          eigenvectors, singular eigenvalues, and left eigenvectors
c       9/ store OEP in a.b.s. representation to JOBARC
c      10/ contract auxiliary AOMEs with the OEP solution, and add the
c          resulting OEP AOMEs to the Hamiltonian
c
c     Memory pointers for Step 2
c      i010 - SVD eigenvalues 
c      i020 - SVD left eigenvectors
c      i030 - SVD right eigenvectors
c      i040 - SVD scratch array
c     
      nsvd = naux
      if (nsvd .gt. nocc*(nbas-nocc)) nsvd = nocc * (nbas-nocc)
c
      oepau = coretop
      coretop = oepau + iintfp * naux
      i010 = coretop
      i020 = i010 + iintfp * nsvd
      i030 = i020 + iintfp * nsvd * nocc * (nbas-nocc)
      i040 = i030 + iintfp * nsvd * naux
      i050 = i040 + iintfp * nsvd
      itop = i050 + iintfp * (5*nsvd + max (naux, nocc*(nbas-nocc)))
c     
      if (itop .gt. maxmem) call insmem (
     $     'OEPFIAITR-F: Singular Value Decomp. of auxiliary MOMEs',
     $     itop+usdmem, maxmem+usdmem)
c
      call dgesvd (
     $     's', 's', nocc*(nbas-nocc), naux, 
     $     icore (auxhp), nocc*(nbas-nocc), icore (i010),
     $     icore (i020), nocc*(nbas-nocc), icore (i030), nsvd,
     $     icore (i040), 5*nsvd + max (naux, nocc*(nbas-nocc)),
     $     itop, ierr)
c     
      if (ierr.ne.0) then
         if (ierr.gt.0) then
            write (6, 1005)
         else
            write (6, 1010) abs (ierr)
         endif
         call errex
      endif
c     
      ndep = 0
      do i=1, nsvd
         call scopy (1, icore (i010 + iintfp *(i-1)), 1, dsvd, 1)
         if ( dabs (dsvd) .gt. svdthr ) then
            dsvd = 1.d0 / dsvd
         else
            dsvd = 0.d0
            ndep = ndep + 1
         endif
         call scopy (1, dsvd, 1, icore (i010 + iintfp *(i-1)), 1)
      enddo
c
      if (ndep.gt.0) write (6, 1020) ndep
c     
      call xgemm(
     &     't', 'n', nsvd, 1, nocc*(nbas-nocc),
     &     1.d0, icore (i020), nocc*(nbas-nocc),
     &     icore (hfxhp), nocc*(nbas-nocc),
     &     0.d0, icore (i050), nsvd)
c     
      do i=1, nsvd
         dsvd = ddot (
     $        1, icore (i050+iintfp*(i-1)), 1,
     $        icore (i010+iintfp*(i-1)), 1)
         call scopy (1, dsvd, 1, icore (i050+iintfp*(i-1)), 1)
      enddo
c      
      call xgemm(
     &     't', 'n', naux, 1, nsvd,
     &     1.d0, icore (i030), nsvd,
     &     icore (i050), nsvd,
     &     0.d0, icore (oepau), naux)
c
      call putrec (
     &     20, 'JOBARC', 'OEPAU', iintfp*naux, icore (oepau))
c     
      call xgemm(
     &     'n', 'n', itriln (nirrep+1), 1, naux,
     &     1.d0, auxao, itriln (nirrep+1),
     &     icore (oepau), naux,
     &     1.d0, scfham, itriln (nirrep+1))
c
      return
      end
      
