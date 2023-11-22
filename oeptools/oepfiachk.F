      subroutine oepfiachk (
     &     nbas, nocc, naux,
     &     icore, maxmem, usdmem,
     &     auxao, scfevl, scfevc, hfxao)
c
c This routine check how well f_ia = 0 would approximate the exchange
c OEP equation.
c
c The arguments are identical to those of OEPITER. The memory is dealt
c with exactly the same way it is done is OEPITER, see comments therein.
c
c Note, however, that HFXAO will be modified and will contain the
c difference between non-local and local exchange.
c
c Input: nbas, nocc, naux, maxmem, usdmem, auxao, scfevl, scfevc, hfxao
c Scratch: icore
c Output: hfxao is modified

c Igor Schweigert, Aug-Sep 2003 $Id: oepfiachk.FPP,v 1.1.1.1 2008/05/29 19:35:39 taube Exp $
c
      implicit none
c
c     Common block
c
#include <symm2.com>
#include <machsp.com>
c
c     The input parameters
c
      integer
     &     nbas, nocc, naux, 
     &     icore (*), maxmem, usdmem
c
      double precision
     &     auxao (itriln (nirrep+1), naux),
     &     scfevl (nbas), scfevc (nbas, nbas),
     &     hfxao (itriln (nirrep+1))
c
c     Local variables
c
      integer
     &     x, i, 
     &     i010, i020, i030, i040, i050, i060, i070, i080, i090,
     &     itop, coretop,
     &     auxhp, oepau, ierr, nsvd, ndep, idamax
c
      double precision
     &     svdthr, dsvd, fiachk (20), ddot, dasum
c
c     Messages to the user
c     
 1005 format (
     &     t3, '@OEPFIACHK-F: DGESVD failed to converge')
 1010 format (
     &     t3, '@OEPFIACHK-F: Illegal argument to DGESVD, argument #',
     &     i2)
 1020 format (
     &     t3, '@OEPFIACHK-I: ',i2,' eigenvalues of the auxiliary ',/,
     &     t3, 'hole-particle MEs are below the threshold of ',
     &     d7.1, /)
 1025 format (
     &     t3,'@OEPFIACHK-I: The smallest eigenvalues of the ',/,
     &     t3,'auxiliary hole-particle MEs is',
     &     f20.10, /)
 1050 format (
     &     t3,'@OEPFIACHK-I: The largest difference between non-',/, 
     &     t3,'local and local exchange ME:', /,
     &     5(t6,a2,1x,f20.10,/))
 1053 format (
     &     t3,'@OEPFIACHK-I: The average difference between non-',/, 
     &     t3,'local and local exchange ME:', /,
     &     5(t6,a2,1x,f20.10,/))
 1055 format (
     &     t3,'@OEPFIACHK-I: The largest difference between non-',/,
     &     t3,'local and local exchange weighted with the energy',/,
     &     t3,'denominator:',/
     &     3(t6,a2,1x,f20.10,/))
 1058 format (
     &     t3,'@OEPFIACHK-I: The average difference between non-',/,
     &     t3,'local and local exchange weighted with the energy',/,
     &     t3,'denominator:',/
     &     3(t6,a2,1x,f20.10,/))
 1060 format (
     &     t3,'@OEPFIACHK-I: The contribution of singles to the ',/,
     &     t3,'second-order correlation energy is',
     &     f20.10, /)
c
c     Set the SVD threshold (later it will be set by a keyword)
c
      svdthr = 1.d-8
c
c     Initialize the top of the core at the begging of the ICORE array
c     (it has been already been shifted in VSCF.)
c
      coretop = 1
c
c     Transform auxiliary AOMEs to MOMEs:
c       1/ allocate memory for MO MEs and a scratch array;
c       2/ check if there is enough memory available;
c       3/ transform auxilary AO MEs to occ.-virt. MEs.
c
      auxhp = coretop
      coretop = auxhp + iintfp * nocc * (nbas - nocc) * naux
      i010 = coretop
      i020 = i010 + iintfp * nbas * nbas
      itop = i020 + iintfp * nocc * nbas
c
      if (itop .gt. maxmem) call insmem (
     &     'OEPFIACHK-F: Transformation of auxiliary b.s. AO MEs',
     &     itop+usdmem, maxmem+usdmem)
c
      do x=1, naux
         call oepao2mo (
     &        nbas, nocc, 'IA  ', scfevc, auxao (1, x), 
     &        icore (i010), icore (i020),
     &        icore (auxhp+iintfp*nocc*(nbas-nocc)*(x-1)))
      enddo
c
c     Perform the SVD procedure on auxiliary MOMEs and find the smallest eigenvalue.
c     This part is very similar to what is done in OEPNSITER:
c      1/ figure out what is the least dimension of the auxiliary MOMEs;
c      2/ allocate memory for the local potential auxiliary expansion, SVD
c         eigenvalues, left and right eigenvectors, and a scratch array;
c      3/ check if there is enough memory available;
c      4/ find SVD of the auxiliary MOMEs;
c      5/ check if SVD converged
c      6/ find the smallest eigenvalue and print it
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
     &     'OEPFIACHK-F: Singular Value Decomp. of auxiliary MOMEs',
     &     itop+usdmem, maxmem+usdmem)
c
      call dgesvd (
     &     's', 's', nocc*(nbas-nocc), naux, 
     &     icore (auxhp), nocc*(nbas-nocc), icore (i010),
     &     icore (i020), nocc*(nbas-nocc), icore (i030), nsvd,
     &     icore (i040), 5*nsvd + max (naux, nocc*(nbas-nocc)),
     &     itop, ierr)
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
      fiachk (1)  = 1.d99
      do i=1, nsvd
         call scopy (1, icore (i010 + iintfp *(i-1)), 1, dsvd, 1)
         if ( dabs (dsvd) .lt. svdthr ) ndep = ndep + 1
         if ( dabs (dsvd) .lt. fiachk (1)) fiachk (1) = dsvd
      enddo
c
      write (6, 1020) ndep, svdthr
      write (6, 1025) fiachk (1)
c
c     Calculate the differences between n.-l. exch. and local exchange
c     MOMEs, and print the largest differences:
c       1/ allocate memory for the OEP solutions and non-local/local
c          difference MOMEs (various block);
c       2/ check if there is enough memory available;
c       3/ read the OEP solution from JOBARC;
c       4/ contract it with auxiliar AOMEs and substract from n.-l.
c          exch. AOMEs;
c       5/ transform these difference to MOs (block by block);
c       6/ copy the diagonals of HH and PP blocks to separate arrays
c          and zero them out to exclude from comparison;
c       7/ find the largest difference in all the block and diagonals
c          and print them;
c       8/ scale differences by energy denominators, block by block,
c          excepts the diagonals;
c       9/ find the largest weighted difference and print them;
c      10/ recalculate occ.-virt. MOMEs and contract them with scaled
c          ones to get  the singles contribution to the second-order
c          energy;
c$$$   11/ since we substracted local exchange AOMEs from n.-l. exch.
c$$$       AOMEs, and the latter would be used later in the post-SCF
c$$$       OEP block, add them back in.
c
c     Memory pointers:
c       i010 - OEP solution in a.b.s. representation;
c       i020 - difference MOMEs (occ.-virt. part);
c       i030 - difference MOMEs (occ.-occ. part);
c       i040 - difference MOMEs (virt.-virt. part);
c       i050, i060 - diagonals of i030, i040
c       i070, i080 - scratch arrays need for AO to MO transformation
c       i040 will be used again to calculate E^(2)_S
c
      i010 = coretop
      i020 = i010 + iintfp * naux
      i030 = i020 + iintfp * nocc * (nbas-nocc)
      i040 = i030 + iintfp * nocc * nocc
      i050 = i040 + iintfp * (nbas-nocc) * (nbas-nocc)
      i060 = i050 + iintfp * nocc
      i070 = i060 + iintfp * (nbas-nocc)
      i080 = i070 + iintfp * nbas * nbas
      itop = i080 + iintfp * nbas * nbas
c     
      if (itop .gt. maxmem) call insmem (
     &     'OEPFIACHK-F: Calculating non-local minus local MOMEs',
     &     itop+usdmem, maxmem+usdmem)
c
      call getrec (
     &     20, 'JOBARC', 'OEPAU', iintfp*naux, icore (i010))
c     
      call xgemm(
     &     'n', 'n', itriln (nirrep+1), 1, naux,
     &     -1.d0, auxao, itriln (nirrep+1),
     &     icore (i010), naux,
     &     1.d0, hfxao, itriln (nirrep+1))
c
      call oepao2mo (
     &     nbas, nocc, 'IA  ',
     &     scfevc, hfxao, icore (i070), icore (i080), icore (i020))
      call oepao2mo (
     &     nbas, nocc, 'IJ  ',
     &     scfevc, hfxao, icore (i070), icore (i080), icore (i030))
      call oepao2mo (
     &     nbas, nocc, 'AB  ',
     &     scfevc, hfxao, icore (i070), icore (i080), icore (i040))
c
      call zero (icore (i050), nocc)
      call zero (icore (i060), nbas-nocc)
      call dswap (
     &     nocc,
     &     icore (i030), nocc+1,
     &     icore (i050), 1)
      call dswap (
     &     nbas-nocc,
     &     icore (i040), nbas-nocc+1,
     &     icore (i060), 1)
c     
      call scopy (
     &     1, icore (i020 - iintfp + 
     &     iintfp*idamax (nocc*(nbas-nocc), icore (i020), 1)), 1,
     &     fiachk(1), 1)
      call scopy (
     &     1, icore (i030 - iintfp + 
     &     iintfp*idamax (nocc**2, icore (i030), 1)), 1,
     &     fiachk(2), 1)
      call scopy (
     &     1, icore (i040 - iintfp + 
     &     iintfp*idamax ((nbas-nocc)**2, icore (i040), 1)), 1,
     &     fiachk(3), 1)
      call scopy (
     &     1, icore (i050 - iintfp + 
     &     iintfp*idamax (nocc, icore (i050), 1)), 1,
     &     fiachk(4), 1)
      call scopy (
     &     1, icore (i060 - iintfp + 
     &     iintfp*idamax ((nbas-nocc), icore (i060), 1)), 1,
     &     fiachk(5), 1)
c     
      write (6, 1050)
     &     'IA', fiachk (1), 'IJ', fiachk (2), 'AB', fiachk (3),
     &     'II', fiachk (4), 'AA', fiachk (5)
c     
      fiachk (1) =
     &     dasum (nocc*(nbas-nocc), icore (i020), 1) /
     &     dble (nocc*(nbas-nocc))
      fiachk (2) =
     &     dasum (nocc*nocc, icore (i030), 1) /
     &     dble (nocc*(nocc-1))
      fiachk (3) =
     &     dasum ((nbas-nocc)**2, icore (i040), 1) /
     &     dble ((nbas-nocc)*(nbas-nocc-1))
      fiachk (4) =
     &     dasum (nocc, icore (i050), 1) /
     &     dble (nocc)
      fiachk (5) =
     &     dasum (nbas-nocc, icore (i060), 1) /
     &     dble (nbas-nocc)
      write (6,1053)
     &     'IA', fiachk (1), 'IJ', fiachk (2), 'AB', fiachk (3),
     &     'II', fiachk (4), 'AA', fiachk (5)
c     
      call oepenrscl (
     &     nbas, nocc, 'IA  ', scfevl, icore (i020))
      call oepenrscl (
     &     nbas, nocc, 'IJ  ', scfevl, icore (i030))
      call oepenrscl (
     &     nbas, nocc, 'AB  ', scfevl, icore (i040))
      call scopy (
     &     1, icore (i020 - iintfp + 
     &     iintfp*idamax (nocc*(nbas-nocc), icore (i020), 1)), 1,
     &     fiachk(1), 1)
      call scopy (
     &     1, icore (i030 - iintfp + 
     &     iintfp*idamax (nocc**2, icore (i030), 1)), 1,
     &     fiachk(2), 1)
      call scopy (
     &     1, icore (i040 - iintfp + 
     &     iintfp*idamax ((nbas-nocc)**2, icore (i040), 1)), 1,
     &     fiachk(3), 1)
      write (6, 1055)
     &     'IA', fiachk (1), 'IJ', fiachk (2), 'AB', fiachk (3)
c     
      fiachk (1) =
     &     dasum (nocc*(nbas-nocc), icore (i020), 1) /
     &     dble (nocc*(nbas-nocc))
      fiachk (2) =
     &     dasum (nocc*nocc, icore (i030), 1) /
     &     dble (nocc*(nocc-1))
      fiachk (3) =
     &     dasum ((nbas-nocc)**2, icore (i040), 1) /
     &     dble ((nbas-nocc)*(nbas-nocc-1))
      write (6, 1058)
     &     'IA', fiachk (1), 'IJ', fiachk (2), 'AB', fiachk (3)
c
      call oepao2mo (
     &     nbas, nocc, 'IA  ',
     &     scfevc, hfxao, icore (i070), icore (i080), icore (i040))
      fiachk (1) = ddot (
     &     nocc*(nbas-nocc),
     &     icore (i020), 1,
     &     icore (i040), 1)
      write (6, 1060) fiachk (1)
c     
      return
      end
      
