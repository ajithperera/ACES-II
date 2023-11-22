      subroutine readdens (
     &     nbas,
     &     iscf, xocc, evec, eval,
     &     iden, scraxm, scrmxm, dens)
c
c This routine reads SCF eigenvalues and eigenvectors, 
c and SCF (IDEN=0) or relaxed (IDEN>0) occupation numbers
c and density matrix
c
c Note that since VSCF_KS has its own record for occupation
c numbers and eigenvalues (if exact_ex=on), ISCF > 0 is supposed
c to make this routine to read those records, but this needs
c further verification. 
c
c Note also that occupation numbers are converted to real
c for the SCF or correspond to eigenvalues of the relaxed
c density.
c Attention here! occupation numbers are not implemented. 
c
c It is limited to RHF and SYMMETRY=OFF
c
c Igor Schweigert, Apr 2003, $Id: readdens.FPP,v 1.1.1.1 2008/05/29 19:35:40 taube Exp $
c
      implicit none
c
#include <machsp.com>
c
c     .. input parameters ..
c
      integer
     &     nbas, iscf, iden
c
      double precision
     &     xocc (nbas), eval (nbas), evec (nbas, nbas),
     &     scraxm (nbas, nbas), scrmxm (nbas, nbas), dens (nbas, nbas)
c
c     .. local variables ..
c
      integer
     &     i, j, k, iocc (nbas)
c
c     .. information messages to the user                             ..
c
 110  format (
     &     t3,'@READDENS-I, VSCF_KS record is used to retrieve ',
     &     'SCF information.')
 120  format (t3,'@READDENS-I, SCF density matrix has been read.')
 130  format (t3,'@READDENS-I, Relaxed density matrix has been read.')
c
c     .. read the SCF eigenvalues                                     .. 
c     .. attention here! KS eigenvalues may be shifted!               ..
c
      call getrec (
     &     20, 'JOBARC','SCFEVLA0', iintfp*nbas, eval)
c
c     .. read the SCF eigenvectors in VMOL AO representation          ..
c
      call getrec(
     &     20,'JOBARC', 'SCFEVCA0', iintfp*nbas*nbas, evec)
c
c     .. read SCF or relaxed occupation numbers and density matrix    ..
c
      call dzero (dens, nbas*nbas)
      call dzero (xocc, nbas)
c
      if (iden .eq. 0) then
c
c        .. read the SCF occupation numbers and convert them to real  ..
c        .. attention here! this needs further verification           ..
c
         if (iscf.eq.0) then
            call getrec (20,'JOBARC','OCCUPYA0', 1, iocc (1))
            do i=1, nbas
               xocc (i) = 0.d0
               if (i.le.iocc (1)) xocc (i) = 2.d0
            enddo
         else
            call getrec(1,'JOBARC', 'SCFOCCA', nbas, iocc)
            do i=1, nbas
               xocc (i) = 2.d0 * dble (iocc (i))
            enddo
         endif
c
c        .. build the SCF density matrix using the eigenvectors and   ..
c        .. occupation numbers                                        ..
c
         do i=1, nbas
            do j=1, nbas
               do k=1, nbas
                  dens (i,j) = dens (i,j) +
     &                 evec (i,k) * xocc (k) * evec (j,k)
               enddo
            enddo
         enddo
c
         write (6,120)
c
      else
c
c        .. read eigenvalues of the relaxed density matrix as         ..
c        .. occupation numbers                                        ..
c        .. attention here! this is not implemented                   ..
c
         call dzero (xocc, nbas)
c
c        .. read the relaxed density matrix from JOBARC               ..
c        .. (see vprops/prop.f) in the MO representation              ..
c
         call getrec (
     &        20,'JOBARC', 'RELDENSA', iintfp*nbas*nbas, scrmxm)
c     
c        .. transform it to the AO representation using SCF           ..
c        .. coefficients                                              ..
c
         call xgemm (
     &        'n', 'n', nbas, nbas, nbas,
     &        1.d0, evec, nbas, 
     &        scrmxm, nbas,
     &        0.d0, scraxm, nbas)
c     
         call xgemm (
     &        'n', 't', nbas, nbas, nbas,
     &        1.d0, scraxm, nbas,
     &        evec, nbas,
     &        0.d0, dens, nbas)
c
         write (6,130)
c
      endif
c
c     .. done ..
c     
      return
      end
