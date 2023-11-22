      subroutine oep1iter (
     &     nbas, nocc, icore, maxmem, usdmem,
     &     scfh1e, scfevl, scfevc, scfden, scfham, scfenr)
c
c This routine solves the first-order Optimized Effective Potential
c equations to get the local exchange OEP.
c
c The OEP solution will be obtained as an expansion over auxiliary basis
c functions. The OEP local exchange AOMEs will be then calculated and
c added to core-Hamiltonian AOMEs (supplied in SCFH1E) and Coulomb AOMEs
c to form the SCF Hamiltonian (passed back through SCFHAM). The
c corresponding OEP energy through first-order is passed back through
c SCFENR.
c
c The ADDOEP2 flag tells whether the correlation potential (read in from
c JOBARC) will be added to the exchange potential. This is used to
c bypass the updating of the second-order correalation potential during
c the second-order OEP iterations.
c 
c Since we need a lot of different scratch arrays to solve the OEP
c equation, the routine itself allocates and releases memory, using
c available core passed in through ICORE. The amount of memory available
c is given by MAXMEM. That is all the routine needs to know, however, in
c a case of insufficient memory, we need to let user know how much total
c memory is needed, so the amount of memory already used in VSCF is
c given by USDMEM.
c
c In:  nbas, nocc, naux, maxmem, usdmem, addoep2
c      scfh1e, scfevl, scfevc, scfden, scfenr, scfham
c Scr: icore
c Out: scfenr and scfham are modified
c
c Igor Schweigert, Aug 2003 
c $Id: oep1iter.FPP,v 1.1.1.1 2008/05/29 19:35:39 taube Exp $
c
      implicit none
c
c     Common block
c
#include <symm2.com>
#include <machsp.com>
#include "oep.com"
c
c     The input parameters
c
      logical addoep2
c 
      integer
     &     nbas, nocc, 
     &     icore (*), maxmem, usdmem
c
      double precision
     &     scfh1e (itriln (nirrep+1)),
     &     scfevl (nbas), scfevc (nbas, nbas),
     &     scfden (itriln (nirrep+1)),
     &     scfham (itriln (nirrep+1)), scfenr
c
c     Local variables
c
      integer
     &     i010, i020, i030, i040, i050, i060, itop, coretop,
     &     i070, i080,
     &     auxaa, auxph, hfxaa, hfxph, 
     &     resinv, oep1au, oep2au, refevl
c
      integer lnbuf, luint
      parameter (lnbuf=600)
      parameter (luint=10)
c
      double precision homo
c
c     Define the one- and two-electron integrals array sizes.
c     
#define DEFINE_OPERATORS
#include <oepintind.com>
#define DEFINE_VARIABLES
#define N_BAS nbas
#define N_OCC nocc
#include <oepintind.com>
c
c     Initialize the top of the core at the begging of the ICORE array
c     (it has been already been shifted in VSCF.)
c
      coretop = 1
c
c     Calculate all the MO integrals (auxiliary, Coulomb, and non-local
c     exchange) needed to solve the OEP1 equation: <> Allocate memory
c     for MO integrals and the scratch arrays necessary to read and
c     transform the AO integrals. <> Check if there is enough memory
c     available. <> If the semicanonicals orbitals are requested, build
c     the Fock matrix and rotate the orbitals to semicanonical. <>
c     Initialize the SCF Hamiltonian with the core Hamiltonian. <> Read
c     and transform AOMEs.
c
c     Memory pointers legend: i010 - AH intermediate, i050 - integral
c     buffer, i060 - index buffer.
c
c
      auxaa = coretop
      auxph = auxaa + iintfp * n_aa * naux
      hfxaa = auxph + iintfp * n_ph * naux 
      hfxph = hfxaa + iintfp * n_aa
      coretop = hfxph + iintfp * n_ph
      i010  = coretop
      i050  = i010  + iintfp * n_ah
      i060  = i050  + iintfp * lnbuf
      itop  = i060  + lnbuf
      if (OEPH0 (1:4) .eq. 'SEMI') then
         i070 = itop + mod (itop-1, iintfp)
         i080 = i070 + iintfp * nbas * nbas
         itop = i080 + iintfp * nbas * nbas
      endif
c     
      if (itop .gt. maxmem) call insmem (
     &     'OEP1ITER-F: Transformation of the auxiliary and '// 
     &     'two-electron AO integrals',
     &     itop-1+usdmem, maxmem+usdmem)
c     
      if (OEPH0 (1:4) .eq. 'SEMI') then
         call oepsemcan (
     &        nbas, nocc, luint, lnbuf, icore (i060),
     &        icore (i050), scfevc, scfden, scfh1e,
     &        icore (hfxaa), icore (auxaa), icore (auxph),
     &        icore (i070), icore (i080))
      endif
c     
      call dcopy (n_aa, scfh1e, 1, scfham, 1)
c
      call oep1ints (
     &     nbas, nocc, naux, luint, lnbuf, OEPAXP, icore (i060), 
     &     icore (i050), scfevc, scfden, scfham, icore (i010), 
     &     icore (auxaa), icore (auxph), 
     &     icore (hfxaa), icore (hfxph))
c
c     Calculate the eigenvalues of the reference Hamiltonian: <>
c     Allocate memory for the eigenvalues. <> Check if there is enough
c     memory available. <> Calculate the reference Hamiltonian.
c     
      refevl = coretop
      coretop = refevl + iintfp * nbas
c     
      if (itop .gt. maxmem) call insmem (
     &     'OEP1ITER-F: Calculating the reference Hamiltonian.',
     &     itop-1+usdmem, maxmem+usdmem)
c
      call oeprefham (
     &     nbas, nocc, naux, OEPH0,
     &     scfevl, scfevc, scfham,
     &     icore (hfxaa), icore (refevl))
c     
c     Calculate the inverse of the response function: <> Allocate memory
c     for the inverse and the scratch array. <> Check if there is enough
c     memory available. <> Construct the response function and calculate
c     its inverse.
c     
      resinv = coretop
      coretop = resinv + iintfp * naux * naux
      i010 = coretop
      i020 = i010 + iintfp * naux
      itop = i020 + iintfp * naux * naux
c     
      if (itop .gt. maxmem) call insmem (
     &     'OEPITER-F: Calculating the response inverse using SVD.',
     &     Itop-1+usdmem, maxmem+usdmem)
c     
      call oepresinv (
     &     nbas, nocc, naux, OEPSVDTHR,
     &     icore (refevl), icore (auxph), 
     &     icore (i010), icore (i020), icore (resinv))
c
c     Solve the first-order OEP equation: <> Allocate memory for OEP1
c     and the scratch arrays. <> Check if there is enough memory. <>
c     Calculate the first-order energy and write it to JOBARC. <>
c     Construct the local exchange potential by calculating the
c     first-order RHS and contracting it with the response inverse. <>
c     Calculate the violation of the OEPX HOMO-condition.  <> Scale the
c     local exchange by 1 - the fraction of the nonlocal exchange. <>
c     Store the (scaled!) OEPX in the auxiliary representation to
c     JOBARC. <> If requested, read the second-order OEP solution from
c     JOBARC and add it the first-order solution. <> Store the total
c     potential into JOBARC. (Note that this record depends on the
c     fraction of the nonlocal exchange.) <> Calculate the OEP AOMEs and
c     add them to the SCF Hamiltonian. <> Add the fraction of nonlocal
c     exchange AOMEs to the SCF Hamiltonian.
c     
      oep1au = coretop
      oep2au = oep1au + iintfp * naux
      coretop = oep2au + iintfp * naux
      i010 = coretop
      itop = i010 + iintfp * naux
c
      if (itop .gt. maxmem) call insmem (
     &     'OEP1ITER-F: Solving the first-order OEP equation.',
     &     itop-1+usdmem, maxmem+usdmem)
c
      call getrec (20, 'JOBARC', 'NUCREP  ', iintfp, scfenr)
      call oep1enr (nbas, scfden, scfh1e, scfham, icore (hfxaa), scfenr)
      call putrec (20, 'JOBARC', 'OEP1ENER', iintfp, scfenr)
c     
      call oep1rhs (
     &     nbas, nocc, naux,
     &     icore (refevl), icore (auxph), icore (hfxph), icore(i010))
      call xgemm(
     &     'n', 'n', naux, 1, naux,
     &     1.d0, icore (resinv), naux,
     &     icore (i010), naux,
     &     0.d0, icore (oep1au), naux)
c     
      homo = 0.d0
      call oep1homo (
     &     nbas, nocc, naux, OEP1HYBRID,
     &     scfevc, icore (hfxaa), icore (oep1au),
     &     icore (auxaa), homo)
c     
      call putrec (20, 'JOBARC', 'OEPXHOMO', iintfp, homo)
c
      call dscal (naux, 1d0-OEP1HYBRID, icore (oep1au), 1)
c
      call putrec (
     &     20, 'JOBARC', 'OEPXAU  ', iintfp*naux, icore (oep1au))
      call putrec (
     &     20, 'JOBARC', 'OEP1AU  ', iintfp*naux, icore (oep1au))
c
      if (addoep2) then
         call getrec (
     &        20, 'JOBARC', 'OEP2AU  ', iintfp*naux, icore (oep2au))
         call daxpy (naux, 1.d0, icore (oep2au), 1, icore (oep1au), 1)
      endif
c     
      call putrec (
     &     20, 'JOBARC', 'OEPAU   ', iintfp*naux, icore (oep1au))
c     
      call xgemm(
     &     'n', 'n', n_aa, 1, naux,
     &     1.d0, icore (auxaa), n_aa,
     &     icore (oep1au), naux,
     &     1.d0, scfham, n_aa)
c
      call daxpy (
     &     n_aa,
     &     OEP1HYBRID, icore (hfxaa), 1,
     &     scfham, 1)
c
c$$$      call daxpy (
c$$$     &     n_aa,
c$$$     &     1.d0, icore (hfxaa), 1,
c$$$     &     scfham, 1)
c     
      return
      end
