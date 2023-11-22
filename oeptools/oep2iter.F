      subroutine oep2iter (
     &     nbas, nocc,
     &     icore, maxmem, usdmem,
     &     scfh1e, scfevl, scfevc, scfden, scfham, scfenr, dor12)
c
c This routine solves the Optimized Effective Potential equations (in
c the first order and second order) to get the local OEP exchange and
c correlation potentials.
c
c The OEP solution will be obtained as an expansion over auxiliary basis
c functions. The OEP local exchange and/or correlation AOMEs will be
c then calculated and added to core-Hamiltonian AOMEs (supplied in
c SCFH1E) and Coulomb AOMEs to form the SCF Hamiltonian (passed back
c through SCFHAM). The corresponding OEP energy is added to the
c nuclear-nuclear repulsion energy (supplied in SCFENR) and passed back
c through SCFENR.
c
c Since we need a lot of different scratch arrays to solve the OEP
c equation, the routine itself allocates and releases memory, using
c available core passed in through ICORE. The amount of memory available
c is given by MAXMEM. That is all the routine needs to know, however, in
c a case of insufficient memory, we need to let user know how much total
c memory is needed, so the amount of memory already used in VSCF is
c given by USDMEM.
c
c In:  nbas, nocc, naux, maxmem, usdmem,
c      scfh1e, scfevl, scfevc, scfden, scfenr, scfham
c Scr: icore
c Out: scfenr and scfham are modified
c
c Igor Schweigert, Jan 2003 
c $Id: oep2iter.FPP,v 1.3 2008/06/06 18:09:17 taube Exp $
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
      integer
     &     nbas, nocc,
     &     icore (*), maxmem, usdmem
c
      double precision
     &     scfh1e (itriln (nirrep+1)),
     &     scfevl (nbas), scfevc (nbas, nbas),
     &     scfden (itriln (nirrep+1)),
     &     scfham (itriln (nirrep+1)), scfenr
     
      logical 
     &      dor12
c
c     Local variables
c
      logical
     &     diagonal
c     
      integer
     &     x, n, 
     &     i010, i020, i030, i040, i050, i060, itop, coretop,
     &     auxaa, auxpp, auxph, auxhh,
     &     hfxaa, hfxpp, hfxph, hfxhh,
     &     ippph, ipphh, iphph, iphhh,
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
c     Calculate all the MO integrals (auxiliary and two-electron) needed
c     to solve the OEP2 equation: <> Allocate memory for the auxiliary,
c     one- and two-electron MO integrals and the scratch arrays
c     necessary to read and transform the AO integrals. <> Check if
c     there is enough memory available. <> Initialize the SCF
c     Hamiltonian with the core Hamiltonian. <> Read and transform
c     AOMEs.
c
c     Memory pointers legend: i010 - PAAH, i020 - AAAH, i030 - PPAH,
c     i040 - PHAH, i050 - integral buffer, i060 - index buffer.
c
c
      auxaa = coretop
      auxpp = auxaa + iintfp * n_aa * naux
      auxph = auxpp + iintfp * n_pp * naux 
      auxhh = auxph + iintfp * n_ph * naux 
      hfxaa = auxhh + iintfp * n_hh * naux 
      hfxpp = hfxaa + iintfp * n_aa
      hfxph = hfxpp + iintfp * n_pp
      hfxhh = hfxph + iintfp * n_ph
      ippph = hfxhh + iintfp * n_hh
      ipphh = ippph + iintfp * n_ppph
      iphph = ipphh + iintfp * n_pphh
      iphhh = iphph + iintfp * n_phph
      coretop = iphhh + iintfp * n_phhh
      i010  = ippph
      i020  = i010  + iintfp * n_paah
      i030  = i020  
      i040  = i030  + iintfp * n_ppah
      i050  = i020  + iintfp * n_aaah
      i060  = i050  + iintfp * lnbuf
      itop  = i060  + lnbuf
c     
      if (itop .gt. maxmem) call insmem (
     &     'OEP2ITER-F: Transformation of the auxiliary and 
     &     two-electron AO integrals',
     &     itop-1+usdmem, maxmem+usdmem)
c     
      if (OEPH0 (1:4) .eq. 'SEMI') then
         call oepsemcan (
     &        nbas, nocc, luint, lnbuf, icore (i060),
     &        icore (i050), scfevc, scfden, scfh1e,
     &        icore (hfxaa), icore (hfxpp), icore (hfxhh),
     &        icore (i010), icore (i020))
      endif
c     
      call dcopy (n_aa, scfh1e, 1, scfham, 1)
c 
      call oep2ints (
     &     .true., nbas, nocc, naux, luint, lnbuf, OEPAXP, icore (i060), 
     &     icore (i050), scfevc, scfham, 
     &     icore (auxaa), icore (auxpp), icore (auxph), icore (auxhh),
     &     icore (hfxaa), icore (hfxpp), icore (hfxph), icore (hfxhh),
     &     icore (i010 ), icore (i020 ), icore (i030 ), icore (i040),
     &     icore (ippph), icore (ipphh), icore (iphph), icore (iphhh))
c
c$$$      print *, 'OEP2INTS: AUXPH'
c$$$      call output (
c$$$     &     icore (auxph), 1, naux*n_ph, 1, 1, naux*n_ph, 1, 1)
c$$$      print *, 'OEP2INTS: SCFH'
c$$$      call output (
c$$$     &     scfham, 1, n_aa, 1, 1, n_aa, 1, 1)
c$$$      print *, 'OEP2INTS: HFXAA'
c$$$      call output (
c$$$     &     icore (hfxaa), 1, n_aa, 1, 1, n_aa, 1, 1)
c$$$      print *, 'OEP2INTS: HFXPH'
c$$$      call output (
c$$$     &     icore (hfxph), 1, n_ph, 1, 1, n_ph, 1, 1)
c
c     Calculate the reference Hamiltonian: <> Allocate memory for the
c     eigenvalues of the reference Hamiltonian. <> Check if there is
c     enough memory available. <> Calculate the reference Hamiltonian.
c     
      refevl = coretop
      coretop = refevl + iintfp * nbas
c     
      if (itop .gt. maxmem) call insmem (
     &     'OEP2ITER-F: Calculating the reference Hamiltonian.',
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
     &     'OEP2ITER-F: Calculating the response inverse using SVD.',
     &     Itop-1+usdmem, maxmem+usdmem)
c     
      call oepresinv (
     &     nbas, nocc, naux, OEPSVDTHR,
     &     icore (refevl), icore (auxph), 
     &     icore (i010), icore (i020), icore (resinv))
c
c     Solve the first-order OEP equation: <> Allocate memory for the
c     OEP1 solution, and the scratch arrays. <> Check if there is enough
c     memory. <> Calculate the first-order energy and write it to
c     JOBARC. <> If the local exchange potential is used, calculate the
c     first-order RHS and contract it with the response inverse to get
c     the the first-order OEP solution. If the non-local exchange
c     potential is used, set the OEPX solution to zero and add the
c     non-local exchange AOMEs to the SCF Hamiltonian. <> Store the OEPX
c     in the auxiliary representation to JOBARC.
c     
      oep1au = coretop
      coretop = oep1au + iintfp * naux
      i010 = coretop
      itop = i010 + iintfp * naux
c
      if (itop .gt. maxmem) call insmem (
     &     'OEP2ITER-F: Solving the first-order OEP equation.',
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
      call putrec (20, 'JOBARC', 'OEPXHOMO', iintfp, homo)
c
c     Solve the second-order OEP equation: <> Allocate memory. <> Check
c     if there is enough memory available. <> Check if diagonal elements
c     have to be included. <> Calculate the RHS of the second-order OEP
c     equation and the second-order energy. <> Store the second-order
c     OEP energy to JOBARC. <> Contract the response inverse with the
c     RHS to get the second-order solution. <> Store it to JOBARC. <>
c     Add the first-order solution to the second-order solution. <>
c     Store it to JOBARC. <> Calculate the AOMEs of the OEP and add them
c     to the SCF Hamiltonian.
c
c     Note that the non-local exchange MOMEs are overwritten in OEP2RHS!
c
c     Memory pointer legend: i010, i020, i030 - PP, PH, HH intemediates,
c     i040 - the second-order OEP RHS.
c
      oep2au = coretop
      coretop = oep2au + iintfp * naux
      i010 = coretop
      i020 = i010 + iintfp * n_pp
      i030 = i020 + iintfp * n_ph
      i040 = i030 + iintfp * n_hh      
      i050 = i040 + iintfp * naux
      if (dor12) then
        itop = i050 + iintfp * n_ph
      else
        itop = i050
      endif  
c
      if (itop .gt. maxmem) call insmem (
     &     'OEP2ITER-F: Solving the second-order OEP equation.',
     &     itop-1+usdmem, maxmem+usdmem)
c
      diagonal = OEPH0 (1:4) .eq. 'OEP '

      call oep2rhs (
     &     nbas, nocc, naux,
     &     OEPDFF, OEPDFW, OEPDWW,
     &     OEPFDF, OEPWDW, diagonal, OEP1HYBRID,
     &     scfenr, icore (refevl), icore (oep1au),
     &     icore (auxpp), icore (auxph), icore (auxhh),
     &     icore (hfxpp), icore (hfxph), icore (hfxhh),
     &     icore (ippph), icore (ipphh), icore (iphph), icore (iphhh),
     &     icore (i010 ), icore (i020 ), icore (i030 ), icore (i050 ), 
     &     icore (i040 ), dor12)

      call putrec (20, 'JOBARC', 'OEP2ENER', iintfp, scfenr)
c     
      call xgemm(
     &     'n', 'n', naux, 1, naux,
     &     1.d0, icore (resinv), naux,
     &     icore (i040), naux,
     &     0.d0, icore (oep2au), naux)
c     
c$$$      write (6,*) 'OEP2 RHS'
c$$$      call output (
c$$$     &     icore (i040), 1, 9, 1, 1, naux, 1, 1)
c$$$      write (6,*) 'OEP2 RESINV'
c$$$      do n=1, naux
c$$$         write (6,*)
c$$$     &        n,
c$$$     &        deref (1+naux*(n-1), icore (resinv)),
c$$$     &        deref (n, icore (i040))
c$$$      enddo
c$$$      write (6,*) 'OEP2 SOLUTION'
c$$$      call output (
c$$$     &     icore (oep2au), 1, 9, 1, 1, naux, 1, 1)
c
      call dscal (naux, 1d0-OEP1HYBRID, icore (oep1au), 1)
c     
      call putrec (
     &     20, 'JOBARC', 'OEP1AU  ', iintfp*naux, icore (oep1au))
      call putrec (
     &     20, 'JOBARC', 'OEPXAU  ', iintfp*naux, icore (oep1au))
      call putrec (
     &     20, 'JOBARC', 'OEP2AU  ', iintfp*naux, icore (oep2au))
c
      call daxpy (naux, 1.d0, icore (oep1au), 1, icore (oep2au), 1)
c     
      call putrec (
     &     20, 'JOBARC', 'OEPAU   ', iintfp*naux, icore (oep2au))
c     
      call xgemm(
     &     'n', 'n', n_aa, 1, naux,
     &     1.d0, icore (auxaa), n_aa,
     &     icore (oep2au), naux,
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
      
      function deref (n, a)
      implicit none
      integer n
      double precision a (*), deref
      deref = a (n)
      return
      end

      subroutine ao2mo (n, c, x, y, z)
      implicit none
      integer n
      double precision c (*), x (*), y (*), z (*)
      call expnd2 (x, z, n)
      call xgemm (
     &     't', 'n', n, n, n,
     &     1.d0, c, n,
     &     z, n,
     &     0.d0, y, n)
      call xgemm (
     &     'n', 'n', n, n, n,
     &     1.d0, y, n,
     &     c, n,
     &     0.d0, z, n)
      return
      end

c$$$      if (OEPH0 (1:4) .eq. 'SEMI') then
c$$$         i010 = coretop
c$$$         i020 = i010 + iintfp * n_aa
c$$$         i030 = i020 + iintfp * n_ma
c$$$         i040 = i030 + iintfp * n_pp
c$$$         i050 = i040 + iintfp * n_ph
c$$$         i060 = i050 + iintfp * n_hh
c$$$         call dcopy (n_aa, scfham, 1, icore (i010), 1)
c$$$         call daxpy (n_aa, 1.d0, icore (hfxaa), 1, icore (i010), 1)
c$$$         call oepao2mo (
c$$$     &        n_a, n_h, scfevc,
c$$$     &        icore (i010), icore (i020),
c$$$     &        icore (i030), icore (i040), icore (i050))
c$$$         print *, 'OEP2ITER: SEMI PP'
c$$$         call output (icore (i030), 1, n_pp, 1, 1, n_pp, 1, 1)
c$$$         print *, 'OEP2ITER: FOCK PH'
c$$$         call output (icore (i040), 1, n_p, 1, n_h, n_p, n_h, 1)
c$$$         print *, 'OEP2ITER: FOCK HH'
c$$$         call output (icore (i050), 1, n_hh, 1, 1, n_hh, 1, 1)
c$$$      endif
