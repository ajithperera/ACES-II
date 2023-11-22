      subroutine oep_dowork(iter, iuhf, nbas, nocc, icore, maxmem,  
     &                      usdmem, scfh1e, scfevl, scfevc, scfden, 
     &                      scfham, scfenr, dmax, itertype, 
     &                      transfer,slat,Zslater,islaterON)
c
      implicit none
c
#include <oep.com>
#include <machsp.com>
#include <flags.com>
#include <flags2.com>
#include <symm2.com>
c 
      character*3 itertype
c
c     The input parameters
c
      logical addoep2, transfer
c 
      integer
     &     nbas, nocc(16), iuhf,
     &     icore (*) 

      logical slat,Zslater
      integer islaterON
      double precision
     &     scfh1e (itriln (nirrep+1)),
     &     scfevl (nbas), scfevc (nbas, nbas),
     &     scfden (itriln (nirrep+1)),
     &     scfham (itriln (nirrep+1)), scfenr
c
c     Local variables
c
      integer
     &     i000, i010, coretop, usdmem, maxmem, iter, luout,
     &     auxaa, auxph, hfxaa, hfxph, 
     &     resinv, oep1au, oep2au, refevl
c
      integer lnbuf, luint
      parameter (lnbuf=600)
      parameter (luint=10, luout=6)
c
      double precision homo, OEPTOL, au2ev, dmax, tol, etot
c
c     Define the one- and two-electron integrals array sizes.
c     
#define DEFINE_OPERATORS
C#include <oepintind.com>
#define DEFINE_VARIABLES
#define N_BAS nbas
#define N_OCC nocc
C#include <oepintind.com>
c     
c     OEP SCF block, Igor Schweigert, Sep 2003,
c     $Id: vscfoep.FPP,v 1.13 2008/06/06 18:10:43 taube Exp $
c     
c     A SCF iteration starts with constructing the Fock matrix. If the
c     OEP method is used, we will solve the OEP equation, build the OEP
c     Hamiltonian matrix and bypass the standard Fock matrix
c     construction and the total energy calculation.
c     
c     Note that we need some initial guess for eigenvalues and
c     eigenvectors. The INITGES routine does not store the 
c     eigenvectors. So if INITGES is used let's perform the HF SCF
c     iteration once as a guess for eigenvectors.
c     
c     Note also that we bypass the "reaction field" contributions, which
c     is actually the first thing done in SCF iterations. The reason for
c     that is if the OEP block is the first in a SCF iteration then we
c     still have eigenvalues and eigenvectors in ICORE (I050) and ICORE
c     (I060) correspondingly after SCFIT. So we don't need to create
c     extra arrays to store them.
c     
c     Note also that all the incompatabilities have already been checked
c     in the OEP preSCF block. So the fact that we get here means that
c     this is RHF calculation (IUHF=0), the VMOL integral package is
c     used (PETITE_LIST=.FALSE.), and the Fock matrix is build out of
c     core (AOFIL=.TRUE.)
c     
c     The OEP SCF block consists of the following steps: <> Figure out
c     what kind of iteration is the current one (all the games with the
c     convergence should be played here.) <> If it is a OEP1 or OEP2
c     iteration: <<>> Set the top of the core behind the SCF
c     eigenvectors. <<>> Solve the corresponding OEP equations and build
c     the OEP Hamiltonian. <<>> Report the iteration type, total energy,
c     and density difference from the previous iteration. <<>> Remember
c     the type of the iteration just performed. <<>> Proceed with the
c     rest of the SCF iteration by jumping to DMPSCF.
C
      r12on=.false.
      call oeptimer ('VSCF: SCFITER', 'ON')

      if(zslater) then 
        if( dmax .le. 10.d0**(-islaterON) .and. iter .gt. 1) then
          if( slat ) then
             call slater
             slat=.false.
          end if 
       end if
      end if
C
      if( .not. slat) then
      if (OEP1.or.OEP2.and..not.OEP1ON) then
         if (dmax.le.10.d0**(-OEP1ONTHR).and.iter.gt.1)
     &        OEP1ON = .TRUE.
      end if
       if (OEP2.and..not.OEP2ON) then
         if (dmax.le.10.d0**(-OEP2ONTHR).and.iter.gt.1) then
            OEP2ON = .TRUE.
            OEP2ITR = .TRUE.
         endif 
      elseif (OEP2.and.OEP2ON) then
         OEP2ITR = .FALSE.
         if (dmax.le.10.d0**(-OEP2ALWAYSTHR))
     &        OEP2ALWAYS = .TRUE.
         if (dmax.le.10.d0**(-OEP2OFFTHR).and..not.OEP2ALWAYS)
     &        OEP2UPDATE = 9999
         if (OEP2ALWAYS.or.iter - (iter/OEP2UPDATE)*OEP2UPDATE .eq. 0)
     &        OEP2ITR = .TRUE.
         if (dmax.le.0.0000) then
            tol = OEPTOL
            OEP2ITR = .TRUE. 
         endif
      endif
c
      i000 = 1

      if (OEP1ON .or. OEP2ON) then
           
         i010 = i000 + (iuhf+1) * isqrln (nirrep+1) * iintfp
         if (.not.OEP2ITR) then
            call oep1iter (
     &           nbas, nocc (1),  
     &           icore (i000), maxmem, usdmem,  
     &           scfh1e, scfevl, scfevc, scfden, scfham, scfenr)
         else
            call oep2iter (
     &           nbas, nocc (1), 
     &           icore (i000), maxmem, usdmem, 
     &           scfh1e, scfevl, scfevc, scfden, scfham, scfenr,
     &           r12on)
         endif
         etot = scfenr
         if (.not.OEP2ON) then
            itertype = 'x  '
         elseif (OEP2ITR) then
            itertype = 'xc '
         else
            itertype = 'xc'//char (39)
         endif
         write (luout, 5950) iter-1, itertype, etot, dmax 
C
         transfer = .true.
      endif 
      end if        
 5950    format (
     &        t3,i6,1x,a3,1x,f20.10,10x,d20.10)
      
      return
      end

