      subroutine slater
#include <aces.h>
#include <flags.com>
#include <flags2.com>
#include <iuhf.com>
#include <mol.com>
#include <kscore.com>
#include <symm2_ks.com>
#include <timeinfo.com>
#include <int.com>
#include <intgrtflags.com>
#include "machsp.com"
      integer
     &     znull, pnull, zevec, zdens, zeval, zscraxm, 
     &     zscrmxm, zscrpxz, zslat, zprm2ao, setptr, i
c
      double precision ddot
c
c     .. messages to user ..
c      
 1010 format (/,'module SLATER', /, 'development version ',
     $     '$Id: slater.FPP,v 1.1.1.1 2009/06/29 19:35:40 verma Exp $')
 1013 format (
     $     t3, '@SLATER-F, The Slater potential for an unrestricted ',
     $     'reference is not implemented yet.')
 1014 format (
     $     t3, '@SLATER-F, The Slater potential for any symmetry other ',
     $     'than C1 is not implemented yet.')
 1020 format (t3, '@SLATER-I, The Slater potential integrates to',
     $     f20.12) 
 1030 format (t3, '@SLATER-I, AO matrix elements of the Slater ',
     $     'potential has been stored into file SLAT51AO',/)
c$$$ 1030 format (t3, '@SLATER-I, AO matrix elements of the Slater ',
c$$$     $     'potential has been stored in JOBARC',/)
c
c     .. print out the banner ..
c
      write (6,1010)
c
c     .. initialize the INTGRT constants ..
c
      call callstack_init('SLATER')
      int_ks = .false.
c      call initintgrt(kshf,intnumradpts)
      if (iflags(11).ge.1) iuhf=1
      call getrec(-1,'JOBARC','COMPNIRR',1,nirrep)
      call sb_com_symm2_ks 
      call sb_com_mol
      needbas=1

      znull = setptr (1, F_REAL, 1)
      pnull = setptr (1, F_INTEGER, 1)
      zdens = setptr (1, F_REAL, nbastot*nbastot)
      zevec = setptr (1, F_REAL, nbastot*nbastot)
      zeval = setptr (1, F_REAL, nbastot)
      zxocc = setptr (1, F_REAL, nbastot)
      zslat = setptr (1, F_REAL, nbastot*nbastot)
c
      zprm2ao = setptr (1, F_REAL, totprim*nbastot)
      zscrpxz = setptr (1, F_REAL, totprim*naobasfn)
      zscraxm = setptr (1, F_REAL, naobasfn*nbastot)
      zscrmxm = setptr (1, F_REAL, nbastot*nbastot)
c
c     .. setup the grid for numerical integration (INTGRT/setuprd.F) ..
      call setupgrd
c
c     .. primitive coefficients (read in LIBRA3/A3PRIM) are for the   ..
c     .. ZMAT AO's, transform them to symmetry adapted (VMOL,         ..
c     .. "computational") AO primitive coefficients                   ..
      call prm2vmol (
     $     natoms, totprim, naobasfn, nbastot,
     $     maxshell, kscore (pnangatom), kscore (pnshellao), iintfp,
     $     dcore (zprimcoef), dcore (zscrpxz), dcore (zscraxm),
     $     dcore (zprm2ao))

c     .. store the new primitive coefficients in the old array        ..
c     .. this will be fixed, once this code is separated form INTGRT  ..
c
      call dzero (dcore (zprimcoef), totprim*naobasfn)
      call dcopy (
     $     totprim*nbastot,
     $     dcore (zprm2ao), 1,
     $     dcore (zprimcoef), 1)
c
c     .. read the density matrix; 
c     .. note that if iflag2 (253) = SCF_TYPE >0 then KS SCF record   ..
c     .. will be used; if iflag (2) = CALC > 0 then relaxed density   ..
c     .. will be read                                                 ..
c     .. for now we always use VSCF to obtain the Slater potential    ..
c     .. otherwise two zmat's have to be used                         ..
c
      call readdens (
     &     nbastot, 0,
     &     dcore (zxocc), dcore (zevec), dcore (zeval), iflags (2),
     &     dcore (zscraxm), dcore (zscrmxm), dcore (zdens))
c     
c     .. calculate the AO matrix element of the Slater potential ..
c
      call dzero ( dcore (zslat), nbastot*nbastot)
      call numintslat (
     &     dcore (zdens), dcore (zslat))
c
c     .. print out the resulting exchange energy ..
c
      write (6,1020) .5d0 * ddot ( 
     &     nbastot*nbastot,
     &     dcore (zslat), 1,
     &     dcore (zdens), 1) 
c
c     .. pack the Slater AO ME (into the triangular form) and store   ..
c     .. them to the JOBARC file                                      ..
c
c     .. ZDENS is not used anymore, so let's store the packed AO there..
c
      call mat_trans (
     $     MAT_FULL, MAT_TRI, dcore (zslat), dcore (zdens), 0)
c
c      call triang (dcore (zslat), dcore (zdens), nbastot, 'y')
      call putrec (1,'JOBARC',
     &     'SLAT51AO', iintfp*nbastot*(nbastot+1)/2, dcore (zdens))
c
c     .. There are some problems with running VSCF and VSCF_KS        ..
c     .. without cleaning out the working directory in between.       ..
c     .. So instead of writing Slater AO ME to JOBARC let's           ..
c     .. store it into a dedicated file which will be carried         ..
c     .. over from VSCF+VSLATER calculation to VSCF_KS calculation    ..
c
c      open (
cx     $     unit=33, file='SLAT51AB', form='unformatted', status='new')
c
c      write (33) (dcore (zdens+i-1), i=1, nbastot*(nbastot+1)/2)
c
c      close (unit=33, status='keep')
c
      write (6,1030)
c
c     .. done ..
c         call aces_ja_fin
      call callstack_term
      end
