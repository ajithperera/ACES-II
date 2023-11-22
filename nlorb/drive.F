#include <flags.h>
      subroutine drive(icore,jjjcrtop,iuhf)
      implicit none

      integer beta,icore(jjjcrtop), memscr, iiicrtop, iuhf, jjjcrtop

      integer iii00000, iiicoefs, iiidtran, iiiinvtr,
     & iiiovrlp, iiidenmx, iiidscrs, iiiachrg, iiishell,
     & iiinbasa, iiinchoo, iiichoos, iiicoord, iiimomto, iiicbfto,
     & iiinshel, iiifockm, iiiordvc, iiiiodvc, iiievals, Inorm, Iscr,
     & Iscr1, Iend, Iscr2, Length, NBO59

      integer jjj00000, jjjcbfnt, jjjmomnt, jjjcbfto, jjjmomto,
     & jjjsjunk, jjjcntsc, jjjcnts2, jjjiscrs, jjjtrans, jjjrowss,
     & jjjscrsp, jjjpivot

      integer nbas, nbastot, nelec, natoms, nocc(2), nvir(2), nreal
      
      integer ione

      integer bondsize, mxchoose, mxnal, largel, angmax, switch12

      double precision a2u, enuc, maxocc, eelec, det, acesener

      double precision occthresh, virthresh

      logical ryd2hyb, mjump, spheric, onlynao, onlynho,
     & onlynbo, onlynlmo, verprn, putaces, setchoose

      logical creatlist, yesno

      parameter (
     & memscr   = 1000000,
     & angmax   = 20,
     & ryd2hyb  = .false.,
     & mjump    = .true.,
     & spheric  = .true.,
!     & onlynao  = .false.,
!     & onlynho  = .false.,
!     & onlynbo  = .false.,
!     & onlynlmo = .true.,
     & verprn   = .false.,
C     & bondsize = 10,
     & mxchoose = 1,
     & a2u      = 1.88972599D0,
     & putaces  = .true.,
     & setchoose= .false.,
     & creatlist= .false.,
C     & occthresh=0.955D0,
C     & virthresh=0.996D0
     & ione=1)

#include "flags.com"
#include "machsp.com"


	onlynao = .false.
	onlynho = .false.
	onlynbo = .false.
	onlynlmo = .true.



      iiicrtop = jjjcrtop - memscr
      jjj00000 = iiicrtop
      iii00000 = 1

      write(*,*)
      write(*,*)
      write(*,43)
      write(*,44)
      write(*,45)
      write(*,43)
      write(*,*)
      write(*,*)
      call getrec(20,'JOBARC','NREALATM',ione,nreal)
C
C The bondsize is the maximum number of atoms that can form bond. Set the
C default as number of real atom and if one wants to change it give the
C bondsize file option.

      open(unit=100,file='bondsize')
         read(unit=100,fmt=*) bondsize
      close(100)
C      open(unit=200,file='threshold')
C         read(unit=200,fmt=*) occthresh
C         read(unit=200,fmt=*) virthresh
C      close(200)
C      write(*,*) 'occ threshold = ', occthresh
C      write(*,*) 'vir threshold = ', virthresh
 
      call getrec(0,'JOBARC','12SWITCH',nbas,switch12)
      if (nbas.eq.-1) then
         nbas = 0
         go to 950
      end if
      nbas = 0
      call getrec(20,'JOBARC','12SWITCH',1,switch12)
      if (switch12.eq.1) then
         write(*,*)
         write(*,40)
         write(*,41)
         write(*,42)
         write(*,*)
         stop
      end if

 950  continue

      call initdim(nbas,nbastot,nocc,nvir,
     & natoms,nelec)
 
      iiicoord = iii00000
      call ckmem(iiicoord,iiicrtop)
      
      iiiachrg = iiicoord+natoms*3*iintfp
      call ckmem(iiiachrg,iiicrtop)

      jjjcbfnt = jjj00000
      call ckmem(jjjcbfnt,jjjcrtop)

      jjjmomnt = jjjcbfnt+nbas
      call ckmem(jjjmomnt,jjjcrtop)

      call initize(icore(iiicoord),icore(iiiachrg),
     & icore(jjjmomnt),icore(jjjcbfnt),natoms,nbas)

      jjjcbfto = jjjmomnt+nbas
      call ckmem(jjjcbfto,jjjcrtop)

      jjjmomto = jjjcbfto+nbas
      call ckmem(jjjmomto,jjjcrtop)

      call iseccopy(icore(jjjcbfnt),icore(jjjcbfto),nbas)
      call iseccopy(icore(jjjmomnt),icore(jjjmomto),nbas)

      if (nbastot.ne.nbas) then
         call izerosec(icore(jjjcbfto),nbas)
         call izerosec(icore(jjjmomto),nbas)

         jjjcbfto = jjjmomnt+nbas
         call ckmem(jjjcbfto,jjjcrtop)

         jjjmomto = jjjcbfto+nbastot
         call ckmem(jjjmomto,jjjcrtop)

         jjjrowss = jjjmomto+nbastot
         call ckmem(jjjrowss,jjjcrtop)

         jjjtrans = jjjrowss+nbas-nbastot
         call ckmem(jjjtrans,jjjcrtop)

         jjjscrsp = jjjtrans+nbastot*nbas
         call ckmem(jjjscrsp,jjjcrtop)

         call doitran(icore(jjjrowss),icore(jjjtrans),
     & icore(jjjscrsp),icore(jjjcbfnt),icore(jjjmomnt),
     & icore(jjjcbfto),icore(jjjmomto),nbas,nbastot,angmax)
      end if

      iiicbfto = iiiachrg+natoms
      call ckmem(iiicbfto,iiicrtop)

      iiimomto = iiicbfto+nbastot
      call ckmem(iiimomto,iiicrtop)

      call iseccopy(icore(jjjcbfto),icore(iiicbfto),nbastot)
      call iseccopy(icore(jjjmomto),icore(iiimomto),nbastot)

      call prncoor(natoms,icore(iiicoord),a2u,icore(iiiachrg))
      call prnang(largel,nbastot,icore(iiimomto),icore(iiicbfto))
     
      iiinshel = iiimomto+nbastot
      call ckmem(iiinshel,iiicrtop)
 
      iiishell = iiinshel+natoms
      call ckmem(iiishell,iiicrtop)

      iiinbasa = iiishell+(largel+1)*natoms
      call ckmem(iiinbasa,iiicrtop)

      call izerosec(icore(jjj00000),nbas*4+nbas-nbastot+
     & nbastot*nbas+nbas*nbas)

      jjjsjunk = jjj00000
      call ckmem(jjjsjunk,jjjcrtop)

      jjjcntsc = jjjsjunk+nbastot
      call ckmem(jjjcntsc,jjjcrtop)

      jjjcnts2 = jjjcntsc+natoms*angmax
      call ckmem(jjjcnts2,jjjcrtop)

      call sortmom(icore(iiicbfto),icore(iiimomto),icore(iiinshel),
     & icore(iiishell),icore(iiinbasa),largel,mxnal,natoms,nbastot,
     & angmax,icore(jjjsjunk),icore(jjjcntsc),icore(jjjcnts2))

      iiidtran = iiinbasa+(largel+1)*natoms
      call ckmem(iiidtran,iiicrtop)

      iiiinvtr = iiidtran+nbastot*nbastot*iintfp
      call ckmem(iiiinvtr,iiicrtop)

      iiiordvc = iiiinvtr+nbastot*nbastot*iintfp
      call ckmem(iiiordvc,iiicrtop)

      iiiiodvc = iiiordvc+nbastot
      call ckmem(iiiiodvc,iiicrtop)

      call dodtrans(icore(iiidtran),icore(iiiinvtr),
     & icore(iiimomto),icore(iiicbfto),icore(iiinbasa),
     & icore(iiinshel),icore(iiishell),
     & largel,natoms,angmax,nbastot)

      call prnmat('Owntrans            ',icore(iiidtran),nbastot,
     & nbastot,.true.)
      call prnmat('Owntrans Inverse    ',icore(iiiinvtr),nbastot,
     & nbastot,.true.)

      call detord(icore(iiidtran),icore(iiiordvc),nbastot)
      call prnord(icore(iiiordvc),nbastot)
      call detord(icore(iiiinvtr),icore(iiiiodvc),nbastot)
      call prnord(icore(iiiiodvc),nbastot)

      write(*,22) 
 
      iiicoefs = iiiiodvc+nbastot
      call ckmem(iiicoefs,iiicrtop)
 
      iiiovrlp = iiicoefs+nbastot*nbastot*iintfp*(iuhf+1)
      call ckmem(iiiovrlp,iiicrtop)

      iiidenmx = iiiovrlp+nbastot*nbastot*iintfp
      call ckmem(iiidenmx,iiicrtop)

      call getrec(20,'JOBARC','SCFEVCA0',nbastot*nbastot*iintfp,
     & icore(iiicoefs))
      call prnmat('Canonical Coef Alpha',icore(iiicoefs),nbastot,
     & nbastot,.true.)

      call mkdens(icore(iiicoefs),icore(iiidenmx),icore(iiiordvc),
     & nbastot,nocc(1),.false.)
      if (iuhf.eq.0) then
         call scalmat(icore(iiidenmx),nbastot,2.0D0)
      end if
      call prnmat('Density Alpha       ',icore(iiidenmx),nbastot,
     & nbastot,.true.)

      if (iuhf.eq.1) then
         call getrec(20,'JOBARC','SCFEVCB0',nbastot*nbastot*iintfp,
     & icore(iiicoefs+nbastot*nbastot*iintfp))
         call prnmat('Canonical Coef Beta ',icore(iiicoefs+
     & nbastot*nbastot*iintfp),nbastot,nbastot,.true.)

         call mkdens(icore(iiicoefs+nbastot*nbastot*iintfp),
     & icore(iiidenmx+nbastot*nbastot*iintfp),icore(iiiordvc),
     & nbastot,nocc(2),.false.)
         call prnmat('Density Beta        ',icore(iiidenmx+
     & nbastot*nbastot*iintfp),nbastot,nbastot,.true.)
      end if

      call dzerosec(nbastot*nbastot*iintfp,icore(iiiovrlp)) 
      call getrec(20,'JOBARC','AOOVRLAP',nbastot*nbastot*iintfp,
     & icore(iiiovrlp))
      call prnmat('Overlap             ',icore(iiiovrlp),nbastot,
     & nbastot,.true.)

      call dzerosec(nbastot*nbastot*iintfp,icore(iiidenmx))
      call transform2(icore(iiidtran),icore(iiidenmx),
     & icore(iiicoefs),nbastot)
      call prnmat('Canon Coef New Alpha',icore(iiicoefs),nbastot,
     & nbastot,verprn)

      call dzerosec(nbastot*nbastot*iintfp,icore(iiidenmx))
      call mkdens(icore(iiicoefs),icore(iiidenmx),icore(iiiordvc),
     & nbastot,nocc(1),.false.)
      if (iuhf.eq.0) then
         call scalmat(icore(iiidenmx),nbastot,2.0D0)
      end if
      call prnmat('Density New Alpha   ',icore(iiidenmx),nbastot,
     & nbastot,verprn)

      if (iuhf.eq.1) then
         call dzerosec(nbastot*nbastot*iintfp,icore(iiidenmx+
     & nbastot*nbastot*iintfp))
         call transform2(icore(iiidtran),icore(iiidenmx+
     & nbastot*nbastot*iintfp),icore(iiicoefs+nbastot*nbastot*
     & iintfp),nbastot)
         call prnmat('Canon Coef New Beta ',icore(iiicoefs+
     & nbastot*nbastot*iintfp),nbastot,nbastot,verprn)

         call dzerosec(nbastot*nbastot*iintfp,icore(iiidenmx+
     & nbastot*nbastot*iintfp))
         call mkdens(icore(iiicoefs+nbastot*nbastot*iintfp),
     & icore(iiidenmx+nbastot*nbastot*iintfp),icore(iiiordvc),
     & nbastot,nocc(2),.false.)
         call prnmat('Density New Beta    ',icore(iiidenmx+
     & nbastot*nbastot*iintfp),nbastot,nbastot,verprn)
      end if

      call getrec(0,'JOBARC','CORRDENA',nbas,switch12)
      if (nbas.eq.-1) then
         nbas = 0
         switch12 = 0
         go to 951
      end if
      nbas = 0
      switch12 = 0
      write(*,*)
      write(*,49)
      write(*,*)
      call dzerosec(nbastot*nbastot,icore(iiidenmx))
      call getrec(20,'JOBARC','CORRDENA',nbastot*nbastot,
     & icore(iiidenmx))
      call prnmat('CC Response AO densA',icore(iiidenmx),nbastot,
     & nbastot,.true.)
      call dzerosec(nbastot*nbastot,icore(jjj00000))
      call transform(icore(iiidtran),icore(jjj00000),
     & icore(iiidenmx),nbastot)
      call prnmat('New CC Resp AO densA',icore(iiidenmx),nbastot,
     & nbastot,verprn)
      call dzerosec(nbastot*nbastot,icore(jjj00000))
      if (iuhf.eq.1) then
         write(*,*)
         write(*,50)
         write(*,*)
         call dzerosec(nbastot*nbastot,icore(iiidenmx+
     & nbastot*nbastot))
         call getrec(20,'JOBARC','CORRDENB',nbastot*nbastot,
     & icore(iiidenmx+nbastot*nbastot))
         call prnmat('CC Response AO densB',icore(iiidenmx+
     & nbastot*nbastot),nbastot,nbastot,.true.)
         call dzerosec(nbastot*nbastot,icore(jjj00000))
         call transform(icore(iiidtran),icore(jjj00000),
     & icore(iiidenmx+nbastot*nbastot),nbastot)
         call prnmat('New CC Resp AO densB',icore(iiidenmx+
     & nbastot*nbastot),nbastot,nbastot,verprn)
         call dzerosec(nbastot*nbastot,icore(jjj00000))
      end if

 951  continue

      call dzerosec(nbastot*nbastot*iintfp,icore(iiicoefs))
      call transform(icore(iiidtran),icore(iiicoefs),
     & icore(iiiovrlp),nbastot)
      call prnmat('Overlap New         ',icore(iiiovrlp),nbastot,
     & nbastot,verprn)

      iiinchoo = iiidenmx+nbastot*nbastot*iintfp*(iuhf+1)
      call ckmem(iiinchoo,iiicrtop)

      iiichoos= iiinchoo+bondsize
      call ckmem(iiichoos,iiicrtop)

      iiidscrs = iiichoos+bondsize*bondsize*mxchoose
      call ckmem(iiidscrs,iiicrtop)

      call izerosec(icore(jjj00000),nbastot+2*angmax*natoms)

      jjjiscrs = jjj00000
      call ckmem(jjjiscrs,jjjcrtop)

      call dzerosec(nbastot*nbastot*iintfp,icore(iiicoefs))
      if (iuhf.eq.1) then
         call dzerosec(nbastot*nbastot*iintfp,icore(iiicoefs+
     & nbastot*nbastot*iintfp))
      end if

      call initmax(maxocc,iuhf)

      write(*,*)
      write(*,24)
      write(*,*)

      call putrec(20,'JOBARC','IIINSHEL',natoms,icore(iiinshel))
      call putrec(20,'JOBARC','IIISHELL',(largel+1)*natoms,
     & icore(iiishell))
      call putrec(20,'JOBARC','IIINBASA',(largel+1)*natoms,
     & icore(iiinbasa))
      call putrec(20,'JOBARC','IIIACHRG',natoms,icore(iiiachrg))
      call putrec(20,'JOBARC','IIINCHOO',bondsize,icore(iiinchoo))
      call putrec(20,'JOBARC','IIICHOOS',bondsize*mxchoose*bondsize,
     & icore(iiichoos))
      call putrec(20,'JOBARC','MYOVRLP ',nbastot*nbastot*iintfp,
     & icore(iiiovrlp))

      if (iuhf.eq.0) then
         write(*,*)
         write(*,*)
         write(*,34)
         write(*,*)
         write(*,*)
      end if
      if (iuhf.eq.1) then
         write(*,*)
         write(*,*)
         write(*,35)
         write(*,*)
         write(*,*)
      end if
 
      if (setchoose) then
         call setchoo(icore(iiinchoo),icore(iiichoos),mxchoose,
     & bondsize)
      end if

	write (*,*) 'started'
C#ifdef _NOSKIP
      call nlo__gener_nlo_orbitals(1,jjjcrtop-jjjiscrs,
     & iiicrtop-iiidscrs,nbastot,natoms,largel,mxnal,maxocc,
     & bondsize,ryd2hyb,icore(iiinshel),icore(iiishell),
     & icore(iiinbasa),icore(iiiachrg),icore(iiidenmx),
     & icore(iiiovrlp),mjump,spheric,onlynao,onlynho,onlynbo,
     & onlynlmo,mxchoose,icore(iiinchoo),icore(iiichoos),
     & icore(jjjiscrs),icore(iiidscrs),icore(iiicoefs))
C#endif

	write(*,*) 'finished'

      if (iuhf.eq.1) then

         call getrec(20,'JOBARC','IIINSHEL',natoms,icore(iiinshel))
         call getrec(20,'JOBARC','IIISHELL',(largel+1)*natoms,
     & icore(iiishell))
         call getrec(20,'JOBARC','IIINBASA',(largel+1)*natoms,
     & icore(iiinbasa))
         call getrec(20,'JOBARC','IIIACHRG',natoms,icore(iiiachrg))
         call getrec(20,'JOBARC','IIINCHOO',bondsize,icore(iiinchoo))
         call getrec(20,'JOBARC','IIICHOOS',bondsize*mxchoose*bondsize,
     & icore(iiichoos))
         call getrec(20,'JOBARC','MYOVRLP ',nbastot*nbastot*iintfp,
     & icore(iiiovrlp))

         write(*,*)
         write(*,*)
         write(*,37)
         write(*,*)
         write(*,*)
C#ifdef _NOSKIP
         call nlo__gener_nlo_orbitals(0,jjjcrtop-jjjiscrs,
     & iiicrtop-iiidscrs,nbastot,natoms,largel,mxnal,maxocc,
     & bondsize,ryd2hyb,icore(iiinshel),icore(iiishell),
     & icore(iiinbasa),icore(iiiachrg),
     & icore(iiidenmx+nbastot*nbastot*iintfp),icore(iiiovrlp),
     & mjump,spheric,onlynao,onlynho,onlynbo,onlynlmo,mxchoose,
     & icore(iiinchoo),icore(iiichoos),icore(jjjiscrs),
     & icore(iiidscrs),icore(iiicoefs+nbastot*nbastot*iintfp))
C#endif

      end if

      write(*,*)
      write(*,25)
      write(*,*)
!!!!! iiicoefs is location of nlmo transformed coeffs.

      call prnmat('Local New Alpha     ',icore(iiicoefs),nbastot,
     & nbastot,verprn)

      call dzerosec(nbastot*nbastot*iintfp,icore(iiidenmx))
      call mkdens(icore(iiicoefs),icore(iiidenmx),icore(iiiordvc),
     & nbastot,nocc(1),.false.)
      if (iuhf.eq.0) then
         call scalmat(icore(iiidenmx),nbastot,2.0D0)
      end if
      call prnmat('Local Dens New Alpha',icore(iiidenmx),nbastot,
     & nbastot,verprn)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (iuhf.eq.1) then
         call prnmat('Local New Beta      ',icore(iiicoefs+
     & nbastot*nbastot*iintfp),nbastot,nbastot,verprn)

         call dzerosec(nbastot*nbastot*iintfp,icore(iiidenmx+
     & nbastot*nbastot*iintfp))
         call mkdens(icore(iiicoefs+nbastot*nbastot*iintfp),
     & icore(iiidenmx+nbastot*nbastot*iintfp),icore(iiiordvc),
     & nbastot,nocc(2),.false.)
         call prnmat('Local Dens New Beta ',icore(iiidenmx+
     & nbastot*nbastot*iintfp),nbastot,nbastot,verprn)
      end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call dzerosec(nbastot*nbastot*iintfp,icore(iiiovrlp))
      call getrec(20,'JOBARC','MYOVRLP ',nbastot*nbastot*iintfp,
     & icore(iiiovrlp))
      call prnmat('Overlap New         ',icore(iiiovrlp),nbastot,
     & nbastot,verprn)

      call dzerosec(nbastot*nbastot*iintfp,icore(iiidenmx))
      call transform2(icore(iiiinvtr),icore(iiidenmx),
     & icore(iiicoefs),nbastot)
      call chksgn(icore(iiicoefs),nbastot) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!
      inquire(file='orbsfile',exist=yesno)
      if (yesno) then
         call dzerosec(nbastot*nbastot*iintfp,icore(iiicoefs))
         call readorbfile(nbastot,icore(iiicoefs))
      end if

      call prnmat('Localized Alpha     ',icore(iiicoefs),nbastot,
     & nbastot,.true.)
!!! iiidenmx seems to be useless index here as its overwritten
      call dzerosec(nbastot*nbastot*iintfp,icore(iiidenmx))
      call getrec(20,'JOBARC','SCFEVCA0',nbastot*nbastot*iintfp,
     & icore(iiidenmx))
!!!! iiidenmx now holds the old-HF coeffs/evecs
      write(*,31) putaces
      if (.not.putaces) then
         write(*,36)
      end if
!!!! puts NLMO coeffs into JOBARC
      if (putaces) then
         write(*,*)
         write(*,*)
         write(*,30)
         write(*,*)
         write(*,*)
         call putrec(20,'JOBARC','SCFEVCA0',nbastot*nbastot*iintfp,
     & icore(iiicoefs))
C
#ifdef _DEBUG_LVL0
        iscr1 = iiidscrs
        inorm = iScr1 + nbastot*nbastot*iintfp
        call ckmem(inorm,iiicrtop)
        iscr2 = inorm  + nbastot*nbastot*iintfp
        call ckmem(iscr2,iiicrtop)
        iend  = iscr2 + nbastot*nbastot*iintfp
        call ckmem(iend,iiicrtop)
C
C Check and see whether Frank-Weinhold's NBO59 has run, if so 
C read the NLMO vectors from NBO59.
C
        LENGTH = -1
        CALL GETREC(0, "JOBARC", "NBO59_WI", LENGTH, NBO59)
        IF (LENGTH .GT. 0) Then
           Write(6,"(a)")"The NLMOS read from F. Weinhold NBO59 Program"
           CALL GETREC(20, "JOBARC", "NLMOS_59",nbastot*nbastot*iintfp,
     &                 icore(iiicoefs))
           Write(6,"(a)") "The NLMOS from NBO59"
           call output(icore(iiicoefs), 1, nbastot, 1, nbastot, nbastot,
     &                 nbastot, 1)
        ENDIF

        call getrec(20,'JOBARC','AOOVRLAP',nbastot*nbastot*iintfp,
     &              icore(iiiovrlp))

        call xgemm("T", "N", nbastot, nbastot, nbastot, 1.0D0, 
     &              icore(iiicoefs), nbastot, icore(iiiovrlp), 
     &              nbastot, 0.0D0, icore(iscr1), nbastot)

        call xgemm("N", "N", nbastot, nbastot, nbastot, 1.0D0, 
     &              icore(iscr1), nbastot, icore(iiicoefs), nbastot,
     &              0.0D0, icore(inorm), nbastot)
C  
        Write(6,"(a)") "The norm check: C^tSC = 1" 
        call output(icore(inorm), 1, nbastot, 1, nbastot, nbastot, 
     &             nbastot, 1)

       Write(6,"(a)") "The original CMOs"
       call output(icore(iiidenmx), 1, nbastot, 1, nbastot, nbastot, 
     &              nbastot, 1)

        call dcopy(nbastot*nbastot,icore(iiidenmx),1,icore(Iscr1),1)
        call minv(icore(iscr1), nbastot, nbastot, icore(iscr2),
     &            DET,1.0D-8, 0,0)

        call xgemm("N", "N", nbastot, nbastot, nbastot, 1.0D0,
     &              icore(Iscr1), nbastot, icore(iiidenmx), nbastot,
     &              0.0D0, icore(inorm), nbastot)
C
C CC^-1 = 1        
C
        Write(6,"(a)") "Check the inverse (C^-1)"
        call output(icore(inorm), 1, nbastot, 1, nbastot, nbastot, 
     &              nbastot, 1)

C Form (C'C^-1)

        call xgemm("N", "N", nbastot, nbastot, nbastot, 1.0D0, 
     &              icore(Iscr1), nbastot, icore(IIIcoefs), nbastot,
     &              0.0D0, icore(inorm), nbastot)

        call output(icore(inorm), 1, nbastot, 1, nbastot, nbastot, 
     &             nbastot, 1)

        Write(6,"(a)") "Check for the Unitary nature:(C'C^-1)^t(C'C^-1)"

        call xgemm("T", "N", nbastot, nbastot, nbastot, 1.0D0, 
     &              icore(inorm), nbastot, icore(inorm), nbastot,
     &              0.0D0, icore(iscr2), nbastot)

        call output(icore(iscr2), 1, nbastot, 1, nbastot, nbastot, 
     &             nbastot, 1)

#endif

         call getrec(0,'JOBARC','CORRDENA',nbas,switch12)
         if (nbas.eq.-1) then
            nbas = 0
            switch12 = 0
            go to 953
         end if
         nbas = 0
         switch12 = 0
         if (creatlist) then
            call orbfile(nbastot,icore(iiicoefs))
         end if
 953  continue

!!!! Puts old HF coeffs into below record
         call putrec(20,'JOBARC','SCFCANA0',nbastot*nbastot*iintfp,
     & icore(iiidenmx))
      end if
      call putrec(20,'JOBARC','SCFLOCA0',nbastot*nbastot*iintfp,
     & icore(iiicoefs))

      jjjpivot = jjj00000
      call ckmem(jjjpivot,jjjcrtop)

      call izerosec(icore(jjjpivot),nbastot)
      call minv(icore(iiidenmx),nbastot,nbastot,icore(jjjpivot),
     & det,1.0D-8,0,0)
      call dzerosec(nbastot*nbastot*iintfp,icore(iiidscrs))

!!!! I dont know what this does. Looks like we multiply old HF coeffs by
!!!! NLMO coeffs, and call it "Can Non Trans Alpha"
      call xgemm('N','N',nbastot,nbastot,nbastot,1.0D0,icore(iiidenmx),
     & nbastot,icore(iiicoefs),nbastot,0.0D0,icore(iiidscrs),nbastot)
      call prnmat('Can Non Trans Alpha ',icore(iiidscrs),nbastot,
     & nbastot,.true.)
      call putrec(20,'JOBARC','CANNONAL',nbastot*nbastot*iintfp,
     & icore(iiidscrs))

      call getrec(0,'JOBARC','CORRDENA',nbas,switch12)
      if (nbas.eq.-1) then
         nbas = 0
         switch12 = 0
         go to 952
      end if
      nbas = 0
      switch12 = 0
      if (creatlist) then
         call dzerosec(nbastot*nbastot*iintfp,icore(iiidenmx))
         call droplist(nocc(1),nvir(1),nbastot,icore(iiidscrs),
     & occthresh,virthresh,icore(iiidenmx),icore(iiidenmx+nbastot))
      end if

 952  continue

      call dzerosec(nbastot*nbastot*iintfp,icore(iiidenmx))
      call izerosec(icore(jjjpivot),nbastot)
!!! iiidscrs is old HF coeff x NLMO coeff; and we are moving this to index
!!! iiidenmx & transposing
      call dtranspose(icore(iiidscrs),icore(iiidenmx),nbastot)
      call xgemm('N','N',nbastot,nbastot,nbastot,1.0D0,icore(iiidenmx),
     & nbastot,icore(iiidscrs),nbastot,0.0D0,icore(jjjpivot),nbastot)
      call prnmat('Trans(T)Trans Alpha ',icore(jjjpivot),nbastot,
     & nbastot,.true.)

      call dzerosec(nbastot*nbastot*iintfp,icore(iiidenmx))
!!! mkdens creates density given NLMO coeff. iiiiodvc contains the ordering of the
!!! NLMOs and the boolean determines whether or not resorting occurs.
      call mkdens(icore(iiicoefs),icore(iiidenmx),icore(iiiiodvc),
     & nbastot,nocc(1),.false.)
      if (iuhf.eq.0) then
         call scalmat(icore(iiidenmx),nbastot,2.0D0)
      end if
      call prnmat('Localized Dens Alpha',icore(iiidenmx),nbastot,
     & nbastot,verprn)

      if (iuhf.eq.1) then
         call dzerosec(nbastot*nbastot*iintfp,
     & icore(iiidenmx+nbastot*nbastot*iintfp))
         call transform2(icore(iiiinvtr),
     & icore(iiidenmx+nbastot*nbastot*iintfp),
     & icore(iiicoefs+nbastot*nbastot*iintfp),nbastot)
         call chksgn(icore(iiicoefs+nbastot*nbastot*iintfp),nbastot)
         call prnmat('Localized Beta      ',icore(iiicoefs+nbastot*
     & nbastot*iintfp),nbastot,nbastot,.true.)

         call dzerosec(nbastot*nbastot*iintfp,icore(iiidenmx+
     & nbastot*nbastot*iintfp))
         call getrec(20,'JOBARC','SCFEVCB0',nbastot*nbastot*iintfp,
     & icore(iiidenmx+nbastot*nbastot*iintfp))
         write(*,31) putaces
         if (.not.putaces) then
            write(*,36)
         end if
         if (putaces) then
            write(*,*)
            write(*,*)
            write(*,30)
            write(*,*)
            write(*,*)
            call putrec(20,'JOBARC','SCFEVCB0',nbastot*nbastot*iintfp,
     & icore(iiicoefs+nbastot*nbastot*iintfp))
            call putrec(20,'JOBARC','SCFCANB0',nbastot*nbastot*iintfp,
     & icore(iiidenmx+nbastot*nbastot*iintfp))
         end if
      call putrec(20,'JOBARC','SCFLOCB0',nbastot*nbastot*iintfp,
     & icore(iiicoefs+nbastot*nbastot*iintfp))

         call izerosec(icore(jjjpivot),nbastot)
         call minv(icore(iiidenmx+nbastot*nbastot*iintfp),
     & nbastot,nbastot,icore(jjjpivot),det,1.0D-8,0,0)
         call dzerosec(nbastot*nbastot*iintfp,icore(iiidscrs))
         call xgemm('N','N',nbastot,nbastot,nbastot,1.0D0,
     & icore(iiidenmx+nbastot*nbastot*iintfp),nbastot,
     & icore(iiicoefs+nbastot*nbastot*iintfp),nbastot,
     & 0.0D0,icore(iiidscrs),nbastot)
         call prnmat('Can Non Trans Beta  ',icore(iiidscrs),nbastot,
     & nbastot,.true.)
         call putrec(20,'JOBARC','CANNONBE',nbastot*nbastot*iintfp,
     & icore(iiidscrs))

         call dzerosec(nbastot*nbastot*iintfp,icore(iiidenmx+
     & nbastot*nbastot*iintfp))
         call izerosec(icore(jjjpivot),nbastot)
         call dtranspose(icore(iiidscrs),icore(iiidenmx+
     & nbastot*nbastot*iintfp),nbastot)
         call xgemm('N','N',nbastot,nbastot,nbastot,1.0D0,
     & icore(iiidenmx+nbastot*nbastot*iintfp),nbastot,
     & icore(iiidscrs),nbastot,0.0D0,icore(jjjpivot),nbastot)
         call prnmat('Trans(T)Trans Beta  ',icore(jjjpivot),nbastot,
     & nbastot,.true.)

         call dzerosec(nbastot*nbastot*iintfp,
     & icore(iiidenmx+nbastot*nbastot*iintfp))
         call mkdens(icore(iiicoefs+nbastot*nbastot*iintfp),
     & icore(iiidenmx+nbastot*nbastot*iintfp),icore(iiiiodvc),
     & nbastot,nocc(2),.false.)
         call prnmat('Localized Dens Beta ',icore(iiidenmx+
     & nbastot*nbastot*iintfp),nbastot,nbastot,verprn)
      end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call dzerosec(nbastot*nbastot*iintfp,icore(iiidscrs))
      call transform(icore(iiiinvtr),icore(iiidscrs),
     & icore(iiiovrlp),nbastot)
      call prnmat('Overlap             ',icore(iiiovrlp),nbastot,
     & nbastot,.true.)

      call dzerosec(nbastot*nbastot*iintfp,icore(iiidscrs))
      call ckorth(icore(iiiovrlp),icore(iiidscrs),icore(iiicoefs),
     & nbastot)
      call prnmat('Overlap in MOS Alpha ',icore(iiidscrs),nbastot,
     & nbastot,.true.)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (iuhf.eq.1) then
         call dzerosec(nbastot*nbastot*iintfp,icore(iiidscrs))
         call ckorth(icore(iiiovrlp),icore(iiidscrs),
     & icore(iiicoefs+nbastot*nbastot*iintfp),nbastot)
         call prnmat('Overlap in MOS Beta  ',icore(iiidscrs),
     & nbastot,nbastot,.true.)
      end if

      iiifockm = iiidscrs+nbastot*nbastot*iintfp
      call ckmem(iiifockm,iiicrtop)

      iiievals = iiifockm+nbastot*nbastot*iintfp*(iuhf+1)
      call ckmem(iiievals,iiicrtop)

      call dzerosec(nbastot*nbastot*iintfp,icore(iiidscrs))
      call dzerosec(nbastot*nbastot*iintfp,icore(iiifockm))

      natoms = 0
      call getrec(0,'JOBARC','ONEHAO  ',natoms,
     & icore(iiifockm))
      call getrec(20,'JOBARC','ONEHAO  ',natoms,
     & icore(iiifockm))
      natoms = 0
      call unpkoneh(icore(iiifockm),icore(iiidscrs),nbastot)
      call prnmat('Core Hamiltonian    ',icore(iiidscrs),
     & nbastot,nbastot,verprn)

      call dzerosec(nbastot*nbastot*iintfp,icore(iiifockm))
      if (iuhf.eq.1) then
         call dzerosec(nbastot*nbastot*iintfp,icore(iiifockm+
     & nbastot*nbastot*iintfp))
      end if
	print*, 'iuhf is: ', iuhf
      if (iuhf.eq.0) then
         print*,'ZWW RHF'
         call getrec(20,'JOBARC','FOCKA   ',nbastot*nbastot*iintfp,
     & icore(iiifockm))
         call prnmat('Fock Alpha          ',icore(iiifockm),
     & nbastot,nbastot,.true.)
!!!!! iiidscrs is currently onehao

         call eneval(icore(iiidenmx),icore(iiidscrs),
     & icore(iiifockm),nbastot,eelec)
         
         call dzerosec(nbastot*nbastot*iintfp,icore(iiidenmx))
!!!!!! localizes core Hamiltonian onto index iiidenmx
         call orthock(icore(iiidenmx),icore(iiidscrs),
     & icore(iiicoefs),nbastot)
         call prnmat('Loc Core Hamil. in MO Alpha',icore(iiidenmx),
     & nbastot,nbastot,.true.)
         call putrec(20,'JOBARC','LOCCOREA',nbastot*nbastot*iintfp,
     & icore(iiidenmx))

         call dzerosec(nbastot*nbastot*iintfp,icore(iiidscrs))
!!!! localizes Fock  matrix 
         call orthock(icore(iiidscrs),icore(iiifockm),
     & icore(iiicoefs),nbastot)
         call prnmat('Loc Fock in MO Alpha',icore(iiidscrs),
     & nbastot,nbastot,.true.)
         call putrec(20,'JOBARC','LOCFOCKA',nbastot*nbastot*iintfp,
     & icore(iiidscrs))
	print*,'ZWW main call'
!!!! try custom resort: use eigenvalues of LOCFOCKA
!!!!  to determine col ordering. 
!!!! then load SCFEVCA0 and reorder cols
	beta=1
	call resort(icore(iiidscrs),nbastot,iintfp,beta)

!!! CALL A CHECK TO COMPUT C^Y F C =E IN LOCAL NLMO BASIS TO SEE EVALUE
!!!	ORDERING
!	 call checkFock(icore(iiidscrs),icore(iiicoefs),
!     &			icore(iiiovrlp),nbastot)
!!! ***********************************************************
!!!
         call dzerosec(2*nocc(1)*iintfp,icore(jjj00000))
         call getrec(20,'JOBARC','NUCREP  ',1*iintfp,enuc)
         call hfperbond(icore(iiidenmx),icore(iiidscrs),
     & nbastot,nocc(1),icore(jjj00000),icore(jjj00000+nocc(1)),enuc)
         call putrec(20,'JOBARC','MFBONDEN',nocc(1)*iintfp,
     & icore(jjj00000))

         call dzerosec(nbastot*iintfp,icore(iiievals))
         call doeval(icore(iiidscrs),icore(iiievals),nbastot)
         call prnmat('Loc Fock Diag Alpha ',icore(iiievals),
     & nbastot,1,.true.)
!!! NOW RESORTING EVALS
         call sortEvals(icore(iiievals),nbastot,iintfp,beta)
         acesener = 0.0D0
         call getfia(nbastot,nocc(1),nvir(1),icore(iiidscrs),
     & acesener)
         acesener = 0.0D0

         call getrec(20,'JOBARC','NUCREP  ',1*iintfp,enuc)
         write(*,*)
         write(*,19) eelec
         write(*,20) enuc
         write(*,21) enuc + eelec
         write(*,*)
         call getrec(20,'JOBARC','SCFENEG ',iintfp,acesener)
         if (abs(acesener-(enuc+eelec)).gt.1.0D-6) then
            write(*,48) acesener
            write(*,*)
            write(*,46)
C            stop
         end if
      end if

      if (iuhf.eq.1) then
        print*,'inside nlorb NLSCC'
         call dzerosec(nbastot*nbastot*iintfp,icore(iiiovrlp))
         call totalden(icore(iiidenmx),icore(iiidenmx+
     & nbastot*nbastot*iintfp),icore(iiiovrlp),nbastot)

         call getrec(20,'JOBARC','FOCKA   ',nbastot*nbastot*iintfp,
     & icore(iiifockm))
         call prnmat('Fock Alpha          ',icore(iiifockm),
     & nbastot,nbastot,verprn)

         call getrec(20,'JOBARC','FOCKB   ',nbastot*nbastot*iintfp,
     & icore(iiifockm+nbastot*nbastot*iintfp))
         call prnmat('Fock Beta           ',icore(iiifockm+
     & nbastot*nbastot*iintfp),nbastot,nbastot,verprn)

         call eneval2(icore(iiiovrlp),icore(iiidenmx),
     & icore(iiidenmx+nbastot*nbastot*iintfp),icore(iiidscrs),
     & icore(iiifockm),icore(iiifockm+nbastot*nbastot*iintfp),
     & nbastot,eelec)

         call dzerosec(nbastot*nbastot*iintfp,icore(iiidscrs))
         call orthock(icore(iiidscrs),icore(iiifockm),
     & icore(iiicoefs),nbastot)
         call prnmat('Loc Fock in MO Alpha',icore(iiidscrs),
     & nbastot,nbastot,.true.)
         call putrec(20,'JOBARC','LOCFOCKA',nbastot*nbastot*iintfp,
     & icore(iiidscrs))
!!! FOR UHF, RESORT COLS OF EVECS
	beta=1
        print*,'inside nlorb calling resort:'
        call resort(icore(iiidscrs),nbastot,iintfp,beta)
         call dzerosec(nbastot*iintfp,icore(iiievals))
         call doeval(icore(iiidscrs),icore(iiievals),nbastot)
         call prnmat('Loc Fock Diag Alpha ',icore(iiievals),
     & nbastot,1,.true.)
!!! RESORTING ALPHA ORBITAL EVALS
        call sortEvals(icore(iiievals),nbastot,iintfp,beta)
!!!
         acesener = 0.0D0
         call getfia(nbastot,nocc(1),nvir(1),icore(iiidscrs),
     & acesener)
         acesener = 0.0D0

         call dzerosec(nbastot*nbastot*iintfp,icore(iiidscrs))
         call orthock(icore(iiidscrs),icore(iiifockm+
     & nbastot*nbastot*iintfp),icore(iiicoefs+nbastot*nbastot*iintfp),
     & nbastot)
         call prnmat('Loc Fock in MO Beta ',icore(iiidscrs),
     & nbastot,nbastot,.true.)
         call putrec(20,'JOBARC','LOCFOCKB',nbastot*nbastot*iintfp,
     & icore(iiidscrs))
!!! FOR UHF, RESORT COLS ON EVECS BETA
	beta=300
        call resort(icore(iiidscrs),nbastot,iintfp,beta)
         call dzerosec(nbastot*iintfp,icore(iiievals))
         call doeval(icore(iiidscrs),icore(iiievals),nbastot)
         call prnmat('Loc Fock Diag Beta  ',icore(iiievals),
     & nbastot,1,.true.)
!!! RESORTING BETA ORBITAL EVALS
        call sortEvals(icore(iiievals),nbastot,iintfp,beta)
!!!
         acesener = 0.0D0
         call getfia(nbastot,nocc(2),nvir(2),icore(iiidscrs),
     & acesener)
         acesener = 0.0D0

         call getrec(20,'JOBARC','NUCREP  ',1*iintfp,enuc)
         write(*,*)
         write(*,19) eelec
         write(*,20) enuc
         write(*,21) enuc + eelec
         write(*,*)
         call getrec(20,'JOBARC','SCFENEG ',iintfp,acesener)
         if (abs(acesener-(enuc+eelec)).gt.1.0D-6) then
            write(*,48) acesener
            write(*,*)
            write(*,47)
C            stop
         end if
      end if

      write(*,*)
      write(*,23) 
      write(*,*)

 16   format(A20, ' Matrix')
 19   format('Electronic Energy                 = ', F18.8)
 20   format('Nuclear Energy                    = ', F18.8)
 21   format('Total Energy                      = ', F18.8)
 22   format('ACES: Initializing NBO-type Localization')
 23   format('ACES: Successful NBO-type Localization - Exiting')
 24   format('Entering NBO-type Localization written by Norbert Flocke')
 25   format('Exiting NBO-type Localization')
 30   format('Putting non-canonical NLMOS into JOBARC')
 36   format('Non-canonical NLMOS will _NOT_ be put into JOBARC')
 31   format('Logical variable putaces = ',L5)
 32   format('--- RHF REFERENCE ---')
 33   format('--- UHF REFERENCE ---')
 34   format('--- LOCALIZING WITH RHF REFERENCE ---')
 35   format('--- LOCALIZING WITH UHF REFERENCE ALPHA ---')
 37   format('--- LOCALIZING WITH UHF REFERENCE BETA  ---')
 40   format('THE FIRST TWO ATOMS OF ZMAT HAVE BEEN SWAPPED')
 41   format('FIX ZMAT ORDERING')
 42   format('--- Localization will exit ---')
 43   format('========================================')
 44   format('Natural Localized Orbital ACES Interface')
 45   format('      written by Thomas F. Hughes       ')
 46   format('Localized Orbitals May NOT be Non-Canonical RHF')
 47   format('Localized Orbitals May NOT be Non-Canonical UHF')
 48   format('HF energy from ACES = ',F20.15)
 49   format('Using coupled cluster response density alpha')
 50   format('Using coupled cluster response density beta')

      call aces_fin

      end

