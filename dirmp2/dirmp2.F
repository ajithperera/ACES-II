c direct RHF/UHF/ROHF MP2, see TCA 95,13,1997 (Table 2, ver.3)

      implicit double precision (a-h,o-z)
      COMMON /MACHSP/ IINTLN, IFLTLN, IINTFP, IALONE, IBITWD
      common // icore(1)
#include "mpif.h"
#include "parallel_aces.com"
      character*3 str

      call c_gtod(is0,ius0)
      call aces_init_rte

c ----------------------------------------------------------------------

c   o initialize the MPI environment
      call MPW_Init(iErr)
      if (iErr.ne.MPI_SUCCESS) then
         print *, '@DIRMP2: unable to initialize the MPI environment'
         call aces_exit(1)
      end if

c   o initialize the parallel_aces common block
      call aces_com_parallel_aces

#ifdef _HAVE_MPI
c   o change to the temp directory
      call cd_noderank(irank,iErr)
      if (iErr.ne.0) then
         print *, '@DIRMP2: unable to cd to the temp directory'
         print *, '         attempting to run in current dir'
      end if
#endif

c ----------------------------------------------------------------------

      call nr2asc (irank,str)
      open (96,file='OUT.'//str)
      open (97,file='DUMP.'//str)
      write (96,'(a,2i3)') 'nprocs,irank:',nprocs,irank

      call CRAPSI (icore(1),iuhf,-1)
      
      if (iuhf.eq.0) then
         write (96,'(a)') 'CLOSED SHELL DIRECT MP2 CALCULATION'
         call mp2_rhf
      else
         write (96,'(a)') 'OPEN SHELL DIRECT MP2 CALCULATION'
         call mp2_uhf
      end if
      call aces_ja_fin

      call MPW_Finalize(iErr)
      if (iErr.ne.MPI_SUCCESS) then
         print *, '@DIRMP2: unable to finalize the MPI environment'
         call aces_exit(1)
      end if

      call c_gtod(is1,ius1)
      print *, '@DIRMP2: Total wallclock time was ',
     &         (is1-is0)*1.d0+(ius1-ius0)*1.d-6,' seconds.'

      end
c------------------------------------------------------------------------
      subroutine mp2_rhf
      implicit double precision (a-h,o-z)
      common // icore(1)
#undef _PARALLEL_ACES_COM_
#include "parallel_aces.com"

      integer        iflags(100)
      common /flags/ iflags
#include "../gamess/gamess.par"
      common /GMSorder/ coeint(MXAO),iorder(MXAO),
     >               nAfrst(MXSH),nAstep(MXSH),sphrcl
      logical sphrcl, kandl

      COMMON /TIMEINFO/ TIMEIN, TIMENOW, TIMETOT, TIMENEW
      COMMON /MACHSP/ IINTLN, IFLTLN, IINTFP, IALONE, IBITWD

c! commons for communication with 'phase' routines
      common /mp2com/ N, nofj,nofk,nofl, kandl,nbaj,nbak,nbal,nrr
      common /rhfcom/ NO,NV,NO1,NOx
      dimension nbaj(15),nbak(15),nbal(15)

      dimension Notab(2)
      character*3 str
      data zer,two /0.d0, 2.d0/

      call timer(0)

      call getrec (1,'','NAOBASFN',1,n)
      write (96,*) 'no. of AOs:', n
      write (97,'(a8,i8)') 'NAOBASFN', 1
      write (97,*) n

      call getrec (1,'','NBASTOT ',1,nmo)
      write (96,*) 'no. of MOs:', nmo
      write (97,'(a8,i8)') 'NBASTOT ', 1
      write (97,*) nmo

      call getrec (1,'','NOCCORB ',2,NOtab)
      NOa=NOtab(1)
      NOb=NOtab(2)
      NO=NOa
      NV=NMO-NO
      write (96,'(a,2i4)') 'no. of occupied:', NOa,NOb
      write (97,'(a8,i8)') 'NOCCORB ', 2
      write (97,*) NOa
      write (97,*) NOb

c---------- initialize GAMESS integrals -----------

      open (98,file='MOL')
      if (irank.eq.0) open (99,file='MP2.LOG')
      call GETINF (98,.true.,99)
      close (98)
      write (96,*)
      write (96,*) 'MP2 GAMESS initialization complete'

      call xBASCHK (Lmax)
      inshmx=((Lmax+1)*(Lmax+2))/2

c------ memory allocation ------

      mem1 = 0
      mem1 = mem1 +nmo         ! Eorb
      mem1 = mem1 +nmo*nmo     ! coeS
      mem1 = mem1 +N*nmo       ! coe 
      mem1 = mem1 +N*nmo       ! coeG 
      mem1 = mem1 +N*nmo       ! Atab  
      mem1 = mem1 +N*inshmx*inshmx*inshmx   ! buf

      mem2 = 0
      mem2 = mem2 +N*inshmx*inshmx   ! tab1 (per one O)
      mem2 = mem2 +NV*inshmx*inshmx  ! tab2 (per one O)
      mem2 = mem2 +NV*NO*N           ! tab3 (per one O)
      mem2 = mem2 +NV*NO*NV          ! tab4 (per one O)

      MEMMAX = iflags(36)/iintfp
      mem = mem1 +NO*mem2
      memmin = mem1 +mem2
      write (96,'(a,i10)') 'memory available (words): ', MEMMAX
      write (96,'(a,i10)') 'memory to avoid splitting:', mem
      write (96,'(a,i10)') 'minimum memory to run:    ', memmin
      if (memmin.gt.MEMMAX) then
         call errex
      end if
      NOlim=NO
      if (mem.gt.MEMMAX) NOlim = (MEMMAX-mem1) / mem2
      npass = NO/NOlim
      memalo = mem1 +NOlim*mem2
      if (npass*NOlim.ne.NO) npass=npass+1
      write (96,'(a,i10)') 'chunk size:               ',NOlim
      write (96,'(a,i10)') 'No. of passes:            ',npass
      write (96,'(a,i10)') 'memory to allocate:       ',memalo
      write (96,*)

      memalo=memalo*IINTFP

      call ACES_MALLOC (memalo,icore,istart)
      if (icore(1).eq.0) then
         write (96,'(a)') 'unable to allocate...'
         call errex
      end if

      i00 = istart                    ! Eorb 
      i01 = i00 + nmo*IINTFP          ! coeS
      i02 = i01 + nmo*nmo*IINTFP      ! coe 
      i03 = i02 + N*nmo*IINTFP        ! coeG
      i04 = i03 + N*nmo*IINTFP        ! Atab
      i05 = i04 + N*nmo*IINTFP        ! buf
      i06 = i05 + N*inshmx*inshmx*inshmx*IINTFP  ! tab1
      i07 = i06 + NOlim*N*inshmx*inshmx*IINTFP   ! tab2
      i08 = i07 + NOlim*NV*inshmx*inshmx*IINTFP  ! tab3
      i09 = i08 + NOlim*NV*NO*N*IINTFP           ! tab4
      i10 = i09 + NOlim*NV*NO*NV*IINTFP          ! end of core

c------------ read in the rest of info --------- 

      call getrec (1,'','SCFEVLA0',iintfp*nmo, icore(i00))
      call getrec (1,'','SCFEVCA0',iintfp*nmo*nmo, icore(i01))
      write (97,'(a8,i8)') 'SCFEVLA0', iintfp*nmo
      call DTdump (97, icore(i00), nmo)
      write (97,'(a8,i8)') 'SCFEVCA0', iintfp*nmo*nmo
      call DTdump (97, icore(i01), nmo*nmo)
      close (97)

c-----------------------------

c transform MO's from SAO's to CAO's
      if (sphrcl) then
         call S2CAO (icore(i04),n,nmo)
         call XGEMM ('N','N', n,nmo,nmo, 1.d0,icore(i04),n, 
     >          icore(i01),nmo, 0.d0,icore(i02),n)
      else
         call dcopy (nmo*nmo, icore(i01), 1, icore(i02), 1)
      end if   

c transform MO's from ACES to GAMESS order:  coe --> coeG;
c this is necessary because the indices of integrals in GHONDO array
c are not consecutive numbers if ACES order is used, which would slow
c down the transformation and prevent BLASing
      call A2G (icore(i02),icore(i03),N,nmo,iorder)

      call timer(1)
      write (96,'(a,f10.2)') 'init time:      ',timenew

c----------------- begin transformation -------------------

      E2aa=zer
      E2ab=zer
      NO1=1-NOlim
      do 500 ipass=1,npass
         NO1=NO1+NOlim
         NOx=min (NOlim, NO-NO1+1)
         call pass_rhf (icore(i00),icore(i03),icore(i05),
     >        icore(i06),icore(i07),icore(i08),icore(i09),E2xaa,E2xab)
         E2aa=E2aa+E2xaa
         E2ab=E2ab+E2xab
 500  continue 
      E2=two*E2aa+E2ab
      call getrec(1,'JOBARC','SCFENEG ',iintfp,ESCF)
      ETOT=ESCF+E2
      call putrec(1,'JOBARC','TOTENERG',iintfp,ETOT)

      if (irank.eq.0) then
      write (96,*)
      write (96,'(a,f25.12)') 'E2(aa)  =',E2aa
      write (96,'(a,f25.12)') 'E2(ab)  =',E2ab
      write (96,'(a,f25.12)') 'E2(tot) =',E2
      write (*,*) ' Total MBPT(2) energy = ',ETOT
      end if
      end
c------------------------------------------------------------------------
      subroutine mp2_uhf
      implicit double precision (a-h,o-z)
      common // icore(1)
#undef _PARALLEL_ACES_COM_
#include "parallel_aces.com"

#include "../gamess/gamess.par"
      common /GMSorder/ coeint(MXAO),iorder(MXAO),
     >               nAfrst(MXSH),nAstep(MXSH),sphrcl
      logical sphrcl, kandl, rohf

      COMMON /TIMEINFO/ TIMEIN, TIMENOW, TIMETOT, TIMENEW
      COMMON /MACHSP/ IINTLN, IFLTLN, IINTFP, IALONE, IBITWD

c! commons for communication with 'phase' routines
      common /mp2com/ N, nofj,nofk,nofl, kandl,nbaj,nbak,nbal,nrr
      common /uhfcom/ NOa,NVa,NOb,NVb,NO1a,NOxa,NO1b,NOxb
      dimension nbaj(15),nbak(15),nbal(15)

      integer        iflags(100)
      common /flags/ iflags

      dimension Notab(2)
      character*3 str
      data zer,half /0.d0, 0.5d0/
      data rohf /.true./

      call timer(0)

      call getrec (1,'','NAOBASFN',1,n)
      write (96,*) 'no. of AOs:', n
      write (97,'(a8,i8)') 'NAOBASFN', 1
      write (97,*) n

      call getrec (1,'','NBASTOT ',1,nmo)
      write (96,*) 'no. of MOs:', nmo
      write (97,'(a8,i8)') 'NBASTOT ', 1
      write (97,*) nmo

      call getrec (1,'','NOCCORB ',2,NOtab)
      NOa=NOtab(1)
      NOb=NOtab(2)
      NVa=NMO-NOa
      NVb=NMO-NOb
      write (96,'(a,2i4)') 'no. of occupied:', NOa,NOb
      write (97,'(a8,i8)') 'NOCCORB ', 2
      write (97,*) NOa
      write (97,*) NOb

c---------- initialize GAMESS integrals -----------

      open (98,file='MOL')
      if (irank.eq.0) open (99,file='MP2.LOG')
      call GETINF (98,.true.,99)
      close (98)
      if (irank.eq.0) then
         write (96,*)
         write (96,*) 'MP2 GAMESS initialization complete'
      end if   

c get the maximum ang.momentum in the whole basis
      call xBASCHK (Lmax)
      inshmx=((Lmax+1)*(Lmax+2))/2

c------ memory allocation ------

      mem1 = 0
      mem1 = mem1 +nmo         ! Eorba
      mem1 = mem1 +nmo         ! Eorbb
      mem1 = mem1 +nmo*nmo     ! coeSa
      mem1 = mem1 +nmo*nmo     ! coeSb
      mem1 = mem1 +N*nmo       ! coea 
      mem1 = mem1 +N*nmo       ! coeb 
      mem1 = mem1 +N*nmo       ! coeGa 
      mem1 = mem1 +N*nmo       ! coeGb 
      mem1 = mem1 +N*nmo       ! Atab  
      mem1 = mem1 +N*inshmx*inshmx*inshmx   ! buf
      if (rohf) then
         mem1 = mem1 +nmo*nmo  ! Fa
         mem1 = mem1 +nmo*nmo  ! Fb
         mem1 = mem1 +nmo*NOa  ! F1a
         mem1 = mem1 +nmo*NOb  ! F1b
         mem1 = mem1 +NVa*NOa  ! F2a
         mem1 = mem1 +NVb*NOb  ! F2b
      end if

      mem2a = 0
      mem2a = mem2a +N*inshmx*inshmx     ! tab1aa (per one O)
      mem2a = mem2a +NVa*inshmx*inshmx   ! tab2aa (per one O)
      mem2a = mem2a +NVa*NOa*N           ! tab3aa (per one O)
      mem2a = mem2a +NVa*NOb*N           ! tab3ab (per one O)
      mem2a = mem2a +NVa*NOa*NVa         ! tab4aa (per one O)
      mem2a = mem2a +NVa*NOb*NVb         ! tab4ab (per one O)

      mem2b = 0
      mem2b = mem2b +N*inshmx*inshmx     ! tab1bb (per one O)
      mem2b = mem2b +NVb*inshmx*inshmx   ! tab2bb (per one O)
      mem2b = mem2b +NVb*NOb*N           ! tab3bb (per one O)
      mem2b = mem2b +NVb*NOb*NVb         ! tab4bb (per one O)

      MEMMAX = iflags(36)/iintfp
      NOmax=max(NOa,NOb)
      mem = mem1 +NOmax*(mem2a+mem2b)
      memmin = mem1 +mem2a +mem2b

      write (96,'(a,i10)') 'memory available (words): ', MEMMAX
      write (96,'(a,i10)') 'memory to avoid splitting:', mem
      write (96,'(a,i10)') 'minimum memory to run:    ', memmin
      if (memmin.gt.MEMMAX) then
         call errex
      end if
      NOlim=NOmax
      if (mem.gt.MEMMAX) NOlim = (MEMMAX-mem1) / (mem2a+mem2b)
      npass=NOmax/NOlim
      memalo = mem1 +NOlim*(mem2a+mem2b)
      if (npass*NOlim.ne.NOmax) npass=npass+1
      write (96,'(a,i10)') 'chunk size:               ',NOlim
      write (96,'(a,i10)') 'No. of passes:            ',npass
      write (96,'(a,i10)') 'memory to allocate:       ',memalo
      write (96,*)

      memalo=memalo*IINTFP

      call ACES_MALLOC (memalo,icore,istart)
      if (icore(1).eq.0) then
         write (96,'(a)') 'unable to allocate...'
         call errex
      end if

      i00 = istart                    ! Eorba 
      j00 = i00 + nmo*IINTFP          ! Eorbb 
      i01 = j00 + nmo*IINTFP          ! coeSa
      j01 = i01 + nmo*nmo*IINTFP      ! coeSb
      i02 = j01 + nmo*nmo*IINTFP      ! coea 
      j02 = i02 + N*nmo*IINTFP        ! coeb 
      i03 = j02 + N*nmo*IINTFP        ! coeGa
      j03 = i03 + N*nmo*IINTFP        ! coeGb
      i04 = j03 + N*nmo*IINTFP        ! Atab
      i05 = i04 + N*nmo*IINTFP        ! buf
      iend1 = i05 + N*inshmx*inshmx*inshmx*IINTFP
      
      if (rohf) then
         i51 = iend1                  ! Fa
         j51 = i51 +nmo*nmo*IINTFP    ! Fb
         i52 = j51 +nmo*nmo*IINTFP    ! F1a
         j52 = i52 +nmo*NOa*IINTFP    ! F1b
         i53 = j52 +nmo*NOb*IINTFP    ! F2a
         j53 = i53 +NVa*NOa*IINTFP    ! F2b
         iend1 = j53 +NVb*NOb*IINTFP
      end if

      i06 = iend1                                ! tab1aa
      j06 = i06 + NOlim*N*inshmx*inshmx*IINTFP   ! tab1bb
      i07 = j06 + NOlim*N*inshmx*inshmx*IINTFP   ! tab2aa
      j07 = i07 + NOlim*NVa*inshmx*inshmx*IINTFP ! tab2bb
      i08 = j07 + NOlim*NVb*inshmx*inshmx*IINTFP ! tab3aa
      j08 = i08 + NOlim*NVa*NOa*N*IINTFP         ! tab3bb
      k08 = j08 + NOlim*NVb*NOb*N*IINTFP         ! tab3ab
      i09 = k08 + NOlim*NVa*NOb*N*IINTFP         ! tab4aa
      j09 = i09 + NOlim*NVa*NOa*NVa*IINTFP       ! tab4bb
      k09 = j09 + NOlim*NVb*NOb*NVb*IINTFP       ! tab4ab
      i10 = k09 + NOlim*NVa*NOb*NVb*IINTFP       ! end of core

c------------ read in the rest of info --------- 

      call getrec (1,'','SCFEVLA0',iintfp*nmo, icore(i00))
      call getrec (1,'','SCFEVLB0',iintfp*nmo, icore(j00))
      call getrec (1,'','SCFEVCA0',iintfp*nmo*nmo, icore(i01))
      call getrec (1,'','SCFEVCB0',iintfp*nmo*nmo, icore(j01))

      write (97,'(a8,i8)') 'SCFEVLA0', iintfp*nmo
      call DTdump (97, icore(i00), nmo)
      write (97,'(a8,i8)') 'SCFEVLB0', iintfp*nmo
      call DTdump (97, icore(j00), nmo)
      write (97,'(a8,i8)') 'SCFEVCA0', iintfp*nmo*nmo
      call DTdump (97, icore(i01), nmo*nmo)
      write (97,'(a8,i8)') 'SCFEVCB0', iintfp*nmo*nmo
      call DTdump (97, icore(j01), nmo*nmo)

c ROHF path assumes that orbitals in JOBARC are already in the
c semi-canonical basis; 

      if (rohf) then
         call getrec (1,'','FOCKA   ',iintfp*nmo*nmo, icore(i51))
         call getrec (1,'','FOCKB   ',iintfp*nmo*nmo, icore(j51))
         write (97,'(a8,i8)') 'FOCKA   ', iintfp*nmo*nmo
         call DTdump (97, icore(i51), nmo*nmo)
         write (97,'(a8,i8)') 'FOCKB   ', iintfp*nmo*nmo
         call DTdump (97, icore(j51), nmo*nmo)
c transform the Fock matrices from SAO to MO basis
         call FockMO (icore(i51),icore(i52),icore(i53),icore(i01),
     >                nmo,NOa,NVa)
         call FockMO (icore(j51),icore(j52),icore(j53),icore(j01),
     >                nmo,NOb,NVb)
      end if

      close (97)

c-----------------------------

c transform MO's from SAO's to CAO's:  coeS --> coe;
c since integral calculation dominates the overall time, 
c we avoid the expensive Cartesian-->spherical integral
c transformation and work entirely in the Cartesian basis

      if (sphrcl) then
         call S2CAO (icore(i04),n,nmo)
         call XGEMM ('N','N', n,nmo,nmo, 1.d0,icore(i04),n, 
     >          icore(i01),nmo, 0.d0,icore(i02),n)
         call XGEMM ('N','N', n,nmo,nmo, 1.d0,icore(i04),n, 
     >          icore(j01),nmo, 0.d0,icore(j02),n)
      else
         call dcopy (nmo*nmo, icore(i01), 1, icore(i02), 1)
         call dcopy (nmo*nmo, icore(j01), 1, icore(j02), 1)
      end if   

c transform MO's from ACES to GAMESS order:  coe --> coeG;
c this is necessary because the indices of integrals in GHONDO array
c are not consecutive numbers if ACES order is used, which would slow
c down the transformation and prevent BLASing
      call A2G (icore(i02),icore(i03),N,nmo,iorder)
      call A2G (icore(j02),icore(j03),N,nmo,iorder)

      call timer(1)
      write (96,'(a,f10.2)') 'init time:      ',timenew

c----------------- begin transformation -------------------

      E2aa=zer
      E2bb=zer
      E2ab=zer
      NO1a=1-NOlim
      NO1b=1-NOlim
c in each pass a different chunk of the first occupied index is transformed
c NO1a (NO1b) is the first index calculated in a given pass
c NOxa (NOxb) is the number of indices in a given pass
c NOxb can be < 1 in the last chunks if NOa > NOb; the calculation of
c tab1bb, tab2bb, tab3bb, and tab4bb arrays is then automatically omitted
      do 500 ipass=1,npass
         NO1a=NO1a+NOlim
         NOxa=min (NOlim, NOa-NO1a+1)
         NO1b=NO1b+NOlim
         NOxb=min (NOlim, NOb-NO1b+1)
         call pass_uhf (icore(i00),icore(j00), icore(i03),icore(j03),
     >       icore(i05), icore(i06),icore(j06), icore(i07),icore(j07),
     >       icore(i08),icore(j08),icore(k08), 
     >       icore(i09),icore(j09),icore(k09), E2xaa,E2xbb,E2xab)
c         write (96,'(a,f25.12)') 'E2x=',E2x
         E2aa=E2aa+E2xaa
         E2bb=E2bb+E2xbb
         E2ab=E2ab+E2xab
 500  continue 

      E2 = E2aa+E2bb+E2ab

      if (rohf) then
         call E2sing (icore(i00),icore(i53),NOa,NVa, E2sa)
         call E2sing (icore(j00),icore(j53),NOb,NVb, E2sb)
         E2 = E2 +E2sa+E2sb
      end if

      call getrec(1,'JOBARC','SCFENEG ',iintfp,ESCF)
      ETOT=ESCF+E2
      call putrec(1,'JOBARC','TOTENERG',iintfp,ETOT)

      if (irank.eq.0) then
      write (96,*)
      if (rohf) then
         write (96,'(a,f25.12)') 'E2s(a)  =',E2sa
         write (96,'(a,f25.12)') 'E2s(b)  =',E2sb
      end if
      write (96,'(a,f25.12)') 'E2(aa)  =',E2aa
      write (96,'(a,f25.12)') 'E2(bb)  =',E2bb
      write (96,'(a,f25.12)') 'E2(ab)  =',E2ab
      write (96,'(a,f25.12)') 'E2(tot) =',E2
      write (*,*) ' Total MBPT(2) energy = ',ETOT
      end if
      end
c--------------------------------------------------------------------------
      subroutine pass_rhf(Eorb,coeG,buf,tab1,tab2,tab3,tab4,E2xaa,E2xab)
      implicit double precision (a-h,o-z)
      dimension Eorb(*),coeG(*)
      dimension buf(*),tab1(*),tab2(*),tab3(*),tab4(*)

#undef _PARALLEL_ACES_COM_
#include "parallel_aces.com"

      COMMON /TIMEINFO/ TIMEIN, TIMENOW, TIMETOT, TIMENEW
      common /mp2com/ N, nofj,nofk,nofl, kandl,nbaj,nbak,nbal,nrr
      common /rhfcom/ NO,NV,NO1,NOx

      dimension nbai(15),nbaj(15),nbak(15),nbal(15)
      logical iandj,kandl,same,done,pople,sphrcl
      integer pp,qq,rr,ss
      logical nosame,noflip

#include "../gamess/gamess.par"
      PARAMETER (MXDD=16*MXGSH*MXGSH)  
      COMMON /xNSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     >                CF(MXGTOT),CG(MXGTOT),
     >                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     >                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      COMMON /xSHLT  / TOL,CUTOFF,OUT,ICOUNT
      dimension ghondo(50625),ddij(MXDD)
      dimension ijind(225),klind(225) 
      common /GMSorder/ coeint(MXAO),iorder(MXAO),
     >               nAfrst(MXSH),nAstep(MXSH),sphrcl

c switch off the (ij) <-> (kl) symmetry in the integrals
      nosame=.true.
c switch off the ish<->jsh and ksh<->lsh shell flipping
      noflip=.true.

      call timer(0)
      iexch=1
      qq4=1.d0
      cutint=CUTOFF
      pople=.false.
      call zero (tab3, NOx*NV*NO*N)

      nloop=-1
      do 10 ss=1,nshell   ! ll
         nss=kmax(ss)-kmin(ss)+1
      do 20 rr=ss,nshell  ! kk
      nloop=nloop+1
      nloop=mod(nloop,nprocs)
      if (nloop.ne.irank) goto 20
         nrr=kmax(rr)-kmin(rr)+1
         call zero (tab1, NOx*N*nrr*nss)
      do 30 qq=1,nshell   ! jj
         nqq=kmax(qq)-kmin(qq)+1
         call zero (buf, N*nqq*nrr*nss)
      do 40 pp=qq,nshell  ! ii
         
         ii=pp
         jj=qq
         kk=rr
         ll=ss
         call GAMINT (ii,jj,kk,ll,iexch,qq4, done,
     >                ghondo,ddij,pople, nosame,noflip,
     >                iandj,kandl,same,nofi,nofj,nofk,nofl,
     >                ijind,klind,nbai,nbaj,nbak,nbal)
         if (.not.done) then
c            write (96,'(a,i4)') 'GAMINT not done:',pp,qq,rr,ss
            goto 40
         end if

c------------ loops over basis functions in four shells 
         kln=0
         lmax=nofl
         do 360 k=1,nofk
           coek=coeint(nbak(k))
           if (kandl) lmax=k
           do 340 l=1,lmax
             coekl=coek*coeint(nbal(l))
             kln=kln+1
             if (pople) then
               stop
             else
               n1=klind(kln)
             end if

             ijn=0
             jmax=nofj
             do 320 i=1,nofi
                ibas=nbai(i)
                coeikl=coekl*coeint(ibas)
                if (iandj) jmax=i
                do 300 j=1,jmax
                   ijn=ijn+1
                   if (same.and.kln.gt.ijn) goto 340
                   if (pople) then
                      stop
                   else
                      nn=n1+ijind(ijn)
                      val=ghondo(nn) *coeikl*coeint(nbaj(j))
                   end if
c        if (abs(val).lt.cutint) goto 300

c       write (99,'(a,4i4)') 
c     >  'basis functions:',ibas,nbaj(j),nbak(k),nbal(l)
c       write (99,'(a,g25.16)') 'val=',ghondo(nn)

       ind=ibas +(j-1)*N +(k-1)*N*nofj +(l-1)*N*nofj*nofk

       buf(ind)=val            
    
 300  continue  
 320  continue 
 340  continue       
 360  continue       
c------------ loops over basis functions in four shells (END)

 40   continue  ! ---------- end loop over pp shells

      call phase1a (coeG,buf,tab1, NO1,NOx)
      call phase1b (coeG,buf,tab1, NO1,NOx)

 30   continue  ! ---------- end loop over qq shells

      call phase2 (coeG,tab1,tab2, NOx,NO,NV)
      call phase3a (coeG,tab2,tab3, NOx,NO,NV)
      call phase3b (coeG,tab2,tab3, NOx,NO,NV)

 20   continue  ! ---------- end loop over rr shells

 10   continue  ! ---------- end loop over ss shells

      call timer(1)
      write (96,'(a,f10.2)') 'phases 1-3 time:',timenew

      call timer(0)
      call phase4 (coeG,tab3,tab4, NOx,NO,NV,NV)
      call timer(1)
      write (96,'(a,f10.2)') 'phase 4 time:   ',timenew

      call timer(0)
      call MPIAR (tab4,tab3,NOx*NV*NO*NV)
      call timer(1)
      write (96,'(a,f10.2)') 'MPIAR time:     ',timenew

      call timer(0)
      call sum_rhf (Eorb,tab4,NO1,NOx,NO,NV, E2xaa,E2xab)
      call timer(1)
      write (96,'(a,f10.2)') 'sum time:       ',timenew
      end
c--------------------------------------------------------------------------
      subroutine pass_uhf (Eorba,Eorbb,coeGa,coeGb,buf,tab1aa,tab1bb,
     >   tab2aa,tab2bb,tab3aa,tab3bb,tab3ab,tab4aa,tab4bb,tab4ab,
     >   E2xaa,E2xbb,E2xab)
      implicit double precision (a-h,o-z)
      dimension Eorba(*),coeGa(*)
      dimension Eorbb(*),coeGb(*)
      dimension buf(*)
      dimension tab1aa(*),tab2aa(*),tab3aa(*),tab4aa(*)
      dimension tab1bb(*),tab2bb(*),tab3bb(*),tab4bb(*)
      dimension tab3ab(*),tab4ab(*)

#undef _PARALLEL_ACES_COM_
#include "parallel_aces.com"

      COMMON /TIMEINFO/ TIMEIN, TIMENOW, TIMETOT, TIMENEW
      common /mp2com/ N, nofj,nofk,nofl, kandl,nbaj,nbak,nbal,nrr
      common /uhfcom/ NOa,NVa,NOb,NVb,NO1a,NOxa,NO1b,NOxb
      dimension nbai(15),nbaj(15),nbak(15),nbal(15)
      logical iandj,kandl,same,done,pople,sphrcl
      integer pp,qq,rr,ss
      logical nosame,noflip

#include "../gamess/gamess.par"
      PARAMETER (MXDD=16*MXGSH*MXGSH)  
      COMMON /xNSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     >                CF(MXGTOT),CG(MXGTOT),
     >                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     >                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      COMMON /xSHLT  / TOL,CUTOFF,OUT,ICOUNT
      dimension ghondo(50625),ddij(MXDD)
      dimension ijind(225),klind(225) 
      common /GMSorder/ coeint(MXAO),iorder(MXAO),
     >               nAfrst(MXSH),nAstep(MXSH),sphrcl

c switch off the (ij) <-> (kl) symmetry in the integrals
      nosame=.true.
c switch off the ish<->jsh and ksh<->lsh shell flipping
      noflip=.true.

      call timer(0)
      iexch=1
      qq4=1.d0
      cutint=CUTOFF
      pople=.false.
      call zero (tab3aa, NOxa*NVa*NOa*N)
      call zero (tab3bb, NOxb*NVb*NOb*N)
      call zero (tab3ab, NOxa*NVa*NOb*N)

      nloop=-1
      do 10 ss=1,nshell   ! ll
         nss=kmax(ss)-kmin(ss)+1
      do 20 rr=ss,nshell  ! kk
      nloop=nloop+1
      nloop=mod(nloop,nprocs)
      if (nloop.ne.irank) goto 20
         nrr=kmax(rr)-kmin(rr)+1
         call zero (tab1aa, NOxa*N*nrr*nss)
         call zero (tab1bb, NOxb*N*nrr*nss)
      do 30 qq=1,nshell   ! jj
         nqq=kmax(qq)-kmin(qq)+1
         call zero (buf, N*nqq*nrr*nss)
      do 40 pp=qq,nshell  ! ii
         
         ii=pp
         jj=qq
         kk=rr
         ll=ss
         call GAMINT (ii,jj,kk,ll,iexch,qq4, done,
     >                ghondo,ddij,pople, nosame,noflip,
     >                iandj,kandl,same,nofi,nofj,nofk,nofl,
     >                ijind,klind,nbai,nbaj,nbak,nbal)
         if (.not.done) then
c            write (96,'(a,i4)') 'GAMINT not done:',pp,qq,rr,ss
            goto 40
         end if

c------------ loops over basis functions in four shells 
         kln=0
         lmax=nofl
         do 360 k=1,nofk
           coek=coeint(nbak(k))
           if (kandl) lmax=k
           do 340 l=1,lmax
             coekl=coek*coeint(nbal(l))
             kln=kln+1
             if (pople) then
               stop
             else
               n1=klind(kln)
             end if

             ijn=0
             jmax=nofj
             do 320 i=1,nofi
                ibas=nbai(i)
                coeikl=coekl*coeint(ibas)
                if (iandj) jmax=i
                do 300 j=1,jmax
                   ijn=ijn+1
                   if (same.and.kln.gt.ijn) goto 340
                   if (pople) then
                      stop
                   else
                      nn=n1+ijind(ijn)
                      val=ghondo(nn) *coeikl*coeint(nbaj(j))
                   end if
c        if (abs(val).lt.cutint) goto 300

c       write (99,'(a,4i4)') 
c     >  'basis functions:',ibas,nbaj(j),nbak(k),nbal(l)
c       write (99,'(a,g25.16)') 'val=',ghondo(nn)

       ind=ibas +(j-1)*N +(k-1)*N*nofj +(l-1)*N*nofj*nofk

       buf(ind)=val            
    
 300  continue  
 320  continue 
 340  continue       
 360  continue       
c------------ loops over basis functions in four shells (END)

 40   continue  ! ---------- end loop over pp shells

      call phase1a (coeGa,buf,tab1aa, NO1a,NOxa)
      call phase1a (coeGb,buf,tab1bb, NO1b,NOxb)

      call phase1b (coeGa,buf,tab1aa, NO1a,NOxa)
      call phase1b (coeGb,buf,tab1bb, NO1b,NOxb)

 30   continue  ! ---------- end loop over qq shells

      call phase2 (coeGa,tab1aa,tab2aa, NOxa,NOa,NVa)
      call phase2 (coeGb,tab1bb,tab2bb, NOxb,NOb,NVb)

      call phase3a (coeGa,tab2aa,tab3aa, NOxa,NOa,NVa)
      call phase3a (coeGb,tab2bb,tab3bb, NOxb,NOb,NVb)
      call phase3a (coeGb,tab2aa,tab3ab, NOxa,NOb,NVa)

      call phase3b (coeGa,tab2aa,tab3aa, NOxa,NOa,NVa)
      call phase3b (coeGb,tab2bb,tab3bb, NOxb,NOb,NVb)
      call phase3b (coeGb,tab2aa,tab3ab, NOxa,NOb,NVa)

 20   continue  ! ---------- end loop over rr shells

 10   continue  ! ---------- end loop over ss shells

      call timer(1)
      write (96,'(a,f10.2)') 'phases 1-3 time:',timenew

      call timer(0)
      call phase4 (coeGa,tab3aa,tab4aa, NOxa,NOa,NVa,NVa)
      call phase4 (coeGb,tab3bb,tab4bb, NOxb,NOb,NVb,NVb)
      call phase4 (coeGb,tab3ab,tab4ab, NOxa,NOb,NVa,NVb)
      call timer(1)
      write (96,'(a,f10.2)') 'phase 4 time:   ',timenew

      call timer(0)
      if (NOxa.gt.0) then
         call MPIAR (tab4aa,tab3aa,NOxa*NVa*NOa*NVa)
         call MPIAR (tab4ab,tab3ab,NOxa*NVa*NOb*NVb)
      end if
      if (NOxb.gt.0) then
         call MPIAR (tab4bb,tab3bb,NOxb*NVb*NOb*NVb)
      end if
      call timer(1)
      write (96,'(a,f10.2)') 'MPIAR time:     ',timenew

      call timer(0)
      call sum_uhf (Eorba,Eorba,tab4aa, 
     >          NO1a,NOxa,NOa,NOa,NVa,NVa, .false., E2xaa)
      call sum_uhf (Eorbb,Eorbb,tab4bb, 
     >          NO1b,NOxb,NOb,NOb,NVb,NVb, .false., E2xbb)
      call sum_uhf (Eorba,Eorbb,tab4ab, 
     >          NO1a,NOxa,NOa,NOb,NVa,NVb, .true. , E2xab)
      call timer(1)
      write (96,'(a,f10.2)') 'sum time:       ',timenew
      end
c--------------------------------------------------------------------------
      subroutine S2CAO (A,n,nmo)
c compute transformation A from n Cartesian AO's to nmo spherical AO's;
c uses Gamess->Aces mappings established previously to avoid another analysis
c of the MOL file;
c the order of both CAO's and SAO's is the internal ACES order
      implicit double precision (a-h,o-z)
      dimension A(n,nmo)
#include "../gamess/gamess.par"

      COMMON /xNSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     >                CF(MXGTOT),CG(MXGTOT),
     >                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     >                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      common /GMSorder/ coeint(MXAO),iorder(MXAO),
     >               nAfrst(MXSH),nAstep(MXSH),sphrcl
      logical sphrcl
      dimension tabP(3,3),tabD(6,5),tabF(10,7),tabG(15,9)
      data tabP / 1.d0, 0.d0, 0.d0,
     >            0.d0, 1.d0, 0.d0,
     >            0.d0, 0.d0, 1.d0/
      data tabD /-1.d0, 0.d0, 0.d0,-1.d0, 0.d0, 2.d0, ! d0
     >            0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 0.d0, ! d-2           
     >            0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, ! d+1
     >            1.d0, 0.d0, 0.d0,-1.d0, 0.d0, 0.d0, ! d+2
     >            0.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0/ ! d-1
      data tabF /-1.d0, 0.d0, 0.d0,-1.d0, 0.d0, 4.d0, 
     >                        0.d0, 0.d0, 0.d0, 0.d0, ! f+1
     >            0.d0,-1.d0, 0.d0, 0.d0, 0.d0, 0.d0, 
     >                       -1.d0, 0.d0, 4.d0, 0.d0, ! f-1
     >            0.d0, 0.d0,-3.d0, 0.d0, 0.d0, 0.d0, 
     >                        0.d0,-3.d0, 0.d0, 2.d0, ! f0
     >            1.d0, 0.d0, 0.d0,-3.d0, 0.d0, 0.d0, 
     >                        0.d0, 0.d0, 0.d0, 0.d0, ! f+3
     >            0.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 
     >                        0.d0, 0.d0, 0.d0, 0.d0, ! f-2
     >            0.d0, 3.d0, 0.d0, 0.d0, 0.d0, 0.d0,  
     >                       -1.d0, 0.d0, 0.d0, 0.d0, ! f-3
     >            0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 
     >                        0.d0,-1.d0, 0.d0, 0.d0/ ! f+2
      data tabG / 3.d0, 0.d0, 0.d0, 6.d0, 0.d0,-24.d0,0.d0, 0.d0,
     >                  0.d0, 0.d0, 3.d0, 0.d0,-24.d0,0.d0, 8.d0, ! g0    
     >            0.d0,-1.d0, 0.d0, 0.d0, 0.d0, 0.d0,-1.d0, 0.d0,
     >                  6.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, ! g-2
     >            0.d0, 0.d0,-3.d0, 0.d0, 0.d0, 0.d0, 0.d0,-3.d0,
     >                  0.d0, 4.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, ! g+1
     >            1.d0, 0.d0, 0.d0,-6.d0, 0.d0, 0.d0, 0.d0, 0.d0,
     >                  0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 0.d0, ! g+4
     >            0.d0, 0.d0, 0.d0, 0.d0, 3.d0, 0.d0, 0.d0, 0.d0,
     >                  0.d0, 0.d0, 0.d0,-1.d0, 0.d0, 0.d0, 0.d0, ! g-3
     >           -1.d0, 0.d0, 0.d0, 0.d0, 0.d0, 6.d0, 0.d0, 0.d0,
     >                  0.d0, 0.d0, 1.d0, 0.d0,-6.d0, 0.d0, 0.d0, ! g+2
     >            0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 0.d0,-1.d0, 0.d0,
     >                  0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, ! g-4
     >            0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 0.d0,-3.d0,
     >                  0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, ! g+3
     >            0.d0, 0.d0, 0.d0, 0.d0,-3.d0, 0.d0, 0.d0, 0.d0,
     >                  0.d0, 0.d0, 0.d0,-3.d0, 0.d0, 4.d0, 0.d0/ ! g-1

      call zero (A,n*nmo)

      do 10 ish=1,NSHELL
         L = ktype(ish)-1
         nSinL = 2*L+1
         nCinL = ((L+1)*(L+2))/2
         nCAO1 = iorder( KLOC(ish) )
         do 20 i=1,nSinL
            nSAO = nAfrst(ish) +(i-1)*nAstep(ish)
            do 30 j=1,nCinL
               nCAO = nCAO1 + (j-1)*nAstep(ish)
               if (L.eq.0) then
                  A(nCAO,nSAO) = 1.d0
               else if (L.eq.1) then   
                  A(nCAO,nSAO) = tabP(j,i)
               else if (L.eq.2) then   
                  A(nCAO,nSAO) = tabD(j,i)
               else if (L.eq.3) then   
                  A(nCAO,nSAO) = tabF(j,i)
               else if (L.eq.4) then
                  A(nCAO,nSAO) = tabG(j,i)
               else
                  write (96,'(a,i6)') 'S2CAO: unsupported function:',L
                  stop
               end if
 30         continue
 20      continue
 10   continue
      end
c--------------------------------------------------------------------------
      subroutine A2G (coe,coeG,N,nmo,iorder)
      implicit double precision (a-h,o-z)
      dimension coe(*),coeG(*),iorder(*)

      do 5 iG=1,N
      iA=iorder(iG)
      ind=iG
      do 5 j=1,nmo
      coeG(ind)=coe(iA)
      iA=iA+N
      ind=ind+N
 5    continue
      end
c--------------------------------------------------------------------------
      subroutine FockMO (F,F1,F2,coeS, nmo,NO,NV)
c transforms Fock matrix from AO basis, F(nmo,nmo), first into F1(nmo,NO)
c and then into MO basis: F2(NV,NO);
c coeS(nmo,nmo) are MO's in the AO basis, ACES order
      implicit double precision (a-h,o-z)
      dimension F(*),F1(*),F2(*),coeS(*)

      call XGEMM ('N','N',nmo,NO,nmo, 1.d0,F,nmo, coeS,nmo, 0.d0,F1,nmo)
      call XGEMM ('T','N',NV,NO,nmo,  1.d0,coeS(nmo*NO+1),nmo, 
     >                        F1,nmo, 0.d0,F2,NV)

      end
c--------------------------------------------------------------------------
      subroutine E2sing (Eorb,F2,NO,NV, E2s)
c calculates singles contribution to E2 for ROHF reference
      implicit double precision (a-h,o-z)
      dimension Eorb(*),F2(*)
      integer a

      sum=0.d0
      ind=0
      do 10 i=1,NO
      Ei=Eorb(i)
      do 10 a=1,NV
      Ea=Eorb(NO+a)
      ind=ind+1
      sum=sum +F2(ind)*F2(ind)/ (Ei-Ea)
 10   continue
      E2s=sum
      end
c--------------------------------------------------------------------------
      subroutine DTdump (io,tab,n)
      double precision tab(*)
      do 10 i=1,n
 10   write (io,'(g23.16)') tab(i)
      end
c--------------------------------------------------------------------------
      subroutine phase1a (coeG,buf,tab1, NO1,NOx)
      implicit double precision (a-h,o-z)
      double precision coeG(*),buf(*),tab1(*)
      common /mp2com/ N, nofj,nofk,nofl, kandl,nbaj,nbak,nbal,nrr
      dimension nbaj(15),nbak(15),nbal(15)
      logical kandl
      integer q

      if (NOx.lt.1) return

      do 50 nj=1,nofj
      q=nbaj(nj)
      indc= q +(NO1-1)*N
      indb= q +(nj-1)*N 
      ind1= 1 +(q-1)*NOx

      call XGEMM ('T','N', NOx,nofk*nofl,N-q+1, 1.d0,coeG(indc),N,
     >      buf(indb),N*nofj, 1.d0,tab1(ind1),NOx*N)
 50   continue
      end


c      lmax=nofl
c      do 50 nk=1,nofk
c      if (kandl) lmax=nk
c      do 50 nl=1,lmax
c      do 50 nj=1,nofj
c      q=nbaj(nj)

c      ind1= 1 +(q-1)*NOx +(nk-1)*NOx*N +(nl-1)*NOx*N*nrr
c      do 54 i=1,NOx

c      sum=0.d0
c      indc= q +(NO1+i-2)*N
c      indb= q +(nj-1)*N +(nk-1)*N*nofj +(nl-1)*N*nofj*nofk
c      do 55 p=q,N
c      sum=sum +coeG(indc)*buf(indb)
c      indc=indc+1
c      indb=indb+1
c 55   continue
c      tab1(ind1)=tab1(ind1) +sumc

c      ind1=ind1+ 1
c 54   continue

c 50   continue
c      end
c--------------------------------------------------------------------------
      subroutine phase1b (coeG,buf,tab1, NO1,NOx)
      implicit double precision (a-h,o-z)
      double precision coeG(*),buf(*),tab1(*)
      common /mp2com/ N, nofj,nofk,nofl, kandl,nbaj,nbak,nbal,nrr
      dimension nbaj(15),nbak(15),nbal(15)
      logical kandl
      integer p,q
       
      lmax=nofl
      do 50 nk=1,nofk
      if (kandl) lmax=nk
      do 50 nl=1,lmax
      do 50 nj=1,nofj
      q=nbaj(nj)

      indc= q +(NO1-1)*N
      do 54 i=1,NOx
      coex=coeG(indc)

      ind1= i +q*NOx +(nk-1)*NOx*N +(nl-1)*NOx*N*nrr
      indb= (q+1) +(nj-1)*N +(nk-1)*N*nofj +(nl-1)*N*nofj*nofk
      do 56 p=q+1,N
      tab1(ind1)=tab1(ind1) +coex*buf(indb)
      ind1=ind1 +NOx
      indb=indb +1
 56   continue

      indc=indc+N
 54   continue

 50   continue
      end
c--------------------------------------------------------------------------
      subroutine phase2 (coeG,tab1,tab2, NOx,NO,NV)
      implicit double precision (a-h,o-z)
      double precision coeG(*),tab1(*),tab2(*)
      common /mp2com/ N, nofj,nofk,nofl, kandl,nbaj,nbak,nbal,nrr
      dimension nbaj(15),nbak(15),nbal(15)
      logical kandl
      integer a
 
      if (NOx.lt.1) return

      indc = NO*N +1
      lmax=nofl
      do 60 nk=1,nofk
      if (kandl) lmax=nk
      do 60 nl=1,lmax
      ind1= 1 +(nk-1)*NOx*N +(nl-1)*NOx*N*nrr
      ind2= 1 +(nk-1)*NOx*NV +(nl-1)*NOx*NV*nrr

      call XGEMM ('N','N',NOx,NV,N, 1.d0,tab1(ind1),NOx, coeG(indc),N,
     >             0.d0,tab2(ind2),NOx) 
 60   continue
      end

c      lmax=nofl
c      do 60 nk=1,nofk
c      if (kandl) lmax=nk
c      do 60 nl=1,lmax
c      do 60 a=1,NV
c      NBa=NO+a   

c      ind2= 1 +(a-1)*NOx +(nk-1)*NOx*NV +(nl-1)*NOx*NV*nrr
c      do 64 i=1,NOx

c      sum=0.d0
c      indc= 1 +(NBa-1)*N
c      ind1= i +(nk-1)*NOx*N +(nl-1)*NOx*N*nrr
c      do 65 q=1,N
c      sum=sum +coeG(indc)*tab1(ind1)
c      indc=indc+1
c      ind1=ind1+NOx
c 65   continue
c      tab2(ind2)=sum

c      ind2=ind2+1
c 64   continue

c 60   continue

c      end
c--------------------------------------------------------------------------
      subroutine phase3a (coeG,tab2,tab3, NOx,NO,NV)
      implicit double precision (a-h,o-z)
      double precision coeG(*),tab2(*),tab3(*)
      common /mp2com/ N, nofj,nofk,nofl, kandl,nbaj,nbak,nbal,nrr
      dimension nbaj(15),nbak(15),nbal(15)
      logical kandl
      integer a,r,s

      do 70 nl=1,nofl
      s=nbal(nl)   
      kmin=1
      if (kandl) kmin=nl
      do 70 a=1,NV
      do 70 i=1,NOx

      ind3= i +(a-1)*NOx +(s-1)*NOx*NV*NO
      do 71 j=1,NO

      sum=0.d0
      r= nbak(kmin)
      indc= r +(j-1)*N
      ind2= i +(a-1)*NOx +(kmin-1)*NOx*NV +(nl-1)*NOx*NV*nrr
      do 75 nk=kmin,nofk
      r=nbak(nk)   
      sum=sum +coeG(indc)*tab2(ind2)
      indc=indc+1
      ind2=ind2+NOx*NV
 75   continue
      tab3(ind3)=tab3(ind3) +sum

      ind3=ind3+NOx*NV
 71   continue

 70   continue

      end
c--------------------------------------------------------------------------
      subroutine phase3b (coeG,tab2,tab3, NOx,NO,NV)
      implicit double precision (a-h,o-z)
      double precision coeG(*),tab2(*),tab3(*)
      common /mp2com/ N, nofj,nofk,nofl, kandl,nbaj,nbak,nbal,nrr
      dimension nbaj(15),nbak(15),nbal(15)
      logical kandl
      integer a,r,s

      do 70 nl=1,nofl
      s=nbal(nl)   
      kmin=1
      if (kandl) kmin=nl+1
      do 70 a=1,NV
      do 70 i=1,NOx

      indc= s
      do 71 j=1,NO

      coex=coeG(indc)
      ind2= i +(a-1)*NOx +(kmin-1)*NOx*NV +(nl-1)*NOx*NV*nrr
      r=nbak(kmin)
      ind3= i +(a-1)*NOx +(j-1)*NOx*NV +(r-1)*NOx*NV*NO
      do 76 nk=kmin,nofk
      r=nbak(nk)   
      tab3(ind3)=tab3(ind3) +coex*tab2(ind2)
      ind2=ind2+NOx*NV
      ind3=ind3+NOx*NV*NO
 76   continue

      indc=indc+N
 71   continue

 70   continue
      end
c--------------------------------------------------------------------------
      subroutine phase4 (coeG,tab3,tab4, NOx,NO,NVl,NVr)
      implicit double precision (a-h,o-z)
      double precision coeG(*),tab3(*),tab4(*)
      common /mp2com/ N, nofj,nofk,nofl, kandl,nbaj,nbak,nbal,nrr
      dimension nbaj(15),nbak(15),nbal(15)
      logical kandl
      integer a,b,s

      if (NOx.lt.1) return
      
      indc = NO*N +1
      call XGEMM ('N','N', NOx*NVl*NO,NVr,N, 1.d0,tab3,NOx*NVl*NO, 
     >             coeG(indc),N, 0.d0,tab4,NOx*NVl*NO)  
      end

c      do 80 b=1,NVr
c      NBb=NO+b
c      do 80 j=1,NO
c      do 80 a=1,NVl

c      ind4= 1 +(a-1)*NOx +(j-1)*NOx*NVl +(b-1)*NOx*NVl*NO
c      do 81 i=1,NOx

c      sum=0.d0
c      indc= 1 +(NBb-1)*N
c      ind3= i +(a-1)*NOx +(j-1)*NOx*NVl
c      do 85 s=1,N
c      sum=sum +coeG(indc)*tab3(ind3)
c      indc=indc+1
c      ind3=ind3+NOx*NVl*NO
c 85   continue
c      tab4(ind4)=sum

c      ind4=ind4+1
c 81   continue

c 80   continue

c      end
c--------------------------------------------------------------------------
      subroutine sum_rhf (Eorb,tab4,NO1,NOx,NO,NV, E2aa,E2ab)
      implicit double precision (a-h,o-z)
      double precision Eorb(*),tab4(*)
      integer a,b
      data zer,one,fourth /0.d0, 1.d0, 0.25d0/

      E2aa=zer
      E2ab=zer

      do b=1,NV
      Eb=Eorb(NO+b)
      do j=1,NO
      Ej=Eorb(j)
      do a=1,NV
      Ea=Eorb(NO+a)
      ind4a= 1 +(a-1)*NOx +(j-1)*NOx*NV +(b-1)*NOx*NV*NO
      ind4b= 1 +(b-1)*NOx +(j-1)*NOx*NV +(a-1)*NOx*NV*NO
      do i=1,NOx
      Ei=Eorb(NO1+i-1)

      xiajb=tab4(ind4a)
      xibja=tab4(ind4b)
      dif=xiajb-xibja
      den=one/(Ei+Ej-Ea-Eb)

      E2aa=E2aa +dif*dif*den
      E2ab=E2ab +xiajb*xiajb*den

c      E2=E2 + (xiajb*xiajb +0.5d0*dif*dif) /den
c      E2=E2 + ( 2.d0*xiajb*xiajb - xiajb*xibja ) /den

      ind4a=ind4a+1
      ind4b=ind4b+1
      end do
      end do
      end do
      end do

      E2aa=fourth*E2aa
      end
c--------------------------------------------------------------------------
      subroutine sum_uhf 
     >          (Eorbl,Eorbr,tab4, NO1,NOx,NOl,NOr,NVl,NVr, ab, E2)
      implicit double precision (a-h,o-z)
      double precision Eorbl(*),Eorbr(*),tab4(*)
      logical ab
      integer a,b
      data zer,fourth /0.d0, 0.25d0/

      E2=zer
      do b=1,NVr
      Eb=Eorbr(NOr+b)
      do j=1,NOr
      Ej=Eorbr(j)
      do a=1,NVl
      Ea=Eorbl(NOl+a)
      ind4a= 1 +(a-1)*NOx +(j-1)*NOx*NVl +(b-1)*NOx*NVl*NOr
      ind4b= 1 +(b-1)*NOx +(j-1)*NOx*NVl +(a-1)*NOx*NVl*NOr
      do i=1,NOx
      Ei=Eorbl(NO1+i-1)

      xiajb=tab4(ind4a)
      if (ab) then
         dif=xiajb
      else
         xibja=tab4(ind4b)
         dif=xiajb-xibja
      end if
      den=Ei+Ej-Ea-Eb

      E2=E2 + dif*dif/den

      ind4a=ind4a+1
      ind4b=ind4b+1
      end do
      end do
      end do
      end do
      if (.not.ab) E2=fourth*E2

      end
c-------------------------------------------------------------
