
C THE PRIMARY PURPOSE OF THE INTProc MODULE IS TO PRODUCE THE
C  ORDERED MO INTEGRAL LISTS WHICH ARE EVENTUALLY WRITTEN TO
C  THE MOINTS FILE.  ALSO, A NUMBER OF TWO INDEX QUANTITIES
C  SUCH AS FOCK MATRIX LISTS ARE FORMED AS WELL.
C
C CODED BY J.F. STANTON AND J. GAUSS
cjp MR-BW-CC extension added by Jiri Pittner (1998-2000)

#include "flags.h"

      PROGRAM INTPROC
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DIRPRD,USEC,ABCDTYPE
      LOGICAL DOPHPH,ABIJ,DOALL,NOABCD
      LOGICAL COMPRESS,SYMMETRIC 
      LOGICAL MACRO_ITER 
      COMMON / / ICORE(1)
      COMMON /ISTART/ I0,ICRSIZ
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /FLAGS2/ IFLAGS2(500)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPY(255,2),DIRPRD(8,8)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /INCORE/ ICREC(2),USEC,IXT(2),IMOD(2)
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      COMMON /ABCD/ ABCDTYPE
      COMMON /DOINTS/ DOALL,ABIJ,NOABCD
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYMLOC/ ISYMOFF(8,8,25)
C
CJDW  KKB stuff.
C
      COMMON /SYMPOP2/ IRP_PD(8,22)
      COMMON /SHIFT/   ISHIFT,NDRGEO
CJDW END
cjp
cjp switches to be later implemented as aces namelist options
cjp
#include "bwcc.com"
cjp
C Parametrized CC parameters. 
#include "parcc.com"
c
      call bwread(.false.,iuhf)
      call parread(iuhf)
C
C INITIALIZE MACHINE DEPENDENT PARAMETERS
C
      if (isbwcc) iref=1
      CALL CRAPSI(ICORE(1),IUHF,0)
      if (isbwcc) call bwprep(nocco,nvrto,iuhf)
C
      CALL SETMET
      MAXCOR = ICRSIZ
      ITWO   = 2
cjp prepare for reopenmo
      do iref=1,nref
         call storemoio
      end do


      do 1234 iref = 1, nref
      if (iref.gt.nproc) call reopenmo

      totmaxdenom=0.

c----------------------------------------------------------------
c----    To handle the drop-mo in energy gradient  --  KB -------
c----------------------------------------------------------------
      call getrec(20,'JOBARC','NDROPGEO',1,NDRGEO)
      if (ndrgeo.eq.2) then
         ISHIFT = 300
      else
         ISHIFT = 0
      end if
c----------------------------------------------------------------

      NMO = NOCCO(1) + NVRTO(1)
cYAU      I0IRRBET=I0+MAXCOR-NMO
cYAU      I0IRRALP=I0IRRBET-NMO
cYAU      CALL GETREC(20,'JOBARC','IRREPALP',NMO,ICORE(I0IRRALP))
cYAU      CALL GETREC(20,'JOBARC','IRREPBET',NMO,ICORE(I0IRRBET))
cYAU      IPASS=I0IRRALP
cYAU      MAXCOR = MAXCOR + ( -2 * NMO )
      MAXCOR = MAXCOR - NMO - NMO
      IPASS = I0 + MAXCOR
      CALL GETREC(20,'JOBARC','IRREPALP',NMO,ICORE(IPASS      ))
      CALL GETREC(20,'JOBARC','IRREPBET',NMO,ICORE(IPASS + NMO))

      ABIJ   = .FALSE.
      NOABCD = .FALSE.
      DOALL  = .TRUE.
      IF (IFLAGS(83).EQ.0) THEN
         IF (IFLAGS(2).LE.1) DOALL=.FALSE.
         IF (IFLAGS(2).EQ.0) NOABCD=.TRUE.
         IF (IFLAGS(2).EQ.1.AND.
     &       IFLAGS(3).LE.1) NOABCD=.TRUE.
         IF (IFLAGS(2).EQ.1.AND.
     &       IFLAGS(18).GE.3.AND.IFLAGS(18).LT.9) THEN
            NOABCD=.TRUE.
         END IF
CMN
         IF (IFLAGS(2).EQ.1.AND.
     &       IFLAGS(3).EQ.0.AND.IFLAGS(74).EQ.0.AND.IFLAGS(22).EQ.0.AND.
     &       IFLAGS(87).EQ.0) THEN
            ABIJ=.TRUE.
         END IF
CMN END
cSG
c 4/25/96 For P-EOM-MBPT(2) we only need <AB|CD> for gradients
         IF (IFLAGS(87).GT.0) THEN
            IF ((IFLAGS(87).EQ.7.OR.IFLAGS(87).EQ.8).AND.
     &          IFLAGS(3).EQ.0.AND.IFLAGS(2).EQ.1) THEN
               NOABCD = .TRUE.
               DOALL  = .FALSE.
            ELSE
               NOABCD = .FALSE.
               DOALL  = .TRUE.
            END IF
         END IF
C
      ELSE IF (IFLAGS(83).EQ.1) THEN
         DOALL  = .TRUE.
         NOABCD = .FALSE.
      ELSE IF(IFLAGS(83).EQ.2) THEN
         DOALL  = .FALSE.
         NOABCD = .TRUE.
      END IF
      IF (IFLAGS(93).EQ.2) NOABCD=.TRUE.
CJDW
c---------------------------------------------------------------------
c---  To handle the drop-mo in energy gradient    --------  KB -------
c---  For Drop-Mo cases : without drop any virtual mo   --------------
c---  MOABCD is not recalculated in the second run of intprc   -------
c---------------------------------------------------------------------
      CALL GETREC(20,'JOBARC','NDROTVRT',1,NDROTVRT)
c     CALL GETREC(20,'JOBARC','IFLDMOGD',1,IFLDMGD)
      IF (.NOT.NOABCD .AND. NDRGEO.EQ.2) THEN
         IF (NDROTVRT.EQ.0) THEN
c        IF (NDROTVRT.EQ.0 .AND. IFLDMGD.NE.0 ) THEN
c           IFLDMGD = 2
            NOABCD  = .TRUE.
c           CALL PUTREC(20,'JOBARC','IFLDMOGD',1,IFLDMGD)
         END IF
      END IF
c---------------------------------------------------------------------
CJDW END
C
      CALL GETREC(20,'JOBARC','NOCCORB ',ITWO,NOCCO)
      CALL GETREC(20,'JOBARC','NVRTORB ',ITWO,NVRTO)
      IF (IUHF.EQ.0) THEN
         NOCC = NOCCO(1)
         NVRT = NVRTO(1)
      END IF
C
C  DO INITIAL SORTING OF INTEGRALS.  FORM PPPP, PPPH, ETC. FILES.
C
      call aces_io_reset
      call aces_cache_reset

      NSTO=NOCCO(1)+NVRTO(1)
C
CJDW  KKB stuff.
C     CALL DGMOI(ICORE(I0),MAXCOR,NSTO,NIRREP,IUHF)
C
      IF (NDRGEO.EQ.1) THEN
         CALL DGMOID(ICORE(I0),MAXCOR,NSTO,NIRREP,IUHF)
      ELSE
         CALL DGMOI(ICORE(I0),MAXCOR,NSTO,NIRREP,IUHF)
      END IF
CJDW END
C
      CALL CLMOIO(ICORE(I0),MAXCOR,NOCCO,NVRTO,IUHF)
C
C PROCESS HHHH INTEGRALS.  THIS IS DONE IN CORE.
C
      IF (.NOT.ABIJ) THEN
         IF (IUHF.NE.0) THEN
            CALL DS16AA(ICORE(I0),MAXCOR,1,IUHF,1,
     &                  ICANT,NMO,ICORE(IPASS))
            CALL DS16AA(ICORE(I0),MAXCOR,1,IUHF,2,
     &                  ICANT,NMO,ICORE(IPASS))
         END IF
         CALL DS16AB(ICORE(I0),MAXCOR,1,IUHF,ICANT,NMO,ICORE(IPASS))
C
C PROCESS PHHH INTEGRALS.  ALSO DONE IN CORE.  DONE FOR
C  E(2) GRADIENT AND HIGHER ENERGY CALCULATIONS.
C
            CALL DS25AA(ICORE(I0),MAXCOR,2,IUHF,1,IOUT,NMO,ICORE(IPASS))
         IF (IUHF.EQ.1) THEN
            CALL DS25AA(ICORE(I0),MAXCOR,2,IUHF,2,IOUT,NMO,ICORE(IPASS))
            CALL DS25AB(ICORE(I0),MAXCOR,2,IUHF,       NMO,ICORE(IPASS))
         ELSE
         IF (IOUT.EQ.0) THEN
            CALL DS25AB(ICORE(I0),MAXCOR,2,IUHF,       NMO,ICORE(IPASS))
         END IF
         END IF
      END IF
C
C PROCESS PPHH AND PHPH INTEGRALS.  THIS IS DONE IN CORE.
C
      DOPHPH = .NOT.ABIJ
      DO I = 1+IUHF, 1, -1
         CALL DS3AA(ICORE(I0),MAXCOR,3,IUHF,I,NMO,ICORE(IPASS))
         IF (DOPHPH) THEN
            CALL DS4AA(ICORE(I0),MAXCOR,4,IUHF,I,NMO,ICORE(IPASS))
         END IF
      END DO
      IF (IUHF.NE.0) THEN
         CALL DS3AB(ICORE(I0),MAXCOR,3,IUHF,3,NMO,ICORE(IPASS))
         IF (DOPHPH) THEN
            CALL DS4AB(ICORE(I0),MAXCOR,4,IUHF,3,NMO,ICORE(IPASS))
         END IF
      END IF
C
C Giving the users to choose whether they want no REDUNDANT lists
C option. This may save some  (i doubt that it is siginificant)
C disk space with a modest (i hope) increase in CPU time. This
C option was included to address some of Anthony's concerns about
C REDUNDANT storage despite these are not the disk intensive
C terms. Ajith Perera 07/2002

      IF (iFlags2(h_IFLAGS2_single_store).eq.0) THEN
         CALL RNABIJ(ICORE(I0),MAXCOR,IUHF,'W')
      END IF
C
C  DO PPPH INTEGRALS.
C
      IF (.NOT.ABIJ) THEN
         IF (IUHF.NE.0) THEN
            CALL DS25AA(ICORE(I0),MAXCOR,5,IUHF,1,IOUT,NMO,ICORE(IPASS))
            CALL DS25AA(ICORE(I0),MAXCOR,5,IUHF,2,IOUT,NMO,ICORE(IPASS))
         END IF
         CALL DS25AB(ICORE(I0),MAXCOR,5,IUHF,NMO,ICORE(IPASS))
      END IF
C
C NOW DO THE PPPP INTEGRALS.
C
      IF (.NOT.NOABCD) THEN
c         IMODE3 = 1
         call aces_io_remove(52,'MOABCD')
         IMODE3 = 0
         IF (IUHF.NE.0) THEN
            CALL INIPCK(1,1,1,231,IMODE3,0,1)
            CALL INIPCK(1,2,2,232,IMODE3,0,1)
            CALL INIPCK(1,13,13,233,IMODE3,0,1)
         ELSE
C
CMN  DETERMINE IF ABCD INTEGRALS ARE TO BE COMPRESSED OR NOT
C
            COMPRESS = IFLAGS2(107) .EQ. 2
            IF (.NOT. COMPRESS) THEN
               CALL INIPCK(1,13,13,233,IMODE3,0,1)
            ELSE
               CALL INIPCK(1,5,13,233,IMODE3,0,1)
            END IF
         END IF
         IF (IUHF.NE.0) THEN
            CALL DS16AA(ICORE(I0),MAXCOR,6,IUHF,1,
     &                  ICANT,NMO,ICORE(IPASS))
            CALL DS16AA(ICORE(I0),MAXCOR,6,IUHF,2,
     &                  ICANT,NMO,ICORE(IPASS))
         END IF
         CALL DS16AB(ICORE(I0),MAXCOR,6,IUHF,ICANT,NMO,ICORE(IPASS))
      END IF
C
C PROCESS FOCK MATRIX
C
      CALL DFOCK(ICORE(I0),MAXCOR,IUHF)
CJDW
c-------------------------------------------------------------------
c-----    end of w2 processing   --------- May 94,  KB  ------------
c-------------------------------------------------------------------
      if (ishift.ne.0) then
         CALL PUTSTF(ICORE(I0),MAXCOR,IUHF)
         call aces_fin
         stop
      end if
c-------------------------------------------------------------------
CJDW END
C
C WRITE DENOMINATOR ARRAYS AND FORM INITIAL T2 VECTOR.
C
      SYMMETRIC = (IFLAGS(h_IFLAGS_calclevel) .EQ. 48 .OR.
     +             IFLAGS(h_IFLAGS_calclevel) .EQ. 49)
#ifdef _DEBUG_LVLM
CSSS        Call getrec(0,"JOBARC","ORBOPITR",Ncycle,junk)
CSSS        If (Ncycle .lt. 0) Macro_iter = .True.
#endif 

      CALL DDIJAB(ICORE(I0),MAXCOR,IUHF)
       CALL DFRMT2(ICORE(I0),MAXCOR,IUHF)
      IF (IUHF.NE.0) THEN
         CALL DMOSMT(ICORE(I0),MAXCOR,IUHF)
         CALL S2PROJ(ICORE(I0),MAXCOR,IUHF,.FALSE.)
      END IF
C
C PRINT OUT DOMINANT AMPLITUDES FROM FIRST-ORDER WAVEFUNCTION.
C
c     CALL AMPSUM(ICORE(I0),MAXCOR,IUHF,0,.FALSE.,'T')
      IF (.NOT. SYMMETRIC) THEN
         CALL AMPSUM(ICORE(I0),MAXCOR,IUHF,0,isbwcc,'T')
      ELSE
         WRITE(6,*)
         WRITE(6,"(16x,a)") " ----------WARNING--------------"
         WRITE(6,"(a,a)") " For rCCD or drCCD the second-order",
     &                    " energy printed here is incorrect!" 
      ENDIF 
      CALL PUTSTF(ICORE(I0),MAXCOR,IUHF)
C
C PERFORM HARTREE-FOCK STABILITY ANALYSIS
C
      IF (IFLAGS(74).NE.0) CALL STABLE(ICORE(I0),MAXCOR,IUHF)
      IF (IFLAGS(1).GE.5) CALL aces_io_summary
      IF (IFLAGS(1).GE.1) THEN
         CALL MEMUSE
         CALL FLOPUSE(IUHF)
      END IF

      if (isbwcc) then
         write(*,*) '@INTPRC: iref=',iref,
     &              ' TOTAL MAX 1/denomin: ',totmaxdenom
      end if

      if (nref.gt.1) then
         call aces_cache_flush
         call aces_cache_reset
         call storemoio
      end if

CSSS      call check_ints(icore(i0),Maxcor,iuhf,.false.)

 1234 continue

      call aces_fin

      stop
      end

