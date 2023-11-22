
      SUBROUTINE EXPORT_GUESS_AO_ORBS(CZAO, EVAL, OCCNUM, REFCOORD,
     &                                COORD, CSO, ANGTYP, ATMASS,
     &                                CNEW, DCORE, ZAOSO, SCR1,
     &                                SCRN, SCRA, SCRB, SCRTMP,
     &                                IUHF, NBAS, NBASX, MAXDCOR,
     &                                LDIM2, NATOMS, IOLDOCC)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      INTEGER ANGTYP, OCCNUM
      LOGICAL MOSYMOK,MOSYMOK2,AOBASMOS_FILE, OLDAOMOS_FILE
C
#include "symm2.com"
#include "machsp.com"
#include "mxatms.par"
C
      DIMENSION CZAO(NBASX*NBAS,(IUHF+1)),REFCOORD(3,NATOMS),
     &          OCCNUM(NBAS*(IUHF+1)),EVAL((IUHF+1)*NBAS),
     &          ANGTYP(NBASX),COORD(3,NATOMS),ATMASS(NATOMS),
     &          WTRAN(3,3),IJUNK(MXATMS),JJUNK(MXATMS),CNEW(NBASX),
     &          DCORE(MAXDCOR),CSO(NBAS*NBAS*(IUHF+1)),
     &          ZAOSO(NBAS*NBASX),SCR1(NBAS*NBAS),SCRN(NBAS*NBAS),
     &          SCRB(LDIM2),SCRTMP(LDIM2),IOLDOCC(*)
C
      COMMON /FILES/LUOUT,MOINTS
      COMMON /FLAGS/IFLAGS(100)
      COMMON /POPUL/NOCC(16)
C
C      print*, "Entered the GUESS_READ_AOS", "IUHF", IUHF
C
C See the comments in initges.F for an explanation.
C
      AOBASMOS_FILE = .FALSE.  
      INQUIRE(FILE='AOBASMOS',EXIST=AOBASMOS_FILE)
      OLDAOMOS_FILE = .FALSE. 
      INQUIRE(FILE='OLDAOMOS',EXIST=OLDAOMOS_FILE)

      IF (OLDAOMOS_FILE) THEN
         OPEN(71,FILE='OLDAOMOS',STATUS='OLD',ACCESS='SEQUENTIAL',
     &        FORM='FORMATTED')
      ELSE IF (AOBASMOS_FILE) THEN
         OPEN(71,FILE='AOBASMOS',STATUS='OLD',ACCESS='SEQUENTIAL', 
     &        FORM='FORMATTED') 
      ENDIF

      NCOLS = 4
      NFULL = NBAS/NCOLS
      NLEFT = NBAS - NFULL*NCOLS
C
C Read first line of OLDAOMOS and check dimensions are as expected.
C
cYAU - what about the nbfirr array? (cf. evcao2)
      READ(71,'(11i7)') MBAS,MBASX,MIUHF
      IF(MBAS .NE. NBAS .OR. MBASX .NE. NBASX .OR. MIUHF .NE. IUHF)THEN
         WRITE(6,*) '  @INITGES-F, Mismatch in dimensions of OLDAOMOS '
         WRITE(6,*) '    Expected and found NBAS  : ', NBAS ,MBAS
         WRITE(6,*) '    Expected and found NBASX : ', NBASX,MBASX
         WRITE(6,*) '    Expected and found IUHF  : ', IUHF ,MIUHF
         CALL ERREX
      ENDIF
C
      DO 560 ISPIN=1,IUHF+1
C
       DO 520 IFULL=1,NFULL
       DO 510 IBAS =1,NBASX
       READ(71,'(4F20.10)')
     & (CZAO((IFULL-1)*NCOLS*NBASX + (J-1)*NBASX + IBAS,ISPIN),J=1,4)
  510  CONTINUE
  520  CONTINUE
C
       IF(NLEFT.GT.0)THEN
        DO 530 IBAS=1,NBASX
        IF(NLEFT.EQ.1) READ(71,'(1F20.10)')
     &   (CZAO(NFULL*NCOLS*NBASX + (J-1)*NBASX + IBAS,ISPIN),J=1,1)
        IF(NLEFT.EQ.2) READ(71,'(2F20.10)')
     &   (CZAO(NFULL*NCOLS*NBASX + (J-1)*NBASX + IBAS,ISPIN),J=1,2)
        IF(NLEFT.EQ.3) READ(71,'(3F20.10)')
     &   (CZAO(NFULL*NCOLS*NBASX + (J-1)*NBASX + IBAS,ISPIN),J=1,3)
  530   CONTINUE
       ENDIF
C
       DO 540 IBAS=1,NBAS
       READ(71,'(I3,F20.10)') OCCNUM((ISPIN-1)*NBAS+IBAS),
     &                        EVAL((ISPIN-1)*NBAS+IBAS)
  540  CONTINUE
C
       DO 550 IBAS=1,NBASX
       READ(71,*) ANGTYP(IBAS)
  550  CONTINUE
  560  CONTINUE
C
       DO 570 IATOM=1,NATOMS
       READ(71,*) (REFCOORD(I,IATOM),I=1,3)
  570  CONTINUE
C
C There is no reordering in SEWARD/ALASKA. The order of the atoms
C is the order that they typed in SEWARD.INP (ZMAT order).
C
       IF(IFLAGS(1).GE.10.AND.IFLAGS(56).NE.4)THEN
        write(6,*) ' @INITGES-I, MAP2ZMAT arrays '
        CALL GETREC(20,'JOBARC','MAP2ZMAT',NATOMS,JJUNK)
        DO 580 IATOM=1,NATOMS
        READ(71,*) IJUNK(IATOM)
        write(6,*) ijunk(IATOM),jjunk(IATOM)
  580   CONTINUE
       ENDIF
       CLOSE(71,STATUS='KEEP')
C
C Calculate matrix WTRAN relating current and reference coordinates.
C
       CALL GETREC(20,'JOBARC','COORD',NATOMS*3*IINTFP,COORD)
       CALL GETREC(20,'JOBARC','ATOMMASS',NATOMS*IINTFP,ATMASS)
       IF(IFLAGS(1) .GE. 10)THEN
        write(6,*) ' @INITGES-I, Current coordinates '
        do 590 iatom=1,natoms
        write(6,'(3F20.10)') (coord(i,iatom),i=1,3)
  590   continue
       ENDIF
       CALL Q2QPRIME(COORD,REFCOORD,ATMASS,SCR1,WTRAN,NATOMS)
C
C Put AO basis MOs in current orientation.
C
       CALL XFORMMOS(CZAO,CNEW,WTRAN,ANGTYP,EVAL,NBASX,NBAS,IUHF)
C
       IF(IFLAGS(1) .GE. 10)THEN
        WRITE(6,*) ' @INITGES-I, Alpha AO vectors after XFORMMOS '
        CALL OUTPUT(CZAO(1,1),1,NBASX,1,NBAS,NBASX,NBAS,1)
        IF(IUHF.GT.0)THEN
         WRITE(6,*) ' @INITGES-I, Beta AO vectors after XFORMMOS '
         CALL OUTPUT(CZAO(1,2),1,NBASX,1,NBAS,NBASX,NBAS,1)
        ENDIF
       ENDIF
C
C Check symmetry of AO basis MOs in current orientation.
C
       CALL MIRREPS(DCORE,MAXDCOR*IINTFP,NBAS,IUHF,EVAL,MOSYMOK)
c       Print*, MOSYMOK
       IF (.NOT.MOSYMOK) THEN
          CALL MODMOS(CZAO,SCR1,EVAL,NBASX,NBAS,IUHF)
          MOSYMOK2 = .TRUE.
          CALL MIRREPS(DCORE,MAXDCOR*IINTFP,NBAS,IUHF,EVAL,MOSYMOK2)
          IF (MOSYMOK2) THEN
             WRITE(*,*) ' after modification, symmetry is OK '
          ELSE
             WRITE(*,*) ' after modification, symmetry is still bad '
             CALL ERREX
          END IF
       END IF
C
C Transform the MOs to the SO basis.
C
       CALL GETREC(20,'JOBARC','ZMAT2CMP',NBAS*NBASX*IINTFP,ZAOSO)
       CALL XGEMM('N','N',NBAS,NBAS,NBASX,1.0D+00,
     &            ZAOSO,NBAS,CZAO(1,1),NBASX,0.0D+00,CSO,NBAS)
       IF(IFLAGS(1).GE.10)THEN
        WRITE(6,*) ' @INITGES-I, Alpha SO vectors '
        CALL OUTPUT(CSO,1,NBAS,1,NBAS,NBAS,NBAS,1)
       ENDIF
       IF(IUHF.GT.0)THEN
        CALL XGEMM('N','N',NBAS,NBAS,NBASX,1.0D+00,
     &             ZAOSO,NBAS,CZAO(1,2),NBASX,0.0D+00,
     &             CSO(NBAS*NBAS+1),NBAS)
        IF(IFLAGS(1).GE.10)THEN
         WRITE(6,*) ' @INITGES-I, Beta  SO vectors '
         CALL OUTPUT(CSO(NBAS*NBAS+1),1,NBAS,1,NBAS,NBAS,NBAS,1)
        ENDIF
       ENDIF
C
C Determine symmetry of SO MOs, block them and separate occupied
C and virtual so that we can calculate an initial density matrix.
C
       CALL GSSORBSYM(CSO,SCR1,OCCNUM,SCRN,NBFIRR,
     &                NBAS,NIRREP,IUHF,IINTFP,NOCC)
       IF(IUHF .EQ. 0)THEN
        CALL ICOPY(8,NOCC,1,NOCC(9),1)
       ENDIF
       CALL ICOPY(16,NOCC,1,IOLDOCC,1)
       IFLAGS(41) = 1
C
C Now we have ordered set of SOs in CSO and we know the occupation.
C Orthogonalize and write to JOBARC.
C
       CALL SORTHO(CSO,SCRN,SCRA,SCRB,SCRTMP,LDIM2,NBAS)
       CALL PUTREC(20,'JOBARC','SCFEVCA0',NBAS*NBAS*IINTFP,CSO)
       CALL PUTREC(20,'JOBARC','GUESSVA0',NBAS*NBAS*IINTFP,CSO)

       IF(IUHF.EQ.1) THEN
         CALL SORTHO(CSO(NBAS*NBAS+1),SCRN,SCRA,SCRB,SCRTMP,LDIM2,NBAS)
         CALL PUTREC(20,'JOBARC','SCFEVCB0',NBAS*NBAS*IINTFP,
     &               CSO(NBAS*NBAS+1))
         CALL PUTREC(20,'JOBARC','GUESSCB0',NBAS*NBAS*IINTFP,
     &               CSO(NBAS*NBAS+1))
       ENDIF
C
      RETURN
      END

