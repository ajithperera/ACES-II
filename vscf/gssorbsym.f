      SUBROUTINE GSSORBSYM(CSO,SCR1,OCCNUM,ORBSYM,NBFIRR,
     &                     NBAS,NIRREP,IUHF,IINTFP,NOCC)
      IMPLICIT NONE
      DOUBLE PRECISION CSO,SCR1
      INTEGER OCCNUM,ORBSYM,NBFIRR,NBAS,NIRREP,IUHF,IINTFP
      INTEGER IFLAGS
      DOUBLE PRECISION TOL,CNORM,SNRM2
      INTEGER IBAS,IRREP,NUM,ISPIN,IOFF,IRRMAX,IOCC,NOCC,
     &        IOFFOCC,IOFFVRT,IOFFSYM,IPOS
      DIMENSION CSO(NBAS,NBAS,IUHF+1),SCR1(NBAS,NBAS)
      DIMENSION OCCNUM(NBAS,2),ORBSYM(NBAS,2),NBFIRR(8)
      DIMENSION CNORM(8)
      DIMENSION IOCC(8),NOCC(8,2),IOFFOCC(8,2),
     &          IOFFVRT(8,2),
     &          IOFFSYM(8)
      COMMON /FLAGS/ IFLAGS(100)
C
      DATA TOL / 1.0D-06 /
C
C     Go through a set of MOs in SO basis and determine computational
C     symmetry. Based on OCCNUM, make occupied-virtual separation.
C
      IF(IFLAGS(1).GE.10)THEN
       write(6,*) ' @GSSORBSYM-I, OCCNUM on input '
       write(6,*) occnum
       write(6,*) ' @GSSORBSYM-I, CSO on input '
       CALL OUTPUT(CSO,1,NBAS,1,NBAS,NBAS,NBAS,1)
       IF(IUHF.GT.0)THEN
        CALL OUTPUT(CSO(1,1,2),1,NBAS,1,NBAS,NBAS,NBAS,1)
       ENDIF
      ENDIF
C
      DO 100 ISPIN=1,IUHF+1
c YAU : old
c     CALL ICOPY(NBAS*NBAS*IINTFP,CSO(1,1,ISPIN),1,SCR1,1)
c YAU : new
      CALL DCOPY(NBAS*NBAS,CSO(1,1,ISPIN),1,SCR1,1)
c YAU : end
C
      IF(IFLAGS(1).GE.10)THEN
       call output(cso(1,1,ispin),1,nbas,1,nbas,nbas,nbas,1)
      ENDIF
C
      DO 40 IBAS =1,NBAS
      CALL ZERO(CNORM,8)
      IOFF = 0
      DO 10 IRREP=1,NIRREP
      IF(NBFIRR(IRREP).GT.0)THEN
       CNORM(IRREP) = SNRM2(NBFIRR(IRREP),CSO(1 + IOFF,IBAS,ISPIN),1)
      ENDIF
      IOFF = IOFF + NBFIRR(IRREP)
   10 CONTINUE
C
C     Check that we have one and only one nonzero norm.
C
      NUM = 0
      DO 20 IRREP=1,NIRREP
      IF(CNORM(IRREP).GT.TOL) NUM = NUM + 1
   20 CONTINUE
C
      IF(NUM.LT.0)THEN
       WRITE(6,*) ' @GSSORBSYM-F, Impossible ! Negative NUM ',NUM
       CALL ERREX
      ENDIF
C
      IF(NUM.EQ.0)THEN
       WRITE(6,*) ' @GSSORBSYM-F, Zero orbital. '
       WRITE(6,*) ' Orbital, spin ',IBAS,ISPIN
       CALL ERREX
      ENDIF
C
      IF(NUM.GT.1)THEN
       WRITE(6,*) ' @GSSORBSYM-F, Ambiguity. NUM > 1 '
       WRITE(6,*) ' CNORM, orbital, spin ',CNORM,IBAS,ISPIN
       CALL ERREX
C
C     This may be dangerous, but we are trying it since the D2 group
C     does not seem to be handled properly (noise in CSO).
C     [ it probably is not even a part of the solution of the D2 problem ]
C
C     Try zeroing out offending part.
C
C     First find maximum component of CNORM.
C
       IRRMAX = 1
       DO 21 IRREP=1,NIRREP
       IF(CNORM(IRREP).GT.CNORM(IRRMAX))THEN
        IRRMAX = IRREP
       ENDIF
   21  CONTINUE
C
C     Now zero out other parts.
C
       IOFF = 0
       DO 22 IRREP=1,NIRREP
C
       IF(NBFIRR(IRREP).GT.0)THEN
        IF(IRREP.NE.IRRMAX)THEN
         CALL ZERO(CSO(1 + IOFF,IBAS,ISPIN),NBFIRR(IRREP))
        ENDIF
       ENDIF
       IOFF = IOFF + NBFIRR(IRREP)
   22  CONTINUE

      ENDIF
C
C     Now that we know we have only one nonzero norm, find the
C     symmetry.
C
      DO 30 IRREP=1,NIRREP
      IF(CNORM(IRREP).GT.TOL) ORBSYM(IBAS,ISPIN) = IRREP
   30 CONTINUE
   40 CONTINUE
C
C     Find number of orbitals of each symmetry. Check this matches
C     with NBFIRR.
C
      CALL IZERO(IOCC,8)
      DO 50 IBAS=1,NBAS
      IRREP = ORBSYM(IBAS,ISPIN)
      IOCC(IRREP) = IOCC(IRREP) + 1
   50 CONTINUE
C
      DO 60 IRREP=1,NIRREP
      IF(IOCC(IRREP).NE.NBFIRR(IRREP))THEN
       WRITE(6,*) ' @GSSORBSYM-F, Mismatch ! IRREP, NBFIRR, IOCC '
       WRITE(6,*) IRREP,NBFIRR(IRREP),IOCC(IRREP)
       CALL ERREX
      ENDIF
   60 CONTINUE
C
C     Find occupation number for each symmetry.
C
      IF(IFLAGS(1) .GE. 10)THEN
       write(6,*) ' @GSSORBSYM-I, orbsym '
       write(6,*) orbsym
      ENDIF
      CALL IZERO(NOCC(1,ISPIN),8)
      DO 70 IBAS=1,NBAS
      IRREP = ORBSYM(IBAS,ISPIN)
      IF(OCCNUM(IBAS,ISPIN).GT.0)THEN
       NOCC(IRREP,ISPIN) = NOCC(IRREP,ISPIN) + 1
      ENDIF
   70 CONTINUE
C      

  100 CONTINUE
C
      IF(IFLAGS(1) .GE. 10)THEN
       WRITE(6,*) ' @GSSORBSYM-I, Occupation array '
       WRITE(6,*) (NOCC(IRREP,1     ),IRREP=1,NIRREP)
       WRITE(6,*) (NOCC(IRREP,1+IUHF),IRREP=1,NIRREP)
      ENDIF
C
      IOFFSYM(1) = 0
      IF(NIRREP.GE.2)THEN
       DO 150 IRREP=2,NIRREP
       IOFFSYM(IRREP) = IOFFSYM(IRREP-1) + NBFIRR(IRREP-1)
  150  CONTINUE
      ENDIF
C
      DO 200 ISPIN=1,IUHF+1
c YAU : old
c     CALL ICOPY(NBAS*NBAS*IINTFP,CSO(1,1,ISPIN),1,SCR1,1)
c YAU : new
      CALL DCOPY(NBAS*NBAS,CSO(1,1,ISPIN),1,SCR1,1)
c YAU : end
C
      IF(IFLAGS(1).GE.10)THEN
       write(6,*) ' @GSSORBSYM-I, 200 loop. CSO for ispin ',ispin
       CALL OUTPUT(CSO(1,1,ISPIN),1,NBAS,1,NBAS,NBAS,NBAS,1)
      ENDIF
C
cSSS Bug fix 10/2000 Ajith      CALL IZERO(IOFFOCC,8)
C
      CALL IZERO(IOFFOCC(1, ISPIN), 8)
      CALL ICOPY(8,NOCC(1,ISPIN),1,IOFFVRT(1,ISPIN),1)
C
      DO 190 IBAS=1,NBAS
C
      DO 180 IRREP=1,NIRREP
C
      IF(ORBSYM(IBAS,ISPIN).EQ.IRREP)THEN
C
       IF(OCCNUM(IBAS,ISPIN).GT.0)THEN
        IOFFOCC(IRREP,ISPIN) = IOFFOCC(IRREP,ISPIN) + 1
        IPOS = IOFFSYM(IRREP) + IOFFOCC(IRREP,ISPIN)
c YAU : old
c       CALL ICOPY(IINTFP*NBAS,SCR1(1,IBAS),1,CSO(1,IPOS,ISPIN),1)
c YAU : new
        CALL DCOPY(NBAS,SCR1(1,IBAS),1,CSO(1,IPOS,ISPIN),1)
c YAU : end
       ELSE
        IOFFVRT(IRREP,ISPIN) = IOFFVRT(IRREP,ISPIN) + 1
        IPOS = IOFFSYM(IRREP) + IOFFVRT(IRREP,ISPIN)
c YAU : old
c       CALL ICOPY(IINTFP*NBAS,SCR1(1,IBAS),1,CSO(1,IPOS,ISPIN),1)
c YAU : new
        CALL DCOPY(NBAS,SCR1(1,IBAS),1,CSO(1,IPOS,ISPIN),1)
c YAU : end
       ENDIF
      ENDIF
  180 CONTINUE
  190 CONTINUE
      IF(IFLAGS(1) .GE. 10)THEN
       write(6,*) ' occnum ',occnum
       write(6,*) ' ioffsym ',ioffsym
       write(6,*) ' ioffocc ',ioffocc
       write(6,*) ' ioffvrt ',ioffvrt
      ENDIF
C
      IF(IFLAGS(1).GE.10)THEN
       WRITE(6,*) ' @GSSORBSYM-I, Reordered SO MOs ',ISPIN
       CALL OUTPUT(CSO(1,1,ISPIN),1,NBAS,1,NBAS,NBAS,NBAS,1)
      ENDIF
  200 CONTINUE
C
      RETURN
      END
