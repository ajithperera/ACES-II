      SUBROUTINE SETMET(IUHF,ICORE)
      IMPLICIT INTEGER (A-Z)
      LOGICAL NONSTD, YESNO, CCSD, MBPT, PARTEOM, NODAVID
      LOGICAL SS, SD, DS, DD, USESPINF, REALFREQ, FREQ_EXIST
      LOGICAL JSC_ALL
      DOUBLE PRECISION FREQ
      CHARACTER*8 PERTSTR,STRING
      DIMENSION ICORE(*)
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /FLAGS2/ IFLAGS2(500)
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /PERT1/ NOPERT
      COMMON /PERTTYP/ PERTSTR(100)
      COMMON/FILES/LUOUT,MOINTS
      COMMON/EOMINFO/CCSD, MBPT, PARTEOM, NODAVID
      COMMON/DRVHBAR/SS, SD, DS, DD
      COMMON /SPINF/ T1F, T2FD, T2FX, T1FI, T2FDI, T2FXI
C
C PROCESS INPUT
C
C
C        DETERMINE THE TYPE OF EOM CALCULATION THAT IS DONE
C
        CCSD = .FALSE.
        PARTEOM = .FALSE.
        IF (IFLAGS2(161) .EQ. 1) PARTEOM =.TRUE.
        IF (IFLAGS2(117) .EQ. 1) CCSD = .TRUE.
        MBPT = .NOT. CCSD
C
      WRITE(6,*)
      CALL BANNERJS(LUOUT)

cYAU : IRREPX used before assigned; set to 1
      IRREPX = 1
      CALL GETAOINF(IUHF,IRREPX)
C
C SEE IF SPECIAL INPUT EXISTS
C
      NONSTD=.FALSE.
      OPEN(UNIT=30,FILE='ZMAT',STATUS='OLD',FORM='FORMATTED')
1     READ(30,'(A)',END=102)STRING
      IF(INDEX(STRING,'CCEH').NE.0.OR.
     &   INDEX(STRING,'cceh').NE.0)THEN
       NONSTD=.TRUE.
      ELSE
       GOTO 1
      ENDIF
102   CONTINUE
C
      IF(NONSTD)THEN
         WRITE(6,1001)
         READ(30,*)NOPERT
         DO 103 IPERT=1,NOPERT
            READ(30,'(A8)')PERTSTR(IPERT)
 103     CONTINUE
      ELSE
         if (iflags2(162) .eq. 0) then
         NOPERT=3
         PERTSTR(1)='DIPOLE_X'
         PERTSTR(2)='DIPOLE_Y'
         PERTSTR(3)='DIPOLE_Z'
         elseif (iflags2(162) .eq. 1) then
         NOPERT=9
         PERTSTR(1)='DIPOLE_X'
         PERTSTR(2)='DIPOLE_Y'
         PERTSTR(3)='DIPOLE_Z'
c
         PERTSTR(4)='QUAD_XX '
         PERTSTR(5)='QUAD_YY'
         PERTSTR(6)='QUAD_ZZ '
         PERTSTR(7)='QUAD_XY '
         PERTSTR(8)='QUAD_XZ '
         PERTSTR(9)='QUAD_YZ '
c     PERTSTR(9)='2NDMO_YZ'
         elseif (iflags2(162) .eq. 2) then
         NOPERT=19
         PERTSTR(1)='DIPOLE_X'
         PERTSTR(2)='DIPOLE_Y'
         PERTSTR(3)='DIPOLE_Z'
c
         PERTSTR(4)='QUAD_XX '
         PERTSTR(5)='QUAD_YY'
         PERTSTR(6)='QUAD_ZZ '
         PERTSTR(7)='QUAD_XY '
         PERTSTR(8)='QUAD_XZ '
         PERTSTR(9)='QUAD_YZ '
c     PERTSTR(9)='2NDMO_YZ'
c
         PERTSTR(10)='OCTUPXXX'
         PERTSTR(11)='OCTUPYYY'
         PERTSTR(12)='OCTUPZZZ'
         PERTSTR(13)='OCTUPXXY'
         PERTSTR(14)='OCTUPXXZ'
         PERTSTR(15)='OCTUPXYY'
         PERTSTR(16)='OCTUPYYZ'
         PERTSTR(17)='OCTUPXZZ'
         PERTSTR(18)='OCTUPYZZ'
         PERTSTR(19)='OCTUPXYZ'
         endif
      ENDIF
C
      IF (IFLAGS(18) .EQ. 8 .OR. IFLAGS(18) .EQ. 9 .OR. 
     &    IFLAGS(18) .EQ. 10 .OR. IFLAGS(18) .EQ. 13) THEN
          YESNO = .TRUE.
       ENDIF
C
       IF (.NOT. YESNO) THEN
C
          WRITE(6,1002)NOPERT
          WRITE(6,*)
C
          DO 104 IPERT=1,NOPERT
             WRITE(6,1003)IPERT,PERTSTR(IPERT)
 104      CONTINUE
          WRITE(6, *)
C
        ENDIF
C
C
C  LOGIC TO DETERMINE COMMON BLOCK /DRVHBAR/
C
      SS = .TRUE.
      SD = .TRUE.
      DS = .TRUE.
      IF (IFLAGS2(161) .EQ. 1) THEN
        DD = .FALSE.
      ELSE
        DD = .TRUE.
      ENDIF
C
C  DETERMINE NODAVID. FIRST DETERMINE REALFREQ
C
      REALFREQ = .TRUE.
      IF (IFLAGS(18) .EQ. 11) THEN
C
        INQUIRE(FILE='frequency',EXIST=FREQ_EXIST)
        IF (FREQ_EXIST) THEN
           OPEN(UNIT=10, FILE='frequency', FORM='FORMATTED', 
     &          STATUS='OLD')
           IJUNK = 0
  100      READ(10, *, END=800) FREQ
           READ(10, *, END=999) IJUNK
  999      REALFREQ = (REALFREQ .AND. IJUNK .EQ. 0)
           IF (REALFREQ) GOTO 100
  800      CLOSE(UNIT=10, STATUS='KEEP')
       ELSE
            FREQ = 0.0D0
       ENDIF 
C
        IF (REALFREQ) THEN
        write(6,*) ' @SETMET: All Frequencies to be considered are real'
        ELSE
          write(6,*) ' @SETMET: At least one imaginary frequency '
        ENDIF
C
        ENDIF
C
C  NOW DETERMINE NODAVID
C
        NODAVID = PARTEOM .AND. REALFREQ
C
        IF (NODAVID) THEN
       write(6,*) ' @SETMET: Linear equations solved using partitioning'
        ELSE
          write(6,*) ' @SETMET: Linear equations solved using Davidson '
        ENDIF
C
 1001  FORMAT(T3,'@SETMET-I, User supplied perturbations will be used.')
 1002    FORMAT(T3,'There are ',I2,' perturbations:.')
 1003    FORMAT(T5,'Perturbation ',I2,' is ',A8,'.')
C
      RETURN
      END 
