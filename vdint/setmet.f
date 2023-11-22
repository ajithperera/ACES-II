      SUBROUTINE SETMET
C
C THIS ROUTINE FILLS THE COMMON BLOCK /METHOD/ WHICH 
C CONTAINS INFORMATION ABOUT THE METHODS FOR WHICH 
C GRADIENTS ARE CALCULATED
C
C/METHOD/
C
C  IUHF ..... RHF/UHF FLAG
C  SCF ...... SCF/CORRELATED METHOD FLAG
C  NONHF .... NONHF FLAG FOR CORRELATED METHODS
C
C/OPTION/
C
C  GRAD .....
C  INTWRIT .. TRUE IF DERIVATIVE INTEGRALS ARE WRITTEN TO DISK
C  FOCK ..... TRUE IF FOCK MATRIX DERIVATIVES ARE FORMED
C
C/ROHF/
C
C  ROHF ..... TRUE FOR ROHF REFERENCE FUNCTIONS
C  SEMI ..... TRUE IF SEMICANONICAL ORBITALS ARE USED IN ROHF CALCULATIONS
C
C/DFGH/
C 
C  IDFGH .... TRUE IF SPHERICAL GAUSSIANS ARE USED
C
C/TREATP
C
C  IFIRST ... FIRST ENTRY TO VDINT
C  XCOMP .... TRUE, IF IN THIS PASS CONSIDERED
C  YCOMP .... TRUE, IF IN THIS PASS CONSIDERED
C  ZCOMP .... TRUE, IF IN THIS PASS CONSIDERED
C  BF, BL ... RANGE FOR LOOP OVER PERTURBATIONS
C
CEND
      IMPLICIT INTEGER(A-Z)
      DOUBLE PRECISION AZERO,GX,GY,GZ 
c&line mod
      LOGICAL SCF,NONHF,ROHF,SEMI
      LOGICAL ELECT,MAGNET,SPIN,GEOM,GRAD,INTWRIT,FOCK,IDFGH
      LOGICAL MAGNET2,GIAO,YESNO,IFIRST,XCOMP,YCOMP,ZCOMP,IALL
      LOGICAL JSO,JFC,JSD
      LOGICAL GEERTSEN,CCSDEH
      INTEGER DIRPRD,POPRHF,VRTRHF,POPDOC,VRTDOC
      CHARACTER*80 FNAME
      CHARACTER*6 STRING
      COMMON/METHOD/IUHF,SCF,NONHF
c&line del      COMMON/QRHF/QRHFP
      COMMON/OPENSH/ROHF,SEMI
      COMMON/SYMINF/NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON/PROP/ELECT,MAGNET,GEOM,MAGNET2,GIAO
      COMMON/OPTION/GRAD,INTWRIT,FOCK
      COMMON/FLAGS/IFLAGS(100)
      COMMON/DFGH/IDFGH
      COMMON/QRHFINF/POPRHF(8),VRTRHF(8),NOSH1(8),NOSH2(8),
     &               POPDOC(8),VRTDOC(8),NAI,N1I,N2A,
c&line mod
     &               NUMISCF,NUMASCF,ISPINP,ISPINM,IQRHF
      COMMON/TREATP/IFIRST,IALL,XCOMP,YCOMP,ZCOMP,IBF,IBL
      COMMON/JNMR/JSO,JFC,JSD
      COMMON/GORIGIN/GX,GY,GZ
      COMMON/DIAMAG/GEERTSEN
      COMMON/EH/CCSDEH
      DATA AZERO/0.0D0/
C
      GX=AZERO
      GY=AZERO
      GZ=AZERO
C
      JFC=.FALSE.
      JSD=.FALSE.
      JSO=.FALSE.
      ROHF=.FALSE.
      SEMI=.FALSE.
      GRAD=.TRUE.
      INTWRIT=.FALSE.
      FOCK=.FALSE.
      IDFGH=.FALSE.
      CCSDEH=.FALSE.
C
      itmp=iflags(18)
      IF(ITMP.GE.100) CCSDEH=.TRUE.
      IF(ITMP.GT.11) CCSDEH=.TRUE.
c      IF(IFLAGS(91).NE.0) CCSDEH=.TRUE.
      iflags(18)=mod(iflags(18),100)
C
      GIAO=.TRUE.
      GEERTSEN=.FALSE.
C
      IFIRST=.TRUE.
      ELECT=.FALSE.
      MAGNET=.FALSE.
      MAGNET2=.FALSE.
      SPIN=.FALSE.
      GEOM=.TRUE.
      IF(IFLAGS(18).EQ.1) THEN
       ELECT=.TRUE.
      ENDIF
      IF(IFLAGS(18).EQ.3.OR.IFLAGS(18).EQ.4.OR.
     &   IFLAGS(18).EQ.5.OR.IFLAGS(18).EQ.12) THEN
C
C CALCULATE NMR SHIFTS, IF IFLAGS(18)=4, CALCULATE IN 
C ADDITION MAGNETIZABILITIES
C
       MAGNET=.TRUE.
       IF(IFLAGS(18).EQ.5.OR.IFLAGS(18).EQ.12) THEN
        MAGNET2=.TRUE.
        GIAO=.FALSE.
        IF(IFLAGS(18).EQ.12) GEERTSEN=.TRUE.
        IF(IFLAGS(18).EQ.5) JSO=.TRUE.
       ENDIF
       GEOM=.FALSE.
      ENDIF
      IF(IFLAGS(18).EQ.4) THEN
       SPIN=.TRUE.
      ENDIF
      IF(IFLAGS(62).NE.0) THEN
       IDFGH=.TRUE.
      ENDIF
C
C FLAGS FOR CALCULATION OF NMR COUPLING CONSTANTS
C
      IF(IFLAGS(18).EQ.8) THEN
       MAGNET=.TRUE.
       JSO=.TRUE.
       GEOM=.FALSE.
      ENDIF
      IF(IFLAGS(18).EQ.9.OR.IFLAGS(18).EQ.10) THEN
       GEOM=.FALSE.
       IF(IFLAGS(18).EQ.9) JFC=.TRUE.
       IF(IFLAGS(18).EQ.10) JSD=.TRUE.
      ENDIF
C
C  DERIVATIVE LEVEL
C
      DERIV=IFLAGS(3)
      write(6,100)
C
C  SET IUHF FLAG
C
      IUHF=IFLAGS(11)
C
C  SET SCF FLAG
C
      SCF=IFLAGS(2).EQ.0.AND.IFLAGS(87).NE.1
C
C  SET NONHF FLAG
C
      NONHF=(IFLAGS(77).NE.0)
c&lines  del and modif
      IF(NONHF) CALL QRHFSET
C
      IF(SCF) THEN
       IF(IUHF.EQ.0) THEN
        IF(DERIV.EQ.1) write(6,1000)
        IF(DERIV.EQ.2.AND.(.NOT.MAGNET).AND.(.NOT.JFC).AND.(.NOT.JSD))
     &        write(6,1010)
        IF(DERIV.EQ.2.AND.MAGNET) THEN
         write(6,1011)
         IF(GIAO) THEN
          write(6,6200)
         ELSE
          write(6,6210) GX,GY,GZ
         ENDIF
        ENDIF
        IF(DERIV.EQ.2.AND.JFC) THEN
         write(6,6211)
        ENDIF
        IF(DERIV.EQ.2.AND.JSD) THEN
         write(6,6212)
        ENDIF
       ELSE IF(IUHF.EQ.1) THEN
        IF(DERIV.EQ.1) write(6,2000)
        IF(DERIV.EQ.2) write(6,2010)
       ELSE IF(IUHF.EQ.2) THEN
        ROHF=.TRUE.
        IF(DERIV.EQ.1) write(6,2100)
        IF(DERIV.EQ.2) write(6,2110)
        IF(IFLAGS(39).NE.0) THEN
         SEMI=.TRUE.
         write(6,2111)
        ELSE
         write(6,2112)
        ENDIF
        IUHF=1
       ELSE
        STOP 
       ENDIF 
       IF(DERIV.GE.2) FOCK=.TRUE.
      ELSE
       IF(NONHF) THEN
         IF(DERIV.EQ.1)write(6,3000)
         IF(DERIV.EQ.2)write(6,3001)
       ELSE
        IF(IUHF.EQ.0) THEN
         IF(DERIV.EQ.1)write(6,4000)
         IF(DERIV.EQ.2.AND.(.NOT.MAGNET))write(6,4001)
         IF(DERIV.EQ.2.AND.MAGNET) THEN
          write(6,4011)
          IF(GIAO) THEN
           write(6,6200)
          ELSE
           write(6,6210)
          ENDIF
         ENDIF
        ELSE IF(IUHF.EQ.1) THEN
         IF(DERIV.EQ.1)write(6,5000)
         IF(DERIV.EQ.2)write(6,5001)
        ELSE IF(IUHF.EQ.2) THEN
         IF(DERIV.EQ.1)write(6,6000)
         IF(DERIV.EQ.2)write(6,6001)
         IUHF=1
        ENDIF
        IF(DERIV.EQ.2) THEN
         IF(.NOT.MAGNET) THEN
          CALL GFNAME('GAMLAM  ',FNAME,ILENGTH)
          INQUIRE(FILE=FNAME(1:ILENGTH),EXIST=YESNO)
          IF(YESNO) THEN
           GRAD=.TRUE.
           FOCK=.TRUE.
           INTWRIT=.TRUE.
          ENDIF
         ELSE
          FOCK=.TRUE.
          INTWRIT=.TRUE.
          GRAD=.TRUE.
         ENDIF
        ENDIF
       ENDIF
      ENDIF
      IF(IDFGH) THEN 
       write(*,7001)
      ELSE
       write(*,7002)
      ENDIF
C
C DEAL WITH THE CASE WHERE MAGNETIC PERTURBATIONS ARE CONSIDERED
C SEPARATELY
C
      IF(MAGNET.AND.(.NOT.SCF)) THEN
C
       IF(IFLAGS(90).EQ.0) THEN
        XCOMP=.TRUE.
        YCOMP=.TRUE.
        ZCOMP=.TRUE.
        IALL=.TRUE.
        IBF=1
        IBL=3
       ELSE
        XCOMP=.FALSE.
        YCOMP=.FALSE.
        ZCOMP=.FALSE.
        IALL=.FALSE.
        IONE=1
        CALL GETREC(0,'JOBARC','TREATPER',IONE,IPERNUM)
        IF(IONE.LE.0) THEN
         IFIRST=.TRUE. 
         IPERNUM=0
         IONE=1
        ELSE
         IFIRST=.FALSE.
         CALL GETREC(20,'JOBARC','TREATPER',IONE,IPERNUM)
         IF(IPERNUM.EQ.0) THEN
          IFIRST=.TRUE.
         ENDIF
        ENDIF
        IPERNUM=IPERNUM+1
        IF(IPERNUM.EQ.1) THEN
         XCOMP=.TRUE.
         write(*,6110)
        ELSE IF(IPERNUM.EQ.2) THEN
         YCOMP=.TRUE.
         write(*,6120)
        ELSE IF(IPERNUM.EQ.3) THEN
         ZCOMP=.TRUE.
         write(*,6130)
        ENDIF
        IBF=IPERNUM
        IBL=IPERNUM
        CALL PUTREC(20,'JOBARC','TREATPER',IONE,IPERNUM)
       ENDIF
      ENDIF
      IF(MAGNET.AND..NOT.GIAO) THEN
       OPEN(UNIT=30,FILE='ZMAT',FORM='FORMATTED')
1      READ(30,'(A)',END=102)STRING
       IF(STRING.NE.'%GAUGE'.AND.
     &    STRING.NE.'%gauge') GO TO 1
C
        READ(30,*) GX,GY,GZ
C
102    CONTINUE
       WRITE(*,6213) GX,GY,GZ
      ENDIF
C
100   FORMAT(' One- and two-electron integral derivatives',
     &       ' are calculated')
1000  FORMAT(' for RHF gradients and dipole moments.')
1010  FORMAT(' for RHF hessians and dipole derivatives.')
1011  FORMAT(' for RHF NMR chemical shifts.')
2000  FORMAT(' for UHF gradients and dipole moments.')
2010  FORMAT(' for UHF hessians and dipole derivatives.')
2100  FORMAT(' for ROHF gradients and dipole moments.')
2110  FORMAT(' for ROHF hessians and dipole',
     &       ' derivatives.')
2111  FORMAT(' Semicanonical orbitals are used.')
2112  FORMAT(' Standard orbitals are used.')
3000  FORMAT(' for QRHF-CCSD gradients',
     &       ' and dipole moments.')
4000  FORMAT(' for RHF-CC/MBPT gradients',
     &       ' and dipole moments.')
5000  FORMAT(' for UHF-CC/MBPT gradients',
     &       ' and dipole moments.')
6000  FORMAT(' for ROHF-CC/MBPT gradients',
     &       ' and dipole moments.')
3001  FORMAT(' for QRHF-CCSD hessians',
     &       ' and dipole derivatives.')
4001  FORMAT(' for RHF-CC/MBPT hessians',
     &       ' and dipole derivatives.')
4011  FORMAT(' for RHF-CC/MBPT NMR',
     &       ' chemical shifts.')
5001  FORMAT(' for UHF-CC/MBPT hessians',
     &       ' and dipole derivatives.')
6001  FORMAT(' for ROHF-CC/MBPT hessians',
     &       ' and dipole moments.')
7001  FORMAT(' Spherical gaussians are used.')
7002  FORMAT(' Cartesian gaussians are used.')
6110  FORMAT(' Derivatives with respect to Bx are calculated.')
6120  FORMAT(' Derivatives with respect to By are calculated.')
6130  FORMAT(' Derivatives with respect to Bz are calculated.')
6200  FORMAT(' Gauge including atomic orbitals (GIAOs) are used.')
6209  FORMAT(' Field-independent atomic orbitals are used.',
     &       ' The gauge-origin is the center of the chosen coordinate',
     &       ' system.')
6210  FORMAT(' Field-independent atomic orbitals are used.')
6213  FORMAT(' The gauge-origin',3F10.7)
6211  FORMAT(' for the Fermi-contact(FC)',
     &       ' contribution to J.')
6212  FORMAT(' for the spin-dipolar(SD)', 
     &       ' contribution to J.')
C
      iflags(18)=itmp
C
      RETURN
      END
