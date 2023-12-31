      SUBROUTINE RN_INTS_NEW(NCNTR,NATOMS,NFCT,NBASP,NBAS,LNP1,IPOPF,
     &  IUATMS,NAPTS,NRADPT,CUTOFF,COORD,TOTWT,ALPHA,PCOEF,XGRDPT,
     &  YGRDPT,ZGRDPT,NUNQSHL,NSHL,NANGMOMSHL,NCONFUNSHL,NPRIMFUNSHL,
     &  NANGMOMTSHL,NOFFSETATMP,NOFFSETATMC,NCONFUNTSHL,ISHL2CNTR_MAP,
     &  NPRMFUNTSHL,NOFFSETPRM,NOFFSETCON,NOFFSETSHL,REPLINT,PRDUTINT,
     &  KINTINT,TMP1,TMP2,TMP3,SCFDENT,RELDENT,SCFDENSD,RELDEND,ANGWT,
     &  RADWT,ATOMCHRG,WORK,ILEFT,NTOTSHL,IFACT,ONEP_ENERGY,
     &  TWO_POTNL_R)
C     
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (MAXCENT = 20, MAXGRD = 100)
      PARAMETER (MAXPAIRS = 190)
      DOUBLE PRECISION KINTINT(NBAS,NBAS), KINEP_R, KINE_POTNL_R
C
      INTEGER ATOMCHRG(NATOMS),NPOINT_PAIR(MAXPAIRS)
C
      CHARACTER*8 PGRP, FILENAME1, FILENAME2, FILENAME4, FILENAME5, 
     &            FILENAME6, FILENAME7
      CHARACTER*14 FILENAME3
C
      LOGICAL FILE_EXIST, COR_EXCH_HOLE, OVERLAP_POTN, TWOPE_POT, 
     &        ZMP_POT, ANALYTIC_INT, DO_PROPS, HALT,SPHERICAL
C
      COMMON /PAR/ PI,PIX4
      COMMON /IPAR/ LUOUT, LUNITA, LUNITB
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /FLAGS2/ IFLAGS2(500)
C
      DATA BTOA, HALF, ONE, TWO, FOUR, TEN, ZLICH /0.529177249D0,
     &      0.5D0, 1.0D0, 2.00D0, 4.00D0, 10.0D00, 0.0D00/
C
      DIMENSION ALPHA(NFCT),XGRDPT(NATOMS,NRADPT,194),IPOPF(IUATMS),
     &  YGRDPT(NATOMS,NRADPT,194),ZGRDPT(NATOMS,NRADPT,194),
     &  COORD(NATOMS*3),NSHL(IUATMS),NANGMOMSHL(IUATMS,NUNQSHL),
     &  NCONFUNSHL(IUATMS,NUNQSHL),NPRIMFUNSHL(IUATMS,NUNQSHL),
     &  PCOEF(NFCT,NBAS),REPLINT(NBAS,NBAS),PRDUTINT(NBAS,NBAS),
     &  TMP1(LNP1,LNP1),TMP2(LNP1,LNP1),TMP3(LNP1,LNP1),WORK(ILEFT),
     &  SCFDENT(NBAS,NBAS),RELDENT(NBAS,NBAS),SCFDENSD(NBAS,NBAS),
     &  RELDEND(NBAS,NBAS),NCONFUNTSHL(NTOTSHL),NANGMOMTSHL(NTOTSHL),
     &  NOFFSETPRM(IUATMS,NUNQSHL),NOFFSETCON(IUATMS,NUNQSHL),
     &  NOFFSETSHL(NTOTSHL),NOFFSETATMP(IUATMS),NOFFSETATMC(IUATMS),
     &  TOTWT(NRADPT,194),UVALUE(50),EXC_COR_HOLE(50),ANGWT(194),
     &  RADWT(50),DFT_EXC_COR_HOLE(50),ISHL2CNTR_MAP(NTOTSHL),
     &  NPRMFUNTSHL(NTOTSHL),DENSIN(NBASP,NBASP),ONER1INT(NBAS,NBAS),
     &  ONER2INT(NBAS,NBAS),ONER3INT(NBAS,NBAS),RTO2INT(NBAS,NBAS),
     &  RTO3INT(NBAS,NBAS),RTO1INT(NBAS,NBAS),RTO0INT(NBAS,NBAS),
     &  SPH2CART(NBAS,NBASP),TMP(NBAS,NBAS),
     &  DENS(NBAS*NBAS)

C
      DIMENSION CC_ARRAY(MAXCENT,MAXGRD*MAXPAIRS),
     &          IRREP_PERT(8), CC_TENSOR_ATR(MAXCENT), 
     &          GRID_POINTS(3, MAXPAIRS*MAXGRD), 
     &          CC_TENSOR(MAXPAIRS, MAXPAIRS*MAXGRD)
C     
C Initialize integral values and symmetry multiplication factor
C
      SPHERICAL = (IFLAGS(62) .EQ. 1)
      FACT   = REAL(IFACT)
      TWOPE_POT = .TRUE.
C
      OPEN(80,FILE='plot_powers_r.csv')
      WRITE(80,*) ' NRADPT R <r^-2> <r^-1> <r^0> <r^1> <r^2> <r^3>'
 19   CLOSE (IUNIT3)

      if(SPHERICAL) then
      CALL ZERO(SPH2CART,NBAS*NBASP)
      CALL ZERO(TMP,NBAS*NBAS)
      CALL ZERO(DENS,NBAS*NBAS)
      CALL GETREC(20,'JOBARC','CMP2CART',NBAS*NBASP,SPH2CART)
      CALL GETREC(-20,'JOBARC','AODENSAB',NBASP*NBASP,DENSIN)
      CALL TRANS_SPH2CART(DENSIN,DENS,SPH2CART,TMP,NBAS,NBASP)
      else
      CALL GETREC(-20,'JOBARC','AODENSAB',NBAS*NBAS,DENS)
      endif

      IF (TWOPE_POT) THEN
C
C I've had trouble with arrays picking up garbage when not zeroed 
C before being used, so all the arrays to be used here are zeroed before 
C entering the loops
C
      CALL ZERO(PRDUTINT,   NBAS*NBAS)
      CALL ZERO(ONER1INT,   NBAS*NBAS)
      CALL ZERO(ONER2INT,   NBAS*NBAS)
      CALL ZERO(ONER3INT,   NBAS*NBAS)
      CALL ZERO(RTO0INT,    NBAS*NBAS)  
      CALL ZERO(RTO1INT,    NBAS*NBAS)  
      CALL ZERO(RTO2INT,    NBAS*NBAS)  
      CALL ZERO(RTO3INT,    NBAS*NBAS)  
c
      IATMOS = 1
      DO IRADPT = 1, NRADPT
         DO IANGPT = 1, NAPTS

            XPOINT = XGRDPT(NCNTR, IRADPT, IANGPT)  + 
     &               COORD((IATMOS-1)*3 + 1)
            YPOINT = YGRDPT(NCNTR, IRADPT, IANGPT)  + 
     &               COORD((IATMOS-1)*3 + 2)
            ZPOINT = ZGRDPT(NCNTR, IRADPT, IANGPT)  +
     &               COORD((IATMOS-1)*3 + 3)
C
            R = DSQRT(XPOINT**2 + YPOINT**2 + ZPOINT**2)
C
            CALL PRDCTINT(IUATMS, NATOMS, NUNQSHL, LNP1, IPOPF,
     &                    NSHL, NFCT, NBAS, NOFFSETATMP,
     &                    NOFFSETATMC, NOFFSETPRM, NOFFSETCON,
     &                    NANGMOMSHL, NPRIMFUNSHL, NCONFUNSHL,
     &                    XPOINT, YPOINT, ZPOINT, ALPHA, PCOEF,
     &                    COORD, REPLINT, PRDUTINT, TMP1, TMP2,
     &                    TMP3)
C
C The PRDUTINT array is copied in to different arrays to calculate the 
C different powers of R
C ONE_R  = 1/R
C ONE_R2 = 1/R^2
C ONE_R3 = 1/R^3
C RTO2   = R^2
C RTO3   = R^3
C
             ONER1INT = PRDUTINT
             ONER2INT = PRDUTINT
             ONER3INT = PRDUTINT
             RTO0INT  = PRDUTINT
             RTO1INT  = PRDUTINT
             RTO2INT  = PRDUTINT
             RTO3INT  = PRDUTINT
C
C Now they are scaled accordingly
C
             CALL DSCAL(NBAS*NBAS,1.d0/R,ONER1INT, 1)
             CALL DSCAL(NBAS*NBAS,1.d0/R**2,ONER2INT, 1)
             CALL DSCAL(NBAS*NBAS,1.d0/R**3,ONER3INT, 1)
             CALL DSCAL(NBAS*NBAS,1.d0,RTO0INT, 1)
             CALL DSCAL(NBAS*NBAS,R,RTO1INT, 1)
             CALL DSCAL(NBAS*NBAS,R**2,RTO2INT, 1)
             CALL DSCAL(NBAS*NBAS,R**3,RTO3INT, 1)
C
C Calculating the value of r^n for a certain r
C ONER1 = 1/r
C ONER2 = 1/r^2
C ONER3 = 1/r^3
C RTO1  = r
C RTO2  = r^2
C RTO3  = r^3
              ONER1 = SDOT(NBAS*NBAS, SCFDENT, 1, ONER1INT, 1)
              ONER2 = SDOT(NBAS*NBAS, SCFDENT, 1, ONER2INT, 1)
              ONER3 = SDOT(NBAS*NBAS, SCFDENT, 1, ONER3INT, 1)
              RTO0  = SDOT(NBAS*NBAS, SCFDENT, 1, RTO0INT, 1)
              RTO1  = SDOT(NBAS*NBAS, SCFDENT, 1, RTO1INT, 1)
              RTO2  = SDOT(NBAS*NBAS, SCFDENT, 1, RTO2INT, 1)
              RTO3  = SDOT(NBAS*NBAS, SCFDENT, 1, RTO3INT, 1)
              TOT_DENSITY = RTO0
              IF (TOT_DENSITY .LT. CUTOFF) GOTO 10

C Calculating the expectation values
C ONER1EXP  = <1/r>
C ONER2EXP  = <1/r^2>
C ONER3EXP  = <1/r^3>
C R1EXP     = <r>
C R2EXP     = <r^2>
C R3EXP     = <r^3>
 
              ONER1EXP  = ONER1EXP  + TOTWT(IRADPT,IANGPT)*ONER1
              ONER2EXP  = ONER2EXP  + TOTWT(IRADPT,IANGPT)*ONER2
              ONER3EXP  = ONER3EXP  + TOTWT(IRADPT,IANGPT)*ONER3
              R0EXP     = R0EXP     + TOTWT(IRADPT,IANGPT)*RTO0
              R1EXP     = R1EXP     + TOTWT(IRADPT,IANGPT)*RTO1
              R2EXP     = R2EXP     + TOTWT(IRADPT,IANGPT)*RTO2
              R3EXP     = R3EXP     + TOTWT(IRADPT,IANGPT)*RTO3
C
C Loop over atoms end here!
C
            ENDDO
C
C The loop over radial points end here!
C
              ONER1PIR2  = FOUR*PI*(R**2)*ONER1
              ONER2PIR2  = FOUR*PI*(R**2)*ONER2
              ONER3PIR2  = FOUR*PI*(R**2)*ONER3
              R0PIR2     = FOUR*PI*(R**2)*RTO0
              R1PIR2     = FOUR*PI*(R**2)*RTO1
              R2PIR2     = FOUR*PI*(R**2)*RTO2
              R3PIR2     = FOUR*PI*(R**2)*RTO3

      WRITE(80,210)IRADPT,R,ONER2PIR2,ONER1PIR2,R0PIR2,R1PIR2,
     &             R2PIR2,R3PIR2

         ENDDO
C
 210   FORMAT(1x,I2,1x,F15.7,3x,F15.7,3x,F15.7,3x,F15.7,3x,F15.7,
     &    3x,F15.7,3x,F15.7)
 10   CONTINUE
C     
      ENDIF
C
      OPEN(100,FILE='powers_r.txt')
      WRITE(100,*) ' n            <r^n>'
      WRITE(100,200) -3,ONER3EXP
      WRITE(100,200) -2,ONER2EXP 
      WRITE(100,200) -1,ONER1EXP 
      WRITE(100,200) 0,R0EXP
      WRITE(100,200) 1,R1EXP  
      WRITE(100,200) 2,R2EXP  
      WRITE(100,200) 3,R3EXP  
 200  FORMAT(1x,I2,8x,F15.7)
      CLOSE(100)
C        
      ONER1 = SDOT(NBAS*NBAS, SCFDENST, 1, TOTONER1INT, 1)
c        print*,oner1
c        do i=1,nbas
c          do j=1,nbas
c            print*,i,j,totoner1int(i,j),dens(i,j)
c          enddo
c        enddo
C
      RETURN
      END
 
