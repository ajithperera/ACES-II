      SUBROUTINE INTGRT_TWODEN(NCNTR,NATOMS,NFCT,NBASP,NBAS,LNP1,IPOPF,
     &  IUATMS,NAPTS,NRADPT,CUTOFF,COORD,TOTWT,ALPHA,PCOEF,XGRDPT,
     &  YGRDPT,ZGRDPT,NUNQSHL,NSHL,NANGMOMSHL,NCONFUNSHL,NPRIMFUNSHL,
     &  NANGMOMTSHL,NOFFSETATMP,NOFFSETATMC,NCONFUNTSHL,ISHL2CNTR_MAP,
     &  NPRMFUNTSHL,NOFFSETPRM,NOFFSETCON,NOFFSETSHL,REPLINT,PRDUTINT,
     &  KINTINT,TMP1,TMP2,TMP3,SCFDENT,RELDENT,SCFDEND,RELDEND,ANGWT,
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
     &        ZMP_POT, ANALYTIC_INT, DO_PROPS, HALT
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
     &  SCFDENT(NBAS,NBAS),RELDENT(NBAS,NBAS),SCFDEND(NBAS,NBAS),
     &  RELDEND(NBAS,NBAS),NCONFUNTSHL(NTOTSHL),NANGMOMTSHL(NTOTSHL),
     &  NOFFSETPRM(IUATMS,NUNQSHL),NOFFSETCON(IUATMS,NUNQSHL),
     &  NOFFSETSHL(NTOTSHL),NOFFSETATMP(IUATMS),NOFFSETATMC(IUATMS),
     &  TOTWT(NRADPT,194),UVALUE(50),EXC_COR_HOLE(50),ANGWT(194),
     &  RADWT(50),DFT_EXC_COR_HOLE(50),ISHL2CNTR_MAP(NTOTSHL),
     &  NPRMFUNTSHL(NTOTSHL)
C
      DIMENSION CC_ARRAY(MAXCENT,MAXGRD*MAXPAIRS),
     &          IRREP_PERT(8), CC_TENSOR_ATR(MAXCENT), 
     &          GRID_POINTS(3, MAXPAIRS*MAXGRD), 
     &          CC_TENSOR(MAXPAIRS, MAXPAIRS*MAXGRD)
C     
C Initialize integral values and symmetry multiplication factor
C
      FACT   = REAL(IFACT)
      TWO_ENERGY = 0.0D0
      ONH_ENERGY = 0.0D0
      TOL        = 1.0E-5
      COR_EXCH_HOLE = .FALSE.
      OVERLAP_POTN  = .FALSE.
      TWOPE_POT     = .FALSE.
      ZMP_POT       = .FALSE.
      ANALYTIC_INT  = .FALSE.
      DO_PROPS      = .FALSE.
C
      FILENAME2 = "RDENSITY"
      FILENAME3 = "pick_potential"
C
      IUNIT1 = 11 + NCNTR
      IUNIT2 = 22 + NCNTR
      IUNIT3 = 33 
      IUNIT4 = 44 + NCNTR
      IUNIT5 = 55 + NCNTR
      IUNIT6 = 66 + NCNTR
      IUNIT7 = 77 + NCNTR
C
      OPEN (UNIT=IUNIT2, FILE=FILENAME2, FORM="FORMATTED")
C
      INQUIRE(FILE="pick_potential",EXIST=FILE_EXIST)
      IF (FILE_EXIST) THEN
         OPEN (UNIT=IUNIT3, FILE=FILENAME3, FORM="FORMATTED")
         READ (IUNIT3,*,END=19) REPLSN_FACT, EXCHNG_FACT, 
     &                          CORLN_FACT_O, CORLN_FACT_T
         READ (IUNIT3,*,END=19) IFLAG1, IFLAG2, IFLAG3, IFLAG4
      ELSe
         REPLSN_FACT  = 1.00D0
         EXCHNG_FACT  = 1.00D0
         CORLN_FACT_O = 1.00D0
         CORLN_FACT_T = 1.00D0
         IFLAG1 = 0
         IFLAG2 = 1
         IFLAG3 = 0
         IFLAG4 = 0
      ENDIF
C
      IF (IFLAG1 .EQ. 1) THEN
         COR_EXCH_HOLE = .TRUE. 
         FILENAME1 = "COR_HOLE"        
         OPEN (UNIT=IUNIT1, FILE=FILENAME1, FORM="FORMATTED")
      ELSE IF (IFLAG2 .EQ. 1) THEN
         TWOPE_POT = .TRUE.
         FILENAME1 = "VCONDPOT"
         FILENAME4 = "ONEHLPOT"
         FILENAME5 = "VKINTPOT"
         FILENAME6 = "VRPLTPOT"
         FILENAME7 = "CORLNPOT"
         OPEN (UNIT=IUNIT1, FILE=FILENAME1, FORM="FORMATTED")
         OPEN (UNIT=IUNIT4, FILE=FILENAME4, FORM="FORMATTED")
         OPEN (UNIT=IUNIT5, FILE=FILENAME5, FORM="FORMATTED")
         OPEN (UNIT=IUNIT6, FILE=FILENAME6, FORM="FORMATTED")
         OPEN (UNIT=IUNIT7, FILE=FILENAME7, FORM="FORMATTED")
         OVERLAP_POTN = .FALSE.
      ELSE IF (IFLAG3 .EQ. 1) THEN
         ZMP_POT = .TRUE.
      ELSE IF (IFLAG4 .EQ. 1) THEN
             DO_PROPS = .TRUE.
         OVERLAP_POTN = .FALSE.
      ENDIF

      WRITE(6,*) 
      WRITE(6,"(a,a)")" Keys and their values:TWOPE_POT,ZMP_POT",
     &                ",DO_PROPS,OVERLAP_POTN:"
      WRITE(6,"(4(1x,l))") TWOPE_POT,ZMP_POT,DO_PROPS,OVERLAP_POTN

      IF (TWOPE_POT) THEN
         WRITE(6,"(2a)")" The scalling factors for repulsion,exchange",
     &                  ",one-particle correlation and"
         WRITE(6,"(a,4F5.2)") " two-particle correlation:",REPLSN_FACT,
     &                         EXCHNG_FACT,CORLN_FACT_O,CORLN_FACT_T
      ENDIF 
#ifdef _DEBUG_LVLM1
        Write(6,*) 
        Write(6,"(a)") "The weight array"
        Write(6,"(4F15.5)") ((TOTWT(i, j),i=1,NRADPT), j=1,194)
#endif
C
 19   CLOSE (IUNIT3)
C
      IF (COR_EXCH_HOLE .AND. .NOT. ANALYTIC_INT) THEN
         IATMOS = 1
         DO JRADPT = 1, NRADPT

                EXC_COR_HOLE(JRADPT) = 0.0D0
            DFT_EXC_COR_HOLE(JRADPT) = 0.0D0
C
            DO IRADPT = 1, NRADPT
               DO IANGPT = 1, NAPTS

                  XPOINT = XGRDPT(NCNTR, IRADPT, IANGPT)  + 
     &                     COORD((IATMOS-1)*3 + 1)
                  YPOINT = YGRDPT(NCNTR, IRADPT, IANGPT)  + 
     &                     COORD((IATMOS-1)*3 + 2)
                  ZPOINT = ZGRDPT(NCNTR, IRADPT, IANGPT)  +
     &                     COORD((IATMOS-1)*3 + 3)
C
                  R = DSQRT(XPOINT**2 + YPOINT**2 + ZPOINT**2)
C
                  CALL ZERO(REPLINT, NBAS*NBAS)
                  POTNL_R = 0.0D0
C
                  DO JANGPT = 1, NAPTS
C
                    IF (JANGPT .EQ. 1) THEN
C     
                        APOINT = XGRDPT(NCNTR, JRADPT, JANGPT)  + 
     &                           COORD((IATMOS-1)*3 + 1) - XPOINT
                        BPOINT = YGRDPT(NCNTR, JRADPT, JANGPT)  + 
     &                           COORD((IATMOS-1)*3 + 2) - YPOINT
                        CPOINT = ZGRDPT(NCNTR, JRADPT, JANGPT)  +
     &                           COORD((IATMOS-1)*3 + 3) - ZPOINT

                        U = DSQRT(APOINT**2 + BPOINT**2 + CPOINT**2)
C
                        APOINT = XGRDPT(NCNTR, JRADPT, JANGPT)  + 
     &                           COORD((IATMOS-1)*3 + 1) 
                        BPOINT = YGRDPT(NCNTR, JRADPT, JANGPT)  + 
     &                           COORD((IATMOS-1)*3 + 2) 
                        CPOINT = ZGRDPT(NCNTR, JRADPT, JANGPT)  +
     &                           COORD((IATMOS-1)*v3 + 3) 
C     
                        CALL PRDCTINT(IUATMS, NATOMS, NUNQSHL, LNP1,
     &                                IPOPF, NSHL, NFCT, NBAS, 
     &                                NOFFSETATMP, NOFFSETATMC, 
     &                                NOFFSETPRM, NOFFSETCON,
     &                                NANGMOMSHL, NPRIMFUNSHL,
     &                                NCONFUNSHL, APOINT, BPOINT,
     &                                CPOINT, ALPHA, PCOEF, COORD, 
     &                                REPLINT, PRDUTINT, TMP1, TMP2, 
     &                                TMP3)
                    ENDIF
C
                     IF (U .GE. TEN) GOTO 50
                     IF (IRADPT .EQ. 1 .AND. IANGPT .EQ.1 .AND. 
     &                   JANGPT .EQ. 1) UVALUE(JRADPT) = U                 
C
                     ANGWEIGHT = ANGWT(JANGPT)
                     CALL DAXPY(NBAS*NBAS, ANGWEIGHT, PRDUTINT, 1,
     &                          REPLINT, 1)
C
                  ENDDO

 50               CONTINUE
C
                  CALL PRDCTINT(IUATMS, NATOMS, NUNQSHL, LNP1, IPOPF,
     &                          NSHL, NFCT, NBAS, NOFFSETATMP,
     &                          NOFFSETATMC, NOFFSETPRM, NOFFSETCON,
     &                          NANGMOMSHL, NPRIMFUNSHL, NCONFUNSHL,
     &                          XPOINT, YPOINT, ZPOINT, ALPHA, PCOEF,
     &                          COORD, REPLINT, PRDUTINT, TMP1, TMP2, 
     &                          TMP3)
C
                  CALL BLTCORLNPOT(IUHF, NTOTSHL, NBAS, NANGMOMTSHL,
     &                             NCONFUNTSHL, NOFFSETSHL, ILEFT,
     &                             REPLINT, PRDUTINT, SCFDENT,
     &                             RELDENT, SCFDEND, RELDEND, WORK,
     &                             POTNL_R, REPLSN_FACT, EXCHNG_FACT,
     &                             CORLN_FACT_O, CORLN_FACT_T)
C
                  EXC_COR_HOLE(JRADPT) = EXC_COR_HOLE(JRADPT) +  
     &                                   POTNL_R*TOTWT(IRADPT, IANGPT)
C
C To get exchange-correlation plots for traditional density functional potentials
C one need to calcualte the following expression. 
C
C     K_F      = [3*PI**2Rho(R)]^{1/3}
C     X        = K_F*U
C     DFT_HOLE = -3/PI{Int[Rho(R)*K_F*{SIN(X)/X^2 - COS(X)/X}^2]}
C
C where Rho(R) is the density at R (calculated using DFT orbitals), the
C variables have there usual meaning. These formulas were given to 
C me by Thomas Henderson on 12/2001. 
C
                  DENSITY = SDOT(NBAS*NBAS, SCFDENT, 1, PRDUTINT, 1)
                  DK_F    = (3*(PI**2)*DENSITY)**(1.0D0/3.0D0)
                  X       = DK_F*UVALUE(JRADPT)
                 
C$$                  DFT_EXC_COR_HOLE(JRADPT)=DFT_EXC_COR_HOLE(JRADPT)   +
C$$     &                                     (DENSITY*DK_F*(DSIN(X)/X**2-
C$$     &                                     DCOS(X)/X)**2)*
C$$     &                                     TOTWT(IRADPT, IANGPT)
C$$$                  DFT_EXC_COR_HOLE(JRADPT) = DFT_EXC_COR_HOLE(JRADPT)*
C$$$     &                                       (-3.0D0/PI*UVALUE(JRADPT))

               ENDDO
            ENDDO
C
         ENDDO
C
      ENDIF
C 
      IF (TWOPE_POT) THEN

C
         DO IRADPT = 1, NRADPT
CSSS         DO IRADPT = 1, 1
C
            POTNL_R = 0.0D+00
C
             DO IANGPT = 1, NAPTS
CSSS             DO IANGPT = 1, 1
C
C The following two subroutine calls calculate the integral 
C
C   V(R) = Sum GAMM(u,n,ls)*[Chi_u(R)*Chi_l(R)]/Rho(R)]*[Int{Chi_n(r)*Chi_s(r)/|r - R|dr}]
C          u,n,l,s
C
C which is the  correlation potential at point R. Here GAMMA is the correlated
C two particle density matrix and Chi's are Gaussian basis function. For convenince,
C we designate [Int{Chi_n(r)*Chi_s(r)/|r - R|dr}] as A-Term [Chi_u(R)*Chi_l(R)]/Rho(R)] 
C as B-Term. Both terms are calculated in  REPPRDINT and stored in the arrays
C PRDUTINT and REPLINT. Implemented in 01/2001 by Ajith Perera. 
C 
C The following loop goes over all the atoms. Note that the grid points are
C (XGRDPT, YGRDPT, ZGRDPT arrays) centered on symmetry uniquie atoms 
C in full point group (NCNTR). 
C
C$$$   DO IATMOS = 1, NATOMS
C               
C Let's get an arbitrary point on space for this unique center! Note that for
C symmetry redundent atoms the final integral has the same magnitude and not
C explicitly calcualted. At the end of the loop approprate scaling is done 
C based on the number of symmetry redundent atoms for this particular symm.
C unique atom!!
C
C$$$  XPOINT = 0.000
C$$$  YPOINT = 0.0D0
C$$$  ZPOINT = 0.0D0
C$$$
               IATMOS = 1 
C     
               IF (IANGPT .EQ. 1) THEN
C     
                  XPOINT = XGRDPT(NCNTR, IRADPT, IANGPT)  + 
     &                     COORD((IATMOS-1)*3 + 1)
                  YPOINT = YGRDPT(NCNTR, IRADPT, IANGPT)  + 
     &                     COORD((IATMOS-1)*3 + 2)
                  ZPOINT = ZGRDPT(NCNTR, IRADPT, IANGPT)  +
     &                    COORD((IATMOS-1)*3 + 3)
C
                  R = DSQRT(XPOINT**2 + YPOINT**2 + ZPOINT**2)

CSSS                  Print*, "Entering REPPRINDT"                
CSSS                 XPOINT = -4.000
CSSS                 YPOINT = -4.0D0
CSSS                 ZPOINT = 0.0D0
C$$$
#ifdef _DEBUG_LVLM
              If(IRADPT .le. 10 .and. IANGPT .le. 10)
     &        Write(6,"(F10.5)")   XGRDPT(NCNTR, IRADPT, IANGPT),
     &                             YGRDPT(NCNTR, IRADPT, IANGPT),
     &                             ZGRDPT(NCNTR, IRADPT, IANGPT)
#endif
                  CALL REPPRDINT(IUATMS, NATOMS, NUNQSHL, LNP1, IPOPF,
     &                           NSHL, NFCT, NBAS, NOFFSETATMP,
     &                           NOFFSETATMC, NOFFSETPRM, NOFFSETCON, 
     &                           NANGMOMSHL, NPRIMFUNSHL, NCONFUNSHL, 
     &                           XPOINT, YPOINT, ZPOINT, ALPHA, PCOEF, 
     &                           COORD, REPLINT, PRDUTINT, KINTINT, 
     &                           TMP1, TMP2, TMP3, OVERLAP_POTN)

                  CALL BLTCORLNPOT(IUHF, NTOTSHL, NBAS, NANGMOMTSHL,
     &                             NCONFUNTSHL, NOFFSETSHL, ILEFT,
     &                             REPLINT, PRDUTINT, SCFDENT,
     &                             RELDENT, SCFDEND, RELDEND, WORK,
     &                             POTNL_R, REPLSN_FACT, EXCHNG_FACT,
     &                             CORLN_FACT_O, CORLN_FACT_T)
C
C Let's do something smart here. The PRDUTINT array contains the
C value of a Gaussian product at a given point. So, we should
C be able to calculate the density at the point by taking the 
C dot product of the PRDUTINT with the total density(RELDENT). 
C 
#ifdef _DEBUG_LVLM
                  Write(6,*) "The SCF density matrix"
                  CALL OUTPUT(SCFDENT, 1, NBAS, 1, NBAS, NBAS, NBAS, 1)
                  Write(6,*) "The total (SCF + CORR) density matrix"
                  CALL OUTPUT(RELDENT, 1, NBAS, 1, NBAS, NBAS, NBAS, 1)
#endif 
                  TOT_DENSITY = SDOT(NBAS*NBAS, RELDENT,1, PRDUTINT, 1)
                  SCF_DENSITY = SDOT(NBAS*NBAS, SCFDENT,1, PRDUTINT, 1)
                  COR_DENSITY = TOT_DENSITY - SCF_DENSITY
C
C Let's calculate the one-electron Hamiltonian energy (Kinetic and Nuclear
C repulsion) at a given point. First let's get the correlated density 
C matrix.
                  CALL SAXPY(NBAS*NBAS, -1.0D0, SCFDENT, 1, RELDENT, 1)
C 
C Now we can very easily calculate the Kinetic energy and nuclear repulsion
C energy at a given point by simply taking dot products with KINTINT and
C PRDUTINT with the correlated density matrix RELDENT. 
C 
                  CALL DSCAL(NBAS*NBAS, HALF, KINTINT, 1)
                  CALL DSCAL(NBAS*NBAS, ATOMCHRG(IATMOS)/R, PRDUTINT, 1)
                  KINEP_R = SDOT(NBAS*NBAS, RELDENT, 1, KINTINT, 1)
                  REPLP_R = SDOT(NBAS*NBAS, RELDENT, 1, PRDUTINT, 1)
                  ONEHM_R = (KINEP_R - REPLP_R)
C 
C Let's turn the density back on where it should be
C
                  CALL SAXPY(NBAS*NBAS, 1.0D0, SCFDENT, 1, RELDENT, 1)
C
#ifdef _DEBUG_LVLM
                  Write(6,*) "The SCF Dens @twoden = ", SCF_DENSITY
#endif 
                 IF (TOT_DENSITY .LT. CUTOFF) GOTO 10
C
C Evaluate the correlation potential for this grid point on this 
C particular symm. unique center. We need to make sure that this 
C division is sane. 
C     
C$$$                  IF (KINEP_R .GE. 1.0D-4) THEN
C$$$                     KINE_POTNL_R =  KINEP_R/TOT_DENSITY
C$$$                  ELSE
C$$$                     KINE_POTNL_R = 0.0D0
C$$$                  ENDIF
                  KINE_POTNL_R =  KINEP_R/TOT_DENSITY
                  REPL_POTNL_R = -REPLP_R/TOT_DENSITY
                  ONEH_POTNL_R =  ONEHM_R/TOT_DENSITY
                  VCON_POTNL_R =  POTNL_R/(TWO*TOT_DENSITY)
                  VTOT_POTNL_R =  (VCON_POTNL_R +
     &                            ONEH_POTNL_R)
C
C$$$                  Write(6,*) "The Correlation Potential and Density= ", 
C$$$     &            ONEHM_R, TOT_DENSITY
C
C Let's write out a formatted file of correlation potential for 
C each point. I am not sure yet how to get a correlation potential
C for a molecule. But, I know how to do it for each symmetry 
C unique atmos (There is however no problems for atom)!!! 
C
                  COR_ENERGY_R =  TOTWT(IRADPT,IANGPT)*POTNL_R
C     
                  WRITE(IUNIT1, 9999) R, VCON_POTNL_R
                  WRITE(IUNIT4, 9999) R, ONEH_POTNL_R 
                  WRITE(IUNIT5, 9999) R, KINE_POTNL_R 
                  WRITE(IUNIT6, 9999) R, REPL_POTNL_R 
                  WRITE(IUNIT7, 9999) R, VTOT_POTNL_R 
                  WRITE(IUNIT2, 9999) R, FOUR*PI*(R**2)*TOT_DENSITY
C     
               END IF
C
               TWO_ENERGY  = TWO_ENERGY + TOTWT(IRADPT,IANGPT)*POTNL_R
               ONH_ENERGY  = ONH_ENERGY + TOTWT(IRADPT,IANGPT)*ONEHM_R 
C
C Loop over atoms end here!
C
CSSS  ENDDO
C
C The loop over angular points end here!
C
            ENDDO
C
C The loop over radial points end here!
C
         ENDDO
C
 10   CONTINUE
C     
      ENDIF
C
      IF (ZMP_POT) THEN

         IATMOS = 1
C
         DO IRADPT = 1, NRADPT
           
            DO IANGPT = 1, NAPTS 
C     
               XPOINT = XGRDPT(NCNTR, IRADPT, IANGPT)  + 
     &                  COORD((IATMOS-1)*3 + 1)
               YPOINT = YGRDPT(NCNTR, IRADPT, IANGPT)  + 
     &                  COORD((IATMOS-1)*3 + 2)
               ZPOINT = ZGRDPT(NCNTR, IRADPT, IANGPT)  +
     &                  COORD((IATMOS-1)*3 + 3)
C
               R = DSQRT(XPOINT**2 + YPOINT**2 + ZPOINT**2)
               
               CALL REPPRDINT(IUATMS, NATOMS, NUNQSHL, LNP1, IPOPF,
     &                        NSHL, NFCT, NBAS, NOFFSETATMP,
     &                        NOFFSETATMC, NOFFSETPRM, NOFFSETCON, 
     &                        NANGMOMSHL, NPRIMFUNSHL, NCONFUNSHL, 
     &                        XPOINT, YPOINT, ZPOINT, ALPHA, PCOEF, 
     &                        COORD, REPLINT, PRDUTINT, KINTINT, 
     &                        TMP1, TMP2, TMP3, OVERLAP_POTN)
C
               CALL  BUILT_ZMP_MATRX(PRDUTINT, WORK, ILEFT, NBAS)
               
               ZMP_POTNTIAL = SDOT(NBAS*NBAS, PRDUTINT, 1, REPLINT, 1)

            ENDDO
         ENDDO
C
      ELSE IF (COR_EXCH_HOLE .AND. ANALYTIC_INT) THEN
         IATMOS = 1
#ifdef _DEBUG_LVLM1
          Write(6,*)
          Write(6,*) "Grid points in twoden"
#endif
 
         DO IRADPT = 1, NRADPT
            DO IANGPT = 1, NAPTS
C
#ifdef _DEBUG_LVL0
              If(IRADPT .le. 10 .and. IANGPT .le. 10)
     &           Write(6,"(3F10.5)")XGRDPT(NCNTR, IRADPT, IANGPT),
     &                              YGRDPT(NCNTR, IRADPT, IANGPT),
     &                              ZGRDPT(NCNTR, IRADPT, IANGPT)
#endif
C
               IF (IANGPT .EQ. 1) THEN
                  XPOINT = XGRDPT(NCNTR, IRADPT, IANGPT)  + 
     &                     COORD((IATMOS-1)*3 + 1)
                  YPOINT = YGRDPT(NCNTR, IRADPT, IANGPT)  + 
     &                     COORD((IATMOS-1)*3 + 2)
                  ZPOINT = ZGRDPT(NCNTR, IRADPT, IANGPT)  +
     &                     COORD((IATMOS-1)*3 + 3)

                  UVALUE(IRADPT) = DSQRT(XPOINT**2 + YPOINT**2 + 
     &                                   ZPOINT**2)                  
C
CSSS                  CALL BLTCOR_HOLE(IUATMS, NATOMS, NUNQSHL, IUHF,
CSSS     &                             NTOTSHL, NBAS, NANGMOMTSHL, 
CSSS     &                             NCONFUNTSHL, ISHL2CNTR_MAP, 
CSSS     &                             NOFFSETSHL, ILEFT, ALPHA, PCOEF, 
CSSS     &                             COORD, NANGMOMSHL, NPRMFUNTSHL, 
CSSS     &                             NCONFUNSHL, NTOTPRIM, 
CSSS     &                             UVALUE, SCFDENT, RELDENT, SCFDEND, 
CSSS     &                             RELDEND, WORK, POTNL_R, 
CSSS     &                             REPLSN_FACT, EXCHNG_FACT, 
CSSS     &                             CORLN_FACT_O, CORLN_FACT_T)
               ENDIF
C
            ENDDO
         ENDDO
 
      END IF
C
      IF (DO_PROPS) THEN
C
        CALL ZERO(CC_TENSOR, MAXCENT*MAXCEN)
C
        CALL GETREC(0, 'JOBARC', 'NTOTPERT', IRECLEN, NTPERT)
        IF (IRECLEN .GT. 0) THEN
           CALL GETREC(20, 'JOBARC', 'NTOTPERT', 1, NTPERT)
        ELSE
           NTPERT = 0
           CALL IZERO(IRREP_PERT, 8)
           CALL GETREC (20, 'JOBARC', 'NPERTIRR', 8, IRREP_PERT)
           DO IREPS = 1, 8
              NTPERT = NTPERT + IRREP_PERT(IREPS)
           END DO
        ENDIF
        CALL A2GET_DEN(WORK, ILEFT, NTPERT, NBAS, 
     &                 IBEGIN_P_DENS, IBEGIN_P_OPRT, IBEGIN_AO_OVRLP,
     &                 IBEGIN_MO_OVRLP, IBEGIN_MO_VECTS,
     &                 IBEGIN_MEM_LEFT)
C
C There is no integration for FC couplings.
C
        IPOINT = 0
        IJPAIR = 0

        CALL IZERO(NPOINT_4PAIR, NATOMS*(NATOMS-1)/2)
C
        DO IATMOS = 1, NATOMS

           DO JATMOS = 1, NATOMS
C
#ifdef _DEBUG_LVLM1
              If(IRADPT .le. 10 .and. IANGPT .le. 10)
     &           Write(6,"(3F10.5)")(XGRDPT(NCNTR, IRADPT, IANGPT) +
     &                              COORD((IATMOS-1)*3 + 1)),
     &                              (YGRDPT(NCNTR, IRADPT, IANGPT) +
     &                              COORD((IATMOS-1)*3 + 2)),
     &                              (ZGRDPT(NCNTR, IRADPT, IANGPT)  +
                                    COORD((IATMOS-1)*3 + 3))
#endif                    
              IF (IATMOS .LT. JATMOS) THEN
c
                  IJPAIR   = IJPAIR + 1

                  X_IPOINT = COORD((IATMOS-1)*3 + 1)
                  Y_IPOINT = COORD((IATMOS-1)*3 + 2)
                  Z_IPOINT = COORD((IATMOS-1)*3 + 3)
C
                  X_JPOINT = COORD((JATMOS-1)*3 + 1)
                  Y_JPOINT = COORD((JATMOS-1)*3 + 2)
                  Z_JPOINT = COORD((JATMOS-1)*3 + 3)
C
C Vector along the IATOM and JATOM
C            
                  X_IJPOINT = X_JPOINT - X_IPOINT
                  Y_IJPOINT = Y_JPOINT - Y_IPOINT
                  Z_IJPOINT = Z_JPOINT - Z_IPOINT
C
                  DSCALE = 0.0D0
                  HALT   = .FALSE.
#ifdef _DEBUG_LVL0
                Write(6,*)
                Write(6, "(a,2I5)") "The Atom Pair", IATMOS, JATMOS
                Write(6,*)
#endif
                  DO While (.NOT. HALT)   
   
                     X_INCREMENT = (X_IJPOINT)*DSCALE
                     Y_INCREMENT = (Y_IJPOINT)*DSCALE
                     Z_INCREMENT = (Z_IJPOINT)*DSCALE

                     XDIFF = DABS(X_INCREMENT - X_IJPOINT)
                     YDIFF = DABS(Y_INCREMENT - Y_IJPOINT)
                     ZDIFF = DABS(Z_INCREMENT - Z_IJPOINT)
                     HALT  = (XDIFF.LT.TOL.AND.YDIFF.LT.TOL.AND.
     &                        ZDIFf .LT. TOL) 
 
                     XPOINT = X_IPOINT + X_INCREMENT
                     YPOINT = Y_IPOINT + Y_INCREMENT
                     ZPOINT = Z_IPOINT + Z_INCREMENT
C
#ifdef _DEBUG_LVL0
             Write(6,"(a,3F10.5)") "Perturbation is centered on ",
     &                                XPOINT, YPOINT, ZPOINT 
             Write(6,"(6F10.5)") X_INCREMENT, Y_INCREMENT, 
     &                             Z_INCREMENT,X_IJPOINT,
     &                             Y_IJPOINT,Z_IJPOINT
             Write(6, "(i2)") IPOINT
#endif
                     R = DSQRT(XPOINT**2 + YPOINT**2 + ZPOINT**2)
 
                     CALL REPPRDINT(IUATMS, NATOMS, NUNQSHL, LNP1,
     &                              IPOPF, NSHL, NFCT, NBAS, 
     &                              NOFFSETATMP, NOFFSETATMC, 
     &                              NOFFSETPRM, NOFFSETCON,
     &                              NANGMOMSHL, NPRIMFUNSHL, 
     &                              NCONFUNSHL, XPOINT, YPOINT, 
     &                              ZPOINT, ALPHA, PCOEF, COORD, 
     &                              REPLINT, PRDUTINT, KINTINT,
     &                              TMP1, TMP2, TMP3, OVERLAP_POTN)
C
                    CALL A2BUILD_TENSOR(WORK, MAXCOR, PRDUTINT, 
     &                                  CC_TENSOR_ATR, NTPERT, NBAS, 
     &                                  NATOMS, IBEGIN_P_DENS, 
     &                                  IBEGIN_P_OPRT, 
     &                                  IBEGIN_AO_OVRLP, 
     &                                  IBEGIN_MO_OVRLP,
     &                                  IBEGIN_MO_VECTS,
     &                                  IBEGIN_MEM_LEFT, IATMOS)
C
                    NPOINT_PAIR(IJPAIR) = NPOINT_PAIR(IJPAIR)  + 1
                                  IPOINT = IPOINT + 1
C
                    GRID_POINTS(1,IPOINT) = XPOINT - X_IPOINT
                    GRID_POINTS(2,IPOINT) = YPOINT - Y_IPOINT
                    GRID_POINTS(3,IPOINT) = ZPOINT - Z_IPOINT
C
                    CALL DCOPY(NTPERT, CC_TENSOR_ATR, 1, 
     &                         CC_ARRAY(1, IPOINT), 1)
                    DSCALE = DSCALE + 0.01D0
C
CSSS                    CALL SYMMET2(CC_ARRAY(1, IPOINT), NTPERT)
CSSS                    CALL DSCAL(NTPERT*NTPERT, 2.0D0, 
CSSS     &                         CC_ARRAY(1,IPOINT), 1)
CSSS
CSSS                    CALL PRNT_SPNSPN_TNSR(CC_ARRAY(1, IPOINT), 
CSSS     &                                    NTPERT, .TRUE., .FALSE.,
CSSS     &                                    .FALSE., .FALSE., 1)
                  ENDDO
C
              ENDIF
C
           ENDDO
        ENDDO
C
      CALL PROCESS_TENSOR(CC_ARRAY, CC_TENSOR,GRID_POINTS, NTPERT, 
     &                    IPOINT, IJPAIR, MAXCENT, MAXGRD, MAXPAIRS,
     &                    NATOMS, NPOINT_PAIR)
#ifdef _DEBUG_LVLM1
            Write(6,*)
            Write(6, "(a)") "@-INTGRT-TWODEN, CC Tensor"
            Index = 0
            DO I =1, IJPAIR
               Write(6,*) "New Pair"
C               Index = (I-1)*NPOINT_PAIR(I)
               DO K = 1, NPOINT_PAIR(I)
                  Index = INdex + 1
                  Write(6,*) Index
                  CALL OUTPUT(CC_ARRAY(1,Index), 1, NTPERT, 1, 1,
     &                        NTPERT, 1, 1)
               ENDDO
C               Index = 0
            ENDDO
#endif
C
      ENDIF 
C
      IF (COR_EXCH_HOLE) THEN 

         EXC_COR_ENERGY = 0.0D0
         EXC_COR_LIMIT  = 0.0D0
C
         DO JRADPT = 1, NRADPT
C     
            IF (UVALUE(JRADPT) .NE. ZLICH) THEN
               IF (IFLAGS2(153) .EQ. 0) THEN
                  WRITE(IUNIT1, 9999) UVALUE(JRADPT),
     &                                EXC_COR_HOLE(JRADPT)*
     &                                UVALUE(JRADPT)*HALF
               ELSE
                  WRITE(IUNIT1, 9999) UVALUE(JRADPT),
     &                                DFT_EXC_COR_HOLE(JRADPT)
               ENDIF
C
               EXC_COR_ENERGY = EXC_COR_ENERGY + EXC_COR_HOLE(JRADPT)*
     &                          RADWT(JRADPT)/UVALUE(JRADPT)
               EXC_COR_LIMIT  = EXC_COR_LIMIT + EXC_COR_HOLE(JRADPT)*
     &                          RADWT(JRADPT)
            ENDIF
        ENDDO
C
C Write the appropriate limits as an external check. In the case of exchange
C operator the limits are exchange energy and -1, and in the case of correlation
C energy the lmits are correlation energy and 0. 
C
        IF (COR_EXCH_HOLE) THEN

           WRITE(6,*) 
           WRITE(LUOUT,9)   ONEP_ENERGY
           WRITE(LUOUT,89)  EXC_COR_LIMIT*HALF
           WRITE(LUOUT,999) EXC_COR_ENERGY*HALF
        ENDIF
C
      ELSE IF (.NOT. DO_PROPS) THEN
C
C Write the appropriate contribution (depending on the what you choose
C in pick_potential the quantity printied here can be HF replusion, 
C HF exchange and two particle correlation contributions to the to
C total energy of the system.
C
         WRITE(6,*) 
         WRITE(LUOUT,9)    ONEP_ENERGY
         WRITE(LUOUT,99)   ONH_ENERGY
         WRITE(LUOUT,999)  TWO_ENERGY*FACT*HALF
         WRITE(LUOUT,888)  ONH_ENERGY + TWO_ENERGY*FACT*HALF

      ENDIF
C
 9    FORMAT('The One Particle Cont. Correlation Energy =',F14.9 )
 99   FORMAT('The One Hamilto. Cont. Correlation Energy =',F14.9 )
 89   FORMAT('The Repulsion/Exchange/Correlation Limit  =',F14.9)
 999  FORMAT('The Repulsion/Exchange/Correlation Energy =',F14.9)
 888  FORMAT('The total Correlation Energy              =',F14.9)
 9999 FORMAT(1X, F14.9, 3X, F14.9)
C
      IF (IFLAG2. EQ. 1) THEN
         CLOSE(IUNIT1)
         CLOSE(IUNIT2)
         CLOSE(IUNIT4)
         CLOSE(IUNIT5)
         CLOSE(IUNIT6)
         CLOSE(IUNIT7)
      ENDIF
C     
      RETURN
      END
 
