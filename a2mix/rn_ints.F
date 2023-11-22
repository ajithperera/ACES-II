      SUBROUTINE RN_INTS(NCNTR,NATOMS,NFCT,NTANGM,IANGX,IANGY,
     $   IANGZ,TOTWT,XCDNT,YCDNT,ZCDNT,RSQRD,COORD,XGRDPT,
     $   YGRDPT,ZGRDPT,ALPHA,NANGMOM,NMOMFCT,
     $   FCTDEN,FCTGRX,FCTGRY,FCTGRZ,IPOPF,IUATMS,IFACT,IUHF,NRADPT,
     $   XALPHA,BECKE,dlda,
     &   ELYP,CUTOFF,NAPTS,NCNTO,NMOA,NMOB,NBASP,
     $   PCOEFFA,PCOEFFB,VALMO,XGVAL,YGVAL,ZGVAL,ANGXVL,ANGYVL,
     $   ANGZVL,XOCC,XNAT,MAXANG,RRR,PGRP,NELEC_COMP,RMOM_T,RMOM_D)
C     
C     This routine does a numerical integration over spherical 
C     coordinates.  The radial quadrature formula is taken from P. 
C     M. W. Gill, B. G. Johnson, J. A. Pople and M. J. Frisch, 
C     Chem. Phys. Lett. 197, 499 (1992).
C     The angular quadrature formula is taken from V. I. Lebedev, 
C     Zh. Vychisl. Mat. Fiz. 15, 48 (1975) and 16, 293 (1976), 
C     (english translation in U.S.S.R. Comput. Math and Math Phys).
C     
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     
      CHARACTER*8 PGRP
      LOGICAL SPIN_DENS 
      CHARACTER*5 FILENAME1,FILENAME2,FILENAME3,FILENAME4,
     &            FILENAME5,FILENAME6,FILENAME7
C
      COMMON /PAR/ PI,PIX4
      COMMON /IPAR/ LUOUT
      COMMON /CNST/ ALPH,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11,X12,X13,X14,X15
      COMMON /FRAC/ F1,F2,F3,F4,F5,F6,F7,F8,F9,F10,F11,F12,F13,F14,F15
C     
      DIMENSION IANGX(56),IANGY(56),IANGZ(56),
     &   TOTWT(NRADPT,194),IPOPF(IUATMS),
     &   XGRDPT(NATOMS,NRADPT,194),YGRDPT(NATOMS,NRADPT,194),
     &   ZGRDPT(NATOMS,NRADPT,194),
     &   XCDNT(NATOMS,NRADPT,194),YCDNT(NATOMS,NRADPT,194),
     &   ZCDNT(NATOMS,NRADPT,194),RSQRD(NATOMS,NRADPT,194),
     &   ALPHA(NFCT),NANGMOM(NATOMS),NMOMFCT(NATOMS*NTANGM),
     &   FCTDEN(NFCT),FCTGRX(NFCT),FCTGRY(NFCT),FCTGRZ(NFCT),
     $   PCOEFFA(NBASP*NFCT),PCOEFFB(NBASP*NFCT),VALMO(NBASP),
     $   XGVAL(NMOA),YGVAL(NMOA),ZGVAL(NMOA),ANGXVL(8),
     $   ANGYVL(8),ANGZVL(8),RRR(NATOMS,NRADPT,194),
     $   RMOM_T(7),RMOM_D(7),XOCC(NBASP,2),XNAT(NBASP),
     $   COORD(NATOMS*3)

      DATA THRES/1.0D-3/
C
C     EULER-MACLAURIN SCHEME, RADIAL INTEGRATION
C     GAUSS-MARKOV QUADRATURE, ANGULAR INTEGRATION
C     
C     Initialize integral values
C
      TOTELE=0.0D+00
      CUTOFF=1.0D-10
      PI    =DATAN(1.0D+00)*4.D+00

#ifdef _DEBUG_LVLM1
        Write(6,*)
        Write(6,"(a)") "The weight array in one-den"
C        Write(6,"(4F15.5)") ((TOTWT(j, i),j=1,NRADPT), i=1,194)
#endif
      IF(PGRP.NE.'C1')THEN
        FACT=FACT*REAL(IPOPF(NCNTO))
      ENDIF
      CALL DZERO(RMOM_T,7)
      CALL DZERO(RMOM_D,7)
      SPIN_DENS = .FALSE.
      IOFF = 1

C Compute <r^-3>,<r^-2>,<r^-1>,<r^0>,<r^1>,<r^2>,<r^3>

      FILENAME1 = "R3INV"
      FILENAME2 = "R2INV"
      FILENAME3 = "R1INV"
      FILENAME4 = "R0PWR"
      FILENAME5 = "R1PWR"
      FILENAME6 = "R2PWR"
      FILENAME7 = "R3PWR"

      IUNIT1 = 10
      IUNIT2 = 11
      IUNIT3 = 12
      IUNIT4 = 13
      IUNIT5 = 14
      IUNIT6 = 15
      IUNIT7 = 16

      OPEN (UNIT=IUNIT1, FILE=FILENAME1, FORM="FORMATTED")
      OPEN (UNIT=IUNIT2, FILE=FILENAME2, FORM="FORMATTED")
      OPEN (UNIT=IUNIT3, FILE=FILENAME3, FORM="FORMATTED")
      OPEN (UNIT=IUNIT4, FILE=FILENAME4, FORM="FORMATTED")
      OPEN (UNIT=IUNIT5, FILE=FILENAME5, FORM="FORMATTED")
      OPEN (UNIT=IUNIT6, FILE=FILENAME6, FORM="FORMATTED")
      OPEN (UNIT=IUNIT7, FILE=FILENAME7, FORM="FORMATTED")

      DO IMOM = -3, 3, 1
C     Loop over angular points
#ifdef _DEBUG_LVL0
         IF (IMOM .EQ.0) SUM=0.0D0
#endif 
C
      DO 10 IRADPT=1,NRADPT
            RMOM_P = 0.0D0 
C     Loop over radial points
         DO 20 IANGPT=1,NAPTS

C         IF (IANGPT .EQ. 1) THEN
C
C     Set function counter
            IFCT=1

C     Loop over atoms
            DO 30 IATM=1,NATOMS
C
               XPOINT = XGRDPT(NCNTR, IRADPT, IANGPT)  +
     &                  COORD((IATM-1)*3 + 1)
               YPOINT = YGRDPT(NCNTR, IRADPT, IANGPT)  +
     &                  COORD((IATM-1)*3 + 2)
               ZPOINT = ZGRDPT(NCNTR, IRADPT, IANGPT)  +
     &                  COORD((IATM-1)*3 + 3)

               R2 = XPOINT**2 + YPOINT**2 + ZPOINT**2
               R  = DSQRT(R2)

C     ANGULAR MOMENTUM
               ANGXVL(1)=1.D+00
               ANGYVL(1)=1.D+00
               ANGZVL(1)=1.D+00
               ANGXVL(2)=1.D+00
               ANGYVL(2)=1.D+00
               ANGZVL(2)=1.D+00
               ANGXVL(3)=XCDNT(IATM,IRADPT,IANGPT)
               ANGYVL(3)=YCDNT(IATM,IRADPT,IANGPT)
               ANGZVL(3)=ZCDNT(IATM,IRADPT,IANGPT)
#ifdef _DEBUG_LVLM
              If(IRADPT .le. 3 .and. IANGPT .le. 3) then
c             Write(6,"(3(1x,F15.7))")  XCDNT(IATM, IRADPT, IANGPT),
c     &                            YCDNT(IATM, IRADPT, IANGPT),
c     &                            ZCDNT(IATM, IRADPT, IANGPT)
	     write(6,"(2(1x,F15.7))") RSQRD(IATM,IRADPT,IANGPT),R2
CSSS             write(6,"(3(1x,F15.7))")XPOINT,YPOINT,ZPOINT
             endif
#endif 
C
               DO 35 IANG=3,MAXANG+2
                  ANGXVL(IANG+1)=ANGXVL(IANG)*XCDNT(IATM,IRADPT,IANGPT)
                  ANGYVL(IANG+1)=ANGYVL(IANG)*YCDNT(IATM,IRADPT,IANGPT)
                  ANGZVL(IANG+1)=ANGZVL(IANG)*ZCDNT(IATM,IRADPT,IANGPT)
  35           CONTINUE
C
C     Loop over angular momentum for each atom 
               DO 40 IANGMOM=1,NANGMOM(IATM) 
C     Angular intermediates for the density and gradient
C     
C     Density
                  ANGDEN=ANGXVL(IANGX(IANGMOM)+2)
     $                   *ANGYVL(IANGY(IANGMOM)+2)
     $                   *ANGZVL(IANGZ(IANGMOM)+2)

C     Loop over number of functions for each angular momentum
                  DO 50 IMOMFCT=1,NMOMFCT((IATM-1)*NTANGM+IANGMOM)
                     TEXP=DEXP(-ALPHA(IFCT)*R2)
                     TWOALP=2.D+00*ALPHA(IFCT)
                     RM =R**(IMOM)
                     FCTDEN(IFCT)=ANGDEN*TEXP*RM
                     IFCT=IFCT+1
 50               CONTINUE
 40            CONTINUE
 30         CONTINUE
C
C     Apha spin
C
C     Contract functions with alpha coefficients and add to get 
C     alpha MO values at this point.
            CALL XGEMM('T','N',NBASP,1,NFCT,1.0D+00,PCOEFFA,NFCT,
     $                  FCTDEN,NFCT,0.0D+00,VALMO,NBASP)
            CALL VECPRD(VALMO,VALMO,XNAT,NBASP) 
            ROA = DDOT(NBASP,XNAT,1,XOCC(1,1),1)

            IF(IUHF.EQ.1)THEN
C
C     Beta Spin
C
C     Contract functions with beta coefficients and add to get beta
C     MO values at this point.
               CALL XGEMM('T','N',NBASP,1,NFCT,1.0D+00,PCOEFFB,NFCT,
     $                     FCTDEN,NFCT,0.0D+00,VALMO,NBASP)
               CALL VECPRD(VALMO,VALMO,XNAT,NBASP) 
               ROB = DDOT(NBASP,XNAT,1,XOCC(1,2),1)
            ELSE
               ROB=ROA
            ENDIF
C
C      Determine total and spin-density properties. 
C
            RMOM_T(IOFF) = RMOM_T(IOFF) + TOTWT(IRADPT,IANGPT)*
     $                                    (ROA+ROB)
            IF ((ROA-ROB) .GT. 0.0D0) THEN
                RMOM_D(IOFF) = RMOM_D(IOFF) +
     $                         TOTWT(IRADPT,IANGPT)*(ROA-ROB)
                RMOM_P = RMOM_P +
     $                         TOTWT(IRADPT,IANGPT)*(ROA-ROB)
            ENDIF
C
C     remaining radial points and go to the next angular point
C
            IF(ROA+ROB.LT.CUTOFF) GOTO 10

C      ENDIF 
C        Next angular point
 20      CONTINUE

         IF (IMOM .EQ. -3) WRITE(IUNIT1,9999) R, RMOM_P
         IF (IMOM .EQ. -2) WRITE(IUNIT2,9999) R, RMOM_P
         IF (IMOM .EQ. -1) WRITE(IUNIT3,9999) R, RMOM_P
         IF (IMOM .EQ.  0) WRITE(IUNIT4,9999) R, RMOM_P
         IF (IMOM .EQ.  1) WRITE(IUNIT5,9999) R, RMOM_P
         IF (IMOM .EQ.  2) WRITE(IUNIT6,9999) R, RMOM_P
         IF (IMOM .EQ.  3) WRITE(IUNIT7,9999) R, RMOM_P
#ifdef _DEBUG_LVL0
         IF (IMOM .EQ. 0) Sum = Sum + RMOM_P
#endif
C     Next radial point
 10   CONTINUE
      IOFF = IOFF + 1
      ENDDO 
 9999 FORMAT(1X,F14.9,3X,F14.9)

#ifdef _DEBUG_LVL0
        Write(6,"(2a,F10.5)") "The total number of unpaired",
     &                        " electrons: ", Sum
#endif
      NELEC_COM = NELEC_COMP + RMOM_T(4)
#ifdef _NOSKIP
      WRITE(6,*)
      WRITE(6,"(9X,2A)") '------------------------------',
     &                       '--------------------------------'
      WRITE(6,"(9X,2A)") " The expectation value of moments",
     &                   " operators from total density"
      WRITE(6,"(2X,2A)") '------------------------------',
     &                 '-----------------------------------------------'
      WRITE(6,"(2x,2A)") "    <r^-3>     <r^-2>     <r^-1>     <r^0>",
     &                   "      <r^1>      <r^2>      <r^3> "

      WRITE(6,"(2x,7(3X,ES8.2E2))") (RMOM_T(I)*FACT, I=1,7)

      WRITE(6,"(2X,2A)") '------------------------------',
     &                 '-----------------------------------------------'

      IF (RMOM_D(4) .GT. THRES) SPIN_DENS=.TRUE.
      WRITE(6,*) 
      IF (SPIN_DENS) THEN
      WRITE(6,"(9X,2A)") '------------------------------',
     &                       '--------------------------------'
      WRITE(6,"(9X,2A)") " The expectation value of moments",
     &                   " operators from spin density"
      WRITE(6,"(2X,2A)") '------------------------------',
     &                 '-----------------------------------------------'
      WRITE(6,"(2x,2A)") "    <r^-3>     <r^-2>     <r^-1>     <r^0>",
     &                   "      <r^1>      <r^2>      <r^3> "

      WRITE(6,"(2x,7(3X,ES8.2E2))") (RMOM_D(I)*FACT, I=1,7)

      WRITE(6,"(2X,2A)") '------------------------------',
     &                 '-----------------------------------------------'
      ENDIF 
#endif
C     
      RETURN
      END
