      SUBROUTINE OCT(ITMP,XDAT,YDAT,ZDAT,XCDNT,YCDNT,ZCDNT,
     & ATMXVEC,ATMYVEC,ATMZVEC,RSQRD,RI,RIJ,NATOMS,AIJ,WTINTR,
     & TOTWT,WGHT,I1,I2,IOFF1,IOFF2,NITR,NCNTR,RAD,NRADPT,NUC,
     & RRR, XGRDPT,YGRDPT,ZGRDPT,ANGWT,RADWT_OCT)
C
C This routine sets up the numerical grid for the xyz octant
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      COMMON /PAR/ PI,PIX4
      COMMON /IPAR/ LUOUT
C
      DIMENSION WGHT(98),XDAT(194),YDAT(194),ZDAT(194),
     & ATMXVEC(NATOMS,NATOMS),ATMYVEC(NATOMS,NATOMS),
     & ATMZVEC(NATOMS,NATOMS),RI(NATOMS),RIJ(NATOMS,NATOMS),
     & WTINTR(NATOMS),TOTWT(NRADPT,194),
     & XCDNT(NATOMS,NRADPT,194),YCDNT(NATOMS,NRADPT,194),
     & ZCDNT(NATOMS,NRADPT,194),RSQRD(NATOMS,NRADPT,194),
     & AIJ(NATOMS,NATOMS),NUC(NATOMS),RRR(NATOMS,NRADPT,194),
     & XGRDPT(NATOMS,NRADPT,194),YGRDPT(NATOMS,NRADPT,194),
     & ZGRDPT(NATOMS,NRADPT,194),ANGWT(194),RADWT_OCT(50)

C Integrate over one octant of angular points
C
C     Number of radial points
      XTMP=REAL(ITMP)
      XNP1=XTMP+1.D+00
         DO 600 IRAD=1,ITMP
C     Determine the radial points and weights
         XI=REAL(IRAD)
C        Determine radial point
         R1=(RAD*XI*XI)/((XNP1-XI)*(XNP1-XI))
C        Calculate weight for radial quadrature
         DEN=(XNP1-XI)*(XNP1-XI)*(XNP1-XI)*(XNP1-XI)
     &                *(XNP1-XI)*(XNP1-XI)*(XNP1-XI)
         RADWT=(2.0*RAD*RAD*RAD*XNP1*XI*XI*XI*XI*XI)/DEN
C
         DO 605 I=I1,I2
            XD=R1*XDAT(I+IOFF1)
            YD=R1*YDAT(I+IOFF1)
            ZD=R1*ZDAT(I+IOFF1)
C
C The correlation potential work need just the grid points for 
C each symmetry unique center. 01/2001, Ajith Perera                   
C
            XGRDPT(NCNTR, IRAD, I) = XD
            YGRDPT(NCNTR, IRAD, I) = YD
            ZGRDPT(NCNTR, IRAD, I) = ZD

         DO 610 IATM=1,NATOMS
            XCDNT(IATM,IRAD,I)=ATMXVEC(NCNTR,IATM)+XD
            YCDNT(IATM,IRAD,I)=ATMYVEC(NCNTR,IATM)+YD
            ZCDNT(IATM,IRAD,I)=ATMZVEC(NCNTR,IATM)+ZD
           RSQRD(IATM,IRAD,I)=XCDNT(IATM,IRAD,I)**2.D+00+
     &        YCDNT(IATM,IRAD,I)**2.D+00+ZCDNT(IATM,IRAD,I)**2.D+00
            RI(IATM)=DSQRT(RSQRD(IATM,IRAD,I))
            RRR(IATM,IRAD,I)=RI(IATM)
  610    CONTINUE
         WTTOT=0.D+00
         DO 611 IATM=1,NATOMS
            IF(NUC(IATM).EQ.110)THEN
               WTINTR(IATM)=0.D+00
               GOTO 611
            ELSE
               WTINTR(IATM)=1.D+00
            ENDIF
            DO 612 JATM=1,IATM-1
              IF(NUC(JATM).EQ.110)GOTO 612
              ZMUIJ=(RI(IATM)-RI(JATM))/RIJ(IATM,JATM)
              XMUIJ=ZMUIJ
     &           +AIJ(IATM,JATM)*(1.D+00-AIJ(IATM,JATM)*AIJ(IATM,JATM))
              F4=FUZZY(NITR,XMUIJ)
              CUTIJ=0.5D+00*(1.D+00-F4)
              WTINTR(IATM)=WTINTR(IATM)*CUTIJ
  612       CONTINUE
            DO 613 JATM=IATM+1,NATOMS
              IF(NUC(JATM).EQ.110)GOTO 613
              ZMUIJ=(RI(IATM)-RI(JATM))/RIJ(IATM,JATM)
              XMUIJ=ZMUIJ
     &           +AIJ(IATM,JATM)*(1.D+00-AIJ(IATM,JATM)*AIJ(IATM,JATM))
              F4=FUZZY(NITR,XMUIJ)
              CUTIJ=0.5D+00*(1.D+00-F4)
              WTINTR(IATM)=WTINTR(IATM)*CUTIJ
  613       CONTINUE
         WTTOT=WTTOT+WTINTR(IATM)
  611    CONTINUE
c         ATMWT=1
         ATMWT=WTINTR(NCNTR)/WTTOT
         TOTWT(IRAD,I)=PIX4*ATMWT*RADWT*WGHT(I-IOFF2)
         IF (IRAD .EQ. 1) ANGWT(I) = WGHT(I-IOFF2)*PIX4
         RADWT_OCT(IRAD) = RADWT
  605    CONTINUE
  600    CONTINUE
C
      RETURN
      END


