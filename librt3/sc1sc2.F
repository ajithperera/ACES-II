
      SUBROUTINE SC1SC2(T3AAA,T3AAB,LEN,IADT3,IRPIJK)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

cYAU - This is very worrisome. For one of our test cases, LEN is 619,
c      but further down, when we start looping over indices, we address
c      1158 elements. The input file is included for what it's worth:
c
c TEST.016 RHF-CCSD[T]/DZP GEOMETRY OPTIMIZATION H2O
c O
c H 1 R*
c H 1 R* 2 A*
c 
c R=0.943056
c A=105.96879
c 
c *ACES2(CALC=CCSD[T],BASIS=DZP,MEMORY=4000000)
c
c    Anyway, the GNU optimizer has trouble with some of the loops in
c this routine. A few were rewritten while the others were left alone.

      INTEGER LEN
      DOUBLE PRECISION T3AAA(LEN), T3AAB(LEN)
      INTEGER A,B,C,AB,AC,BC,ABC,ACB,BCA
      INTEGER POP,VRT,DIRPRD
      DIMENSION IADT3(8)
#include "maxbasfn.par"
      DIMENSION SIGN(MAXBASFN),IADAC(MAXBASFN),IADBC(MAXBASFN)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     1                DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYM/    POP(8,2),VRT(8,2),NTAA,NTBB,NF1AA,NF2AA,
     1                NF1BB,NF2BB
      COMMON /INFO/   NOCCO(2),NVRTO(2)
C
      COMMON /T3OFF/  IOFFVV(8,8,10),IOFFOO(8,8,10),IOFFVO(8,8,4)
C
      INDEX(I) = ISHFT((I*(I-1)),-1)
C
C     SUBROUTINE THAT IS SUPPOSED TO GENERATE ALL ALPHA TRIPLES FROM
C     CASE 2 (AAB) TRIPLES FOR RHF, ACCORDING TO
C
C      ABC       ABc        BCa       ACb
C     T      =  T     +    T     -   T
C      IJK       IJK        IJK       IJK
C
C     IN SYMMETRY
C
      DO 500 IRPC=1,NIRREP
      IF (VRT(IRPC,1).EQ.0) GOTO 500
      IRPAB = DIRPRD(IRPC,IRPIJK)
      IF (IRPAB.EQ.1) THEN

c ----------------------------------------------------------------------

      DO IRPB=1,NIRREP
      IF (VRT(IRPB,1).GE.2) THEN

c *****************************
c *** BEGIN (IRPAB=1) CASES ***
c *****************************

      MAXC = VRT(IRPC,1)
      MAXB = VRT(IRPB,1)
c     MAXA = nothing, since the upper bound changes per situation

      IRPA  = IRPB
      IRPAC = DIRPRD(IRPA,IRPC)
      IRPBC = DIRPRD(IRPB,IRPC)

c   o Case 1
      IF (IRPB.EQ.IRPC) THEN
         DO C=1,MAXC
            DO I=1,MAXC
               IF (I.LT.C) THEN
                  SIGN(I) = 1.0D0
                  IADAC(I) = IOFFVV(IRPC,IRPAC,1) + INDEX(C-1) + I
               END IF
               IF (I.EQ.C) SIGN(I) = 0.0D0
               IF (I.GT.C) THEN
                  SIGN(I) = -1.0D0
                  IADAC(I) = IOFFVV(IRPA,IRPAC,1) + INDEX(I-1) + C
               END IF
               IADBC(I) = IADAC(I)
            END DO
            DO B=2,MAXB
               MAXA=B-1
               DO A=1,MAXA
      ABC = IADT3(IRPC) + (C*IRPDPD(IRPAB,1)) - IRPDPD(IRPAB,1)
     &    + IOFFVV(IRPB,IRPAB,1) + INDEX(MAXA) + A
      BCA = IADT3(IRPA) + (A*IRPDPD(IRPBC,1)) - IRPDPD(IRPBC,1)
     &    + IADBC(B)
      ACB = IADT3(IRPB) + (B*IRPDPD(IRPAC,1)) - IRPDPD(IRPAC,1)
     &    + IADAC(A)
      IF ((B.NE.C).AND.(A.NE.C)) T3AAA(ABC-1) = T3AAA(ABC-1)
     &                           + SIGN(B) * T3AAB(BCA-1)
     &                           - SIGN(A) * T3AAB(ACB-1)
      IF ((B.EQ.C).AND.(A.NE.C)) T3AAA(ABC-1) = T3AAA(ABC-1)
     &                           - SIGN(A) * T3AAB(ACB-1)
      IF ((B.NE.C).AND.(A.EQ.C)) T3AAA(ABC-1) = T3AAA(ABC-1)
     &                           + SIGN(B) * T3AAB(BCA-1)
               END DO
            END DO
         END DO
c     END IF (IRPB.EQ.IRPC)
      END IF

c   o Case 2
      IF (IRPB.LT.IRPC) THEN
         SIGNBC = 1.0D0
         SIGNAC = 1.0D0
         DO C=1,MAXC
            DO B=2,MAXB
               MAXA=B-1
               DO A=1,MAXA
      ABC = IADT3(IRPC) + (C*IRPDPD(IRPAB,1)) - IRPDPD(IRPAB,1)
     &    + IOFFVV(IRPB,IRPAB,1) + INDEX(MAXA) + A
      BCA = IADT3(IRPA) + (A*IRPDPD(IRPBC,1)) - IRPDPD(IRPBC,1)
     &    + IOFFVV(IRPC,IRPBC,1) + C*MAXB - MAXB + B
      ACB = IADT3(IRPB) + (B*IRPDPD(IRPAC,1)) - IRPDPD(IRPAC,1)
     &    + IOFFVV(IRPC,IRPAC,1) + C*VRT(IRPA,1) - VRT(IRPA,1) + A
      T3AAA(ABC-1) = T3AAA(ABC-1) + SIGNBC * T3AAB(BCA-1)
     &                            - SIGNAC * T3AAB(ACB-1)
               END DO
            END DO
         END DO
c     END IF (IRPB.LT.IRPC)
      END IF

c   o Case 3
      IF (IRPB.GT.IRPC) THEN
         SIGNBC = -1.0D0
         SIGNAC = -1.0D0
         DO C=1,MAXC
            DO B=2,MAXB
               MAXA=B-1
               DO A=1,MAXA
      ABC = IADT3(IRPC) + (C*IRPDPD(IRPAB,1)) - IRPDPD(IRPAB,1)
     &    + IOFFVV(IRPB,IRPAB,1) + INDEX(MAXA) + A
      BCA = IADT3(IRPA) + (A*IRPDPD(IRPBC,1)) - IRPDPD(IRPBC,1)
     &    + IOFFVV(IRPB,IRPBC,1) + MAXA*MAXC + C
      ACB = IADT3(IRPB) + (B*IRPDPD(IRPAC,1)) - IRPDPD(IRPAC,1)
     &    + IOFFVV(IRPA,IRPAC,1) + A*MAXC - MAXC + C
      T3AAA(ABC-1) = T3AAA(ABC-1) + SIGNBC * T3AAB(BCA-1)
     &                            - SIGNAC * T3AAB(ACB-1)
               END DO
            END DO
         END DO
c     END IF (IRPB.GT.IRPC)
      END IF

c ***************************
c *** END (IRPAB=1) CASES ***
c ***************************

c     END IF (VRT(IRPB,1).GE.2)
      END IF
c     END DO IRPB=1,NIRREP
      END DO

c ----------------------------------------------------------------------

c     END IF (IRPAB.EQ.1)
      END IF

      IF(IRPAB.NE.1)THEN
C
      DO  400 IRPB=1,NIRREP
C
      IRPA  = DIRPRD(IRPB,IRPAB)
      IRPAC = DIRPRD(IRPA,IRPC)
      IRPBC = DIRPRD(IRPB,IRPC)
C
      IF(IRPA.GT.IRPB) GOTO 400
      IF(VRT(IRPB,1).EQ.0.OR.VRT(IRPA,1).EQ.0) GOTO 400
C
        IF(IRPA.EQ.IRPC)THEN
C
C       It follows that IRPB > IRPC
C
        SIGNBC = -1.0D+00
        DO  230 C=1,VRT(IRPC,1)
C
        DO  205    I=1,VRT(IRPC,1)
        IF(I.LT.C)THEN
        SIGN(I) =  1.0D+00
        IADAC(I) = IOFFVV(IRPC,IRPAC,1) + INDEX(C-1) + I
        ENDIF
        IF(I.EQ.C) SIGN(I) = 0.0D+00
        IF(I.GT.C)THEN
        SIGN(I) = -1.0D+00
        IADAC(I) = IOFFVV(IRPA,IRPAC,1) + INDEX(I-1) + C
        ENDIF
  205 CONTINUE
C
        DO  220 B=1,VRT(IRPB,1)
        DO  210 A=1,VRT(IRPA,1)
C
        ABC = IADT3(IRPC) + (C-1) * IRPDPD(IRPAB,1) + 
     1        IOFFVV(IRPB,IRPAB,1) + (B-1)*VRT(IRPA,1) + A
C
        BC = IOFFVV(IRPB,IRPBC,1) + (B-1)*VRT(IRPC,1) + C
C        IF(A.LT.C)THEN
C        SIGNAC =  1.0D+00
C        AC     =  IOFFVV(IRPC,IRPAC,1) + INDEX(C-1) + A
C        ELSE
C        SIGNAC = -1.0D+00
C        AC     =  IOFFVV(IRPA,IRPAC,1) + INDEX(A-1) + C
C        ENDIF
C
        BCA = IADT3(IRPA) + (A-1) * IRPDPD(IRPBC,1) + BC
        ACB = IADT3(IRPB) + (B-1) * IRPDPD(IRPAC,1) + IADAC(A)
C
        IF(A.NE.C)THEN
C        T3AAA(ABC-1) = T3AAA(ABC-1) + SIGNBC * T3AAB(BCA-1)
C     1                              - SIGNAC * T3AAB(ACB-1)
        T3AAA(ABC-1) = T3AAA(ABC-1) + SIGNBC * T3AAB(BCA-1)
     1                              - SIGN(A) * T3AAB(ACB-1)
        ELSE
        T3AAA(ABC-1) = T3AAA(ABC-1) + SIGNBC * T3AAB(BCA-1)
        ENDIF
  210   CONTINUE
  220   CONTINUE
  230   CONTINUE
        ENDIF
C
C     o It follows that IRPB > IRPC
        IF (IRPA.GT.IRPC) THEN
           SIGNAC = -1.0D0
           SIGNBC = -1.0D0
           MAXC=VRT(IRPC,1)
           MAXB=VRT(IRPB,1)
           MAXA=VRT(IRPA,1)
           DO C=1,MAXC
              DO B=1,MAXB
                 DO A=1,MAXA
      ABC = IADT3(IRPC) + (C*IRPDPD(IRPAB,1)) - IRPDPD(IRPAB,1)
     &    + IOFFVV(IRPB,IRPAB,1) + B*MAXA - MAXA + A
      BCA = IADT3(IRPA) + (A*IRPDPD(IRPBC,1)) - IRPDPD(IRPBC,1)
     &    + IOFFVV(IRPB,IRPBC,1) + B*MAXC - MAXC + C
      ACB = IADT3(IRPB) + (B*IRPDPD(IRPAC,1)) - IRPDPD(IRPAC,1)
     &    + IOFFVV(IRPA,IRPAC,1) + A*MAXC - MAXC + C
      T3AAA(ABC-1) = T3AAA(ABC-1) + SIGNBC * T3AAB(BCA-1)
     &                            - SIGNAC * T3AAB(ACB-1)
                 END DO
              END DO
           END DO
        END IF
C
        IF(IRPB.EQ.IRPC)THEN
C
C       It follows that IRPC > IRPA
C
        SIGNAC =  1.0D+00
        DO  290 C=1,VRT(IRPC,1)
        DO  280 B=1,VRT(IRPB,1)
        DO  270 A=1,VRT(IRPA,1)
C
        ABC = IADT3(IRPC) + (C-1) * IRPDPD(IRPAB,1) + 
     1        IOFFVV(IRPB,IRPAB,1) + (B-1)*VRT(IRPA,1) + A
C
        AC = IOFFVV(IRPC,IRPAC,1) + (C-1)*VRT(IRPA,1) + A
        IF(B.LT.C)THEN
        SIGNBC =  1.0D+00
        BC     =  IOFFVV(IRPC,IRPBC,1) + INDEX(C-1) + B
        ELSE
        SIGNBC = -1.0D+00
        BC     =  IOFFVV(IRPB,IRPBC,1) + INDEX(B-1) + C
        ENDIF
C
        BCA = IADT3(IRPA) + (A-1) * IRPDPD(IRPBC,1) + BC
        ACB = IADT3(IRPB) + (B-1) * IRPDPD(IRPAC,1) + AC
C
        IF(B.NE.C)THEN
        T3AAA(ABC-1) = T3AAA(ABC-1) + SIGNBC * T3AAB(BCA-1)
     1                              - SIGNAC * T3AAB(ACB-1)
        ELSE
        T3AAA(ABC-1) = T3AAA(ABC-1) - SIGNAC * T3AAB(ACB-1)
        ENDIF
  270   CONTINUE
  280   CONTINUE
  290   CONTINUE
        ENDIF
C
        IF(IRPB.LT.IRPC)THEN
C
C       It follows that IRPC > IRPA
C
        SIGNBC =  1.0D+00
        SIGNAC =  1.0D+00
        DO  320 C=1,VRT(IRPC,1)
        DO  310 B=1,VRT(IRPB,1)
        DO  300 A=1,VRT(IRPA,1)
C
        ABC = IADT3(IRPC) + (C-1) * IRPDPD(IRPAB,1) + 
     1        IOFFVV(IRPB,IRPAB,1) + (B-1)*VRT(IRPA,1) + A
C
        AC = IOFFVV(IRPC,IRPAC,1) + (C-1)*VRT(IRPA,1) + A
        BC = IOFFVV(IRPC,IRPBC,1) + (C-1)*VRT(IRPB,1) + B
C
        BCA = IADT3(IRPA) + (A-1) * IRPDPD(IRPBC,1) + BC
        ACB = IADT3(IRPB) + (B-1) * IRPDPD(IRPAC,1) + AC
C
        T3AAA(ABC-1) = T3AAA(ABC-1) + SIGNBC * T3AAB(BCA-1)
     1                              - SIGNAC * T3AAB(ACB-1)
  300   CONTINUE
  310   CONTINUE
  320   CONTINUE
        ENDIF

        IF(IRPB.GT.IRPC.AND.IRPA.LT.IRPC)THEN
C
        SIGNBC = -1.0D+00
        SIGNAC =  1.0D+00
        DO  350 C=1,VRT(IRPC,1)
        DO  340 B=1,VRT(IRPB,1)
        DO  330 A=1,VRT(IRPA,1)
C
        ABC = IADT3(IRPC) + (C-1) * IRPDPD(IRPAB,1) + 
     1        IOFFVV(IRPB,IRPAB,1) + (B-1)*VRT(IRPA,1) + A
C
        AC = IOFFVV(IRPC,IRPAC,1) + (C-1)*VRT(IRPA,1) + A
        BC = IOFFVV(IRPB,IRPBC,1) + (B-1)*VRT(IRPC,1) + C
C
        BCA = IADT3(IRPA) + (A-1) * IRPDPD(IRPBC,1) + BC
        ACB = IADT3(IRPB) + (B-1) * IRPDPD(IRPAC,1) + AC
C
        T3AAA(ABC-1) = T3AAA(ABC-1) + SIGNBC * T3AAB(BCA-1)
     1                              - SIGNAC * T3AAB(ACB-1)
  330   CONTINUE
  340   CONTINUE
  350   CONTINUE
        ENDIF
C
  400 CONTINUE
C
      ENDIF
  500 CONTINUE
      RETURN
      END
