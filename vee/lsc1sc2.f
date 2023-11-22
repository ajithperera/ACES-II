      SUBROUTINE LSC1SC2(T3AAA,T3AAB,IADT3,IRPIJK)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION T3AAA(1),T3AAB(1)
      INTEGER A,B,C,AC,BC,ABC,ACB,BCA
      INTEGER POP,VRT,DIRPRD
      DIMENSION IADT3(8)
      DIMENSION SIGN(600),IADAC(600),IADBC(600)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     1                DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYM/    POP(8,2),VRT(8,2),NTAA,NTBB,NF1AA,NF2AA,
     1                NF1BB,NF2BB
      COMMON /T3OFF/  IOFFVV(8,8,10),IOFFOO(8,8,10),IOFFVO(8,8,4)
C
      INDEX(I) = I*(I-1)/2
C
c      write(6,*) ioffvv
C     SUBROUTINE THAT IS SUPPOSED TO GENERATE ALL ALPHA TRIPLES FROM
C     CASE 2 (AAB) TRIPLES FOR RHF, ACCORDING TO
C
C      ABC       ABc        BCa       ACb
C     T      =  T     +    T     -   T
C      IJK       IJK        IJK       IJK
C
C     IN SYMMETRY
C
C
      DO  500 IRPC=1,NIRREP
C
c      write(6,*) ' @sc1sc2-i, write 1 '
c      write(6,*) irpc,nirrep,vrt,pop
C
      IF(VRT(IRPC,1).EQ.0) GOTO 500
C
      IRPAB = DIRPRD(IRPC,IRPIJK)
C
      IF(IRPAB.EQ.1)THEN
C
c      write(6,*) ' @sc1sc2-i, write 2 '
C
      DO  100 IRPB=1,NIRREP
c      write(6,*) ' @sc1sc2-i, write 3 '
C
      IF(VRT(IRPB,1).LT.2) GOTO 100
c      write(6,*) ' @sc1sc2-i, write 4 '
C
      IRPA  = IRPB
      IRPAC = DIRPRD(IRPA,IRPC)
      IRPBC = DIRPRD(IRPB,IRPC)
C
        IF(IRPB.EQ.IRPC)THEN
C
C       All virtual irreps are equal. Signs horrid.
C
c      write(6,*) ' @sc1sc2-i, write 5 '
        DO   30    C=1,VRT(IRPC,1)
C
        DO    5    I=1,VRT(IRPC,1)
c      write(6,*) ' @sc1sc2-i, write 6 '
        IF(I.LT.C)THEN
        SIGN(I) =  1.0D+00
        IADAC(I) = IOFFVV(IRPC,IRPAC,1) + INDEX(C-1) + I
C        IADBC(I) = IOFFVV(IRPC,IRPBC,1) + INDEX(C-1) + I
        ENDIF
        IF(I.EQ.C) SIGN(I) = 0.0D+00
        IF(I.GT.C)THEN
        SIGN(I) = -1.0D+00
        IADAC(I) = IOFFVV(IRPA,IRPAC,1) + INDEX(I-1) + C
C        IADBC(I) = IOFFVV(IRPB,IRPBC,1) + INDEX(I-1) + C
        ENDIF
      IADBC(I) = IADAC(I)
c      write(6,*) ' iadac(i),i ', iadac(i),i
c      write(6,*) ' @sc1sc2-i, write 7 '
    5 CONTINUE
c      write(6,*) ' @sc1sc2-i, write 5A '
C
        DO   20    B=2,VRT(IRPB,1)
        DO   10    A=1,B-1
C
c      write(6,*) ' @sc1sc2-i, write 8 '
        ABC = IADT3(IRPC) + (C-1) * IRPDPD(IRPAB,1) + 
     1        IOFFVV(IRPB,IRPAB,1) + INDEX(B-1) + A
C
C        IF(B.LT.C)THEN
C        BC = IOFFVV(IRPC,IRPBC,1) + INDEX(C-1) + B
CC        SIGNBC =  1.0D+00
C        ELSE
C        BC = IOFFVV(IRPB,IRPBC,1) + INDEX(B-1) + C
CC        SIGNBC = -1.0D+00
C        ENDIF
CC
C        IF(A.LT.C)THEN
C        AC = IOFFVV(IRPC,IRPAC,1) + INDEX(C-1) + A
CC        SIGNAC =  1.0D+00
C        ELSE
C        AC = IOFFVV(IRPA,IRPAC,1) + INDEX(A-1) + C
CC        SIGNAC = -1.0D+00
C        ENDIF
C
C        BCA = IADT3(IRPA) + (A-1) * IRPDPD(IRPBC,1) + BC
C        ACB = IADT3(IRPB) + (B-1) * IRPDPD(IRPAC,1) + AC
        BCA = IADT3(IRPA) + (A-1) * IRPDPD(IRPBC,1) + IADBC(B)
c      write(6,*) ' @sc1sc2-i, write 9 '
        ACB = IADT3(IRPB) + (B-1) * IRPDPD(IRPAC,1) + IADAC(A)
c      write(6,*) ' @sc1sc2-i, write 10 '
C
CC        IF(B.NE.C.AND.A.NE.C)THEN
C        T3AAA(ABC-1) = T3AAA(ABC-1) + SIGNBC * T3AAB(BCA-1)
C     1                              - SIGNAC * T3AAB(ACB-1)
c      write(6,*) ' @sc1sc2-i, write 1 '
c      write(6,*) abc,bca,acb,b,a
c       write(6,*) iadt3(irpb),irpa,irpb,a,b,c,irpdpd(irpac,1),iadac(a)
c       write(6,*) sign(a),sign(b)
   
      IF(B.NE.C.AND.A.NE.C)THEN
        T3AAA(ABC-1) = T3AAA(ABC-1) + SIGN(B) * T3AAB(BCA-1)
     1                              - SIGN(A) * T3AAB(ACB-1)
      ENDIF
C
      IF(B.EQ.C.AND.A.NE.C)THEN
        T3AAA(ABC-1) = T3AAA(ABC-1)
     1                              - SIGN(A) * T3AAB(ACB-1)
      ENDIF
C
      IF(B.NE.C.AND.A.EQ.C)THEN
        T3AAA(ABC-1) = T3AAA(ABC-1) + SIGN(B) * T3AAB(BCA-1)
      ENDIF
C
c      write(6,*) ' @sc1sc2-i, write 1 '
CC        ENDIF
C        IF(B.NE.C.AND.A.EQ.C)THEN
CC        T3AAA(ABC-1) = T3AAA(ABC-1) + SIGNBC * T3AAB(BCA-1)
C        T3AAA(ABC-1) = T3AAA(ABC-1) + SIGN(B) * T3AAB(BCA-1)
C        ENDIF
C        IF(B.EQ.C.AND.A.NE.C)THEN
CC        T3AAA(ABC-1) = T3AAA(ABC-1) - SIGNAC * T3AAB(ACB-1)
C        T3AAA(ABC-1) = T3AAA(ABC-1) - SIGN(A) * T3AAB(ACB-1)
C        ENDIF
cc      write(6,*) ' @sc1sc2-i, write 11 '
   10   CONTINUE
c      write(6,*) ' @sc1sc2-i, write 12 '
   20   CONTINUE
c      write(6,*) ' @sc1sc2-i, write 13 '
   30   CONTINUE
c      write(6,*) ' @sc1sc2-i, write 6 '
C
        ENDIF
C
        IF(IRPB.LT.IRPC)THEN
C
C       Both signs as in formula.
C
        SIGNBC =  1.0D+00
        SIGNAC =  1.0D+00
        DO   60    C=1,VRT(IRPC,1)
        DO   50    B=2,VRT(IRPB,1)
        DO   40    A=1,B-1
C
        ABC = IADT3(IRPC) + (C-1) * IRPDPD(IRPAB,1) + 
     1        IOFFVV(IRPB,IRPAB,1) + INDEX(B-1) + A
C
        BC = IOFFVV(IRPC,IRPBC,1) + (C-1)*VRT(IRPB,1) + B
        AC = IOFFVV(IRPC,IRPAC,1) + (C-1)*VRT(IRPA,1) + A
C
        BCA = IADT3(IRPA) + (A-1) * IRPDPD(IRPBC,1) + BC
        ACB = IADT3(IRPB) + (B-1) * IRPDPD(IRPAC,1) + AC
C
        T3AAA(ABC-1) = T3AAA(ABC-1) + SIGNBC * T3AAB(BCA-1)
     1                              - SIGNAC * T3AAB(ACB-1)
c      write(6,*) ' @sc1sc2-i, write 14 '
   40   CONTINUE
   50   CONTINUE
   60   CONTINUE
C
        ENDIF
C
        IF(IRPB.GT.IRPC)THEN
C
C       Both signs reversed.
C
        SIGNBC = -1.0D+00
        SIGNAC = -1.0D+00
        DO   90    C=1,VRT(IRPC,1)
        DO   80    B=2,VRT(IRPB,1)
        DO   70    A=1,B-1
C
        ABC = IADT3(IRPC) + (C-1) * IRPDPD(IRPAB,1) + 
     1        IOFFVV(IRPB,IRPAB,1) + INDEX(B-1) + A
C
        BC = IOFFVV(IRPB,IRPBC,1) + (B-1)*VRT(IRPC,1) + C
        AC = IOFFVV(IRPA,IRPAC,1) + (A-1)*VRT(IRPC,1) + C
C
        BCA = IADT3(IRPA) + (A-1) * IRPDPD(IRPBC,1) + BC
        ACB = IADT3(IRPB) + (B-1) * IRPDPD(IRPAC,1) + AC
C
        T3AAA(ABC-1) = T3AAA(ABC-1) + SIGNBC * T3AAB(BCA-1)
     1                              - SIGNAC * T3AAB(ACB-1)
c      write(6,*) ' @sc1sc2-i, write 15 '
   70   CONTINUE
   80   CONTINUE
   90   CONTINUE
C
        ENDIF
C
  100 CONTINUE
C
      ENDIF
C
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
        IF(IRPA.GT.IRPC)THEN
C
C       It follows that IRPB > IRPC
C
        SIGNAC = -1.0D+00
        SIGNBC = -1.0D+00
        DO  260 C=1,VRT(IRPC,1)
        DO  250 B=1,VRT(IRPB,1)
        DO  240 A=1,VRT(IRPA,1)
C
        ABC = IADT3(IRPC) + (C-1) * IRPDPD(IRPAB,1) + 
     1        IOFFVV(IRPB,IRPAB,1) + (B-1)*VRT(IRPA,1) + A
C
        BC = IOFFVV(IRPB,IRPBC,1) + (B-1)*VRT(IRPC,1) + C
        AC = IOFFVV(IRPA,IRPAC,1) + (A-1)*VRT(IRPC,1) + C
C
        BCA = IADT3(IRPA) + (A-1) * IRPDPD(IRPBC,1) + BC
        ACB = IADT3(IRPB) + (B-1) * IRPDPD(IRPAC,1) + AC
C
        T3AAA(ABC-1) = T3AAA(ABC-1) + SIGNBC * T3AAB(BCA-1)
     1                              - SIGNAC * T3AAB(ACB-1)
  240   CONTINUE
  250   CONTINUE
  260   CONTINUE
        ENDIF
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
