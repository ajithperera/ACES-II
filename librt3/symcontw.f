      SUBROUTINE SYMCONTW(T3,W,IADW,ISPIN,IRPIJK)
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION T3(1),W(1)
      DIMENSION IADW(8)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     1                DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYM/    POP(8,2),VRT(8,2),NTAA,NTBB,NF1AA,NF2AA,
     1                NF1BB,NF2BB
C
      COMMON /T3OFF/  IOFFVV(8,8,10),IOFFOO(8,8,10),IOFFVO(8,8,4)
C
      INDEX(I) = I*(I-1)/2
C
C     SUBROUTINE THAT IS SUPPOSED TO DO THE EXPANSION
C
C     T3(A<B<C) = W(B<C,A) - W(A<C,B) + W(A<B,C)
C
C     IN SYMMETRY
C
      ISTART = 0
C
      IF(NIRREP.LT.3) GOTO 1010
C
C     ALL IRREPS OF A,B,C DIFFERENT
C
      DO  1000 IRPC=3,NIRREP
      IF(VRT(IRPC,ISPIN).EQ.0) GOTO 1000
      DO   990 IRPB=2,IRPC-1
      IF(VRT(IRPB,ISPIN).EQ.0) GOTO  990
      DO   980 IRPA=1,IRPB-1
      IF(VRT(IRPA,ISPIN).EQ.0) GOTO  980
C
      IRPAB =  DIRPRD(IRPA,IRPB)
      IRPAC =  DIRPRD(IRPA,IRPC)
      IRPBC =  DIRPRD(IRPB,IRPC)
      IRPABC = DIRPRD(IRPA,IRPBC)
C
      IF(IRPABC.NE.IRPIJK) GOTO 980
C
      DO  430 C=1,VRT(IRPC,ISPIN)
      DO  420 B=1,VRT(IRPB,ISPIN)
      DO  410 A=1,VRT(IRPA,ISPIN)
C
      ABC = ISTART + (C-1)*VRT(IRPB,ISPIN)*VRT(IRPA,ISPIN)
     1             + (B-1)*VRT(IRPA,ISPIN) + A
C
      T3(ABC) = T3(ABC) + W(IADW(IRPC) + (C-1)*IRPDPD(IRPAB,ISPIN)
     1                  + IOFFVV(IRPB,IRPAB,ISPIN) 
     1                  + (B-1)*VRT(IRPA,ISPIN)
     1                  + A - 1)
     1                  - W(IADW(IRPB) + (B-1)*IRPDPD(IRPAC,ISPIN)
     1                  + IOFFVV(IRPC,IRPAC,ISPIN) 
     1                  + (C-1)*VRT(IRPA,ISPIN)
     1                  + A - 1)
     1                  + W(IADW(IRPA) + (A-1)*IRPDPD(IRPBC,ISPIN)
     1                  + IOFFVV(IRPC,IRPBC,ISPIN) 
     1                  + (C-1)*VRT(IRPB,ISPIN)
     1                  + B - 1)
C
C     IF THIS IS RIGHT, WHAT WE EFFECTIVELY HAVE IS
C
C     T3(A,B,C) = W1(A,B,C) - W2(A,C,B) + W3(B,C,A)
C
C     (NOTE DIFFERENT DIMENSIONS OF W1, W2, W3). ULTIMATELY REPLACE
C     INNER LOOPS BY SUBROUTINE CALL SINCE THAT WILL ALLOW BETTER
C     VECTORIZATION.
C
  410 CONTINUE
  420 CONTINUE
  430 CONTINUE
C
      ISTART = ISTART + VRT(IRPC,ISPIN)*VRT(IRPB,ISPIN)*VRT(IRPA,ISPIN)
C
  980 CONTINUE
  990 CONTINUE
 1000 CONTINUE
C
 1010 CONTINUE
C
      IF(NIRREP.LT.2) GOTO 3010
C
C     IRREPS OF B,C ARE THE SAME, IRREP OF A IS DIFFERENT
C
      DO  2000 IRPC=2,NIRREP
      IRPB = IRPC
      IF(VRT(IRPC,ISPIN).LT.2) GOTO 2000
      DO  1980 IRPA=1,IRPC-1
      IF(VRT(IRPA,ISPIN).EQ.0) GOTO 1980
C
      IRPAB =  DIRPRD(IRPA,IRPB)
      IRPAC =  DIRPRD(IRPA,IRPC)
      IRPBC =  DIRPRD(IRPB,IRPC)
      IRPABC = DIRPRD(IRPA,IRPBC)
C
      IF(IRPABC.NE.IRPIJK) GOTO 1980
C
      DO  1430 C=2,VRT(IRPC,ISPIN)
      DO  1420 B=1,C-1
      DO  1410 A=1,VRT(IRPA,ISPIN)
C
C      ABC = ISTART + ((VRT(IRPC,ISPIN)-1)*VRT(IRPC,ISPIN)/2 - 1)*
C     1                VRT(IRPA,ISPIN) + A
      ABC = ISTART + (INDEX(C-1) + B - 1) * VRT(IRPA,ISPIN) + A
C
      T3(ABC) = T3(ABC) + W(IADW(IRPC) + (C-1)*IRPDPD(IRPAB,ISPIN)
     1                  + IOFFVV(IRPB,IRPAB,ISPIN) 
     1                  + (B-1)*VRT(IRPA,ISPIN)
     1                  + A - 1)
     1                  - W(IADW(IRPB) + (B-1)*IRPDPD(IRPAC,ISPIN)
     1                  + IOFFVV(IRPC,IRPAC,ISPIN) 
     1                  + (C-1)*VRT(IRPA,ISPIN)
     1                  + A - 1)
     1                  + W(IADW(IRPA) + (A-1)*IRPDPD(IRPBC,ISPIN)
     1                  + IOFFVV(IRPC,IRPBC,ISPIN) 
     1                  + INDEX(C-1) + B - 1)
C
C     IF THIS IS RIGHT, WHAT WE EFFECTIVELY HAVE IS
C
C     T3(A,B,C) = W1(A,B,C) - W1(A,C,B) + W2(B<C,A)
C
 1410 CONTINUE
 1420 CONTINUE
 1430 CONTINUE
C
      ISTART = ISTART + ((VRT(IRPC,ISPIN)-1)*VRT(IRPC,ISPIN)/2)
     1                  *VRT(IRPA,ISPIN)
C
 1980 CONTINUE
 2000 CONTINUE
C
C     IRREPS OF A,B ARE THE SAME, IRREP OF C IS DIFFERENT
C
      DO  3000 IRPC=2,NIRREP
      IF(VRT(IRPC,ISPIN).EQ.0) GOTO 3000
      DO  2990 IRPB=1,IRPC-1
      IRPA = IRPB
      IF(VRT(IRPB,ISPIN).LT.2) GOTO 2990
C
      IRPAB =  DIRPRD(IRPA,IRPB)
      IRPAC =  DIRPRD(IRPA,IRPC)
      IRPBC =  DIRPRD(IRPB,IRPC)
      IRPABC = DIRPRD(IRPA,IRPBC)
C
      IF(IRPABC.NE.IRPIJK) GOTO 2990
C
      DO  2430 C=1,VRT(IRPC,ISPIN)
      DO  2420 B=2,VRT(IRPB,ISPIN)
      DO  2410 A=1,B-1
C
      ABC = ISTART + (C-1)*(VRT(IRPB,ISPIN)-1)*VRT(IRPB,ISPIN)/2
     1               + INDEX(B-1) + A
C
      T3(ABC) = T3(ABC) + W(IADW(IRPC) + (C-1)*IRPDPD(IRPAB,ISPIN)
     1                  + IOFFVV(IRPB,IRPAB,ISPIN) 
     1                  + INDEX(B-1) + A - 1)
     1                  - W(IADW(IRPB) + (B-1)*IRPDPD(IRPAC,ISPIN)
     1                  + IOFFVV(IRPC,IRPAC,ISPIN) 
     1                  + (C-1)*VRT(IRPA,ISPIN)
     1                  + A - 1)
     1                  + W(IADW(IRPA) + (A-1)*IRPDPD(IRPBC,ISPIN)
     1                  + IOFFVV(IRPC,IRPBC,ISPIN) 
     1                  + (C-1)*VRT(IRPB,ISPIN)
     1                  + B - 1)
C
C     IF THIS IS RIGHT, WHAT WE EFFECTIVELY HAVE IS
C
C     T3(A,B,C) = W1(A<B,C) - W2(A,C,B) + W2(B,C,A)
C
 2410 CONTINUE
 2420 CONTINUE
 2430 CONTINUE
C
      ISTART = ISTART + ((VRT(IRPB,ISPIN)-1)*VRT(IRPB,ISPIN)/2)
     1                  *VRT(IRPC,ISPIN)
C
 2990 CONTINUE
 3000 CONTINUE
C
 3010 CONTINUE
C
C     ALL IRREPS OF ARE THE SAME. THIS IS THE ONLY PIECE OF CODE
C     EXECUTED IN C1 JOBS.
C
      DO  4000 IRPC=1,NIRREP
      IRPB = IRPC
      IRPA = IRPC
      IF(VRT(IRPC,ISPIN).LT.3) GOTO 4000
C
      IRPAB =  DIRPRD(IRPA,IRPB)
      IRPAC =  DIRPRD(IRPA,IRPC)
      IRPBC =  DIRPRD(IRPB,IRPC)
      IRPABC = DIRPRD(IRPA,IRPBC)
C
      IF(IRPABC.NE.IRPIJK) GOTO 4000
C
      DO  3430 C=3,VRT(IRPC,ISPIN)
      DO  3420 B=2,C-1
      DO  3410 A=1,B-1
C
      ABC = ISTART + (C-3)*(C-2)*(C-1)/6 + INDEX(B-1) + A
C
      T3(ABC) = T3(ABC) + W(IADW(IRPC) + (C-1)*IRPDPD(IRPAB,ISPIN)
     1                  + IOFFVV(IRPB,IRPAB,ISPIN) 
     1                  + INDEX(B-1) + A - 1)
     1                  - W(IADW(IRPB) + (B-1)*IRPDPD(IRPAC,ISPIN)
     1                  + IOFFVV(IRPC,IRPAC,ISPIN) 
     1                  + INDEX(C-1) + A - 1)
     1                  + W(IADW(IRPA) + (A-1)*IRPDPD(IRPBC,ISPIN)
     1                  + IOFFVV(IRPC,IRPBC,ISPIN) 
     1                  + INDEX(C-1) + B - 1)
C
C     IF THIS IS RIGHT, WHAT WE EFFECTIVELY HAVE IS
C
C     T3(A,B,C) = W1(A<B,C) - W1(A<C,B) + W1(B<C,A)
C
 3410 CONTINUE
 3420 CONTINUE
 3430 CONTINUE
C
      ISTART = ISTART + (VRT(IRPC,ISPIN)-2)*(VRT(IRPC,ISPIN)-1)
     1                  *VRT(IRPC,ISPIN)/6
C
 4000 CONTINUE
C
      RETURN
      END
