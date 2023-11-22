C
      SUBROUTINE FLDNUC(NOC,VLIST,C,NOCMX,NT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION VLIST(NOCMX,*),C(3),SNUC(11), SLASK(10)
      COMMON /PARM/ISL2,ISL3,ISL4,islx,THRESH,MTYP,ILNMCM,INT,INUI,
     & NCENTS,
     1 ID20,NSYM,NAO(12),NMO(12),junk,PTOT(3,10), INTERF, JOBIPH,
     2 RUNMOS
      DO20 I=1,3
20    SNUC(I)=0.0
      DO21 I=1,NOC
      DST=(VLIST(I,1)-C(1))**2+(VLIST(I,2)-C(2))**2+(VLIST(I,3)-C(3))**2
      IF(DST.EQ.0.0) GO TO 21
      DO 23 J=1,3
      SNUC(J)=SNUC(J)+VLIST(I,4)*(VLIST(I,J)-C(J))/(SQRT(DST))**3
23    CONTINUE
21    CONTINUE
      DO22 L=1,3
22    PTOT(1,L) = -SNUC(L)
      NT = 3
      RETURN
C
CRAY 1
      ENTRY NONUC(NOC,VLIST,C,NOCMX,NT)
CDC 1
C      ENTRY NONUC
C     ZERO OUT EXTRA TWO WORDS FOR REL
      PTOT(1,1)=0.0
      PTOT(1,2)=0.0
      PTOT(1,3)=0.0
      NT=1
      RETURN
CRAY 1
      ENTRY POTNUC(NOC,VLIST,C,NOCMX,NT)
CDC 1
C      ENTRY POTNUC
      PNUC=0.0
      DO 12 I=1,NOC
      DST=(VLIST(I,1)-C(1))**2+(VLIST(I,2)-C(2))**2+(VLIST(I,3)-C(3))**2
      IF(DST.EQ.0.0) GO TO 12
      PNUC=PNUC+VLIST(I,4)/SQRT(DST)
12    CONTINUE
      PTOT(1,1)=PNUC
      NT = 1
      RETURN
CRAY 1
      ENTRY GRDNUC(NOC,VLIST,C,NOCMX,NT)
CDC 1
C      ENTRY GRDNUC
      DO 30 I=1,6
30    SNUC(I)=0.0
      DO 31 I=1,NOC
      DST=(VLIST(I,1)-C(1))**2+(VLIST(I,2)-C(2))**2+(VLIST(I,3)-C(3))**2
      IF(DST.EQ.0.0) GO TO 31
      MU=0
      DO 32 J=1,3
      DO 33 K=1,J
      MU=MU+1
      SNUC(MU)=SNUC(MU)+VLIST(I,4)*(VLIST(I,J)-C(J))*(VLIST(I,K)-C(K))
     1/(SQRT(DST))**5
33    CONTINUE
32    CONTINUE
31    CONTINUE
      DST=SNUC(1)+SNUC(3)+SNUC(6)
      PTOT(1,1)=3.*SNUC(1)-DST
      PTOT(1,2)=3.*SNUC(3)-DST
      PTOT(1,3)=3.*SNUC(6)-DST
      PTOT(1,4)=SNUC(2)*3.0
      PTOT(1,5)=SNUC(4)*3.0
      PTOT(1,6)=SNUC(5)*3.0
      NT = 6
      RETURN
CRAY 1
      ENTRY QUDNUC(NOC,VLIST,C,NOCMX,NT)
CDC 1
C      ENTRY QUDNUC
      DO50 I=1,6
50    SNUC(I)=0.0
      DO 51 I=1,NOC
      MU=0
      DO52 J=1,3
      DO53 K=1,J
      MU=MU+1
      SNUC(MU)=SNUC(MU)+VLIST(I,4)*(VLIST(I,J)-C(J))*(VLIST(I,K)-C(K))
53    CONTINUE
52    CONTINUE
51    CONTINUE
      PTOT(1,1)=SNUC(1)
      PTOT(1,2)=SNUC(3)
      PTOT(1,3)=SNUC(6)
      PTOT(1,4)=1.5*SNUC(2)
      PTOT(1,5)=1.5*SNUC(4)
      PTOT(1,6)=1.5*SNUC(5)
      PTOT(1,7)=PTOT(1,1)+PTOT(1,2)+PTOT(1,3)
      PTOT(1,8) = SNUC(1)
      PTOT(1,9) = SNUC(3)
      PTOT(1,10) = SNUC(6)
      PTOT(1,1)=(3.*PTOT(1,1)-PTOT(1,7))/2.
      PTOT(1,2)=(3.*PTOT(1,2)-PTOT(1,7))/2.
      PTOT(1,3) =(3.*PTOT(1,3)-PTOT(1,7))/2.
      NT = 10
      RETURN
CRAY 1
      ENTRY SECNUC(NOC,VLIST,C,NOCMX,NT)
CDC 1
c
c jfs - modified to get rid of nuclear charge contribution
c
C      ENTRY SECNUC
      DO 150 I=1,6
150   SNUC(I)=0.0
c      DO 151 I=1,NOC
c      MU=0
c      DO 152 J=1,3
c      DO 153 K=1,J
c      MU=MU+1
c      SNUC(MU)=SNUC(MU)+VLIST(I,4)*(VLIST(I,J)-C(J))*(VLIST(I,K)-C(K))
c153   CONTINUE
c152   CONTINUE
c151   CONTINUE
      PTOT(1,1) = SNUC(1)
      PTOT(1,2) = SNUC(3)
      PTOT(1,3) = SNUC(6)
      PTOT(1,4) = SNUC(2)
      PTOT(1,5) = SNUC(4)
      PTOT(1,6) = SNUC(5)
      NT = 6
      RETURN
CRAY 1
      ENTRY MAGNUC(NOC,VLIST,C,NOCMX,NT)
CDC 1
C      ENTRY MAGNUC
      DO70 I=1,6
70    SNUC(I)=0.0
      DO 71 I=1,NOC
      MU=0
      DO72 J=1,3
      DO73 K=1,J
      MU=MU+1
      SNUC(MU)=SNUC(MU)+VLIST(I,4)*(VLIST(I,J)-C(J))*(VLIST(I,K)-C(K))
73    CONTINUE
72    CONTINUE
71    CONTINUE
      SLASK(1)=SNUC(1)
      SLASK(2)=SNUC(3)
      SLASK(3)=SNUC(6)
      SLASK(4)=SNUC(2)
      SLASK(5)=SNUC(4)
      SLASK(6)=SNUC(5)
C
C....    NOW TRANSFORM FROM 2ND MOMENTS TO DIAMAGNETIC SUSCEPTIBILITIES
C
      PTOT(1,1) = 1.5*(SLASK(2) + SLASK(3))
      PTOT(1,2) = 1.5*(SLASK(1) + SLASK(3))
      PTOT(1,3) = 1.5*(SLASK(1) + SLASK(2))
      PTOT(1,4) = -1.5*SLASK(4)
      PTOT(1,5) = -1.5*SLASK(5)
      PTOT(1,6) = -1.5*SLASK(6)
      NT = 6
      RETURN
CRAY 1
      ENTRY DIPNUC(NOC,VLIST,C,NOCMX,NT)
CDC 1
C      ENTRY DIPNUC
      DO40 I=1,3
40    SNUC(I)=0.0
      DO41 I=1,NOC
      DO41 J=1,3
41    SNUC(J)=SNUC(J)+VLIST(I,4)*(VLIST(I,J)-C(J))
      DO42 K=1,3
42    PTOT(1,K)=SNUC(K)
      NT = 3
      RETURN
CRAY 1
      ENTRY OCTNUC(NOC,VLIST,C,NOCMX,NT)
CDC 1
C      ENTRY OCTNUC
      DO60 I=1,10
60    SNUC(I)=0.0
      DO61 I=1,NOC
      MU=0
      DO62 J=1,3
      DO63 K=1,J
      DO64 M=1,K
      MU=MU+1
      SNUC(MU)=SNUC(MU)+VLIST(I,4)*(VLIST(I,J)-C(J))*(VLIST(I,K)-C(K))*(
     1VLIST(I,M)-C(M))
64    CONTINUE
63    CONTINUE
62    CONTINUE
61    CONTINUE
      SLASK(1)=SNUC(1)
      SLASK(2)=SNUC(4)
      SLASK(3)=SNUC(10)
      SLASK(4)=SNUC(2)
      SLASK(5)=SNUC(5)
      SLASK(6)=SNUC(3)
      SLASK(7)=SNUC(7)
      SLASK(8)=SNUC(8)
      SLASK(9)=SNUC(9)
      SLASK(10)=SNUC(6)
C
C....    NOW TRANSFORM FROM 3RD MOMENTS TO OCTOPOLE
C
      PTOT(1,1) = SLASK(1) - 1.5*(SLASK(6) + SLASK(8))
      PTOT(1,2) = SLASK(2) - 1.5*(SLASK(4) + SLASK(9))
      PTOT(1,3) = SLASK(3) - 1.5*(SLASK(5) + SLASK(7))
      PTOT(1,4) = 2.*SLASK(4) - 0.5*(SLASK(2) + SLASK(9))
      PTOT(1,5) = 2.*SLASK(5) - 0.5*(SLASK(3) + SLASK(7))
      PTOT(1,6) = 2.*SLASK(6) - 0.5*(SLASK(1) + SLASK(8))
      PTOT(1,7) = 2.*SLASK(7) - 0.5*(SLASK(3) + SLASK(5))
      PTOT(1,8) = 2.*SLASK(8) - 0.5*(SLASK(1) + SLASK(6))
      PTOT(1,9) = 2.*SLASK(9) - 0.5*(SLASK(2) + SLASK(4))
      PTOT(1,10) = SLASK(10)
      NT = 10
      RETURN
CRAY 1
      ENTRY THINUC(NOC,VLIST,C,NOCMX,NT)
CDC 1
C      ENTRY THINUC
      DO 160 I=1,10
160   SNUC(I)=0.0
      DO 161 I=1,NOC
      MU=0
      DO 162 J=1,3
      DO 163 K=1,J
      DO 164 M=1,K
      MU=MU+1
      SNUC(MU)=SNUC(MU)+VLIST(I,4)*(VLIST(I,J)-C(J))*(VLIST(I,K)-C(K))*(
     1VLIST(I,M)-C(M))
164   CONTINUE
163   CONTINUE
162   CONTINUE
161   CONTINUE
      PTOT(1,1) = SNUC(1)
      PTOT(1,2) = SNUC(4)
      PTOT(1,3) = SNUC(10)
      PTOT(1,4) = SNUC(2)
      PTOT(1,5) = SNUC(5)
      PTOT(1,6) = SNUC(3)
      PTOT(1,7) = SNUC(7)
      PTOT(1,8) = SNUC(8)
      PTOT(1,9) = SNUC(9)
      PTOT(1,10) = SNUC(6)
      NT = 10
      RETURN
      END
