C     SUBROUTINE FINDO(NAO,N,BIA,IVRT,IOCC,C,PRS,F)
      SUBROUTINE FINDO(NAO,N,IVRT,IOCC,C,PRS,F)
C
C   F-matrix(Two-electron Part) by INDO METHOD 
C   integrals are passed in eV
C
C                                        Hideo Sekino
C                                        Q.T.P.
C                                        UFL @ GNV
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/PIND/ NAT,N1(500),N2(500),N3(100)
      COMMON/PIND1/ G(100,100),G2(5),G1(5)
      DIMENSION PRS(NAO,NAO),F(NAO,NAO)
     X ,C(NAO,NAO),AP(100),AO2S(100)
C     DIMENSION BIA(N),IVRT(N),IOCC(N)
      DIMENSION IVRT(N),IOCC(N)
C     DIMENSION F2(5),G1(5)
C     DATA F2/0., 4.727,5.961, 7.294,8.593/
C     DATA G1/0., 7.285, 9.416, 11.816, 14.485 /
C
      EAUNIT=27.2113957D0
      C6=1.D0/6.D0
C     CALL DENSP(NAO,N,IVRT,IOCC,C,BIA,PRS)
      DO 180 I=1,NAT
      AP(I)=0.D0
      DO 180 J=1,NAO
      IF(N2(J) .NE.I) GO TO 180
      AP(I)=AP(I)+PRS(J,J)
      IF(N1(J).EQ.2) AO2S(I)=PRS(J,J)
  180 CONTINUE
      DO 300 I=1,NAO
      DO 300 J=1,I
      N1I=N1(I)
      N1J=N1(J)
      N2I=N2(I)
      N2J=N2(J)
      N3I=N3(N2I)
      IF(I.NE.J) GO TO 250
      GO TO(215,220,225,225,225),N1I
  215 F(I,I)= 0.5D0*PRS(I,I)*G(N2I,N2I)
      GO TO 230
  220 F(I,I)= 0.5D0*PRS(I,I)*G(N2I,N2I)+(AP(N2I)
     1       -PRS(I,I))*(G(N2I,N2I)-C6*G1(N3I))
      GO TO 230
  225 F(I,I)= 0.5D0*PRS(I,I)*(G(N2I,N2I)+0.16D0*G2(N3I))
     1       +(AP(N2I)-PRS(I,I))*(G(N2I,N2I)-0.14D0
     2       *G2(N3I))+AO2S(N2I)*(0.14D0*G2(N3I)-C6*G1(N3I))
  230 DO 240 K=1,NAT
      N3K=N3(K)
      IF(N2I.NE.K) F(I,I)=F(I,I)+AP(K)*G(N2I,K)
  240 CONTINUE
      F(I,I)=F(I,I)/EAUNIT
      GO TO 300
  250 IF(N2I.NE.N2J) GO TO 270
      IF(N1I.NE.2.AND. N1J. NE. 2) GO TO 260
      F(I,J)=-0.5D0*PRS(I,J)*(G(N2I,N2I)-G1(N3I))
      GO TO 280
  260 F(I,J)=-0.5D0*PRS(I,J)*(G(N2I,N2I)-0.44D0*G2(N3I))
      GO TO 280
  270 F(I,J)= -0.5D0*PRS(I,J)*G(N2I,N2J)
  280 F(I,J)=F(I,J)/EAUNIT
      F(J,I)=F(I,J)
  300 CONTINUE
      CALL MATMUL(F(1,1),C,PRS,NAO,NAO,NAO,1,0)
      CALL TRANSQ(C,NAO)
      CALL MATMUL(C,PRS,F(1,1),NAO,NAO,NAO,1,0)
      CALL TRANSQ(C,NAO)
      RETURN      
      END       
