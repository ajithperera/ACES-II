C     SUBROUTINE FIMDO(NAO,N,BIA,IVRT,IOCC,C,PRS,F)
      SUBROUTINE FIMDO(NAO,N,IVRT,IOCC,C,PRS,F)
C   
C   Anti-symmetric Fock Matrix
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
     X ,C(NAO,NAO)
C     DIMENSION BIA(N),IVRT(N),IOCC(N)
      DIMENSION IVRT(N),IOCC(N)
C     DIMENSION F2(5),G1(5)
C     DATA F2/0., 4.727,5.961, 7.294,8.593/
C     DATA G1/0., 7.285, 9.416, 11.816, 14.485 /
C
      EAUNIT=27.2113957D0
      C3= 1.D0/3.D0
C     CALL DENSM(NAO,N,IVRT,IOCC,C,BIA,PRS)
      DO 300 I=1,NAO
      DO 300 J=1,I
      IF(I.NE.J) THEN
      N1I=N1(I)
      N1J=N1(J)
      N2I=N2(I)
      N2J=N2(J)
      N3I=N3(N2I)
      IF(N2I.NE.N2J) GO TO 270
      IF(N1I.NE.2.AND. N1J. NE. 2) GO TO 260
      F(I,J)=-0.5D0*PRS(I,J)*(G(N2I,N2I)-C3*G1(N3I))
      GO TO 280
  260 F(I,J)=-0.5D0*PRS(I,J)*(G(N2I,N2I)-0.2*G2(N3I))
      GO TO 280
  270 F(I,J)= -0.5D0*PRS(I,J)*G(N2I,N2J)
  280 F(I,J)=F(I,J)/EAUNIT
      F(J,I)=-F(I,J)
      ELSE
      F(I,J)=0.D0
      END IF
  300 CONTINUE
      CALL MATMUL(F(1,1),C,PRS,NAO,NAO,NAO,1,0)
      CALL TRANSQ(C,NAO)
      CALL MATMUL(C,PRS,F(1,1),NAO,NAO,NAO,1,0)
      CALL TRANSQ(C,NAO)
      RETURN      
      END       
