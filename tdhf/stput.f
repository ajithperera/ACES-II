      SUBROUTINE STPUT(A,B,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(N,N),B(1)
      COMMON /INFSYM/NSYMHF,NSO(8),NOCS(8),NVTS(8),IDPR(8,8),NVOS(8)
     X ,NIJS(8),NIJS2(8),NRDS(8),NIJSR(8)
      IJ=0
      DO 1 I=1,N
      DO 1 J=1,I
      A(I,J)=0.D0
      NSYMN=1
      NSYMX=0
      DO 2 IS=1,NSYMHF
      NSYMX=NSO(IS)+NSYMX
      IF(IS.GT.1) NSYMN=NSYMX-NSO(IS)+1
      IF((I.GE.NSYMN.AND.I.LE.NSYMX).AND.(J.GE.NSYMN.AND.J.LE.NSYMX))
     X THEN
      IJ=IJ+1
      A(I,J)=B(IJ)
      A(J,I)=A(I,J)
      END IF
    2 CONTINUE
    1 CONTINUE
      RETURN
      END