      SUBROUTINE DIAG1I(EIP,FIJ,F,NOCC,NIP,NAO)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z) 
C
      DIMENSION EIP(NIP,NIP)
      DIMENSION FIJ(NOCC,NOCC)
      DIMENSION F(NAO)
C
C EIP = SUM_I -FIJ
C
      DO 20 I1=1, NOCC
         DO 30 I2=1, NOCC
            EIP(I1,I2)=-FIJ(I2,I1)
            IF (I1.EQ.I2) THEN
               EIP(I1,I2)=EIP(I1,I2)-F(I1)
            ENDIF
 30      CONTINUE
 20   CONTINUE
C
      RETURN
      END
