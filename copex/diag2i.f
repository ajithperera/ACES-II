      SUBROUTINE DIAG2I(EIP,WBAR,NVIR,NIP,NOCC,NAO)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z) 
C
      DIMENSION EIP(NIP,NIP)
      DIMENSION WBAR(NAO,NAO,NAO,NAO)
C
C EIP = SUM_I -<IB|KL>
C
      DO 10 IB=1, NVIR
         DO 20 II=1, NOCC
            DO 30 IL=1, NOCC
               DO 40 IK=1, NOCC
                  IF (MOD(IK,NOCC).EQ.0) THEN
                     EIP(NOCC+NOCC+(IL-1)*NOCC+(IB-1)*NOCC*NOCC,II)
     &=-WBAR(II,IB+NOCC,IK,IL) 
                  ELSE
                     EIP(NOCC+MOD(IK,NOCC)+(IL-1)*NOCC+(IB-1)
     &*NOCC*NOCC,II)=-WBAR(II,IB+NOCC,IK,IL) 
                  ENDIF
 40            CONTINUE
 30         CONTINUE
 20      CONTINUE
 10   CONTINUE
C
      RETURN
      END
