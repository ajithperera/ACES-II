      SUBROUTINE DIAG3I(EIP,WBAR,NVIR,NIP,NOCC,NAO)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z) 
C
      DIMENSION EIP(NIP,NIP)
      DIMENSION WBAR(NAO,NAO,NAO,NAO)
C
C EIP = SUM_I,J,A ( <IJ|AK> - 2 <IJ|KA> )
C
      DO 10 IA=1, NVIR
         DO 20 IK=1, NOCC
            DO 30 IJ=1, NOCC
               DO 40 II=1, NOCC
                  IF (MOD(II,NOCC).EQ.0) THEN
                     EIP(IK,NOCC+NOCC+(IJ-1)*NOCC+(IA-1)*NOCC*NOCC)
     &=-2*WBAR(II,IJ,IK,IA+NOCC)+WBAR(IJ,II,IK,IA+NOCC)  
                  ELSE
                     EIP(IK,NOCC+MOD(II,NOCC)+(IJ-1)*NOCC+(IA-1)*NOCC*
     &NOCC)=-2*WBAR(II,IJ,IK,IA+NOCC)+WBAR(IJ,II,IK,IA+NOCC)  
                  ENDIF
 40            CONTINUE
 30         CONTINUE
 20      CONTINUE
 10   CONTINUE
C
      RETURN
      END      

