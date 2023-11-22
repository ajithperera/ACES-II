      SUBROUTINE DIAG4I(EIP,FIA,NVIR,NIP,NOCC,NAO)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z) 
C
      DIMENSION EIP(NIP,NIP)
      DIMENSION FIA(NOCC,NVIR)
C
C EIP = SUM_J,A 2*FJA - SUM_A,I FIA
C
      DO 20 IK=1, NOCC
         DO 10 IA=1, NVIR
            DO 30 IJ=1, NOCC
               DO 40 II=1, NOCC
                  IF (II.EQ.IK) THEN
                     IF (MOD(II,NOCC).EQ.0) THEN
                        EIP(IK,NOCC+NOCC+(IJ-1)*NOCC+
     &(IA-1)*NOCC*NOCC)=EIP(IK,NOCC+NOCC+(IJ-1)*NOCC+
     &(IA-1)*NOCC*NOCC)+2*FIA(IJ,IA)
                     ELSE
                        EIP(IK,NOCC+MOD(II,NOCC)+
     &(IJ-1)*NOCC+(IA-1)*NOCC*NOCC)=EIP(IK,NOCC+MOD(II,NOCC)+
     &(IJ-1)*NOCC+(IA-1)*NOCC*NOCC)+2*FIA(IJ,IA)
                     ENDIF
                  ENDIF
                  IF (IK.EQ.IJ) THEN
                     IF (MOD(II,NOCC).EQ.0) THEN
                        EIP(IK,NOCC+NOCC+(IJ-1)*NOCC+
     &(IA-1)*NOCC*NOCC)=EIP(IK,NOCC+NOCC+(IJ-1)*NOCC+
     &(IA-1)*NOCC*NOCC)-FIA(II,IA)
                     ELSE
                        EIP(IK,NOCC+MOD(II,NOCC)+
     &(IJ-1)*NOCC+(IA-1)*NOCC*NOCC)=EIP(IK,NOCC+MOD(II,NOCC)+
     &(IJ-1)*NOCC+(IA-1)*NOCC*NOCC)-FIA(II,IA)
                     ENDIF
                  ENDIF
 40            CONTINUE
 30         CONTINUE
 10      CONTINUE
 20   CONTINUE
C
      RETURN
      END      

      
