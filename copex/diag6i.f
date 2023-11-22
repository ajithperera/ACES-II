      SUBROUTINE DIAG6I(EIP,WBAR,NVIR,NIP,NOCC,NAO)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z) 
C
      DIMENSION EIP(NIP,NIP)
      DIMENSION WBAR(NAO,NAO,NAO,NAO)
C
C EIP = SUM_I,J <IJ|KL> 
C
      DONE = 1.0D+0
      DONEM = -1.0D+0
      DZERO = 0.D+0
C
      DO 10 IB=1, NVIR
         DO 20 IL=1, NOCC
            DO 30 IK=1, NOCC
               DO 40 IA=1, NVIR
                  DO 50 IJ=1, NOCC
                     DO 60 II=1, NOCC
                        IF ((MOD(IK,NOCC).EQ.0).AND.(MOD(II
     &,NOCC).EQ.0)) THEN
                           IF (IA.EQ.IB) THEN
                              EIP(NOCC+NOCC+(IL-1)*NOCC+(IB-1)*NOCC*NOCC
     &,NOCC+NOCC+(IJ-1)*NOCC+(IA-1)*NOCC*NOCC)=EIP(NOCC+NOCC+(IL-1)*NOCC
     &+(IB-1)*NOCC*NOCC,NOCC+NOCC+(IJ-1)*NOCC+(IA-1)*NOCC*NOCC)+
     &WBAR(II,IJ,IK,IL)
                           ENDIF
                        ELSE
                           IF ((MOD(IK,NOCC).EQ.0).OR.(MOD(II
     &,NOCC).EQ.0)) THEN
                              IF (MOD(IK,NOCC).EQ.0) THEN
                                 IF (IA.EQ.IB) THEN
                                    EIP(NOCC+NOCC+(IL-1)*NOCC+(IB-1)*
     &NOCC*NOCC,NOCC+MOD(II,NOCC)+(IJ-1)*NOCC+(IA-1)*NOCC*NOCC)=EIP(NOCC
     &+NOCC+(IL-1)*NOCC+(IB-1)*NOCC*NOCC,NOCC+MOD(II,NOCC)+(IJ-1)*NOCC+
     &(IA-1)*NOCC*NOCC)+WBAR(II,IJ,IK,IL)
                                 ENDIF
                              ELSE
                                 IF (IA.EQ.IB) THEN
                                    EIP(NOCC+MOD(IK,NOCC)+(IL-1)*NOCC+
     &(IB-1)*NOCC*NOCC,NOCC+NOCC+(IJ-1)*NOCC+(IA-1)*NOCC*NOCC)=EIP(NOCC+
     &MOD(IK,NOCC)+(IL-1)*NOCC+(IB-1)*NOCC*NOCC,NOCC+NOCC+(IJ-1)*NOCC+
     &(IA-1)*NOCC*NOCC)+WBAR(II,IJ,IK,IL)
                                 ENDIF
                              ENDIF
                           ELSE
                              IF (IA.EQ.IB) THEN
                                 EIP(NOCC+MOD(IK,NOCC)+(IL-1)*NOCC+(IB-1
     &)*NOCC*NOCC,NOCC+MOD(II,NOCC)+(IJ-1)*NOCC+(IA-1)*NOCC*NOCC)=EIP
     &(NOCC+MOD(IK,NOCC)+(IL-1)*NOCC+(IB-1)*NOCC*NOCC,NOCC+MOD(II,NOCC)+
     &(IJ-1)*NOCC+(IA-1)*NOCC*NOCC)+WBAR(II,IJ,IK,IL)
                              ENDIF
                           ENDIF                      
                        ENDIF
 60                  CONTINUE
 50               CONTINUE
 40            CONTINUE
 30         CONTINUE
 20      CONTINUE
 10   CONTINUE
C
      RETURN
      END
