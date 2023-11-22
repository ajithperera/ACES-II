      SUBROUTINE DIAG7I(EIP,WBAR,T2,NVIR,NIP,NOCC,NAO)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z) 
C
      DIMENSION EIP(NIP,NIP)
      DIMENSION WBAR(NAO,NAO,NAO,NAO)
      DIMENSION T2(NOCC,NOCC,NVIR,NVIR)
C
C EIP = SUM_I,J,A,C ( <JI|CA> -2*<IJ|CA> ) T(LK,BC)
C
      DO 10 IB=1, NVIR
         DO 20 IL=1, NOCC
            DO 30 IK=1, NOCC
               DO 40 IA=1, NVIR
                  DO 50 IJ=1, NOCC
                     DO 60 II=1, NOCC
                        IF ((MOD(IK,NOCC).EQ.0).AND.(MOD(II
     &,NOCC).EQ.0)) THEN
                           DO 70 IC=1, NVIR
                              EIP(NOCC+NOCC+(IL-1)*NOCC+(IB-1)*NOCC*NOCC
     &,NOCC+NOCC+(IJ-1)*NOCC+(IA-1)*NOCC*NOCC)=EIP(NOCC+NOCC+(IL-1)*NOCC
     &+(IB-1)*NOCC*NOCC,NOCC+NOCC+(IJ-1)*NOCC+(IA-1)*NOCC*NOCC)-
     &(2*WBAR(II,IJ,IC+NOCC,IA+NOCC)-WBAR(IJ,II,IC+NOCC,IA+NOCC))*
     &T2(IK,IL,IC,IB)
 70                        CONTINUE
                        ELSE
                           IF ((MOD(IK,NOCC).EQ.0).OR.(MOD(II
     &,NOCC).EQ.0)) THEN
                              IF (MOD(IK,NOCC).EQ.0) THEN
                                 DO 90 IC=1, NVIR
                                    EIP(NOCC+NOCC+(IL-1)*NOCC+(IB-1)*
     &NOCC*NOCC,NOCC+MOD(II,NOCC)+(IJ-1)*NOCC+(IA-1)*NOCC*NOCC)=EIP(NOCC
     &+NOCC+(IL-1)*NOCC+(IB-1)*NOCC*NOCC,NOCC+MOD(II,NOCC)+(IJ-1)*NOCC+
     &(IA-1)*NOCC*NOCC)-(2*WBAR(II,IJ,IC+NOCC,IA+NOCC)-WBAR(IJ,II,IC+
     &NOCC,IA+NOCC))*T2(IK,IL,IC,IB)
 90                              CONTINUE
                              ELSE
                                 DO 11 IC=1, NVIR
                                    EIP(NOCC+MOD(IK,NOCC)+(IL-1)*NOCC+
     &(IB-1)*NOCC*NOCC,NOCC+NOCC+(IJ-1)*NOCC+(IA-1)*NOCC*NOCC)=EIP(NOCC+
     &MOD(IK,NOCC)+(IL-1)*NOCC+(IB-1)*NOCC*NOCC,NOCC+NOCC+(IJ-1)*NOCC+
     &(IA-1)*NOCC*NOCC)-(2*WBAR(II,IJ,IC+NOCC,IA+NOCC)-WBAR(IJ,II,IC+
     &NOCC,IA+NOCC))*T2(IK,IL,IC,IB)
 11                              CONTINUE
                              ENDIF
                           ELSE
                              DO 12 IC=1, NVIR
                                 EIP(NOCC+MOD(IK,NOCC)+(IL-1)*NOCC+(IB-1
     &)*NOCC*NOCC,NOCC+MOD(II,NOCC)+(IJ-1)*NOCC+(IA-1)*NOCC*NOCC)=EIP(
     &NOCC+MOD(IK,NOCC)+(IL-1)*NOCC+(IB-1)*NOCC*NOCC,NOCC+MOD(II,NOCC)+
     &(IJ-1)*NOCC+(IA-1)*NOCC*NOCC)-(2*WBAR(II,IJ,IC+NOCC,IA+NOCC)-WBAR(
     &IJ,II,IC+NOCC,IA+NOCC))*T2(IK,IL,IC,IB)
 12                           CONTINUE
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
