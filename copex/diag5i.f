      SUBROUTINE DIAG5I(EIP,WBAR,NVIR,NIP,NOCC,NAO)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z) 
C
      DIMENSION EIP(NIP,NIP)
      DIMENSION WBAR(NAO,NAO,NAO,NAO)     
C
C EIP = SUM_J,A (2*<JB|AL> - <BJ|AL>) + SUM_I,A (-<IB|AL> - <IB|KA>)
C
      DO 10 IB=1, NVIR
         DO 20 IL=1, NOCC
            DO 30 IK=1, NOCC
               DO 40 IA=1, NVIR
                  DO 50 IJ=1, NOCC
                     DO 60 II=1, NOCC
                        IF ((MOD(IK,NOCC).EQ.0).AND.(MOD(II
     &,NOCC).EQ.0)) THEN
                           IF (IK.EQ.II) THEN
                              EIP(NOCC+NOCC+(IL-1)*NOCC+(IB-1)*NOCC*NOCC
     &,NOCC+NOCC+(IJ-1)*NOCC+(IA-1)*NOCC*NOCC)=EIP(NOCC+NOCC+(IL-1)*NOCC
     &+(IB-1)*NOCC*NOCC,NOCC+NOCC+(IJ-1)*NOCC+(IA-1)*NOCC*NOCC)+
     &2*WBAR(IJ,IB+NOCC,IA+NOCC,IL)-WBAR(IJ,IB+NOCC,IL,IA+NOCC)
                           ENDIF
                           IF (IK.EQ.IJ) THEN
                              EIP(NOCC+NOCC+(IL-1)*NOCC+(IB-1)*NOCC*NOCC
     &,NOCC+NOCC+(IJ-1)*NOCC+(IA-1)*NOCC*NOCC)=EIP(NOCC+NOCC+(IL-1)*NOCC
     &+(IB-1)*NOCC*NOCC,NOCC+NOCC+(IJ-1)*NOCC+(IA-1)*NOCC*NOCC)-WBAR(II,
     &IB+NOCC,IA+NOCC,IL)
                           ENDIF
                           IF (IL.EQ.IJ) THEN
                              EIP(NOCC+NOCC+(IL-1)*NOCC+(IB-1)*NOCC*NOCC
     &,NOCC+NOCC+(IJ-1)*NOCC+(IA-1)*NOCC*NOCC)=EIP(NOCC+NOCC+(IL-1)*NOCC
     &+(IB-1)*NOCC*NOCC,NOCC+NOCC+(IJ-1)*NOCC+(IA-1)*NOCC*NOCC)-WBAR(II,
     &IB+NOCC,IK,IA+NOCC)
                           ENDIF
                        ELSE
                           IF ((MOD(IK,NOCC).EQ.0).OR.(MOD(II
     &,NOCC).EQ.0)) THEN
                              IF (MOD(IK,NOCC).EQ.0) THEN
                                 IF (IK.EQ.II) THEN
                                    EIP(NOCC+NOCC+(IL-1)*NOCC+(IB-1)*
     &NOCC*NOCC,NOCC+MOD(II,NOCC)+(IJ-1)*NOCC+(IA-1)*NOCC*NOCC)=EIP(NOCC
     &+NOCC+(IL-1)*NOCC+(IB-1)*NOCC*NOCC,NOCC+MOD(II,NOCC)+(IJ-1)*NOCC+
     &(IA-1)*NOCC*NOCC)+2*WBAR(IJ,IB+NOCC,IA+NOCC,IL)-WBAR(IJ,IB+NOCC,IL
     &,IA+NOCC)
                                 ENDIF
                                 IF (IK.EQ.IJ) THEN
                                    EIP(NOCC+NOCC+(IL-1)*NOCC+(IB-1)
     &*NOCC*NOCC,NOCC+MOD(II,NOCC)+(IJ-1)*NOCC+(IA-1)*NOCC*NOCC)=EIP(
     &NOCC+NOCC+(IL-1)*NOCC+(IB-1)*NOCC*NOCC,NOCC+MOD(II,NOCC)+(IJ-1)*
     &NOCC+(IA-1)*NOCC*NOCC)-WBAR(II,IB+NOCC,IA+NOCC,IL)
                                 ENDIF
                                 IF (IL.EQ.IJ) THEN
                                    EIP(NOCC+NOCC+(IL-1)*NOCC+(IB-1)*
     &NOCC*NOCC,NOCC+MOD(II,NOCC)+(IJ-1)*NOCC+(IA-1)*NOCC*NOCC)=EIP(NOCC
     &+NOCC+(IL-1)*NOCC+(IB-1)*NOCC*NOCC,NOCC+MOD(II,NOCC)+(IJ-1)*NOCC+
     &(IA-1)*NOCC*NOCC)-WBAR(II,IB+NOCC,IK,IA+NOCC)
                                 ENDIF
                              ELSE
                                 IF (IK.EQ.II) THEN
                                    EIP(NOCC+MOD(IK,NOCC)+(IL-1)*NOCC+
     &(IB-1)*NOCC*NOCC,NOCC+NOCC+(IJ-1)*NOCC+(IA-1)*NOCC*NOCC)=EIP(NOCC+
     &MOD(IK,NOCC)+(IL-1)*NOCC+(IB-1)*NOCC*NOCC,NOCC+NOCC+(IJ-1)*NOCC+
     &(IA-1)*NOCC*NOCC)+2*WBAR(IJ,IB+NOCC,IA+NOCC,IL)-WBAR(IJ,IB+NOCC,IL
     &,IA+NOCC)
                                 ENDIF
                                 IF (IK.EQ.IJ) THEN
                                    EIP(NOCC+MOD(IK,NOCC)+(IL-1)*NOCC+
     &(IB-1)*NOCC*NOCC,NOCC+NOCC+(IJ-1)*NOCC+(IA-1)*NOCC*NOCC)=EIP(NOCC+
     &MOD(IK,NOCC)+(IL-1)*NOCC+(IB-1)*NOCC*NOCC,NOCC+NOCC+(IJ-1)*NOCC+
     &(IA-1)*NOCC*NOCC)-WBAR(II,IB+NOCC,IA+NOCC,IL)
                                 ENDIF
                                 IF (IL.EQ.IJ) THEN
                                    EIP(NOCC+MOD(IK,NOCC)+(IL-1)*NOCC+
     &(IB-1)*NOCC*NOCC,NOCC+NOCC+(IJ-1)*NOCC+(IA-1)*NOCC*NOCC)=EIP(NOCC+
     &MOD(IK,NOCC)+(IL-1)*NOCC+(IB-1)*NOCC*NOCC,NOCC+NOCC+(IJ-1)*NOCC+
     &(IA-1)*NOCC*NOCC)-WBAR(II,IB+NOCC,IK,IA+NOCC)
                                 ENDIF
                              ENDIF
                           ELSE
                              IF (IK.EQ.II) THEN
                                 EIP(NOCC+MOD(IK,NOCC)+(IL-1)*NOCC+(IB-1
     &)*NOCC*NOCC,NOCC+MOD(II,NOCC)+(IJ-1)*NOCC+(IA-1)*NOCC*NOCC)=EIP(
     &NOCC+MOD(IK,NOCC)+(IL-1)*NOCC+(IB-1)*NOCC*NOCC,NOCC+MOD(II,NOCC)+
     &(IJ-1)*NOCC+(IA-1)*NOCC*NOCC)+2*WBAR(IJ,IB+NOCC,IA+NOCC,IL)-WBAR(
     &IJ,IB+NOCC,IL,IA+NOCC)
                              ENDIF
                              IF (IK.EQ.IJ) THEN
                                 EIP(NOCC+MOD(IK,NOCC)+(IL-1)*NOCC+(IB-1
     &)*NOCC*NOCC,NOCC+MOD(II,NOCC)+(IJ-1)*NOCC+(IA-1)*NOCC*NOCC)=EIP(
     &NOCC+MOD(IK,NOCC)+(IL-1)*NOCC+(IB-1)*NOCC*NOCC,NOCC+MOD(II,NOCC)+
     &(IJ-1)*NOCC+(IA-1)*NOCC*NOCC)-WBAR(II,IB+NOCC,IA+NOCC,IL)
                              ENDIF
                              IF (IL.EQ.IJ) THEN
                                 EIP(NOCC+MOD(IK,NOCC)+(IL-1)*NOCC+(IB-1
     &)*NOCC*NOCC,NOCC+MOD(II,NOCC)+(IJ-1)*NOCC+(IA-1)*NOCC*NOCC)=EIP(
     &NOCC+MOD(IK,NOCC)+(IL-1)*NOCC+(IB-1)*NOCC*NOCC,NOCC+MOD(II,NOCC)+
     &(IJ-1)*NOCC+(IA-1)*NOCC*NOCC)-WBAR(II,IB+NOCC,IK,IA+NOCC)
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
