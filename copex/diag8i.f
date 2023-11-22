      SUBROUTINE DIAG8I(EIP,FIJ,FAB,F,NVIR,NIP,NOCC,NAO)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z) 
C
      DIMENSION EIP(NIP,NIP)
      DIMENSION FIJ(NOCC,NOCC)
      DIMENSION FAB(NVIR,NVIR)
      DIMENSION F(NAO)
C
C EIP = SUM_A FBA - SUM_J FJL - SUM_I FIK 
C
      DO 10 IB=1, NVIR
         DO 20 IL=1, NOCC
            DO 30 IK=1, NOCC
               DO 40 IA=1, NVIR
                  DO 50 IJ=1, NOCC
                     DO 60 II=1, NOCC
                        IF ((MOD(IK,NOCC).EQ.0).AND.(MOD(II
     &,NOCC).EQ.0)) THEN
                           IF ((IK.EQ.II).AND.(IL.EQ.IJ)) THEN
                              EIP(NOCC+NOCC+(IL-1)*NOCC+(IB-1)*NOCC*NOCC
     &,NOCC+NOCC+(IJ-1)*NOCC+(IA-1)*NOCC*NOCC)=EIP(NOCC+NOCC+(IL-1)*NOCC
     &+(IB-1)*NOCC*NOCC,NOCC+NOCC+(IJ-1)*NOCC+(IA-1)*NOCC*NOCC)+
     &FAB(IB,IA)
                              IF (IA.EQ.IB) THEN
                                 EIP(NOCC+NOCC+(IL-1)*NOCC+(IB-1)*NOCC*
     &NOCC,NOCC+NOCC+(IJ-1)*NOCC+(IA-1)*NOCC*NOCC)=EIP(NOCC+NOCC+(IL-1)
     &*NOCC+(IB-1)*NOCC*NOCC,NOCC+NOCC+(IJ-1)*NOCC+(IA-1)*NOCC*NOCC)
     &+F(IA+NOCC)
                              ENDIF
                           ENDIF
                           IF ((IK.EQ.II).AND.(IA.EQ.IB)) THEN
                              EIP(NOCC+NOCC+(IL-1)*NOCC+(IB-1)*NOCC*NOCC
     &,NOCC+NOCC+(IJ-1)*NOCC+(IA-1)*NOCC*NOCC)=EIP(NOCC+NOCC+(IL-1)*NOCC
     &+(IB-1)*NOCC*NOCC,NOCC+NOCC+(IJ-1)*NOCC+(IA-1)*NOCC*NOCC)-
     &FIJ(IJ,IL)
                              IF (IL.EQ.IJ) THEN
                                 EIP(NOCC+NOCC+(IL-1)*NOCC+(IB-1)*NOCC*
     &NOCC,NOCC+NOCC+(IJ-1)*NOCC+(IA-1)*NOCC*NOCC)=EIP(NOCC+NOCC+(IL-1)
     &*NOCC+(IB-1)*NOCC*NOCC,NOCC+NOCC+(IJ-1)*NOCC+(IA-1)*NOCC*NOCC)-
     &F(IL)
                              ENDIF
                           ENDIF
                           IF ((IL.EQ.IJ).AND.(IA.EQ.IB)) THEN
                              EIP(NOCC+NOCC+(IL-1)*NOCC+(IB-1)*NOCC*NOCC
     &,NOCC+NOCC+(IJ-1)*NOCC+(IA-1)*NOCC*NOCC)=EIP(NOCC+NOCC+(IL-1)*NOCC
     &+(IB-1)*NOCC*NOCC,NOCC+NOCC+(IJ-1)*NOCC+(IA-1)*NOCC*NOCC)-
     &FIJ(II,IK)
                              IF (IK.EQ.II) THEN
                                 EIP(NOCC+NOCC+(IL-1)*NOCC+(IB-1)*NOCC*
     &NOCC,NOCC+NOCC+(IJ-1)*NOCC+(IA-1)*NOCC*NOCC)=EIP(NOCC+NOCC+(IL-1)*
     &NOCC+(IB-1)*NOCC*NOCC,NOCC+NOCC+(IJ-1)*NOCC+(IA-1)*NOCC*NOCC)-
     &F(II)
                              ENDIF
                           ENDIF
                        ELSE
                           IF ((MOD(IK,NOCC).EQ.0).OR.(MOD(II
     &,NOCC).EQ.0)) THEN
                              IF (MOD(IK,NOCC).EQ.0) THEN
                                 IF ((IK.EQ.II).AND.(IL.EQ.IJ)) THEN
                                    EIP(NOCC+NOCC+(IL-1)*NOCC+(IB-1)*
     &NOCC*NOCC,NOCC+MOD(II,NOCC)+(IJ-1)*NOCC+(IA-1)*NOCC*NOCC)=EIP(NOCC
     &+NOCC+(IL-1)*NOCC+(IB-1)*NOCC*NOCC,NOCC+MOD(II,NOCC)+(IJ-1)*NOCC+
     &(IA-1)*NOCC*NOCC)+FAB(IB,IA)
                                    IF (IA.EQ.IB) THEN
                                       EIP(NOCC+NOCC+(IL-1)*NOCC+(IB-1)*
     &NOCC*NOCC,NOCC+MOD(II,NOCC)+(IJ-1)*NOCC+(IA-1)*NOCC*NOCC)=EIP(NOCC
     &+NOCC+(IL-1)*NOCC+(IB-1)*NOCC*NOCC,NOCC+MOD(II,NOCC)+(IJ-1)*NOCC+
     &(IA-1)*NOCC*NOCC)+F(IA+NOCC)
                                    ENDIF
                                 ENDIF
                                 IF ((IK.EQ.II).AND.(IA.EQ.IB)) THEN
                                    EIP(NOCC+NOCC+(IL-1)*NOCC+(IB-1)*
     &NOCC*NOCC,NOCC+MOD(II,NOCC)+(IJ-1)*NOCC+(IA-1)*NOCC*NOCC)=EIP(NOCC
     &+NOCC+(IL-1)*NOCC+(IB-1)*NOCC*NOCC,NOCC+MOD(II,NOCC)+(IJ-1)*NOCC+
     &(IA-1)*NOCC*NOCC)-FIJ(IJ,IL)
                                    IF (IL.EQ.IJ) THEN
                                       EIP(NOCC+NOCC+(IL-1)*NOCC+(IB-1)*
     &NOCC*NOCC,NOCC+MOD(II,NOCC)+(IJ-1)*NOCC+(IA-1)*NOCC*NOCC)=EIP(NOCC
     &+NOCC+(IL-1)*NOCC+(IB-1)*NOCC*NOCC,NOCC+MOD(II,NOCC)+(IJ-1)*NOCC+
     &(IA-1)*NOCC*NOCC)-F(IJ)
                                    ENDIF
                                 ENDIF
                                 IF ((IL.EQ.IJ).AND.(IA.EQ.IB)) THEN
                                    EIP(NOCC+NOCC+(IL-1)*NOCC+(IB-1)*
     &NOCC*NOCC,NOCC+MOD(II,NOCC)+(IJ-1)*NOCC+(IA-1)*NOCC*NOCC)=EIP(NOCC
     &+NOCC+(IL-1)*NOCC+(IB-1)*NOCC*NOCC,NOCC+MOD(II,NOCC)+(IJ-1)*NOCC+
     &(IA-1)*NOCC*NOCC)-FIJ(II,IK)
                                    IF (IK.EQ.II) THEN
                                       EIP(NOCC+NOCC+(IL-1)*NOCC+(IB-1)*
     &NOCC*NOCC,NOCC+MOD(II,NOCC)+(IJ-1)*NOCC+(IA-1)*NOCC*NOCC)=EIP(NOCC
     &+NOCC+(IL-1)*NOCC+(IB-1)*NOCC*NOCC,NOCC+MOD(II,NOCC)+(IJ-1)*NOCC+
     &(IA-1)*NOCC*NOCC)-F(II)
                                    ENDIF
                                 ENDIF
                              ELSE
                                 IF ((IK.EQ.II).AND.(IL.EQ.IJ)) THEN
                                    EIP(NOCC+MOD(IK,NOCC)+(IL-1)*NOCC+
     &(IB-1)*NOCC*NOCC,NOCC+NOCC+(IJ-1)*NOCC+(IA-1)*NOCC*NOCC)=EIP(NOCC+
     &MOD(IK,NOCC)+(IL-1)*NOCC+(IB-1)*NOCC*NOCC,NOCC+NOCC+(IJ-1)*NOCC+
     &(IA-1)*NOCC*NOCC)+FAB(IB,IA)
                                    IF (IA.EQ.IB) THEN
                                       EIP(NOCC+MOD(IK,NOCC)+(IL-1)*NOCC
     &+(IB-1)*NOCC*NOCC,NOCC+NOCC+(IJ-1)*NOCC+(IA-1)*NOCC*NOCC)=EIP(NOCC
     &+MOD(IK,NOCC)+(IL-1)*NOCC+(IB-1)*NOCC*NOCC,NOCC+NOCC+(IJ-1)*NOCC+
     &(IA-1)*NOCC*NOCC)+F(IA+NOCC)
                                    ENDIF
                                 ENDIF
                                 IF ((IK.EQ.II).AND.(IA.EQ.IB)) THEN
                                    EIP(NOCC+MOD(IK,NOCC)+(IL-1)*NOCC+
     &(IB-1)*NOCC*NOCC,NOCC+NOCC+(IJ-1)*NOCC+(IA-1)*NOCC*NOCC)=EIP(NOCC+
     &MOD(IK,NOCC)+(IL-1)*NOCC+(IB-1)*NOCC*NOCC,NOCC+NOCC+(IJ-1)*NOCC+
     &(IA-1)*NOCC*NOCC)-FIJ(IJ,IL)
                                    IF (IL.EQ.IJ) THEN
                                       EIP(NOCC+MOD(IK,NOCC)+(IL-1)*NOCC
     &+(IB-1)*NOCC*NOCC,NOCC+NOCC+(IJ-1)*NOCC+(IA-1)*NOCC*NOCC)=EIP(NOCC
     &+MOD(IK,NOCC)+(IL-1)*NOCC+(IB-1)*NOCC*NOCC,NOCC+NOCC+(IJ-1)*NOCC+
     &(IA-1)*NOCC*NOCC)-F(IJ)
                                    ENDIF
                                 ENDIF
                                 IF ((IL.EQ.IJ).AND.(IA.EQ.IB)) THEN
                                    EIP(NOCC+MOD(IK,NOCC)+(IL-1)*NOCC+
     &(IB-1)*NOCC*NOCC,NOCC+NOCC+(IJ-1)*NOCC+(IA-1)*NOCC*NOCC)=EIP(NOCC+
     &MOD(IK,NOCC)+(IL-1)*NOCC+(IB-1)*NOCC*NOCC,NOCC+NOCC+(IJ-1)*NOCC+
     &(IA-1)*NOCC*NOCC)-FIJ(II,IK)
                                    IF (IK.EQ.II) THEN
                                       EIP(NOCC+MOD(IK,NOCC)+(IL-1)*NOCC
     &+(IB-1)*NOCC*NOCC,NOCC+NOCC+(IJ-1)*NOCC+(IA-1)*NOCC*NOCC)=EIP(NOCC
     &+MOD(IK,NOCC)+(IL-1)*NOCC+(IB-1)*NOCC*NOCC,NOCC+NOCC+(IJ-1)*NOCC+
     &(IA-1)*NOCC*NOCC)-F(II)
                                    ENDIF
                                 ENDIF
                              ENDIF
                           ELSE
                              IF ((IK.EQ.II).AND.(IL.EQ.IJ)) THEN
                                 EIP(NOCC+MOD(IK,NOCC)+(IL-1)*NOCC+(IB-1
     &)*NOCC*NOCC,NOCC+MOD(II,NOCC)+(IJ-1)*NOCC+(IA-1)*NOCC*NOCC)=EIP(
     &NOCC+MOD(IK,NOCC)+(IL-1)*NOCC+(IB-1)*NOCC*NOCC,NOCC+MOD(II,NOCC)+
     &(IJ-1)*NOCC+(IA-1)*NOCC*NOCC)+FAB(IB,IA)
                                 IF (IA.EQ.IB) THEN
                                    EIP(NOCC+MOD(IK,NOCC)+(IL-1)*NOCC+
     &(IB-1)*NOCC*NOCC,NOCC+MOD(II,NOCC)+(IJ-1)*NOCC+(IA-1)*NOCC*NOCC)=
     &EIP(NOCC+MOD(IK,NOCC)+(IL-1)*NOCC+(IB-1)*NOCC*NOCC,NOCC+MOD(II,
     &NOCC)+(IJ-1)*NOCC+(IA-1)*NOCC*NOCC)+F(IA+NOCC)
                                 ENDIF
                              ENDIF
                              IF ((IK.EQ.II).AND.(IA.EQ.IB)) THEN
                                 EIP(NOCC+MOD(IK,NOCC)+(IL-1)*NOCC+(IB-1
     &)*NOCC*NOCC,NOCC+MOD(II,NOCC)+(IJ-1)*NOCC+(IA-1)*NOCC*NOCC)=EIP(
     &NOCC+MOD(IK,NOCC)+(IL-1)*NOCC+(IB-1)*NOCC*NOCC,NOCC+MOD(II,NOCC)+
     &(IJ-1)*NOCC+(IA-1)*NOCC*NOCC)-FIJ(IJ,IL)
                                 IF (IL.EQ.IJ) THEN
                                    EIP(NOCC+MOD(IK,NOCC)+(IL-1)*NOCC+
     &(IB-1)*NOCC*NOCC,NOCC+MOD(II,NOCC)+(IJ-1)*NOCC+(IA-1)*NOCC*NOCC)=
     &EIP(NOCC+MOD(IK,NOCC)+(IL-1)*NOCC+(IB-1)*NOCC*NOCC,NOCC+MOD(II,
     &NOCC)+(IJ-1)*NOCC+(IA-1)*NOCC*NOCC)-F(IJ)
                                 ENDIF
                              ENDIF
                              IF ((IL.EQ.IJ).AND.(IA.EQ.IB)) THEN
                                 EIP(NOCC+MOD(IK,NOCC)+(IL-1)*NOCC+(IB-1
     &)*NOCC*NOCC,NOCC+MOD(II,NOCC)+(IJ-1)*NOCC+(IA-1)*NOCC*NOCC)=EIP(
     &NOCC+MOD(IK,NOCC)+(IL-1)*NOCC+(IB-1)*NOCC*NOCC,NOCC+MOD(II,NOCC)+
     &(IJ-1)*NOCC+(IA-1)*NOCC*NOCC)-FIJ(II,IK)
                                 IF (IK.EQ.II) THEN
                                    EIP(NOCC+MOD(IK,NOCC)+(IL-1)*NOCC+
     &(IB-1)*NOCC*NOCC,NOCC+MOD(II,NOCC)+(IJ-1)*NOCC+(IA-1)*NOCC*NOCC)=
     &EIP(NOCC+MOD(IK,NOCC)+(IL-1)*NOCC+(IB-1)*NOCC*NOCC,NOCC+MOD(II,
     &NOCC)+(IJ-1)*NOCC+(IA-1)*NOCC*NOCC)-F(II)
                                 ENDIF
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
