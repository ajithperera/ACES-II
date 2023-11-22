      SUBROUTINE DIAG8E(EEA,FIJ,FAB,F,NVIR,NEA,NOCC,NAO)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z) 
C
      DIMENSION EEA(NEA,NEA)
      DIMENSION FIJ(NOCC,NOCC)
      DIMENSION FAB(NVIR,NVIR)
      DIMENSION F(NAO)
C
C EEA = SUM_A FAB - SUM_J FJL - SUM_I FIK
C
      DO 10 IJ=1, NOCC
         DO 20 IC=1, NVIR
            DO 30 ID=1, NVIR
               DO 40 II=1, NOCC
                  DO 50 IA=1, NVIR
                     DO 60 IB=1, NVIR
                        IF ((MOD(ID,NVIR).EQ.0).AND.(MOD(IB
     &,NVIR).EQ.0)) THEN
                           IF ((ID.EQ.IB).AND.(II.EQ.IJ)) THEN
                              EEA(NVIR+NVIR+(IC-1)*NVIR+(IJ-1)*NVIR*NVIR
     &,NVIR+NVIR+(IA-1)*NVIR+(II-1)*NVIR*NVIR)=EEA(NVIR+NVIR+(IC-1)*NVIR
     &+(IJ-1)*NVIR*NVIR,NVIR+NVIR+(IA-1)*NVIR+(II-1)*NVIR*NVIR)+FAB(IC,
     &IA)
                              IF (IC.EQ.IA) THEN
                                 EEA(NVIR+NVIR+(IC-1)*NVIR+(IJ-1)*NVIR*
     &NVIR,NVIR+NVIR+(IA-1)*NVIR+(II-1)*NVIR*NVIR)=EEA(NVIR+NVIR+(IC-1)*
     &NVIR+(IJ-1)*NVIR*NVIR,NVIR+NVIR+(IA-1)*NVIR+(II-1)*NVIR*NVIR)+
     &F(IA+NOCC)
                              ENDIF
                           ENDIF
                           IF ((ID.EQ.IB).AND.(IC.EQ.IA)) THEN
                              EEA(NVIR+NVIR+(IC-1)*NVIR+(IJ-1)*NVIR*NVIR
     &,NVIR+NVIR+(IA-1)*NVIR+(II-1)*NVIR*NVIR)=EEA(NVIR+NVIR+(IC-1)*NVIR
     &+(IJ-1)*NVIR*NVIR,NVIR+NVIR+(IA-1)*NVIR+(II-1)*NVIR*NVIR)-
     &FIJ(II,IJ)
                              IF (II.EQ.IJ) THEN
                                 EEA(NVIR+NVIR+(IC-1)*NVIR+(IJ-1)*NVIR*
     &NVIR,NVIR+NVIR+(IA-1)*NVIR+(II-1)*NVIR*NVIR)=EEA(NVIR+NVIR+(IC-1)*
     &NVIR+(IJ-1)*NVIR*NVIR,NVIR+NVIR+(IA-1)*NVIR+(II-1)*NVIR*NVIR)-
     &F(II)
                              ENDIF
                           ENDIF
                           IF ((IC.EQ.IA).AND.(II.EQ.IJ)) THEN
                              EEA(NVIR+NVIR+(IC-1)*NVIR+(IJ-1)*NVIR*NVIR
     &,NVIR+NVIR+(IA-1)*NVIR+(II-1)*NVIR*NVIR)=EEA(NVIR+NVIR+(IC-1)*NVIR
     &+(IJ-1)*NVIR*NVIR,NVIR+NVIR+(IA-1)*NVIR+(II-1)*NVIR*NVIR)+
     &FAB(ID,IB)
                              IF (ID.EQ.IB) THEN
                                 EEA(NVIR+NVIR+(IC-1)*NVIR+(IJ-1)*NVIR
     &*NVIR,NVIR+NVIR+(IA-1)*NVIR+(II-1)*NVIR*NVIR)=EEA(NVIR+NVIR+(IC-1)
     &*NVIR+(IJ-1)*NVIR*NVIR,NVIR+NVIR+(IA-1)*NVIR+(II-1)*NVIR*NVIR)+
     &F(IB+NOCC)
                              ENDIF
                           ENDIF
                        ELSE
                           IF ((MOD(ID,NVIR).EQ.0).OR.(MOD(IB
     &,NVIR).EQ.0)) THEN
                              IF (MOD(ID,NVIR).EQ.0) THEN
                                 IF ((ID.EQ.IB).AND.(II.EQ.IJ)) THEN
                                    EEA(NVIR+NVIR+(IC-1)*NVIR+(IJ-1)*
     &NVIR*NVIR,NVIR+MOD(IB,NVIR)+(IA-1)*NVIR+(II-1)*NVIR*NVIR)=EEA(NVIR
     &+NVIR+(IC-1)*NVIR+(IJ-1)*NVIR*NVIR,NVIR+MOD(IB,NVIR)+(IA-1)*NVIR+
     &(II-1)*NVIR*NVIR)+FAB(IC,IA)
                                    IF (IC.EQ.IA) THEN
                                       EEA(NVIR+NVIR+(IC-1)*NVIR+(IJ-1)*
     &NVIR*NVIR,NVIR+MOD(IB,NVIR)+(IA-1)*NVIR+(II-1)*NVIR*NVIR)=EEA(NVIR
     &+NVIR+(IC-1)*NVIR+(IJ-1)*NVIR*NVIR,NVIR+MOD(IB,NVIR)+(IA-1)*NVIR+
     &(II-1)*NVIR*NVIR)+F(IA+NOCC)
                                    ENDIF
                                 ENDIF
                                 IF ((ID.EQ.IB).AND.(IA.EQ.IC)) THEN
                                    EEA(NVIR+NVIR+(IC-1)*NVIR+(IJ-1)*
     &NVIR*NVIR,NVIR+MOD(IB,NVIR)+(IA-1)*NVIR+(II-1)*NVIR*NVIR)=EEA(NVIR
     &+NVIR+(IC-1)*NVIR+(IJ-1)*NVIR*NVIR,NVIR+MOD(IB,NVIR)+(IA-1)*NVIR+
     &(II-1)*NVIR*NVIR)-FIJ(II,IJ)
                                    IF (II.EQ.IJ) THEN
                                       EEA(NVIR+NVIR+(IC-1)*NVIR+(IJ-1)*
     &NVIR*NVIR,NVIR+MOD(IB,NVIR)+(IA-1)*NVIR+(II-1)*NVIR*NVIR)=EEA(NVIR
     &+NVIR+(IC-1)*NVIR+(IJ-1)*NVIR*NVIR,NVIR+MOD(IB,NVIR)+(IA-1)*NVIR+
     &(II-1)*NVIR*NVIR)-F(II)
                                    ENDIF
                                 ENDIF
                                 IF ((IC.EQ.IA).AND.(II.EQ.IJ)) THEN
                                    EEA(NVIR+NVIR+(IC-1)*NVIR+(IJ-1)*
     &NVIR*NVIR,NVIR+MOD(IB,NVIR)+(IA-1)*NVIR+(II-1)*NVIR*NVIR)=EEA(NVIR
     &+NVIR+(IC-1)*NVIR+(IJ-1)*NVIR*NVIR,NVIR+MOD(IB,NVIR)+(IA-1)*NVIR+
     &(II-1)*NVIR*NVIR)+FAB(ID,IB)
                                    IF (ID.EQ.IB) THEN
                                       EEA(NVIR+NVIR+(IC-1)*NVIR+(IJ-1)*
     &NVIR*NVIR,NVIR+MOD(IB,NVIR)+(IA-1)*NVIR+(II-1)*NVIR*NVIR)=EEA(NVIR
     &+NVIR+(IC-1)*NVIR+(IJ-1)*NVIR*NVIR,NVIR+MOD(IB,NVIR)+(IA-1)*NVIR+
     &(II-1)*NVIR*NVIR)+F(IB+NOCC)
                                    ENDIF
                                 ENDIF
                              ELSE
                                 IF ((ID.EQ.IB).AND.(II.EQ.IJ)) THEN
                                    EEA(NVIR+MOD(ID,NVIR)+(IC-1)*NVIR+
     &(IJ-1)*NVIR*NVIR,NVIR+NVIR+(IA-1)*NVIR+(II-1)*NVIR*NVIR)=EEA(NVIR+
     &MOD(ID,NVIR)+(IC-1)*NVIR+(IJ-1)*NVIR*NVIR,NVIR+NVIR+(IA-1)*NVIR+
     &(II-1)*NVIR*NVIR)+FAB(IC,IA)
                                    IF (IC.EQ.IA) THEN
                                       EEA(NVIR+MOD(ID,NVIR)+(IC-1)*NVIR
     &+(IJ-1)*NVIR*NVIR,NVIR+NVIR+(IA-1)*NVIR+(II-1)*NVIR*NVIR)=EEA(NVIR
     &+MOD(ID,NVIR)+(IC-1)*NVIR+(IJ-1)*NVIR*NVIR,NVIR+NVIR+(IA-1)*NVIR+
     &(II-1)*NVIR*NVIR)+F(IA+NOCC)
                                    ENDIF
                                 ENDIF
                                 IF ((ID.EQ.IB).AND.(IA.EQ.IC)) THEN
                                    EEA(NVIR+MOD(ID,NVIR)+(IC-1)*NVIR+
     &(IJ-1)*NVIR*NVIR,NVIR+NVIR+(IA-1)*NVIR+(II-1)*NVIR*NVIR)=EEA(NVIR+
     &MOD(ID,NVIR)+(IC-1)*NVIR+(IJ-1)*NVIR*NVIR,NVIR+NVIR+(IA-1)*NVIR+
     &(II-1)*NVIR*NVIR)-FIJ(II,IJ)
                                    IF (II.EQ.IJ) THEN
                                       EEA(NVIR+MOD(ID,NVIR)+(IC-1)*NVIR
     &+(IJ-1)*NVIR*NVIR,NVIR+NVIR+(IA-1)*NVIR+(II-1)*NVIR*NVIR)=EEA(NVIR
     &+MOD(ID,NVIR)+(IC-1)*NVIR+(II-1)*NVIR*NVIR,NVIR+NVIR+(IA-1)*NVIR+
     &(II-1)*NVIR*NVIR)-F(II)
                                    ENDIF
                                 ENDIF
                                 IF ((IC.EQ.IA).AND.(II.EQ.IJ)) THEN
                                    EEA(NVIR+MOD(ID,NVIR)+(IC-1)*NVIR+
     &(IJ-1)*NVIR*NVIR,NVIR+NVIR+(IA-1)*NVIR+(II-1)*NVIR*NVIR)=EEA(NVIR+
     &MOD(ID,NVIR)+(IC-1)*NVIR+(IJ-1)*NVIR*NVIR,NVIR+NVIR+(IA-1)*NVIR+
     &(II-1)*NVIR*NVIR)+FAB(ID,IB)
                                    IF (ID.EQ.IB) THEN
                                       EEA(NVIR+MOD(ID,NVIR)+(IC-1)*NVIR
     &+(IJ-1)*NVIR*NVIR,NVIR+NVIR+(IA-1)*NVIR+(II-1)*NVIR*NVIR)=EEA(NVIR
     &+MOD(ID,NVIR)+(IC-1)*NVIR+(IJ-1)*NVIR*NVIR,NVIR+NVIR+(IA-1)*NVIR+
     &(II-1)*NVIR*NVIR)+F(IB+NOCC)
                                    ENDIF
                                 ENDIF
                              ENDIF
                           ELSE
                              IF ((ID.EQ.IB).AND.(II.EQ.IJ)) THEN
                                 EEA(NVIR+MOD(ID,NVIR)+(IC-1)*NVIR+(IJ-1
     &)*NVIR*NVIR,NVIR+MOD(IB,NVIR)+(IA-1)*NVIR+(II-1)*NVIR*NVIR)=EEA(
     &NVIR+MOD(ID,NVIR)+(IC-1)*NVIR+(IJ-1)*NVIR*NVIR,NVIR+MOD(IB,NVIR)+
     &(IA-1)*NVIR+(II-1)*NVIR*NVIR)+FAB(IC,IA)
                                 IF (IC.EQ.IA) THEN
                                    EEA(NVIR+MOD(ID,NVIR)+(IC-1)*NVIR+
     &(IJ-1)*NVIR*NVIR,NVIR+MOD(IB,NVIR)+(IA-1)*NVIR+(II-1)*NVIR*NVIR)=
     &EEA(NVIR+MOD(ID,NVIR)+(IC-1)*NVIR+(IJ-1)*NVIR*NVIR,NVIR+MOD(IB,
     &NVIR)+(IA-1)*NVIR+(II-1)*NVIR*NVIR)+F(IA+NOCC)
                                 ENDIF
                              ENDIF
                              IF ((ID.EQ.IB).AND.(IA.EQ.IC)) THEN
                                 EEA(NVIR+MOD(ID,NVIR)+(IC-1)*NVIR+(IJ-1
     &)*NVIR*NVIR,NVIR+MOD(IB,NVIR)+(IA-1)*NVIR+(II-1)*NVIR*NVIR)=EEA(
     &NVIR+MOD(ID,NVIR)+(IC-1)*NVIR+(IJ-1)*NVIR*NVIR,NVIR+MOD(IB,NVIR)+
     &(IA-1)*NVIR+(II-1)*NVIR*NVIR)-FIJ(II,IJ)
                                 IF (II.EQ.IJ) THEN
                                    EEA(NVIR+MOD(ID,NVIR)+(IC-1)*NVIR+
     &(IJ-1)*NVIR*NVIR,NVIR+MOD(IB,NVIR)+(IA-1)*NVIR+(II-1)*NVIR*NVIR)=
     &EEA(NVIR+MOD(ID,NVIR)+(IC-1)*NVIR+(IJ-1)*NVIR*NVIR,NVIR+MOD(IB,
     &NVIR)+(IA-1)*NVIR+(II-1)*NVIR*NVIR)-F(II)
                                 ENDIF
                              ENDIF
                              IF ((IC.EQ.IA).AND.(II.EQ.IJ)) THEN
                                 EEA(NVIR+MOD(ID,NVIR)+(IC-1)*NVIR+(IJ-1
     &)*NVIR*NVIR,NVIR+MOD(IB,NVIR)+(IA-1)*NVIR+(II-1)*NVIR*NVIR)=EEA(
     &NVIR+MOD(ID,NVIR)+(IC-1)*NVIR+(IJ-1)*NVIR*NVIR,NVIR+MOD(IB,NVIR)+
     &(IA-1)*NVIR+(II-1)*NVIR*NVIR)+FAB(ID,IB)
                                 IF (ID.EQ.IB) THEN
                                    EEA(NVIR+MOD(ID,NVIR)+(IC-1)*NVIR+
     &(IJ-1)*NVIR*NVIR,NVIR+MOD(IB,NVIR)+(IA-1)*NVIR+(II-1)*NVIR*NVIR)=
     &EEA(NVIR+MOD(ID,NVIR)+(IC-1)*NVIR+(IJ-1)*NVIR*NVIR,NVIR+MOD(IB,
     &NVIR)+(IA-1)*NVIR+(II-1)*NVIR*NVIR)+F(IB+NOCC)
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
