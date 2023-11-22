      SUBROUTINE DIAG5E(EEA,WBAR,NVIR,NEA,NOCC,NAO)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z) 
C
      DIMENSION EEA(NEA,NEA)
      DIMENSION WBAR(NAO,NAO,NAO,NAO)    
C
C EEA = SUM_J,A (2<JB|AL> - <BJ|AL>) + SUM_I,A (-<IB|AL> - <IB|KA>)
C
      DONE = 1.0D+0
      DONEM = -1.0D+0
      DZERO = 0.D+0
C
      DO 10 IJ=1, NOCC
         DO 20 IC=1, NVIR
            DO 30 ID=1, NVIR
               DO 40 II=1, NOCC
                  DO 50 IA=1, NVIR
                     DO 60 IB=1, NVIR
                        IF ((MOD(ID,NVIR).EQ.0).AND.(MOD(IB
     &,NVIR).EQ.0)) THEN
                           IF (ID.EQ.IB) THEN
                              EEA(NVIR+NVIR+(IC-1)*NVIR+(IJ-1)*NVIR*NVIR
     &,NVIR+NVIR+(IA-1)*NVIR+(II-1)*NVIR*NVIR)=EEA(NVIR+NVIR+(IC-1)*NVIR
     &+(IJ-1)*NVIR*NVIR,NVIR+NVIR+(IA-1)*NVIR+(II-1)*NVIR*NVIR)+
     &2*WBAR(II,IC+NOCC,IA+NOCC,IJ)-WBAR(II,IC+NOCC,IJ,IA+NOCC)
                           ENDIF
                           IF (IA.EQ.ID) THEN
                              EEA(NVIR+NVIR+(IC-1)*NVIR+(IJ-1)*NVIR*NVIR
     &,NVIR+NVIR+(IA-1)*NVIR+(II-1)*NVIR*NVIR)=EEA(NVIR+NVIR+(IC-1)*NVIR
     &+(IJ-1)*NVIR*NVIR,NVIR+NVIR+(IA-1)*NVIR+(II-1)*NVIR*NVIR)-WBAR(II,
     &IC+NOCC,IB+NOCC,IJ)
                           ENDIF
                           IF (IC.EQ.IA) THEN
                              EEA(NVIR+NVIR+(IC-1)*NVIR+(IJ-1)*NVIR*NVIR
     &,NVIR+NVIR+(IA-1)*NVIR+(II-1)*NVIR*NVIR)=EEA(NVIR+NVIR+(IC-1)*NVIR
     &+(IJ-1)*NVIR*NVIR,NVIR+NVIR+(IA-1)*NVIR+(II-1)*NVIR*NVIR)-WBAR(II,
     &ID+NOCC,IJ,IB+NOCC)
                           ENDIF
                        ELSE
                           IF ((MOD(ID,NVIR).EQ.0).OR.(MOD(IB
     &,NVIR).EQ.0)) THEN
                              IF (MOD(ID,NVIR).EQ.0) THEN
                                 IF (ID.EQ.IB) THEN
                                    EEA(NVIR+NVIR+(IC-1)*NVIR+(IJ-1)*
     &NVIR*NVIR,NVIR+MOD(IB,NVIR)+(IA-1)*NVIR+(II-1)*NVIR*NVIR)=EEA(NVIR
     &+NVIR+(IC-1)*NVIR+(IJ-1)*NVIR*NVIR,NVIR+MOD(IB,NVIR)+(IA-1)*NVIR+
     &(II-1)*NVIR*NVIR)+2*WBAR(II,IC+NOCC,IA+NOCC,IJ)-WBAR(II,IC+NOCC,IJ
     &,IA+NOCC)
                                 ENDIF
                                 IF (ID.EQ.IA) THEN
                                    EEA(NVIR+NVIR+(IC-1)*NVIR+(IJ-1)*
     &NVIR*NVIR,NVIR+MOD(IB,NVIR)+(IA-1)*NVIR+(II-1)*NVIR*NVIR)=EEA(NVIR
     &+NVIR+(IC-1)*NVIR+(IJ-1)*NVIR*NVIR,NVIR+MOD(IB,NVIR)+(IA-1)*NVIR+
     &(II-1)*NVIR*NVIR)-WBAR(II,IC+NOCC,IB+NOCC,IJ)
                                 ENDIF
                                 IF (IC.EQ.IA) THEN
                                    EEA(NVIR+NVIR+(IC-1)*NVIR+(IJ-1)
     &*NVIR*NVIR,NVIR+MOD(IB,NVIR)+(IA-1)*NVIR+(II-1)*NVIR*NVIR)=EEA(
     &NVIR+NVIR+(IC-1)*NVIR+(IJ-1)*NVIR*NVIR,NVIR+MOD(IB,NVIR)+(IA-1)
     &*NVIR+(II-1)*NVIR*NVIR)-WBAR(II,ID+NOCC,IJ,IB+NOCC)
                                 ENDIF
                              ELSE
                                 IF (ID.EQ.IB) THEN
                                    EEA(NVIR+MOD(ID,NVIR)+(IC-1)*NVIR+
     &(IJ-1)*NVIR*NVIR,NVIR+NVIR+(IA-1)*NVIR+(II-1)*NVIR*NVIR)=EEA(NVIR+
     &MOD(ID,NVIR)+(IC-1)*NVIR+(IJ-1)*NVIR*NVIR,NVIR+NVIR+(IA-1)*NVIR+
     &(II-1)*NVIR*NVIR)+2*WBAR(II,IC+NOCC,IA+NOCC,IJ)-WBAR(II,IC+NOCC,IJ
     &,IA+NOCC)
                                 ENDIF
                                 IF (ID.EQ.IA) THEN
                                    EEA(NVIR+MOD(ID,NVIR)+(IC-1)*NVIR+
     &(IJ-1)*NVIR*NVIR,NVIR+NVIR+(IA-1)*NVIR+(II-1)*NVIR*NVIR)=EEA(NVIR+
     &MOD(ID,NVIR)+(IC-1)*NVIR+(IJ-1)*NVIR*NVIR,NVIR+NVIR+(IA-1)*NVIR+
     &(II-1)*NVIR*NVIR)-WBAR(II,IC+NOCC,IB+NOCC,IJ)
                                 ENDIF
                                 IF (IC.EQ.IA) THEN
                                    EEA(NVIR+MOD(ID,NVIR)+(IC-1)*NVIR+
     &(IJ-1)*NVIR*NVIR,NVIR+NVIR+(IA-1)*NVIR+(II-1)*NVIR*NVIR)=EEA(NVIR+
     &MOD(ID,NVIR)+(IC-1)*NVIR+(IJ-1)*NVIR*NVIR,NVIR+NVIR+(IA-1)*NVIR+
     &(II-1)*NVIR*NVIR)-WBAR(II,ID+NOCC,IJ,IB+NOCC)
                                 ENDIF
                              ENDIF
                           ELSE
                              IF (ID.EQ.IB) THEN
                                 EEA(NVIR+MOD(ID,NVIR)+(IC-1)*NVIR+
     &(IJ-1)*NVIR*NVIR,NVIR+MOD(IB,NVIR)+(IA-1)*NVIR+(II-1)*NVIR*NVIR)=
     &EEA(NVIR+MOD(ID,NVIR)+(IC-1)*NVIR+(IJ-1)*NVIR*NVIR,NVIR+MOD(IB,
     &NVIR)+(IA-1)*NVIR+(II-1)*NVIR*NVIR)+2*WBAR(II,IC+NOCC,IA+NOCC,IJ)
     &-WBAR(II,IC+NOCC,IJ,IA+NOCC)
                              ENDIF
                              IF (ID.EQ.IA) THEN
                                 EEA(NVIR+MOD(ID,NVIR)+(IC-1)*NVIR+(IJ-
     &1)*NVIR*NVIR,NVIR+MOD(IB,NVIR)+(IA-1)*NVIR+(II-1)*NVIR*NVIR)=EEA(
     &NVIR+MOD(ID,NVIR)+(IC-1)*NVIR+(IJ-1)*NVIR*NVIR,NVIR+MOD(IB,NVIR)+
     &(IA-1)*NVIR+(II-1)*NVIR*NVIR)-WBAR(II,IC+NOCC,IB+NOCC,IJ)
                              ENDIF
                              IF (IC.EQ.IA) THEN
                                 EEA(NVIR+MOD(ID,NVIR)+(IC-1)*NVIR+(IJ-
     &1)*NVIR*NVIR,NVIR+MOD(IB,NVIR)+(IA-1)*NVIR+(II-1)*NVIR*NVIR)=EEA(
     &NVIR+MOD(ID,NVIR)+(IC-1)*NVIR+(IJ-1)*NVIR*NVIR,NVIR+MOD(IB,NVIR)+
     &(IA-1)*NVIR+(II-1)*NVIR*NVIR)-WBAR(II,ID+NOCC,IJ,IB+NOCC)
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
