      SUBROUTINE DIAG6E(EEA,WBAR,NVIR,NEA,NOCC,NAO)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z) 
C
      DIMENSION EEA(NEA,NEA)
      DIMENSION WBAR(NAO,NAO,NAO,NAO) 
C
C EEA = SUM_I,J <IJ|KL>
C
      DO 10 IJ=1, NOCC
         DO 20 IC=1, NVIR
            DO 30 ID=1, NVIR
               DO 40 II=1, NOCC
                  DO 50 IA=1, NVIR
                     DO 60 IB=1, NVIR
                        IF ((MOD(ID,NVIR).EQ.0).AND.(MOD(IB
     &,NVIR).EQ.0)) THEN
                           IF (II.EQ.IJ) THEN
                              EEA(NVIR+NVIR+(IC-1)*NVIR+(IJ-1)*NVIR*NVIR
     &,NVIR+NVIR+(IA-1)*NVIR+(II-1)*NVIR*NVIR)=EEA(NVIR+NVIR+(IC-1)*NVIR
     &+(IJ-1)*NVIR*NVIR,NVIR+NVIR+(IA-1)*NVIR+(II-1)*NVIR*NVIR)+
     &WBAR(ID+NOCC,IC+NOCC,IB+NOCC,IA+NOCC)
                           ENDIF
                        ELSE
                           IF ((MOD(ID,NVIR).EQ.0).OR.(MOD(IB
     &,NVIR).EQ.0)) THEN
                              IF (MOD(ID,NVIR).EQ.0) THEN
                                 IF (II.EQ.IJ) THEN
                                    EEA(NVIR+NVIR+(IC-1)*NVIR+(IJ-1)*
     &NVIR*NVIR,NVIR+MOD(IB,NVIR)+(IA-1)*NVIR+(II-1)*NVIR*NVIR)=EEA(NVIR
     &+NVIR+(IC-1)*NVIR+(IJ-1)*NVIR*NVIR,NVIR+MOD(IB,NVIR)+(IA-1)*NVIR+
     &(II-1)*NVIR*NVIR)+WBAR(ID+NOCC,IC+NOCC,IB+NOCC,IA+NOCC)
                                 ENDIF
                              ELSE
                                 IF (II.EQ.IJ) THEN
                                    EEA(NVIR+MOD(ID,NVIR)+(IC-1)*NVIR+
     &(IJ-1)*NVIR*NVIR,NVIR+NVIR+(IA-1)*NVIR+(II-1)*NVIR*NVIR)=EEA(NVIR+
     &MOD(ID,NVIR)+(IC-1)*NVIR+(IJ-1)*NVIR*NVIR,NVIR+NVIR+(IA-1)*NVIR+
     &(II-1)*NVIR*NVIR)+WBAR(ID+NOCC,IC+NOCC,IB+NOCC,IA+NOCC)
                                 ENDIF
                              ENDIF
                           ELSE
                              IF (II.EQ.IJ) THEN
                                 EEA(NVIR+MOD(ID,NVIR)+(IC-1)*NVIR+(IJ-1
     &)*NVIR*NVIR,NVIR+MOD(IB,NVIR)+(IA-1)*NVIR+(II-1)*NVIR*NVIR)=EEA(
     &NVIR+MOD(ID,NVIR)+(IC-1)*NVIR+(IJ-1)*NVIR*NVIR,NVIR+MOD(IB,NVIR)+
     &(IA-1)*NVIR+(II-1)*NVIR*NVIR)+WBAR(ID+NOCC,IC+NOCC,IB+NOCC,IA+
     &NOCC)
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
