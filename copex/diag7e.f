      SUBROUTINE DIAG7E(EEA,WBAR,T2,NVIR,NEA,NOCC,NAO)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z) 
C
      DIMENSION EEA(NEA,NEA)
      DIMENSION WBAR(NAO,NAO,NAO,NAO)
      DIMENSION T2(NOCC,NOCC,NVIR,NVIR)
C
C EEA = SUM_I,J,A,C (-2<IJ|CA> + <JI|CA>) * T(BC,LK)
C
      DO 10 IJ=1, NOCC
         DO 20 IC=1, NVIR
            DO 30 ID=1, NVIR
               DO 40 II=1, NOCC
                  DO 50 IA=1, NVIR
                     DO 60 IB=1, NVIR
                        IF ((MOD(ID,NVIR).EQ.0).AND.(MOD(IB
     &,NVIR).EQ.0)) THEN
                           DO 70 IK=1, NOCC
                              EEA(NVIR+NVIR+(IC-1)*NVIR+(IJ-1)*NVIR*NVIR
     &,NVIR+NVIR+(IA-1)*NVIR+(II-1)*NVIR*NVIR)=EEA(NVIR+NVIR+(IC-1)*NVIR
     &+(IJ-1)*NVIR*NVIR,NVIR+NVIR+(IA-1)*NVIR+(II-1)*NVIR*NVIR)-
     &(2*WBAR(IK,II,IB+NOCC,IA+NOCC)-WBAR(IK,II,IA+NOCC,IB+NOCC))*
     &T2(IK,IJ,ID,IC)
 70                        CONTINUE
                        ELSE
                           IF ((MOD(ID,NVIR).EQ.0).OR.(MOD(IB
     &,NVIR).EQ.0)) THEN
                              IF (MOD(ID,NVIR).EQ.0) THEN
                                 DO 90 IK=1, NOCC
                                    EEA(NVIR+NVIR+(IC-1)*NVIR+(IJ-1)*
     &NVIR*NVIR,NVIR+MOD(IB,NVIR)+(IA-1)*NVIR+(II-1)*NVIR*NVIR)=EEA(NVIR
     &+NVIR+(IC-1)*NVIR+(IJ-1)*NVIR*NVIR,NVIR+MOD(IB,NVIR)+(IA-1)*NVIR+
     &(II-1)*NVIR*NVIR)-(2*WBAR(IK,II,IB+NOCC,IA+NOCC)-WBAR(IK,II,IA+
     &NOCC,IB+NOCC))*T2(IK,IJ,ID,IC)
 90                              CONTINUE
                              ELSE
                                 DO 11 IK=1, NOCC
                                    EEA(NVIR+MOD(ID,NVIR)+(IC-1)*NVIR+
     &(IJ-1)*NVIR*NVIR,NVIR+NVIR+(IA-1)*NVIR+(II-1)*NVIR*NVIR)=EEA(NVIR+
     &MOD(ID,NVIR)+(IC-1)*NVIR+(IJ-1)*NVIR*NVIR,NVIR+NVIR+(IA-1)*NVIR+
     &(II-1)*NVIR*NVIR)-(2*WBAR(IK,II,IB+NOCC,IA+NOCC)-WBAR(IK,II,IA+
     &NOCC,IB+NOCC))*T2(IK,IJ,ID,IC)
 11                              CONTINUE
                              ENDIF
                           ELSE
                              DO 12 IK=1, NOCC
                                 EEA(NVIR+MOD(ID,NVIR)+(IC-1)*NVIR+(IJ-1
     &)*NVIR*NVIR,NVIR+MOD(IB,NVIR)+(IA-1)*NVIR+(II-1)*NVIR*NVIR)=EEA(
     &NVIR+MOD(ID,NVIR)+(IC-1)*NVIR+(IJ-1)*NVIR*NVIR,NVIR+MOD(IB,NVIR)+
     &(IA-1)*NVIR+(II-1)*NVIR*NVIR)-(2*WBAR(IK,II,IB+NOCC,IA+NOCC)-WBAR
     &(IK,II,IA+NOCC,IB+NOCC))*T2(IK,IJ,ID,IC)
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
