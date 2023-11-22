      SUBROUTINE DIAG3E(EEA,WBAR,NVIR,NEA,NOCC,NAO)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z) 
C
      DIMENSION EEA(NEA,NEA)
      DIMENSION WBAR(NAO,NAO,NAO,NAO)
C
C EEA = SUM_C,D,L <AL|CD>
C
      DO 10 ID=1, NVIR
         DO 20 II=1, NOCC
            DO 30 IA=1, NVIR
               DO 40 IB=1, NVIR
                  IF (MOD(IB,NVIR).EQ.0) THEN
                     EEA(ID,NVIR+NVIR+(IA-1)*NVIR+(II-1)*NVIR*NVIR)
     &=2*WBAR(ID+NOCC,II,IB+NOCC,IA+NOCC)-WBAR(ID+NOCC,II,IA+NOCC,
     &IB+NOCC)  
                  ELSE
                     EEA(ID,NVIR+MOD(IB,NVIR)+(IA-1)*NVIR+(II-1)*NVIR
     &*NVIR)=2*WBAR(ID+NOCC,II,IB+NOCC,IA+NOCC)-WBAR(ID+NOCC,II,IA+NOCC,
     &IB+NOCC)  
                  ENDIF
 40            CONTINUE
 30         CONTINUE
 20      CONTINUE
 10   CONTINUE
C
      RETURN
      END      

