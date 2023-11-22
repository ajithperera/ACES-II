      SUBROUTINE DIAG4E(EEA,FIA,NVIR,NEA,NOCC,NAO)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z) 
C
      DIMENSION EEA(NEA,NEA)
      DIMENSION FIA(NOCC,NVIR)
C
C EEA = SUM_I,A FAI(I,A) 
C
      DO 20 ID=1, NVIR
         DO 10 II=1, NOCC
            DO 30 IA=1, NVIR
               DO 40 IB=1, NVIR
                  IF (ID.EQ.IB) THEN
                     IF (MOD(IB,NVIR).EQ.0) THEN
                        EEA(ID,NVIR+NVIR+(IA-1)*NVIR+
     &(II-1)*NVIR*NVIR)=EEA(ID,NVIR+NVIR+(IA-1)*NVIR+
     &(II-1)*NVIR*NVIR)+2*FIA(II,IA)
                     ELSE
                        EEA(ID,NVIR+MOD(IB,NVIR)+
     &(IA-1)*NVIR+(II-1)*NVIR*NVIR)=EEA(ID,NVIR+MOD(IB,NVIR)+
     &(IA-1)*NVIR+(II-1)*NVIR*NVIR)+2*FIA(II,IA)
                     ENDIF
                  ENDIF
                  IF (IA.EQ.ID) THEN
                     IF (MOD(IB,NVIR).EQ.0) THEN
                        EEA(ID,NVIR+NVIR+(IA-1)*NVIR+
     &(II-1)*NVIR*NVIR)=EEA(ID,NVIR+NVIR+(IA-1)*NVIR+
     &(II-1)*NVIR*NVIR)-FIA(II,IB)
                     ELSE
                        EEA(ID,NVIR+MOD(IB,NVIR)+
     &(IA-1)*NVIR+(II-1)*NVIR*NVIR)=EEA(ID,NVIR+MOD(IB,NVIR)+
     &(IA-1)*NVIR+(II-1)*NVIR*NVIR)-FIA(II,IB)
                     ENDIF
                  ENDIF
 40            CONTINUE
 30         CONTINUE
 10      CONTINUE
 20   CONTINUE
C
      RETURN
      END      

      
