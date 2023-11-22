      SUBROUTINE FBAR(FIJ,FIA,FAB,AUX,NOCC,NVIR,NAO)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION FIJ(NOCC,NOCC)
      DIMENSION FIA(NOCC,NVIR)
      DIMENSION FAB(NVIR,NVIR)
      DIMENSION AUX(NAO*NAO)
C
      COMMON /LISTS/ MOIO(10,500),MOIOWD(10,500),MOIOSZ(10,500)
     &,MOIODS(10,500),MOIOFL(10,500)
C
      CALL GETLST(AUX,1,1,1,1,91)
      K=0
      DO 3 I2=1,NOCC
         DO 4 I1=1,NOCC
            K=K+1
            FIJ(I1,I2)=AUX(K)
#ifdef _DEBUG_LVLM
            WRITE(6,"(a,F12.7") 'FIJ(',I1,I2,')=',FIJ(I1,I2)
#endif 
 4       CONTINUE
 3    CONTINUE 
C
      CALL GETLST(AUX,1,1,1,1,92)
      K=0
      DO 7 I1=1,NVIR
         DO 8 I2=1,NVIR
            K=K+1
            FAB(I1,I2)=AUX(K)
#ifdef _DEBUG_LVLM
            WRITE(6,"(a,F12.7)") 'FAB(',I1,I2,')=',FAB(I1,I2)
#endif 
 8       CONTINUE
 7    CONTINUE       
C
      CALL GETLST(AUX,1,1,1,1,93)
      K=0
      DO 11 I1=1,NOCC
         DO 12 I2=1,NVIR
            K=K+1
            FIA(I1,I2)=AUX(K)
#ifdef _DEBUG_LVLM
            WRITE(6,"(a,F12.7)") 'FIA(',I1,I2,')=',FIA(I1,I2)
#endif 
 12      CONTINUE
 11   CONTINUE       
C
      RETURN
      END
