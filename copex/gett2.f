      SUBROUTINE GETT2(TIJAB,AUX,NOCC,NVIR,NAO)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION TIJAB(NOCC,NOCC,NVIR,NVIR)
      DIMENSION AUX(NAO*NAO*NAO*NAO)
C
      COMMON /LISTS/ MOIO(10,500),MOIOWD(10,500),MOIOSZ(10,500)
     &,MOIODS(10,500),MOIOFL(10,500)    
C
C      WRITE(6,*)
C      WRITE(6,*) 'READING T2 AMPLITUDES'
C      WRITE(6,*)
C
      CALL GETLST(AUX,1,MOIODS(1,46),1,1,46)      
      L=0
      DO 1 IJ=1, NOCC
         DO 2 II=1, NOCC
            DO 3 IB=1, NVIR
               DO 4 IA=1, NVIR
                  L=L+1
                  TIJAB(II,IJ,IA,IB)=AUX(L)
C                  WRITE(6,*) 'TIJAB(',II,IJ,IA,IB,')=',
C     &TIJAB(II,IJ,IA,IB)
 4             CONTINUE
 3          CONTINUE
 2       CONTINUE
 1    CONTINUE
C
C      WRITE(6,*) 'T2 IN GETT2'
C      CALL OUTPUT(TIJAB,1,81,1,1,81,1,1)
      RETURN
      END
