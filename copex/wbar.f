      SUBROUTINE WBAR(AUX,HBAR,NOCC,NVIR,NAO)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION AUX(NAO*NAO*NAO*NAO)
      DIMENSION HBAR(NAO,NAO,NAO,NAO)
C
      COMMON /LISTS/ MOIO(10,500),MOIOWD(10,500),MOIOSZ(10,500)
     &,MOIODS(10,500),MOIOFL(10,500)  
C
C      WRITE(6,*)
C      WRITE(6,*) 'READING WBAR'
C      WRITE(6,*)  
C
      DZERO = 0.D+0  
C
C      WRITE(6,*)
C      WRITE(6,*) 'WBAR OCC OCC OCC VIR (W_IJKA) ***********************'
C      WRITE(6,*)
C
C      WRITE(6,*) 'THE NUMBER OF ELEMENTS OF LIST 10:',MOIODS(1,10)
C      WRITE(6,*)
C
      CALL GETLST(AUX,1,MOIODS(1,10),1,1,10)      
      L=0
      DO 1 IA=1, NVIR 
         DO 2 IK=1, NOCC
            DO 3 IJ=1, NOCC
               DO 4 II=1, NOCC
                  L=L+1
                  NA=IA+NOCC
                  HBAR(II,IJ,IK,NA)=AUX(L)
C                  WRITE(6,*) '<',II,IJ,'|',IK,NA,'>=',HBAR(II,IJ,IK,NA)
                  HBAR(IJ,II,NA,IK)=AUX(L)
C                  WRITE(6,*) '<',IJ,II,'|',NA,IK,'>=',HBAR(IJ,II,NA,IK)
 4             CONTINUE
 3          CONTINUE
 2       CONTINUE
 1    CONTINUE
C
C      WRITE(6,*)
C      WRITE(6,*) 'WBAR OCC OCC OCC OCC (W_IJKL) ******************'
C      WRITE(6,*)
C
C      WRITE(6,*) 'THE NUMBER OF ELEMENTS OF LIST 53 IS:',MOIODS(1,53)
C      WRITE(6,*)
C
      CALL GETLST(AUX,1,MOIODS(1,53),1,1,53)  
      L=0
      DO 5 IJ=1, NOCC
         DO 6 II=1, NOCC
            DO 7 IN=1, NOCC
               DO 8 IM=1, NOCC
                  L=L+1
                  HBAR(IM,IN,II,IJ)=AUX(L)
C                  WRITE(6,*) '<',IM,IN,'|',II,IJ,'>=',HBAR(IM,IN,II,IJ)
 8             CONTINUE
 7          CONTINUE
 6       CONTINUE
 5    CONTINUE
C
C      WRITE(6,*)
C      WRITE(6,*) 'WBAR OCC VIR OCC VIR (W_MBJE) ***********************'
C      WRITE(6,*)
C
C      WRITE(6,*) 'THE NUMBER OF ELEMENTS OF LIST 58 IS:',MOIODS(1,58)
C      WRITE(6,*)
C
      CALL GETLST(AUX,1,MOIODS(1,58),1,1,58)   
      L=0
      DO 29 IJ=1, NOCC
         DO 30 IB=1, NVIR
            DO 31 IM=1, NOCC
               DO 32 IE=1, NVIR
                  L=L+1
                  NB=IB+NOCC
                  NE=IE+NOCC
                  HBAR(IM,NB,IJ,NE)=AUX(L)
C                  WRITE(6,*) '<',IM,NB,'|',IJ,NE,'>=',HBAR(IM,NB,IJ,NE)
 32            CONTINUE
 31         CONTINUE
 30      CONTINUE
 29    CONTINUE
C
C      WRITE(6,*)
C      WRITE(6,*) 'WBAR OCC VIR VIR OCC (W_MBEJ) ***********************'
C      WRITE(6,*)
C
C      WRITE(6,*) 'THE NUMBER OF ELEMENTS OF LIST 56 IS:',MOIODS(1,56)
C      WRITE(6,*)
C
      CALL GETLST(AUX,1,MOIODS(1,56),1,1,56)   
      L=0
      DO 9 IJ=1, NOCC
         DO 10 IB=1, NVIR
            DO 11 IM=1, NOCC
               DO 12 IE=1, NVIR
                  L=L+1
                  NB=IB+NOCC
                  NE=IE+NOCC
                  HBAR(IM,NB,NE,IJ)=AUX(L)
C                  WRITE(6,*) '<',IM,NB,'|',NE,IJ,'>=',HBAR(IM,NB,NE,IJ)
 12            CONTINUE
 11         CONTINUE
 10      CONTINUE
 9    CONTINUE
C
C      WRITE(6,*)
C      WRITE(6,*) 'WBAR OCC OCC VIR VIR (W_IJAB) ***********************'
C      WRITE(6,*)
C
C      WRITE(6,*) 'THE NUMBER OF ELEMENTS OF LIST 21 IS:',MOIODS(1,21)
C      WRITE(6,*)
C
C      CALL GETLST(AUX,1,MOIODS(1,21),1,1,21)
C      L=0
C      DO 33 II=1, NOCC
C         DO 34 IB=1, NVIR
C            DO 35 IJ=1, NOCC
C               DO 36 IA=1, NVIR
C                  L=L+1
C                  NA=IA+NOCC
C                  NB=IB+NOCC
C                  HBAR(II,IJ,NA,NB)=AUX(L)
C                  WRITE(6,*) '<',II,IJ,'|',NA,NB,'>=',HBAR(II,IJ,NA,NB)
C 36            CONTINUE
C 35         CONTINUE
C 34      CONTINUE
C 33   CONTINUE
C
      CALL GETLST(AUX,1,MOIODS(1,16),1,1,16)
      L=0
      DO 33 IJ=1, NOCC
         DO 34 II=1, NOCC
            DO 35 IB=1, NVIR
               DO 36 IA=1, NVIR
                  L=L+1
                  NA=IA+NOCC
                  NB=IB+NOCC
                  HBAR(II,IJ,NA,NB)=AUX(L)
C                  WRITE(6,*) '<',II,IJ,'|',NA,NB,'>=',HBAR(II,IJ,NA,NB)
 36            CONTINUE
 35         CONTINUE
 34      CONTINUE
 33   CONTINUE
C
C      WRITE(6,*)
C      WRITE(6,*) 'WBAR OCC VIR OCC OCC (W_KAIJ) **********************'
C      WRITE(6,*)
C
C      WRITE(6,*)'THE NUMBER OF ELEMENTS OF LIST 110 IS:',MOIODS(1,110)
C      WRITE(6,*)
C
      CALL GETLST(AUX,1,MOIODS(1,110),1,1,110)
      L=0
      DO 41 IA=1, NVIR
         DO 42 IK=1, NOCC
            DO 43 IJ=1, NOCC
               DO 44 II=1, NOCC
                  L=L+1
                  NA=IA+NOCC
                  HBAR(IK,NA,II,IJ)=AUX(L)
C                  WRITE(6,*) '<',IK,NA,'|',II,IJ,'>=',HBAR(IK,NA,II,IJ)
 44           CONTINUE
 43        CONTINUE
 42     CONTINUE
 41   CONTINUE
C
C      WRITE(6,*)
C      WRITE(6,*) 'WBAR VIR OCC VIR VIR (W_CIAB) **********************'
C      WRITE(6,*)
C
C      WRITE(6,*)'LIST 30 IS'
C      WRITE(6,*)
C
      CALL GETLST(AUX,1,MOIODS(1,30),1,1,30)
      L=0
      DO 49 II=1, NOCC
         DO 50 IC=1, NVIR
            DO 51 IB=1, NVIR
               DO 52 IA=1, NVIR
                  L=L+1
                  NC=IC+NOCC
                  NB=IB+NOCC
                  NA=IA+NOCC
                  HBAR(NC,II,NA,NB)=AUX(L)
C                  WRITE(6,*) '<',NC,II,'|',NA,NB,'>=',HBAR(NC,II,NA
C     &               ,NB)
 52           CONTINUE
 51        CONTINUE
 50     CONTINUE
 49   CONTINUE
C
C      WRITE(6,*)
C      WRITE(6,*) 'WBAR VIR VIR VIR OCC (W_ABCI) **********************'
C      WRITE(6,*)
C
C      WRITE(6,*)'LIST 130'
C      WRITE(6,*)
C
      CALL GETLST(AUX,1,MOIODS(1,130),1,1,130)
      L=0
      DO 57 II=1, NOCC
         DO 58 IC=1, NVIR
            DO 59 IB=1, NVIR
              DO 60 IA=1, NVIR
                  L=L+1
                  NC=IC+NOCC
                  NB=IB+NOCC
                  NA=IA+NOCC
                  HBAR(NA,NB,NC,II)=AUX(L)
C                  WRITE(6,*) '<',NA,NB,'|',NC,II,'>=',HBAR(NA,NB,NC,II)
 60           CONTINUE
 59        CONTINUE
 58     CONTINUE
 57   CONTINUE
C
C      WRITE(6,*)
C      WRITE(6,*) 'WBAR VIR VIR VIR VIR (W_ABCD) **********************'
C      WRITE(6,*)
C
C      WRITE(6,*)'THE NUMBER OF ELEMENTS OF LIST 233 IS:',MOIODS(1,233)
C      WRITE(6,*)
C
      CALL GETLST(AUX,1,MOIODS(1,233),1,1,233)
      L=0
      DO 61 ID=1, NVIR
         DO 62 IC=1, NVIR
            DO 63 IB=1, NVIR
               DO 64 IA=1, NVIR
                  L=L+1
                  NC=IC+NOCC
                  NB=IB+NOCC
                  NA=IA+NOCC
                  ND=ID+NOCC
                  HBAR(NA,NB,NC,ND)=AUX(L)
C                  WRITE(6,*) '<',NA,NB,'|',NC,ND,'>=',HBAR(NA,NB,NC,ND)
 64           CONTINUE
 63        CONTINUE
 62     CONTINUE
 61   CONTINUE
C
      RETURN
      END

      
