      SUBROUTINE DRCCOF(EA,HPP,HPQ,WORKEA,EIVELEA,EIVEREA,EIVAEA,WIEA,
     &XQP,CP,CQ,CPIN,LWORKEA,NEA,NEAN,NEASH,NIP)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z) 
C
      DIMENSION EA(NEA,NEA)
      DIMENSION HPP(NEAN,NEAN)
      DIMENSION HPQ(NEAN,NEASH)
      DIMENSION WORKEA(LWORKEA)
      DIMENSION EIVELEA(NEA,NEA)
      DIMENSION EIVEREA(NEA,NEA)
      DIMENSION EIVAEA(NEA)
      DIMENSION WIEA(NEA)
      DIMENSION XQP(NEASH,NEAN)
      DIMENSION CP(NEAN,NEAN)
      DIMENSION CQ(NEASH,NEAN)
      DIMENSION CPIN(NEAN,NEAN)
C
      DONE = 1.0D+0
      DZERO = 0.D+0 
C
      DO 1 II=1, NEAN
         DO 2 IJ=1, NEAN
            HPP(II,IJ) = EA(II,IJ)
 2       CONTINUE
 1    CONTINUE 
C
      DO 3 II=1, NEAN
         DO 4 IJ=1, NEASH
            HPQ(II,IJ) = EA(II,IJ+NEAN)
 4       CONTINUE
 3    CONTINUE
C
C      WRITE(6,*) 
C      WRITE(6,*) 'FULL HBAR'
C      N=NEAN+NEASH
C      CALL OUTPUT(EA,1,N,1,N,N,N,1)
C
      CALL DGEEV('V','V',NEA,EA,NEA,EIVAEA,WIEA,EIVELEA,NEA,EIVEREA,NEA,
     &WORKEA,LWORKEA,INFO)
C
      IF (INFO.NE.0) THEN
         WRITE(6,*) 'ERROR IN DIAGONALIZATION'
         WRITE(6,*) 'INFO=',INFO
      ENDIF
C
C      WRITE(6,*)
C      WRITE(6,*) 'EIGENVALUES'
C      CALL OUTPUT(EIVAEA,1,NEA,1,1,NEA,1,1)
C
C      WRITE(6,*)
C      WRITE(6,*) 'RIGHT EIGENVECTORS'
C      CALL OUTPUT(EIVEREA,1,NEA,1,NEA,NEA,NEA,1)
C
      IF (NIP.EQ.0) THEN
         DO 5 I1=1, NEAN
            CP(I1,1)=EIVEREA(I1,126)
 5       CONTINUE
         DO 6 I1=1, NEAN
            CP(I1,2)=EIVEREA(I1,146)
 6       CONTINUE
         DO 7 I1=1, NEAN
            CP(I1,3)=EIVEREA(I1,147)
 7       CONTINUE
         DO 8 I1=1, NEAN
            CP(I1,4)=EIVEREA(I1,145)
 8       CONTINUE
         DO 9 I1=1, NEAN
            CP(I1,5)=EIVEREA(I1,231)
 9       CONTINUE
         DO 10 I1=1, NEAN
            CP(I1,6)=EIVEREA(I1,232)
 10      CONTINUE
         DO 11 I1=1, NEAN
            CP(I1,7)=EIVEREA(I1,233)
 11      CONTINUE
         DO 12 I1=1, NEAN
            CP(I1,8)=EIVEREA(I1,27)
 12      CONTINUE
C
         DO 15 I1=1, NEASH
            CQ(I1,1)=EIVEREA(I1+NEAN,126)
 15      CONTINUE
         DO 16 I1=1, NEASH
            CQ(I1,2)=EIVEREA(I1+NEAN,146)
 16      CONTINUE
         DO 17 I1=1, NEASH
            CQ(I1,3)=EIVEREA(I1+NEAN,147)
 17      CONTINUE
         DO 18 I1=1, NEASH
            CQ(I1,4)=EIVEREA(I1+NEAN,145)
 18      CONTINUE
         DO 19 I1=1, NEASH
            CQ(I1,5)=EIVEREA(I1+NEAN,231)
 19      CONTINUE
         DO 20 I1=1, NEASH
            CQ(I1,6)=EIVEREA(I1+NEAN,232)
 20      CONTINUE
         DO 21 I1=1, NEASH
            CQ(I1,7)=EIVEREA(I1+NEAN,233)
 21      CONTINUE
         DO 22 I1=1, NEASH
            CQ(I1,8)=EIVEREA(I1+NEAN,27)
 22      CONTINUE
      ELSE
C         write(6,*) 'in alternative branch'
         DO 23 I1=1, NEAN
            CP(I1,1)=EIVEREA(I1,98)
 23      CONTINUE
         DO 24 I1=1, NEAN
            CP(I1,2)=EIVEREA(I1,99)
 24      CONTINUE
         DO 25 I1=1, NEAN
            CP(I1,3)=EIVEREA(I1,100)
 25      CONTINUE
         DO 26 I1=1, NEAN
            CP(I1,4)=EIVEREA(I1,104)
 26      CONTINUE
         DO 27 I1=1, NEAN
            CP(I1,5)=EIVEREA(I1,35)
 27      CONTINUE
C
         DO 28 I1=1, NEASH
            CQ(I1,1)=EIVEREA(I1+NEAN,98)
 28      CONTINUE
         DO 29 I1=1, NEASH
            CQ(I1,2)=EIVEREA(I1+NEAN,99)
 29      CONTINUE
         DO 30 I1=1, NEASH
            CQ(I1,3)=EIVEREA(I1+NEAN,100)
 30      CONTINUE
         DO 31 I1=1, NEASH
            CQ(I1,4)=EIVEREA(I1+NEAN,104)
 31      CONTINUE
         DO 32 I1=1, NEASH
            CQ(I1,5)=EIVEREA(I1+NEAN,35)
 32      CONTINUE
      ENDIF     
C
C      WRITE(6,*)
C      WRITE(6,*) 'CP'
C      CALL OUTPUT(CP,1,NEAN,1,NEAN,NEAN,NEAN,1)
C
C      WRITE(6,*)
C      WRITE(6,*) 'CQ'
C      CALL OUTPUT(CQ,1,NEASH,1,NEAN,NEASH,NEAN,1)
C
      CALL MINV(CP,NEAN,NEAN,WORKEA,DET,0,0)      
C
      CALL XGEMM('N','N',NEASH,NEAN,NEAN,DONE,CQ,NEASH,CP,NEAN,
     &DZERO,XQP,NEASH)
C
C      WRITE(6,*)
C      WRITE(6,*) 'XQP'
C      CALL OUTPUT(XQP,1,NEASH,1,NEAN,NEASH,NEAN,1)
C
      RETURN
      END