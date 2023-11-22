      SUBROUTINE DRCCOF(EA,HPP,HPQ,WORKEA,EIVELEA,EIVEREA,EIVAEA,WIEA,
     &XQP,CP,CQ,LWORKEA,NEA,NEAN,NEASH,NIP)
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
C
C HE AUG-CC-PVTZ
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
      CALL DGEEV('V','V',NEA,EA,NEA,EIVAEA,WIEA,EIVELEA,NEA,EIVEREA,NEA,
     &WORKEA,LWORKEA,INFO)
C
      IF (INFO.NE.0) THEN
         WRITE(6,*) 'ERROR IN DIAGONALIZATION'
         WRITE(6,*) 'INFO=',INFO
      ENDIF
C
      IF (NIP.EQ.0) THEN
         DO 5 I1=1, NEAN
            CP(I1,1)=EIVEREA(I1,36)
 5       CONTINUE
         DO 6 I1=1, NEAN
            CP(I1,2)=EIVEREA(I1,63)
C            WRITE(6,*) CP(I1,2)
 6       CONTINUE
         DO 7 I1=1, NEAN
            CP(I1,3)=EIVEREA(I1,66)
 7       CONTINUE
         DO 8 I1=1, NEAN
            CP(I1,4)=EIVEREA(I1,67)
 8       CONTINUE
         DO 9 I1=1, NEAN
            CP(I1,5)=EIVEREA(I1,62)
 9       CONTINUE
         DO 10 I1=1, NEAN
            CP(I1,6)=EIVEREA(I1,154)
 10      CONTINUE
         DO 11 I1=1, NEAN
            CP(I1,7)=EIVEREA(I1,155)
 11      CONTINUE
         DO 12 I1=1, NEAN
            CP(I1,8)=EIVEREA(I1,158)
 12      CONTINUE
         DO 13 I1=1, NEAN
            CP(I1,9)=EIVEREA(I1,159)
 13      CONTINUE
         DO 14 I1=1, NEAN
            CP(I1,10)=EIVEREA(I1,164)
 14      CONTINUE
         DO 15 I1=1, NEAN
            CP(I1,11)=EIVEREA(I1,141)
 15      CONTINUE
         DO 16 I1=1, NEAN
            CP(I1,12)=EIVEREA(I1,165)
 16      CONTINUE
         DO 17 I1=1, NEAN
            CP(I1,13)=EIVEREA(I1,166)
 17      CONTINUE
         DO 18 I1=1, NEAN
            CP(I1,14)=EIVEREA(I1,167)
 18      CONTINUE
         DO 19 I1=1, NEAN
            CP(I1,15)=EIVEREA(I1,170)
 19      CONTINUE
         DO 20 I1=1, NEAN
            CP(I1,16)=EIVEREA(I1,351)
 20      CONTINUE
         DO 21 I1=1, NEAN
            CP(I1,17)=EIVEREA(I1,352)
 21      CONTINUE
         DO 41 I1=1, NEAN
            CP(I1,18)=EIVEREA(I1,357)
 41      CONTINUE
         DO 42 I1=1, NEAN
            CP(I1,19)=EIVEREA(I1,358)
 42      CONTINUE
         DO 43 I1=1, NEAN
            CP(I1,20)=EIVEREA(I1,359)
 43      CONTINUE
         DO 44 I1=1, NEAN
            CP(I1,21)=EIVEREA(I1,466)
 44      CONTINUE
         DO 45 I1=1, NEAN
            CP(I1,22)=EIVEREA(I1,467)
 45      CONTINUE
         DO 46 I1=1, NEAN
            CP(I1,23)=EIVEREA(I1,468)
 46      CONTINUE
         DO 47 I1=1, NEAN
            CP(I1,24)=EIVEREA(I1,347)
 47      CONTINUE
C
         DO 22 I1=1, NEASH
            CQ(I1,1)=EIVEREA(I1+NEAN,36)
 22      CONTINUE
         DO 23 I1=1, NEASH
            CQ(I1,2)=EIVEREA(I1+NEAN,63)
 23      CONTINUE
         DO 24 I1=1, NEASH
            CQ(I1,3)=EIVEREA(I1+NEAN,66)
 24      CONTINUE
         DO 25 I1=1, NEASH
            CQ(I1,4)=EIVEREA(I1+NEAN,67)
 25      CONTINUE
         DO 26 I1=1, NEASH
            CQ(I1,5)=EIVEREA(I1+NEAN,62)
 26      CONTINUE
         DO 27 I1=1, NEASH
            CQ(I1,6)=EIVEREA(I1+NEAN,154)
 27      CONTINUE         
         DO 28 I1=1, NEASH
            CQ(I1,7)=EIVEREA(I1+NEAN,155)
 28      CONTINUE
         DO 29 I1=1, NEASH
            CQ(I1,8)=EIVEREA(I1+NEAN,158)
 29      CONTINUE
         DO 30 I1=1, NEASH
            CQ(I1,9)=EIVEREA(I1+NEAN,159)
 30      CONTINUE
         DO 31 I1=1, NEASH
            CQ(I1,10)=EIVEREA(I1+NEAN,164)
 31      CONTINUE
         DO 32 I1=1, NEASH
            CQ(I1,11)=EIVEREA(I1+NEAN,141)
 32      CONTINUE
         DO 33 I1=1, NEASH
            CQ(I1,12)=EIVEREA(I1+NEAN,165)
 33      CONTINUE
         DO 34 I1=1, NEASH
            CQ(I1,13)=EIVEREA(I1+NEAN,166)
 34      CONTINUE
         DO 35 I1=1, NEASH
            CQ(I1,14)=EIVEREA(I1+NEAN,167)
 35      CONTINUE         
         DO 36 I1=1, NEASH
            CQ(I1,15)=EIVEREA(I1+NEAN,170)
 36      CONTINUE
         DO 37 I1=1, NEASH
            CQ(I1,16)=EIVEREA(I1+NEAN,351)
 37      CONTINUE
         DO 38 I1=1, NEASH
            CQ(I1,17)=EIVEREA(I1+NEAN,352)
 38      CONTINUE
         DO 48 I1=1, NEASH
            CQ(I1,18)=EIVEREA(I1+NEAN,357)
 48      CONTINUE
         DO 49 I1=1, NEASH
            CQ(I1,19)=EIVEREA(I1+NEAN,358)
 49      CONTINUE
         DO 50 I1=1, NEASH
            CQ(I1,20)=EIVEREA(I1+NEAN,359)
 50      CONTINUE
         DO 51 I1=1, NEASH
            CQ(I1,21)=EIVEREA(I1+NEAN,466)
 51      CONTINUE         
         DO 52 I1=1, NEASH
            CQ(I1,22)=EIVEREA(I1+NEAN,467)
 52      CONTINUE
         DO 53 I1=1, NEASH
            CQ(I1,23)=EIVEREA(I1+NEAN,468)
 53      CONTINUE
         DO 54 I1=1, NEASH
            CQ(I1,24)=EIVEREA(I1+NEAN,347)
 54      CONTINUE
      ELSE
C         write(6,*) 'in alternative branch'
         DO 39 I1=1, NEAN
            CP(I1,1)=EIVEREA(I1,3)
 39      CONTINUE
C
         DO 40 I1=1, NEASH
            CQ(I1,1)=EIVEREA(I1+NEAN,3)
 40      CONTINUE
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
      RETURN
      END
