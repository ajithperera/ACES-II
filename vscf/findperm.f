      SUBROUTINE FINDPERM(W,PERMUT)
      IMPLICIT NONE
      DOUBLE PRECISION W,P1,P2,P3,P4,P5,P6,
     &                 MINDEV,DEV,XSIGN,YSIGN,ZSIGN
      INTEGER IPERM,I,J,JMAXX,JMAXY,JMAXZ
      LOGICAL PERMUT
      DIMENSION W(3,3),P1(3,3),P2(3,3),
     &          P3(3,3),P4(3,3),P5(3,3),P6(3,3)
C
      DATA P1 / 1.0D+00,0.0D+00,0.0D+00,0.0D+00,1.0D+00,0.0D+00,
     &          0.0D+00,0.0D+00,1.0D+00 /
      DATA P2 / 0.0D+00,1.0D+00,0.0D+00,1.0D+00,0.0D+00,0.0D+00,
     &          0.0D+00,0.0D+00,1.0D+00 /
      DATA P3 / 0.0D+00,0.0D+00,1.0D+00,0.0D+00,1.0D+00,0.0D+00,
     &          1.0D+00,0.0D+00,0.0D+00 /
      DATA P4 / 1.0D+00,0.0D+00,0.0D+00,0.0D+00,0.0D+00,1.0D+00,
     &          0.0D+00,1.0D+00,0.0D+00 /
      DATA P5 / 0.0D+00,1.0D+00,0.0D+00,0.0D+00,0.0D+00,1.0D+00,
     &          1.0D+00,0.0D+00,0.0D+00 /
      DATA P6 / 0.0D+00,0.0D+00,1.0D+00,1.0D+00,0.0D+00,0.0D+00,
     &          0.0D+00,1.0D+00,0.0D+00 /
C
      MINDEV = 100.0D+00
      IPERM  =  99
C
      DEV = 0.0D+00
      DO 10 J=1,3
      DO 10 I=1,3
      DEV = DEV + (DABS(W(I,J)) - P1(I,J)) **2
   10 CONTINUE
      IF(DEV.LT.MINDEV)THEN
       MINDEV = DEV
       IPERM  = 1
      ENDIF
C
      DEV = 0.0D+00
      DO 20 J=1,3
      DO 20 I=1,3
      DEV = DEV + (DABS(W(I,J)) - P2(I,J)) **2
   20 CONTINUE
      IF(DEV.LT.MINDEV)THEN
       MINDEV = DEV
       IPERM  = 2
      ENDIF
C
      DEV = 0.0D+00
      DO 30 J=1,3
      DO 30 I=1,3
      DEV = DEV + (DABS(W(I,J)) - P3(I,J)) **2
   30 CONTINUE
      IF(DEV.LT.MINDEV)THEN
       MINDEV = DEV
       IPERM  = 3
      ENDIF
C
      DEV = 0.0D+00
      DO 40 J=1,3
      DO 40 I=1,3
      DEV = DEV + (DABS(W(I,J)) - P4(I,J)) **2
   40 CONTINUE
      IF(DEV.LT.MINDEV)THEN
       MINDEV = DEV
       IPERM  = 4
      ENDIF
C
      DEV = 0.0D+00
      DO 50 J=1,3
      DO 50 I=1,3
      DEV = DEV + (DABS(W(I,J)) - P5(I,J)) **2
   50 CONTINUE
      IF(DEV.LT.MINDEV)THEN
       MINDEV = DEV
       IPERM  = 5
      ENDIF
C
      DEV = 0.0D+00
      DO 60 J=1,3
      DO 60 I=1,3
      DEV = DEV + (DABS(W(I,J)) - P6(I,J)) **2
   60 CONTINUE
      IF(DEV.LT.MINDEV)THEN
       MINDEV = DEV
       IPERM  = 6
      ENDIF
C
c      write(6,*) '  @FINDPERM-I, mindev ',mindev
c      write(6,*) '  @FINDPERM-I, iperm  ',iperm
C
      if(mindev.ge.100.0D+00) call errex
      if(iperm .gt.  6 .or. iperm .lt. 1) call errex
C
c      IF(MINDEV .GT. 0.1D+00)THEN
c       WRITE(6,*) '  @FINDPERM-I, No suitable permutation matrix. '
c       WRITE(6,*) '               Trying other method.            '
c       PERMUT = .FALSE.
c       RETURN
c      ELSE
c       WRITE(6,*) '  @FINDPERM-I, Suitable permutation found. Number ',
c     &            IPERM
c       PERMUT = .TRUE.
c      ENDIF
C
C     Find out signs. Look at rows of W to see if images of x,y, and z
C     are positive or negative.
C
      JMAXX=1
      JMAXY=1
      JMAXZ=1
      DO  70 J=1,3
       IF( DABS(W(1,J)) .GT. DABS(W(1,JMAXX)) ) JMAXX = J
       IF( DABS(W(2,J)) .GT. DABS(W(2,JMAXY)) ) JMAXY = J
       IF( DABS(W(3,J)) .GT. DABS(W(3,JMAXZ)) ) JMAXZ = J
   70 CONTINUE
C
      XSIGN = 1.0D+00
      YSIGN = 1.0D+00
      ZSIGN = 1.0D+00
      IF( W(1,JMAXX) .LT. 0) XSIGN = -1.0D+00
      IF( W(2,JMAXY) .LT. 0) YSIGN = -1.0D+00
      IF( W(3,JMAXZ) .LT. 0) ZSIGN = -1.0D+00
C
      WRITE(6,*) '  @FINDPERM-I, signs ', XSIGN,YSIGN,ZSIGN
C
      CALL ZERO(W,9)
      W(1,JMAXX) = XSIGN
      W(2,JMAXY) = YSIGN
      W(3,JMAXZ) = ZSIGN
C
      RETURN
      END
