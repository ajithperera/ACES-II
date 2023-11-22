      subroutine stofit(J,NP,ZA,NGAUSS,lone,COUT,
     $EOUT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION EXX(NGAUSS),C(NGAUSS)
     $,EOUT(NGAUSS),COUT(NGAUSS)
       pi=3.1415926535898D0
       pi3=31.00627666d0
       
C============================================================
C NP     =PRINCIPAL QUANTUM NUMBER
C NGAUSS =NUMBER OF GAUSSIAN FUNCTIONS STO-NG
C ZA     =STO EXPONENT
C EXX    =PRIMITIVE EXPONENT SCALED AND RETURNED
C C      =NORMALIZED AND SCALED RETURNED COEFS
C J      =ANGULAR MOMENTUM (1=S),(2=P),(3=D)
C============================================================
      IF(NP.EQ.1)THEN
         CALL S1S(EXX,C,NGAUSS)
      ELSE IF(NP.EQ.2)THEN
         IF(J.EQ.1)CALL S2S(EXX,C,NGAUSS)
         IF(J.EQ.2)CALL S2P(EXX,C,NGAUSS)
      ELSE IF(NP.EQ.3)THEN
         IF(J.EQ.1)CALL S3S(EXX,C,NGAUSS)
         IF(J.EQ.2)CALL S3P(EXX,C,NGAUSS)
         IF(J.EQ.3)CALL S3D(EXX,C,NGAUSS)
      ELSE IF(NP.EQ.4)THEN
         IF(J.EQ.1)CALL S4S(EXX,C,NGAUSS)
         IF(J.EQ.2)CALL S4P(EXX,C,NGAUSS)
         IF(J.EQ.3)CALL S4D(EXX,C,NGAUSS)
      ELSE IF(NP.EQ.5)THEN
         IF(J.EQ.1)CALL S5S(EXX,C,NGAUSS)
         IF(J.EQ.2)CALL S5P(EXX,C,NGAUSS)
         IF(J.EQ.3)CALL S5D(EXX,C,NGAUSS)
      ENDIF
C SCALE EXPONENT, NORMALIZE, AND CREATE COEFS.
      ZEE=ZA*ZA
!      XPI=(0.25D0/1.570796327D0)**0.75D0
      DO 100 I=1,NGAUSS
         EXX(I)=EXX(I)*ZEE
         EOUT(I)=EXX(I)
 100      CONTINUE
      SUM=0.0D0
c      XJA=DBLE(J)
       index=(2*j)+1

c      DO 101 I=1,NGAUSS
c         DO 102 JJ=1,I
c            T=C(I)*C(JJ)*
c     X   (2.0D0*DSQRT(EXX(I)*EXX(JJ))/(EXX(I)+EXX(JJ)))**(XJA+0.5D0)
c            SUM=SUM+T
c            IF(I.NE.JJ)THEN
c               SUM=SUM+T
c            ENDIF
c 102     CONTINUE
c 101  CONTINUE

c     a modification to this loop b/c as its written above, by ZERNER, slows the code down
c IMMENSELY due to the call to pow()
            
            DO 101 I=1,NGAUSS
         DO 102 JJ=1,I
            rfac=1.0d0
            T=C(I)*C(JJ)
            factor=2.0D0*DSQRT(EXX(I)*EXX(JJ))/(EXX(I)+EXX(JJ))
             do 800 kk=1,index
             rfac=rfac*factor
 800         continue
              T=T*dsqrt(rfac)
            SUM=SUM+T
            IF(I.NE.JJ)THEN
               SUM=SUM+T
            ENDIF
 102     CONTINUE
 101  CONTINUE




     
      XN=1.0D0/SQRT(SUM)
      if((lone.gt.4).and.(lone.lt.8))then
         do 104 I=1,NGAUSS
        C(I)=XN*C(I)*((2048.0D0*EXX(I)**7)/(9.0D0*pi**3))**(.25)
         COUT(I)=C(I)
 104     continue
         go to 10
         end if

      if(lone.eq.1)then
          do 105 I=1,NGAUSS
c       C(I)=XN*C(I)*((2.0D0*EXX(I))/pi)**(.75)
           factor=((2.0D0*EXX(I))/pi)
           factor=factor*factor*factor
           factor=sqrt(factor)
           factor=sqrt(factor)
        C(I)=XN*C(I)*factor
         COUT(I)=C(I)
 105  continue
      go to 10
      end if
      if((lone.gt.1).and.(lone.lt.5))then
         do 106 I=1,NGAUSS
c       C(I)=XN*C(I)*((128.0D0*EXX(I)**5)/pi**3)**(.25)
            factor=128.0D0*EXX(I)*EXX(I)*EXX(I)*EXX(I)*EXX(I)/pi3
            factor=sqrt(factor)
            factor=sqrt(factor)


            C(I)=XN*C(I)*factor

         COUT(I)=C(I)
 106  continue
      go to 10
      end if
      if((lone.gt.7).and.(lone.lt.11))then
         do 107 I=1,NGAUSS
       C(I)=XN*C(I)*((2048.0D0*EXX(I)**7)/pi**3)**(.25)
         COUT(I)=C(I)
 107  continue
      go to 10
      end if         
 10   return
      END
C---------------------------------------------------1S------------
      SUBROUTINE S1S(EXX,CS,NGAUSS)
C     COPYRIGHT  UF 1993  DAVID BAKER AND MICHAEL ZERNER
C*MODULE BASSTO  *DECK S1S
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C        *****  LEAST SQUARES FITS TO SLATER 1S FUNCTIONS       *****
C        *****  BY GAUSSIAN 1S FUNCTIONS                        *****
      DIMENSION EXX(6),CS(6)
      GO TO (100,120,140,160,180,200),NGAUSS
  100 CONTINUE
C        *****  STO(1S---1(1S))                                 *****
C        *****  ER1 4.3191D-02                                  *****
      EXX(1) = 2.709498091D-01
      CS(1) = 1.000000000D+00
      RETURN
  120 CONTINUE
C        *****  STO(1S---2(1S))                                 *****
C        *****  ER1 3.1606D-03                                  *****
      EXX(1) = 8.518186635D-01
      CS(1) = 4.301284983D-01
      EXX(2) = 1.516232927D-01
      CS(2) = 6.789135305D-01
      RETURN
  140 CONTINUE
C        *****  STO(1S---3(1S))                                 *****
C        *****  ER1 3.3053D-04                                  *****
      EXX(1) = 2.227660584D+00
      CS(1) = 1.543289673D-01
      EXX(2) = 4.057711562D-01
      CS(2) = 5.353281423D-01
      EXX(3) = 1.098175104D-01
      CS(3) = 4.446345422D-01
      RETURN
  160 CONTINUE
C        *****  STO(1S---4(1S))                                 *****
C        *****  ER1 4.3763D-05                                  *****
      EXX(1) = 5.216844534D+00
      CS(1) = 5.675242080D-02
      EXX(2) = 9.546182760D-01
      CS(2) = 2.601413550D-01
      EXX(3) = 2.652034102D-01
      CS(3) = 5.328461143D-01
      EXX(4) = 8.801862774D-02
      CS(4) = 2.916254405D-01
      RETURN
  180 CONTINUE
C        *****  STO(1S---5(1S))                                 *****
C        *****  ER1 6.8840D-06                                  *****
      EXX(1) = 1.130563696D+01
      CS(1) = 2.214055312D-02
      EXX(2) = 2.071728178D+00
      CS(2) = 1.135411520D-01
      EXX(3) = 5.786484833D-01
      CS(3) = 3.318161484D-01
      EXX(4) = 1.975724573D-01
      CS(4) = 4.825700713D-01
      EXX(5) = 7.445271746D-02
      CS(5) = 1.935721966D-01
      RETURN
  200 CONTINUE
C        *****  STO(1S---6(1S))                                 *****
C        *****  ER1 1.2372D-06                                  *****
      EXX(1) = 2.310303149D+01
      CS(1) = 9.163596280D-03
      EXX(2) = 4.235915534D+00
      CS(2) = 4.936149294D-02
      EXX(3) = 1.185056519D+00
      CS(3) = 1.685383049D-01
      EXX(4) = 4.070988982D-01
      CS(4) = 3.705627997D-01
      EXX(5) = 1.580884151D-01
      CS(5) = 4.164915298D-01
      EXX(6) = 6.510953954D-02
      CS(6) = 1.303340841D-01
      RETURN
      END
C---------------------------------------------------2P------------
C*MODULE BASSTO  *DECK S2P
      SUBROUTINE S2P(EXX,CP,NGAUSS)
C     COPYRIGHT  UF 1993  DAVID BAKER AND MICHAEL ZERNER
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C        *****  LEAST SQUARES FITS TO SLATER 2P FUNCTIONS       *****
C        *****  BY GAUSSIAN 2P FUNCTIONS                        *****
      DIMENSION EXX(6),CP(6)
      GO TO (100,120,140,160,180,200),NGAUSS
  100 CONTINUE
C        *****  STO(2P---1(2P))                                 *****
C        *****  ER1 4.8232D-02                                  *****
      EXX(1) = 1.759666885D-01
      CP(1) = 1.000000000D+00
      RETURN
  120 CONTINUE
C        *****  STO(2P---2(2P))                                 *****
C        *****  ER1 3.0947D-03                                  *****
      EXX(1) = 4.323908358D-01
      CP(1) = 4.522627513D-01
      EXX(2) = 1.069439065D-01
      CP(2) = 6.713122642D-01
      RETURN
  140 CONTINUE
C        *****  STO(2P---3(2P))                                 *****
C        *****  ER1 2.6858D-04                                  *****
      EXX(1) = 9.192379002D-01
      CP(1) = 1.623948553D-01
      EXX(2) = 2.359194503D-01
      CP(2) = 5.661708862D-01
      EXX(3) = 8.009805746D-02
      CP(3) = 4.223071752D-01
      RETURN
  160 CONTINUE
C        *****  STO(2P---4(2P))                                 *****
C        *****  ER1 2.9037D-05                                  *****
      EXX(1) = 1.798260992D+00
      CP(1) = 5.713170255D-02
      EXX(2) = 4.662622228D-01
      CP(2) = 2.857455515D-01
      EXX(3) = 1.643718620D-01
      CP(3) = 5.517873105D-01
      EXX(4) = 6.543927065D-02
      CP(4) = 2.632314924D-01
      RETURN
  180 CONTINUE
C        *****  STO(2P---5(2P))                                 *****
C        *****  ER1 3.7171D-06                                  *****
      EXX(1) = 3.320386533D+00
      CP(1) = 2.079051117D-02
      EXX(2) = 8.643257633D-01
      CP(2) = 1.235472099D-01
      EXX(3) = 3.079819284D-01
      CP(3) = 3.667738886D-01
      EXX(4) = 1.273309895D-01
      CP(4) = 4.834930290D-01
      EXX(5) = 5.606243164D-02
      CP(5) = 1.653444074D-01
      RETURN
  200 CONTINUE
C        *****  STO(2P---6(2P))                                 *****
C        *****  ER1 5.4444D-07                                  *****
      EXX(1) = 5.868285913D+00
      CP(1) = 7.924233646D-03
      EXX(2) = 1.530329631D+00
      CP(2) = 5.144104825D-02
      EXX(3) = 5.475665231D-01
      CP(3) = 1.898400060D-01
      EXX(4) = 2.288932733D-01
      CP(4) = 4.049863191D-01
      EXX(5) = 1.046655969D-01
      CP(5) = 4.012362861D-01
      EXX(6) = 4.948220127D-02
      CP(6) = 1.051855189D-01
      RETURN
      END
C---------------------------------------------------2S------------
C*MODULE BASSTO  *DECK S2S
      SUBROUTINE S2S(EXX,CS,NGAUSS)
C     COPYRIGHT  UF 1993  DAVID BAKER AND MICHAEL ZERNER
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C        *****  LEAST SQUARES FITS TO SLATER 2S FUNCTIONS       *****
C        *****  BY GUASSIAN 1S FUNCTIONS                        *****
      DIMENSION EXX(6),CS(6)
      GO TO (100,120,140,160,180,200),NGAUSS
  100 CONTINUE
C        *****  STO(2S---1(1S))                                 *****
C        *****  ER1 6.2923D-03                                  *****
      EXX(1) = 1.012151084D-01
      CS(1) = 1.000000000D+00
      RETURN
  120 CONTINUE
C        *****  STO(2S---2(1S))                                 *****
C        *****  ER1 2.6728D-03                                  *****
      EXX(1) = 1.292278611D-01
      CS(1) = 7.470867124D-01
      EXX(2) = 4.908584205D-02
      CS(2) = 2.855980556D-01
      RETURN
  140 CONTINUE
C        *****  STO(2S---3(1S))                                 *****
C        *****  ER1 6.8661D-05                                  *****
      EXX(1) = 2.581578398D+00
      CS(1) = -5.994474934D-02
      EXX(2) = 1.567622104D-01
      CS(2) = 5.960385398D-01
      EXX(3) = 6.018332272D-02
      CS(3) = 4.581786291D-01
      RETURN
  160 CONTINUE
C        *****  STO(2S---4(1S))                                 *****
C        *****  ER1 2.6974D-05                                  *****
      EXX(1) = 1.161525551D+01
      CS(1) = -1.198411747D-02
      EXX(2) = 2.000243111D+00
      CS(2) = -5.472052539D-02
      EXX(3) = 1.607280687D-01
      CS(3) = 5.805587176D-01
      EXX(4) = 6.125744532D-02
      CS(4) = 4.770079976D-01
      RETURN
  180 CONTINUE
C        *****  STO(2S---5(1S))                                 *****
C        *****  ER1 2.7009D-06                                  *****
      EXX(1) = 8.984956862D+00
      CS(1) = -1.596349096D-02
      EXX(2) = 1.673710636D+00
      CS(2) = -5.685884883D-02
      EXX(3) = 1.944726668D-01
      CS(3) = 3.698265599D-01
      EXX(4) = 8.806345634D-02
      CS(4) = 5.480512593D-01
      EXX(5) = 4.249068522D-02
      CS(5) = 1.472634893D-01
      RETURN
  200 CONTINUE
C        *****  STO(2S---6(1S))                                 *****
C        *****  ER1 3.2726D-07                                  *****
      EXX(1) = 2.768496241D+01
      CS(1) = -4.151277819D-03
      EXX(2) = 5.077140627D+00
      CS(2) = -2.067024148D-02
      EXX(3) = 1.426786050D+00
      CS(3) = -5.150303337D-02
      EXX(4) = 2.040335729D-01
      CS(4) = 3.346271174D-01
      EXX(5) = 9.260298399D-02
      CS(5) = 5.621061301D-01
      EXX(6) = 4.416183978D-02
      CS(6) = 1.712994697D-01
      RETURN
      END
C---------------------------------------------------3D------------
C*MODULE BASSTO  *DECK S3D
      SUBROUTINE S3D(EXX,CS,NGAUSS)
C     COPYRIGHT  UF 1993  DAVID BAKER AND MICHAEL ZERNER
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION EXX(6),CS(6)
      GO TO (100,120,140,160,180,200),NGAUSS
  100 CONTINUE
C*****************      STO(3D---1(3D))
C*********           ER1 5.0728D-02
      EXX(1) = 1.302270363D-01
      CS(1) = 1.000000000D+00
      RETURN
  120 CONTINUE
C*****************      STO(3D---2(3D))
C*********           ER1 2.9834D-03
      EXX(1) = 2.777427345D-01
      CS(1) = 4.666137923D-01
      EXX(2) = 8.336507714D-02
      CS(2) = 6.644706516D-01
      RETURN
  140 CONTINUE
C*****************      STO(3D---3(3D))
C*********           ER1 2.2687D-04
      EXX(1) = 5.229112225D-01
      CS(1) = 1.686596060D-01
      EXX(2) = 1.639595876D-01
      CS(2) = 5.847984817D-01
      EXX(3) = 6.386630021D-02
      CS(3) = 4.056779523D-01
      RETURN
  160 CONTINUE
C*****************      STO(3D---4(3D))
C*********           ER1 2.1116D-05
      EXX(1) = 9.185846715D-01
      CS(1) = 5.799057705D-02
      EXX(2) = 2.920461109D-01
      CS(2) = 3.045581349D-01
      EXX(3) = 1.187568890D-01
      CS(3) = 5.601358038D-01
      EXX(4) = 5.286755896D-02
      CS(4) = 2.432423313D-01
      RETURN
  180 CONTINUE
C*****************      STO(3D---5(3D))
C*********           ER1 2.3115D-06
      EXX(1) = 1.539033958D+00
      CS(1) = 2.020869128D-02
      EXX(2) = 4.922090297D-01
      CS(2) = 1.321157923D-01
      EXX(3) = 2.029756928D-01
      CS(3) = 3.911240346D-01
      EXX(4) = 9.424112917D-02
      CS(4) = 4.779609701D-01
      EXX(5) = 4.569058269D-02
      CS(5) = 1.463662294D-01
      RETURN
  200 CONTINUE
C*****************      STO(3D---6(3D))
C*********           ER1 2.8899D-07
      EXX(1) = 2.488296923D+00
      CS(1) = 7.283828112D-03
      EXX(2) = 7.981487853D-01
      CS(2) = 5.386799363D-02
      EXX(3) = 3.311327490D-01
      CS(3) = 2.072139149D-01
      EXX(4) = 1.559114463D-01
      CS(4) = 4.266269092D-01
      EXX(5) = 7.877734732D-02
      CS(5) = 3.843100204D-01
      EXX(6) = 4.058484363D-02
      CS(6) = 8.902827546D-02
      RETURN
      END
C---------------------------------------------------3P------------
C*MODULE BASSTO  *DECK S3P
      SUBROUTINE S3P(EXX,CP,NGAUSS)
C     COPYRIGHT  UF 1993  DAVID BAKER AND MICHAEL ZERNER
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C        *****  LEAST SQUARES FITS TO SLATER 3P FUNCTIONS       *****
C        *****  WITH GAUSSIAN 2P FUNCTIONS                      *****
      DIMENSION EXX(6),CP(6)
      GO TO (100,120,140,160,180,200),NGAUSS
  100 CONTINUE
C        *****  STO(3P---1(2P))                                 *****
C        *****  ER1 1.2745D-02                                  *****
      EXX(1) = 9.113614253D-02
      CP(1) = 1.000000000D+00
      RETURN
  120 CONTINUE
C        *****  STO(3P---2(2P))                                 *****
C        *****  ER1 2.3771D-04                                  *****
      EXX(1) = 1.458620964D-01
      CP(1) = 5.349653144D-01
      EXX(2) = 5.664210742D-02
      CP(2) = 5.299607212D-01
      RETURN
  140 CONTINUE
C        *****  STO(3P---3(2P))                                 *****
C        *****  ER1 1.3487D-04                                  *****
      EXX(1) = 2.692880368D+00
      CP(1) = -1.061945788D-02
      EXX(2) = 1.489359592D-01
      CP(2) = 5.218564264D-01
      EXX(3) = 5.739585040D-02
      CP(3) = 5.450015143D-01
      RETURN
  160 CONTINUE
C        *****  STO(3P---4(2P))                                 *****
C        *****  ER1 2.9785D-06                                  *****
      EXX(1) = 1.853180239D+00
      CP(1) = -1.434249391D-02
      EXX(2) = 1.915075719D-01
      CP(2) = 2.755177589D-01
      EXX(3) = 8.655487938D-02
      CP(3) = 5.846750879D-01
      EXX(4) = 4.184253862D-02
      CP(4) = 2.144986514D-01
      RETURN
  180 CONTINUE
C        *****  STO(3P---5(2P))                                 *****
C        *****  ER1 1.3387D-06                                  *****
      EXX(1) = 6.466803859D+00
      CP(1) = -2.329023747D-03
      EXX(2) = 1.555914802D+00
      CP(2) = -1.357395221D-02
      EXX(3) = 1.955925255D-01
      CP(3) = 2.632185383D-01
      EXX(4) = 8.809647701D-02
      CP(4) = 5.880427024D-01
      EXX(5) = 4.234835707D-02
      CP(5) = 2.242794445D-01
      RETURN
  200 CONTINUE
C        *****  STO(3P---6(2P))                                 *****
C        *****  ER1 7.9285D-08                                  *****
      EXX(1) = 5.077973607D+00
      CP(1) = -3.329929840D-03
      EXX(2) = 1.340786940D+00
      CP(2) = -1.419488340D-02
      EXX(3) = 2.248434849D-01
      CP(3) = 1.639395770D-01
      EXX(4) = 1.131741848D-01
      CP(4) = 4.485358256D-01
      EXX(5) = 6.076408893D-02
      CP(5) = 3.908813050D-01
      EXX(6) = 3.315424265D-02
      CP(6) = 7.411456232D-02
      RETURN
      END
C---------------------------------------------------3S------------
C*MODULE BASSTO  *DECK S3S
      SUBROUTINE S3S(EXX,CS,NGAUSS)
C     COPYRIGHT  UF 1993  DAVID BAKER AND MICHAEL ZERNER
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C        *****  LEAST SQUARES FITS TO SLATER 3S FUNCTIONS BY    *****
C        *****  GAUSSIAN 1S FUNCTIONS                           *****
      DIMENSION EXX(6),CS(6)
      GO TO (100,120,140,160,180,200),NGAUSS
  100 CONTINUE
C        *****  STO(3S---1(1S))                                 *****
C        *****  ER1 1.6764D-02                                  *****
      EXX(1) = 5.296881757D-02
      CS(1) = 1.000000000D+00
      RETURN
  120 CONTINUE
C        *****  STO(3S---2(1S))                                 *****
C        *****  ER1 7.6424D-04                                  *****
      EXX(1) = 6.694095822D-01
      CS(1) = -1.529645716D-01
      EXX(2) = 5.837135094D-02
      CS(2) = 1.051370110D+00
      RETURN
  140 CONTINUE
C        *****  STO(3S---3(1S))                                 *****
C        *****  ER1 8.0718D-05                                  *****
      EXX(1) = 5.641487709D-01
      CS(1) = -1.782577972D-01
      EXX(2) = 6.924421391D-02
      CS(2) = 8.612761663D-01
      EXX(3) = 3.269529097D-02
      CS(3) = 2.261841969D-01
      RETURN
  160 CONTINUE
C        *****  STO(3S---4(1S))                                 *****
C        *****  ER1 1.7254D-06                                  *****
      EXX(1) = 1.513265591D+00
      CS(1) = -3.295496352D-02
      EXX(2) = 4.262497508D-01
      CS(2) = -1.724516959D-01
      EXX(3) = 7.643320863D-02
      CS(3) = 7.518511194D-01
      EXX(4) = 3.760545063D-02
      CS(4) = 3.589627317D-01
      RETURN
  180 CONTINUE
C        *****  STO(3S---5(1S))                                 *****
C        *****  ER1 7.9816D-07                                  *****
      EXX(1) = 4.275877914D+00
      CS(1) = -3.920358850D-03
      EXX(2) = 1.132409433D+00
      CS(2) = -4.168430506D-02
      EXX(3) = 4.016256968D-01
      CS(3) = -1.637440990D-01
      EXX(4) = 7.732370620D-02
      CS(4) = 7.419373723D-01
      EXX(5) = 3.800708627D-02
      CS(5) = 3.724364929D-01
      RETURN
  200 CONTINUE
C        *****  STO(3S---6(1S))                                 *****
C        *****  ER1 4.0662D-08                                  *****
      EXX(1) = 3.273031938D+00
      CS(1) = -6.775596947D-03
      EXX(2) = 9.200611311D-01
      CS(2) = -5.639325779D-02
      EXX(3) = 3.593349765D-01
      CS(3) = -1.587856086D-01
      EXX(4) = 8.636686991D-02
      CS(4) = 5.534527651D-01
      EXX(5) = 4.797373812D-02
      CS(5) = 5.015351020D-01
      EXX(6) = 2.724741144D-02
      CS(6) = 7.223633674D-02
      RETURN
      END
C---------------------------------------------------4D------------
C*MODULE BASSTO  *DECK S4D
      SUBROUTINE S4D(EXX,CD,NGAUSS)
C     COPYRIGHT  UF 1993  DAVID BAKER AND MICHAEL ZERNER
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION EXX(6),CD(6)
      GO TO (100,120,140,160,180,200),NGAUSS
  100 CONTINUE
C*****************      STO(4D---1(3D))
C*********           ER1 2.0360D-02
      EXX(1) = 7.941656339D-02
      CD(1) = 1.000000000D+00
      RETURN
  120 CONTINUE
C*****************      STO(4D---2(3D))
C*********           ER1 3.7488D-04
      EXX(1) = 1.330958892D-01
      CD(1) = 4.932764167D-01
      EXX(2) = 5.272119659D-02
      CD(2) = 5.918727866D-01
      RETURN
  140 CONTINUE
C*****************      STO(4D---3(3D))
C*********           ER1 1.2729D-05
      EXX(1) = 1.777717219D-01
      CD(1) = 2.308552718D-01
      EXX(2) = 8.040647350D-02
      CD(2) = 6.042409177D-01
      EXX(3) = 3.949855551D-02
      CD(3) = 2.595768926D-01
      RETURN
  160 CONTINUE
C*****************      STO(4D---4(3D))
C*********           ER1 5.7655D-06
      EXX(1) = 1.995825422D+00
      CD(1) = -2.816702620D-03
      EXX(2) = 1.823461280D-01
      CD(2) = 2.177095871D-01
      EXX(3) = 8.197240896D-02
      CD(3) = 6.058047348D-01
      EXX(4) = 4.000634951D-02
      CD(4) = 2.717811257D-01
      RETURN
  180 CONTINUE
C*****************      STO(4D---5(3D))
C*********           ER1 1.5122D-07
      EXX(1) = 1.522122079D+00
      CD(1) = -3.673711876D-03
      EXX(2) = 2.173041823D-01
      CD(2) = 1.167122499D-01
      EXX(3) = 1.084876577D-01
      CD(3) = 4.216476416D-01
      EXX(4) = 5.836797641D-02
      CD(4) = 4.547673415D-01
      EXX(5) = 3.206682246D-02
      CD(5) = 1.037803318D-01
      RETURN
  200 CONTINUE
C*****************      STO(4D---6(3D))
C*********           ER1 8.0888D-08
      EXX(1) = 4.634239420D+00
      CD(1) = -4.749842876D-04
      EXX(2) = 1.341648295D+00
      CD(2) = -3.566777891D-03
      EXX(3) = 2.209593028D-01
      CD(3) = 1.108670481D-01
      EXX(4) = 1.101467943D-01
      CD(4) = 4.159646930D-01
      EXX(5) = 5.904190370D-02
      CD(5) = 4.621672517D-01
      EXX(6) = 3.232628887D-02
      CD(6) = 1.081250196D-01
      RETURN
      END
C---------------------------------------------------4P------------
C*MODULE BASSTO  *DECK S4P
      SUBROUTINE S4P(EXX,CP,NGAUSS)
C     COPYRIGHT  UF 1993  DAVID BAKER AND MICHAEL ZERNER
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION EXX(6),CP(6)
      GO TO (100,120,140,160,180,200),NGAUSS
  100 CONTINUE
C*****************      STO(4P---1(2P))
C*********           ER1 3.9516D-03
      EXX(1) = 5.578350235D-02
      CP(1) = 1.000000000D+00
      RETURN
  120 CONTINUE
C*****************      STO(4P---2(2P))
C*********           ER1 2.1579D-03
      EXX(1) = 6.190052680D-02
      CP(1) = 8.743116767D-01
      EXX(2) = 2.648418407D-02
      CP(2) = 1.513640107D-01
      RETURN
  140 CONTINUE
C*****************      STO(4P---3(2P))
C*********           ER1 1.2621D-05
      EXX(1) = 4.859692220D-01
      CP(1) = -6.147823411D-02
      EXX(2) = 7.430216918D-02
      CP(2) = 6.604172234D-01
      EXX(3) = 3.653340923D-02
      CP(3) = 3.932639495D-01
      RETURN
  160 CONTINUE
C*****************      STO(4P---4(2P))
C*********           ER1 5.4740D-06
      EXX(1) = 1.492607880D+00
      CP(1) = -6.035216774D-03
      EXX(2) = 4.327619272D-01
      CP(2) = -6.013310874D-02
      EXX(3) = 7.553156064D-02
      CP(3) = 6.451518200D-01
      EXX(4) = 3.706272183D-02
      CP(4) = 4.117923820D-01
      RETURN
  180 CONTINUE
C*****************      STO(4P---5(2P))
C*********           ER1 1.3595D-07
      EXX(1) = 1.091977298D+00
      CP(1) = -1.143929558D-02
      EXX(2) = 3.719985051D-01
      CP(2) = -6.322651538D-02
      EXX(3) = 8.590019352D-02
      CP(3) = 4.398907721D-01
      EXX(4) = 4.786503860D-02
      CP(4) = 5.245859166D-01
      EXX(5) = 2.730479990D-02
      CP(5) = 1.017072253D-01
      RETURN
  200 CONTINUE
C*****************      STO(4P---6(2P))
C*********           ER1 1.3897D-08
      EXX(1) = 2.389722618D+00
      CP(1) = -1.665913575D-03
      EXX(2) = 7.960947826D-01
      CP(2) = -1.657464971D-02
      EXX(3) = 3.415541380D-01
      CP(3) = -5.958513378D-02
      EXX(4) = 8.847434525D-02
      CP(4) = 4.053115554D-01
      EXX(5) = 4.958248334D-02
      CP(5) = 5.433958189D-01
      EXX(6) = 2.816929784D-02
      CP(6) = 1.204970491D-01
      RETURN
      END
C---------------------------------------------------4S------------
C*MODULE BASSTO  *DECK S4S
      SUBROUTINE S4S(EXX,CS,NGAUSS)
C     COPYRIGHT  UF 1993  DAVID BAKER AND MICHAEL ZERNER
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION EXX(6),CS(6)
      GO TO (100,120,140,160,180,200),NGAUSS
  100 CONTINUE
C*****************      STO(4S---1(1S))
C*********           ER1 4.1760D-02
      EXX(1) = 3.264600274D-02
      CS(1) = 1.000000000D+00
      RETURN
  120 CONTINUE
C*****************      STO(4S---2(1S))
C*********           ER1 1.7480D-04
      EXX(1) = 2.441785453D-01
      CS(1) = -3.046656896D-01
      EXX(2) = 4.051097664D-02
      CS(2) = 1.146877294D+00
      RETURN
  140 CONTINUE
C*****************      STO(4S---3(1S))
C*********           ER1 1.6867D-05
      EXX(1) = 2.267938753D-01
      CS(1) = -3.349048323D-01
      EXX(2) = 4.448178019D-02
      CS(2) = 1.056744667D+00
      EXX(3) = 2.195294664D-02
      CS(3) = 1.256661680D-01
      RETURN
  160 CONTINUE
C*****************      STO(4S---4(1S))
C*********           ER1 8.4105D-07
      EXX(1) = 3.242212833D-01
      CS(1) = -1.120682822D-01
      EXX(2) = 1.663217177D-01
      CS(2) = -2.845426863D-01
      EXX(3) = 5.081097451D-02
      CS(3) = 8.909873788D-01
      EXX(4) = 2.829066600D-02
      CS(4) = 3.517811205D-01
      RETURN
  180 CONTINUE
C*****************      STO(4S---5(1S))
C*********           ER1 6.0273D-08
      EXX(1) = 2.980263783D+00
      CS(1) = 1.513948997D-03
      EXX(2) = 3.792228833D-01
      CS(2) = -7.316801548D-02
      EXX(3) = 1.789717224D-01
      CS(3) = -3.143703799D-01
      EXX(4) = 5.002110360D-02
      CS(4) = 9.032615169D-01
      EXX(5) = 2.789361681D-02
      CS(5) = 3.294210848D-01
      RETURN
  200 CONTINUE
C*****************      STO(4S---6(1S))
C*********           ER1 5.9700D-09
      EXX(1) = 3.232838646D+00
      CS(1) = 1.374817488D-03
      EXX(2) = 3.605788802D-01
      CS(2) = -8.666390043D-02
      EXX(3) = 1.717905487D-01
      CS(3) = -3.130627309D-01
      EXX(4) = 5.277666487D-02
      CS(4) = 7.812787397D-01
      EXX(5) = 3.163400284D-02
      CS(5) = 4.389247988D-01
      EXX(6) = 1.874093091D-02
      CS(6) = 2.487178756D-02
      RETURN
      END
C---------------------------------------------------5D------------
C*MODULE BASSTO  *DECK S5D
      SUBROUTINE S5D(EXX,CD,NGAUSS)
C     COPYRIGHT  UF 1993  DAVID BAKER AND MICHAEL ZERNER
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION EXX(6),CD(6)
      GO TO (100,120,140,160,180,200),NGAUSS
  100 CONTINUE
C*****************      STO(5D---1(3D))
C*********           ER1 6.7308D-03
      EXX(1) = 5.352200793D-02
      CD(1) = 1.000000000D+00
      RETURN
  120 CONTINUE
C*****************      STO(5D---2(3D))
C*********           ER1 3.1257D-04
      EXX(1) = 6.906014388D-02
      CD(1) = 6.539405185D-01
      EXX(2) = 3.399457777D-02
      CD(2) = 3.948945302D-01
      RETURN
  140 CONTINUE
C*****************      STO(5D---3(3D))
C*********           ER1 2.4664D-05
      EXX(1) = 4.913352950D-01
      CD(1) = -2.010175008D-02
      EXX(2) = 7.329090601D-02
      CD(2) = 5.899370608D-01
      EXX(3) = 3.594209290D-02
      CD(3) = 4.658445960D-01
      RETURN
  160 CONTINUE
C*****************      STO(5D---4(3D))
C*********           ER1 1.2737D-06
      EXX(1) = 4.230617826D-01
      CD(1) = -2.421626009D-02
      EXX(2) = 8.293863702D-02
      CD(2) = 3.937644956D-01
      EXX(3) = 4.590326388D-02
      CD(3) = 5.489520286D-01
      EXX(4) = 2.628744797D-02
      CD(4) = 1.190436963D-01
      RETURN
  180 CONTINUE
C*****************      STO(5D---5(3D))
C*********           ER1 8.6929D-08
      EXX(1) = 9.702946470D-01
      CD(1) = -3.231527611D-03
      EXX(2) = 3.603270196D-01
      CD(2) = -2.434931372D-02
      EXX(3) = 8.668717752D-02
      CD(3) = 3.440817054D-01
      EXX(4) = 4.833708379D-02
      CD(4) = 5.693674376D-01
      EXX(5) = 2.751899341D-02
      CD(5) = 1.511340183D-01
      RETURN
  200 CONTINUE
C*****************      STO(5D---6(3D))
C*********           ER1 1.4086D-08
      EXX(1) = 8.820520428D-01
      CD(1) = -4.097377019D-03
      EXX(2) = 3.410838409D-01
      CD(2) = -2.508271857D-02
      EXX(3) = 9.204308840D-02
      CD(3) = 2.648458555D-01
      EXX(4) = 5.472831774D-02
      CD(4) = 5.097437054D-01
      EXX(5) = 3.391202830D-02
      CD(5) = 2.654483467D-01
      EXX(6) = 2.108227374D-02
      CD(6) = 2.623132212D-02
      RETURN
      END
C---------------------------------------------------5P------------
C*MODULE BASSTO  *DECK S5P
      SUBROUTINE S5P(EXX,CP,NGAUSS)
C     COPYRIGHT  UF 1993  DAVID BAKER AND MICHAEL ZERNER
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION EXX(6),CP(6)
      GO TO (100,120,140,160,180,200),NGAUSS
  100 CONTINUE
C*****************      STO(5P---1(2P))
C*********           ER1 8.0406D-03
      EXX(1) = 3.769845216D-02
      CP(1) = 1.000000000D+00
      RETURN
  120 CONTINUE
C*****************      STO(5P---2(2P))
C*********           ER1 6.8912D-04
      EXX(1) = 2.691294191D-01
      CP(1) = -1.034227010D-01
      EXX(2) = 3.980805011D-02
      CP(2) = 1.033376378D+00
      RETURN
  140 CONTINUE
C*****************      STO(5P---3(2P))
C*********           ER1 6.5766D-06
      EXX(1) = 2.127482317D-01
      CP(1) = -1.389529695D-01
      EXX(2) = 4.729648620D-02
      CP(2) = 8.076691064D-01
      EXX(3) = 2.604865324D-02
      CP(3) = 2.726029342D-01
      RETURN
  160 CONTINUE
C*****************      STO(5P---4(2P))
C*********           ER1 4.4493D-07
      EXX(1) = 3.962838833D-01
      CP(1) = -1.801459207D-02
      EXX(2) = 1.838858552D-01
      CP(2) = -1.360777372D-01
      EXX(3) = 4.943555157D-02
      CP(3) = 7.533973719D-01
      EXX(4) = 2.750222273D-02
      CP(4) = 3.409304859D-01
      RETURN
  180 CONTINUE
C*****************      STO(5P---5(2P))
C*********           ER1 1.3483D-08
      EXX(1) = 3.422168934D-01
      CP(1) = -3.113958289D-02
      EXX(2) = 1.665099900D-01
      CP(2) = -1.374007017D-01
      EXX(3) = 5.443732013D-02
      CP(3) = 5.573881018D-01
      EXX(4) = 3.367775277D-02
      CP(4) = 4.855428100D-01
      EXX(5) = 2.091949042D-02
      CP(5) = 6.605423564D-02
      RETURN
  200 CONTINUE
C*****************      STO(5P---6(2P))
C*********           ER1 3.9301D-09
      EXX(1) = 3.778623374D+00
      CP(1) = 1.163246387D-04
      EXX(2) = 3.499121109D-01
      CP(2) = -2.920771322D-02
      EXX(3) = 1.683175469D-01
      CP(3) = -1.381051233D-01
      EXX(4) = 5.404070736D-02
      CP(4) = 5.706134877D-01
      EXX(5) = 3.328911801D-02
      CP(5) = 4.768808140D-01
      EXX(6) = 2.063815019D-02
      CP(6) = 6.021665516D-02
      RETURN
      END
C---------------------------------------------------5S------------
C*MODULE BASSTO  *DECK S5S
      SUBROUTINE S5S(EXX,CS,NGAUSS)
C     COPYRIGHT  UF 1993  DAVID BAKER AND MICHAEL ZERNER
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION EXX(6),CS(6)
      GO TO (100,120,140,160,180,200),NGAUSS
  100 CONTINUE
C*****************      STO(5S---1(1S))
C*********           ER1 7.1317D-02
      EXX(1) = 2.216912938D-02
      CS(1) = 1.000000000D+00
      RETURN
  120 CONTINUE
C*****************      STO(5S---2(1S))
C*********           ER1 1.4289D-04
      EXX(1) = 1.213425654D-01
      CS(1) = -5.114756049D-01
      EXX(2) = 3.133152144D-02
      CS(2) = 1.307377277D+00
      RETURN
  140 CONTINUE
C*****************      STO(5S---3(1S))
C*********           ER1 4.1124D-05
      EXX(1) = 1.080198458D-01
      CS(1) = -6.617401158D-01
      EXX(2) = 4.408119382D-02
      CS(2) = 7.467595004D-01
      EXX(3) = 2.610811810D-02
      CS(3) = 7.146490945D-01
      RETURN
  160 CONTINUE
C*****************      STO(5S---4(1S))
C*********           ER1 5.4159D-07
      EXX(1) = 8.602284252D-01
      CS(1) = 1.103657561D-02
      EXX(2) = 1.189050200D-01
      CS(2) = -5.606519023D-01
      EXX(3) = 3.446076176D-02
      CS(3) = 1.179429987D+00
      EXX(4) = 1.974798796D-02
      CS(4) = 1.734974376D-01
      RETURN
  180 CONTINUE
C*****************      STO(5S---5(1S))
C*********           ER1 7.0816D-08
      EXX(1) = 7.403763257D-01
      CS(1) = 1.375523371D-02
      EXX(2) = 1.367990863D-01
      CS(2) = -3.097344179D-01
      EXX(3) = 9.135301779D-02
      CS(3) = -3.199192259D-01
      EXX(4) = 3.726907315D-02
      CS(4) = 1.084547038D+00
      EXX(5) = 2.241490836D-02
      CS(5) = 3.345288361D-01
      RETURN
  200 CONTINUE
C*****************      STO(5S---6(1S))
C*********           ER1 7.9988D-09
      EXX(1) = 1.410128298D+00
      CS(1) = 2.695439582D-03
      EXX(2) = 5.077878915D-01
      CS(2) = 1.850157487D-02
      EXX(3) = 1.847926858D-01
      CS(3) = -9.588628125D-02
      EXX(4) = 1.061070594D-01
      CS(4) = -5.200673560D-01
      EXX(5) = 3.669584901D-02
      CS(5) = 1.087619490D+00
      EXX(6) = 2.213558430D-02
      CS(6) = 3.103964343D-01
      RETURN
      END
