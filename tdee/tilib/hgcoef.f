      SUBROUTINE HGCOEF (METH, NQ, EL, TQ, MAXDER)
C
C
C     Description of routine.
C
C     HGCOEF is called by the integrator and sets coefficients used
C     there. The vector EL, of length NQ + 1, determines the basic
C     method. The vector TQ, of length 4, is involved in adjusting
C     the step size in relation to truncation error.  Its values are
C     given by the PERTST array.
C     The vectors EL and TQ depend on METH and NQ.
C     HGCOEF also sets MAXDER, the maximum order of the method available
C     Currently it is 12 for the Adams methods and 5 for the
C     BDF methods.
C     The maximum order used may be reduced simply by decreasing the
C     numbers in statements 1 and/or 2 below.
C
C     The coefficients in PERTST need be given to only about
C     one percent accuracy.  The order in which the groups appear below
C     is:
C
C     coefficients for order NQ - 1,
C
C     coefficients for order NQ,
C
C     coefficients for order NQ + 1.
C
C     Within each group are the coefficients for the Adams methods,
C     followed by those for the BDF methods.
C
C     The EL coefficients should be defined to machine accuracy.
C     For a given order NQ, they can be calculated by use of the
C     generating polynomial L(T), whose coefficients are EL(I):
C     
C       L(T) = EL(1) + EL(2)*T + ... + EL(NQ+1)*T**NQ.
C     
C     for the implicit Adams methods, L(T) is given by
C     
C       DL/DT = (T+1)*(T+2)* ... *(T+NQ-1)/K,    L(-1) = 0,
C     
C     where
C     
C     K = FACTORIAL(NQ-1).
C     
C     For the BDF methods,
C     
C       L(T) = (T+1)*(T+2)* ... *(T+NQ)/K,
C     
C     where
C     
C     K = FACTORIAL(NQ)*(1 + 1/2 + ... + 1/NQ).
C     
C     E(2) is always 1. and is initialized in the calling routine.
C
C     The order in which the groups appear below is:
C
C     implicit Adams methods of orders 1 to 12,
C
C     BDF methods of orders 1 to 5.
C
C     Original code: copyright A. C. Hindmarsh and LLL, 1975.
C     Modifications: copyright Erik Deumens and QTP, 1986, 1992.
C
      IMPLICIT NONE

      DOUBLE PRECISION PERTST(12,2,3),EL(13),TQ(4)
      INTEGER METH, NQ, MAXDER, K
      DATA  PERTST /
     $     1.D0,1.D0,2.D0,1.D0,.3158D0,.07407D0,.01391D0,
     $     .002182D0,.0002945D0,.00003492D0,.000003692D0,.0000003524D0,
     $     1.D0,1.D0,.5D0,.1667D0,.04167D0,1.D0,
     $     1.D0,1.D0,1.D0,1.D0,1.D0,1.D0,
     $     2.D0,12.D0,24.D0,37.89D0,53.33D0,70.08D0,
     $     87.97D0,106.9D0,126.7D0,147.4D0,168.8D0,191.D0,
     $     2.D0,4.5D0,7.333D0,10.42D0,13.7D0,1.D0,
     $     1.D0,1.D0,1.D0,1.D0,1.D0,1.D0,
     $     12.D0,24.D0,37.89D0,53.33D0,70.08D0,87.97D0,
     $     106.9D0,126.7D0,147.4D0,168.8D0,191.D0,1.D0,
     $     3.D0,6.D0,9.167D0,12.5D0,1.D0,1.D0,
     $     1.D0,1.D0,1.D0,1.D0,1.D0,1.D0 /
C
      GO TO (1,2),METH
 1    MAXDER = 12
      GO TO (101,102,103,104,105,106,107,108,109,110,111,112),NQ
 2    MAXDER = 5
      GO TO (201,202,203,204,205),NQ
 101  EL(1) = 1.D0
      GO TO 900
 102  EL(1) = 0.5D0
      EL(3) = 0.5D0
      GO TO 900
 103  EL(1) = 4.1666666666667D-01
      EL(3) = 0.75D0
      EL(4) = 1.6666666666667D-01
      GO TO 900
 104  EL(1) = 0.375D0
      EL(3) = 9.1666666666667D-01
      EL(4) = 3.3333333333333D-01
      EL(5) = 4.1666666666667D-02
      GO TO 900
 105  EL(1) = 3.4861111111111D-01
      EL(3) = 1.0416666666667D0
      EL(4) = 4.8611111111111D-01
      EL(5) = 1.0416666666667D-01
      EL(6) = 8.3333333333333D-03
      GO TO 900
 106  EL(1) = 3.2986111111111D-01
      EL(3) = 1.1416666666667D0
      EL(4) = 0.625D0
      EL(5) = 1.7708333333333D-01
      EL(6) = 0.025D0
      EL(7) = 1.3888888888889D-03
      GO TO 900
 107  EL(1) = 3.1559193121693D-01
      EL(3) = 1.225D0
      EL(4) = 7.5185185185185D-01
      EL(5) = 2.5520833333333D-01
      EL(6) = 4.8611111111111D-02
      EL(7) = 4.8611111111111D-03
      EL(8) = 1.9841269841270D-04
      GO TO 900
 108  EL(1) = 3.0422453703704D-01
      EL(3) = 1.2964285714286D0
      EL(4) = 8.6851851851852D-01
      EL(5) = 3.3576388888889D-01
      EL(6) = 7.7777777777778D-02
      EL(7) = 1.0648148148148D-02
      EL(8) = 7.9365079365079D-04
      EL(9) = 2.4801587301587D-05
      GO TO 900
 109  EL(1) = 2.9486800044092D-01
      EL(3) = 1.3589285714286D0
      EL(4) = 9.7655423280423D-01
      EL(5) = 0.4171875D0
      EL(6) = 1.1135416666667D-01
      EL(7) = 0.01875D0
      EL(8) = 1.9345238095238D-03
      EL(9) = 1.1160714285714D-04
      EL(10)= 2.7557319223986D-06
      GO TO 900
 110  EL(1) = 2.8697544642857D-01
      EL(3) = 1.4144841269841D0
      EL(4) = 1.0772156084656D0
      EL(5) = 4.9856701940035D-01
      EL(6) = 0.1484375D0
      EL(7) = 2.9060570987654D-02
      EL(8) = 3.7202380952381D-03
      EL(9) = 2.9968584656085D-04
      EL(10)= 1.3778659611993D-05
      EL(11)= 2.7557319223986D-07
      GO TO 900
 111  EL(1) = 2.8018959644394D-01
      EL(3) = 1.4644841269841D0
      EL(4) = 1.1715145502646D0
      EL(5) = 5.7935819003527D-01
      EL(6) = 1.8832286155203D-01
      EL(7) = 4.1430362654321D-02
      EL(8) = 6.2111441798942D-03
      EL(9) = 6.2520667989418D-04
      EL(10)= 4.0417401528513D-05
      EL(11)= 1.5156525573192D-06
      EL(12)= 2.5052108385442D-08
      GO TO 900
 112  EL(1) = 2.7426554003160D-01
      EL(3) = 1.5099386724387D0
      EL(4) = 1.2602711640212D0
      EL(5) = 6.5923418209877D-01
      EL(6) = 2.3045800264550D-01
      EL(7) = 5.5697246105232D-02
      EL(8) = 9.4394841269841D-03
      EL(9) = 1.1192749669312D-03
      EL(10)= 9.0939153439153D-05
      EL(11)= 4.8225308641975D-06
      EL(12)= 1.5031265031265D-07
      EL(13)= 2.0876756987868D-09
      GO TO 900
 201  EL(1) = 1.D0
      GO TO 900
 202  EL(1) = 6.6666666666667D-01
      EL(3) = 3.3333333333333D-01
      GO TO 900
 203  EL(1) = 5.4545454545455D-01
      EL(3) = EL(1)
      EL(4) = 9.0909090909091D-02
      GO TO 900
 204  EL(1) = 0.48D0
      EL(3) = 0.7D0
      EL(4) = 0.2D0
      EL(5) = 0.02D0
      GO TO 900
 205  EL(1) = 4.3795620437956D-01
      EL(3) = 8.2116788321168D-01
      EL(4) = 3.1021897810219D-01
      EL(5) = 5.4744525547445D-02
      EL(6) = 3.6496350364964D-03
C
  900 DO 910 K = 1,3
        TQ(K) = PERTST(NQ,METH,K)
  910 CONTINUE
      TQ(4) = .5D0*TQ(2)/FLOAT(NQ+2)
      RETURN
      END
