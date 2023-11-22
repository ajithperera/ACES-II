      SUBROUTINE EHTPARAM(ECOR,EVAL,BETA0)
      IMPLICIT NONE
      DOUBLE PRECISION ECOR,EVAL,BETA0
      DOUBLE PRECISION FACTOR
      DIMENSION ECOR(5,18),EVAL(4,18),BETA0(18)
C
      DATA FACTOR /27.2114D+00/
C
      CALL ZERO(ECOR, 5*18)
      CALL ZERO(EVAL, 4*18)
      CALL ZERO(BETA0,  18)
C
      BETA0( 1) =  -9.0D+00 / FACTOR
C
      BETA0( 3) =  -9.0D+00 / FACTOR
      BETA0( 4) = -13.0D+00 / FACTOR
      BETA0( 5) = -17.0D+00 / FACTOR
      BETA0( 6) = -21.0D+00 / FACTOR
      BETA0( 7) = -25.0D+00 / FACTOR
      BETA0( 8) = -31.0D+00 / FACTOR
      BETA0( 9) = -39.0D+00 / FACTOR
C
      BETA0(11) =  -7.720D+00 / FACTOR
      BETA0(12) =  -9.447D+00 / FACTOR
      BETA0(13) = -11.301D+00 / FACTOR
      BETA0(14) = -13.065D+00 / FACTOR
      BETA0(15) = -15.070D+00 / FACTOR
      BETA0(16) = -18.150D+00 / FACTOR
      BETA0(17) = -22.330D+00 / FACTOR
C
      EVAL(1, 1) =  -7.176D+00 / FACTOR
C
      EVAL(1, 3) =  -3.106D+00 / FACTOR
      EVAL(2, 3) =  -1.258D+00 / FACTOR
      EVAL(3, 3) =  -1.258D+00 / FACTOR
      EVAL(4, 3) =  -1.258D+00 / FACTOR
C
      EVAL(1, 4) =  -5.946D+00 / FACTOR
      EVAL(2, 4) =  -2.563D+00 / FACTOR
      EVAL(3, 4) =  -2.563D+00 / FACTOR
      EVAL(4, 4) =  -2.563D+00 / FACTOR
C
      EVAL(1, 5) =  -9.594D+00 / FACTOR
      EVAL(2, 5) =  -4.001D+00 / FACTOR
      EVAL(3, 5) =  -4.001D+00 / FACTOR
      EVAL(4, 5) =  -4.001D+00 / FACTOR
C
      EVAL(1, 6) = -14.051D+00 / FACTOR
      EVAL(2, 6) =  -5.572D+00 / FACTOR
      EVAL(3, 6) =  -5.572D+00 / FACTOR
      EVAL(4, 6) =  -5.572D+00 / FACTOR
C
      EVAL(1, 7) = -19.316D+00 / FACTOR
      EVAL(2, 7) =  -7.275D+00 / FACTOR
      EVAL(3, 7) =  -7.275D+00 / FACTOR
      EVAL(4, 7) =  -7.275D+00 / FACTOR
C
      EVAL(1, 8) = -25.390D+00 / FACTOR
      EVAL(2, 8) =  -9.111D+00 / FACTOR
      EVAL(3, 8) =  -9.111D+00 / FACTOR
      EVAL(4, 8) =  -9.111D+00 / FACTOR
C
      EVAL(1, 9) = -32.272D+00 / FACTOR
      EVAL(2, 9) = -11.080D+00 / FACTOR
      EVAL(3, 9) = -11.080D+00 / FACTOR
      EVAL(4, 9) = -11.080D+00 / FACTOR
C
      EVAL(1,11) =  -2.804D+00 / FACTOR
      EVAL(2,11) =  -1.302D+00 / FACTOR
      EVAL(3,11) =  -1.302D+00 / FACTOR
      EVAL(4,11) =  -1.302D+00 / FACTOR
C
      EVAL(1,12) =  -5.125D+00 / FACTOR
      EVAL(2,12) =  -2.052D+00 / FACTOR
      EVAL(3,12) =  -2.052D+00 / FACTOR
      EVAL(4,12) =  -2.052D+00 / FACTOR
C
      EVAL(1,13) =  -7.771D+00 / FACTOR
      EVAL(2,13) =  -2.995D+00 / FACTOR
      EVAL(3,13) =  -2.995D+00 / FACTOR
      EVAL(4,13) =  -2.995D+00 / FACTOR
C
      EVAL(1,14) = -10.033D+00 / FACTOR
      EVAL(2,14) =  -4.133D+00 / FACTOR
      EVAL(3,14) =  -4.133D+00 / FACTOR
      EVAL(4,14) =  -4.133D+00 / FACTOR
C
      EVAL(1,15) = -14.033D+00 / FACTOR
      EVAL(2,15) =  -5.464D+00 / FACTOR
      EVAL(3,15) =  -5.464D+00 / FACTOR
      EVAL(4,15) =  -5.464D+00 / FACTOR
C
      EVAL(1,16) = -17.650D+00 / FACTOR
      EVAL(2,16) =  -6.989D+00 / FACTOR
      EVAL(3,16) =  -6.989D+00 / FACTOR
      EVAL(4,16) =  -6.989D+00 / FACTOR
C
      EVAL(1,17) = -21.591D+00 / FACTOR
      EVAL(2,17) =  -8.708D+00 / FACTOR
      EVAL(3,17) =  -8.708D+00 / FACTOR
      EVAL(4,17) =  -8.708D+00 / FACTOR
C
C-----------------------------------------------------------------------
C     Core parameters.
C-----------------------------------------------------------------------
      ECOR(1, 6) =   -300.0D+00 / FACTOR
      ECOR(1, 7) =   -420.0D+00 / FACTOR
      ECOR(1, 8) =   -551.0D+00 / FACTOR
      ECOR(1, 9) =   -705.0D+00 / FACTOR
C
      ECOR(1,14) =  -1850.0D+00 / FACTOR
      ECOR(2,14) =   -162.0D+00 / FACTOR
      ECOR(3,14) =   -108.0D+00 / FACTOR
      ECOR(4,14) =   -108.0D+00 / FACTOR
      ECOR(5,14) =   -108.0D+00 / FACTOR
C
      ECOR(1,15) =  -2142.0D+00 / FACTOR
      ECOR(2,15) =   -193.0D+00 / FACTOR
      ECOR(3,15) =   -134.0D+00 / FACTOR
      ECOR(4,15) =   -134.0D+00 / FACTOR
      ECOR(5,15) =   -134.0D+00 / FACTOR
C
      ECOR(1,16) =  -2469.0D+00 / FACTOR
      ECOR(2,16) =   -235.0D+00 / FACTOR
      ECOR(3,16) =   -180.0D+00 / FACTOR
      ECOR(4,16) =   -180.0D+00 / FACTOR
      ECOR(5,16) =   -180.0D+00 / FACTOR
C
      ECOR(1,17) =  -2822.0D+00 / FACTOR
      ECOR(2,17) =   -283.0D+00 / FACTOR
      ECOR(3,17) =   -213.0D+00 / FACTOR
      ECOR(4,17) =   -213.0D+00 / FACTOR
      ECOR(5,17) =   -213.0D+00 / FACTOR
C
      RETURN
      END
