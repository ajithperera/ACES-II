C
      SUBROUTINE FILDGR (L, M, N, GA, V, NT, D)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
C &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C &  Calculate electric field gradient with out Fermi contact   &
C &  contribution. Hacked up version from vprop.f.              &
C &  08/93 Ajith                                                &
C &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C 
      COMMON /PRTBUG/ IPRTFS
      COMMON /VARS/ FNU(20),DP(20),FN(15),FD(15),PI,DFTR(11)
      DIMENSION V(10), D(10), DD(10)
C
      DATA FQR/1.3333333333333D0/
      DATA SQRPI/1.772453851D0/
C
      NT = 6
      ITOT = L + M + N
      ITOTH = ITOT/2
      PRE = (8.00D+00*PI*GA)*FN(L+1)*FN(M+1)*FN(N+1)
      IF(2*ITOTH - ITOT)103,104,103
103   PRE = -PRE
104   DEL = 0.25D+00/GA
      X = GA*(D(1)**2 + D(2)**2 + D(3)**2)
      XX = X + X
C
      CALL FMC(14, X, EXPMX, FMCH)
C
      FNU(15) = FMCH
      FNU(14) = (EXPMX + XX*FNU(15))/27.00D+00
      FNU(13) = (EXPMX + XX*FNU(14))/25.00D+00
      FNU(12) = (EXPMX + XX*FNU(13))/23.00D+00
      FNU(11) = (EXPMX + XX*FNU(12))/21.00D+00
      FNU(10) = (EXPMX + XX*FNU(11))/19.00D+00
      FNU(9)  = (EXPMX + XX*FNU(10))/17.00D+00
      FNU(8)  = (EXPMX + XX*FNU(9))/15.00D+00
      FNU(7)  = (EXPMX + XX*FNU(8))/13.00D+00
      FNU(6)  = (EXPMX + XX*FNU(7))/11.00D+00
      FNU(5)  = (EXPMX + XX*FNU(6))/9.00D+00
      FNU(4)  = (EXPMX + XX*FNU(5))/7.00D+00
      FNU(3)  = (EXPMX + XX*FNU(4))/5.00D+00
      FNU(2)  = (EXPMX + XX*FNU(3))/3.00D+00
      FNU(1)  = EXPMX  + XX*FNU(2)
C
      DP(1)  = 1.00D+00
      DP(2)  = DEL
      DP(3)  = DEL**2
      DP(4)  = DEL*DP(3)
      DP(5)  = DP(4)*DEL
      DP(6)  = DP(5)*DEL
      DP(7)  = DP(6)*DEL
      DP(8)  = DP(7)*DEL
      DP(9)  = DP(8)*DEL
      DP(10) = DP(9)*DEL
      DP(11) = DP(10)*DEL
C
      V(1) = PRE*AAINER(2,0,0,L,M,N,D)
      V(2) = PRE*AAINER(1,1,0,L,M,N,D)
      V(3) = PRE*AAINER(1,0,1,L,M,N,D)
      V(4) = PRE*AAINER(0,2,0,L,M,N,D)
      V(5) = PRE*AAINER(0,1,1,L,M,N,D)
      V(6) = PRE*AAINER(0,0,2,L,M,N,D)
C
      DDSQ = 0.00D+00
C
      DO 411 I = 1, 3
         DD(I) = -D(I)
         DDSQ  = DDSQ+D(I)**2
 411  CONTINUE
C
      DP(1) = 1.00D+00
      IF(L) 412, 413, 412
 412  DP(1) = DD(1)**L
C
 413  IF(M) 414, 415, 414
 414  DP(1) = DP(1)*DD(2)**M
C
 415  IF(N) 416, 417, 416
 416  DP(1) = DP(1)*DD(3)**N
C
 417  DP(1) = DP(1)*EXP(-GA*DDSQ)
C
      DP(1) = FQR*PI*DP(1)
      V(1)  = V(1) + DP(1)
      V(4)  = V(4) + DP(1)
      V(6)  = V(6) + DP(1)
C
      DO 418 I = 1, 6
      V(I) = -V(I)
 418  CONTINUE
C
      NT=6
C
      RETURN
      END
