C
      SUBROUTINE OPAA0(L,M,N,GA,V,NT,D)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /PRTBUG/ IPRTFS
C
C....    EXTENDED FOR F AND G FOR ALL PROPERTIES
C
      DIMENSION V(10),D(10), VS(10)
      DIMENSION DD(10),DS(10)
      COMMON /VARS/ FNU(20),DP(20),FN(15),FD(15),PI,DFTR(11)
C      EQUIVALENCE (DD(1),FNU(1)), (DS(1),FNU(11))
      DATA FQR/1.3333333333333D0/
      DATA SQRPI/1.772453851D0/
      NT=1
      ITOT=L+M+N
      ITOTH=ITOT/2
      PRE=(2.*PI/GA)*FN(L+1)*FN(M+1)*FN(N+1)
      IF(2*ITOTH-ITOT) 1,2,1
1     PRE=-PRE
2     DEL=.25/GA
      X=GA*(D(1)**2+D(2)**2+D(3)**2)
      XX = X + X
      CALL FMC(12,X,EXPMX,FMCH)
      FNU(13)=FMCH
      FNU(12)=(EXPMX+XX*FNU(13))/23.
      FNU(11)=(EXPMX+XX*FNU(12))/21.
      FNU(10)=(EXPMX+XX*FNU(11))/19.
      FNU(9)=(EXPMX+XX*FNU(10))/17.
      FNU(8)=(EXPMX+XX*FNU(9))/15.
      FNU(7)=(EXPMX+XX*FNU(8))/13.
      FNU(6)=(EXPMX+XX*FNU(7))/11.
      FNU(5)=(EXPMX+XX*FNU(6))/9.
      FNU(4)=(EXPMX+XX*FNU(5))/7.
      FNU(3)=(EXPMX+XX*FNU(4))/5.
      FNU(2)=(EXPMX+XX*FNU(3))/3.
      FNU(1)=EXPMX+XX*FNU(2)
      DP(1)=1.
      DP(2)=DEL
      DP(3)=DEL**2
      DP(4)=DP(3)*DEL
      DP(5)=DP(4)*DEL
      DP(6)=DP(5)*DEL
      DP(7)=DP(6)*DEL
      V(1) = -PRE*AAINER(0,0,0,L,M,N,D)
      RETURN
CRAY 1
      ENTRY OPAA1(L,M,N,GA,V,NT,D)
CDC 1
C      ENTRY OPAA1
      NT=3
      ITOT=L+M+N
      ITOTH=ITOT/2
      PRE=4.*PI*FN(L+1)*FN(M+1)*FN(N+1)
      IF(2*ITOTH-ITOT) 101,102,101
101   PRE=-PRE
102   DEL=.25/GA
      X=GA*(D(1)**2+D(2)**2+D(3)**2)
      XX = X + X
      CALL FMC(13,X,EXPMX,FMCH)
      FNU(14)=FMCH
      FNU(13)=(EXPMX+XX*FNU(14))/25.
      FNU(12)=(EXPMX+XX*FNU(13))/23.
      FNU(11)=(EXPMX+XX*FNU(12))/21.
      FNU(10)=(EXPMX+XX*FNU(11))/19.
      FNU(9)=(EXPMX+XX*FNU(10))/17.
      FNU(8)=(EXPMX+XX*FNU(9))/15.
      FNU(7)=(EXPMX+XX*FNU(8))/13.
      FNU(6)=(EXPMX+XX*FNU(7))/11.
      FNU(5)=(EXPMX+XX*FNU(6))/9.
      FNU(4)=(EXPMX+XX*FNU(5))/7.
      FNU(3)=(EXPMX+XX*FNU(4))/5.
      FNU(2)=(EXPMX+XX*FNU(3))/3.
      FNU(1)=EXPMX+XX*FNU(2)
      DP(1)=1.
      DP(2)=DEL
      DP(3)=DEL**2
      DP(4)=DEL*DP(3)
      DP(5)=DP(4)*DEL
      DP(6)=DP(5)*DEL
      DP(7)=DP(6)*DEL
      DP(8)=DP(7)*DEL
      DP(9)=DP(8)*DEL
      V(1) = PRE*AAINER(1,0,0,L,M,N,D)
      V(2) = PRE*AAINER(0,1,0,L,M,N,D)
      V(3) = PRE*AAINER(0,0,1,L,M,N,D)
      RETURN
CRAY 1
      ENTRY OPAA2(L,M,N,GA,V,NT,D)
CDC 1
C      ENTRY OPAA2
      ITOT=L+M+N
      ITOTH=ITOT/2
      PRE=(8.*PI*GA)*FN(L+1)*FN(M+1)*FN(N+1)
      IF(2*ITOTH-ITOT)103,104,103
103   PRE=-PRE
104   DEL=.25/GA
      X=GA*(D(1)**2+D(2)**2+D(3)**2)
      XX = X + X
      CALL FMC(14,X,EXPMX,FMCH)
      FNU(15)=FMCH
      FNU(14)=(EXPMX+XX*FNU(15))/27.
      FNU(13)=(EXPMX+XX*FNU(14))/25.
      FNU(12)=(EXPMX+XX*FNU(13))/23.
      FNU(11)=(EXPMX+XX*FNU(12))/21.
      FNU(10)=(EXPMX+XX*FNU(11))/19.
      FNU(9)=(EXPMX+XX*FNU(10))/17.
      FNU(8)=(EXPMX+XX*FNU(9))/15.
      FNU(7)=(EXPMX+XX*FNU(8))/13.
      FNU(6)=(EXPMX+XX*FNU(7))/11.
      FNU(5)=(EXPMX+XX*FNU(6))/9.
      FNU(4)=(EXPMX+XX*FNU(5))/7.
      FNU(3)=(EXPMX+XX*FNU(4))/5.
      FNU(2)=(EXPMX+XX*FNU(3))/3.
      FNU(1)=EXPMX+XX*FNU(2)
      DP(1)=1.
      DP(2)=DEL  
      DP(3)=DEL**2
      DP(4)=DEL*DP(3)
      DP(5)=DP(4)*DEL
      DP(6)=DP(5)*DEL
      DP(7)=DP(6)*DEL
      DP(8)=DP(7)*DEL
      DP(9)=DP(8)*DEL
      DP(10)=DP(9)*DEL
      DP(11)=DP(10)*DEL
      V(1) = PRE*AAINER(2,0,0,L,M,N,D)
      V(2) = PRE*AAINER(0,2,0,L,M,N,D)
      V(3) = PRE*AAINER(0,0,2,L,M,N,D)
      V(4) = PRE*AAINER(1,1,0,L,M,N,D)
      V(5) = PRE*AAINER(1,0,1,L,M,N,D)
      V(6) = PRE*AAINER(0,1,1,L,M,N,D)
C      IF (IPRTFS .EQ. 1) WRITE(6,87880) L, M, N, PRE, DEL, V
C87880 FORMAT(3I6,E16.8,E16.8,/,(3(4X,E16.8)))
C     CALL OPAC3(L,M,N,GA,DP(1),NT,D)
      DDSQ=0.
      DO 411 I=1,3
      DD(I)=-D(I)
      DDSQ=DDSQ+D(I)**2
411   CONTINUE
      DP(1) =1.
      IF(L) 412,413,412
412   DP(1)=DD(1)**L
413   IF(M) 414,415,414
414   DP(1)=DP(1)*DD(2)**M
415   IF(N) 416,417,416
416   DP(1)=DP(1)*DD(3)**N
417   DP(1) = DP(1)*EXP(-GA*DDSQ)
C     RETURN
      DP(1) = FQR*PI*DP(1)
      V(1)=V(1)+DP(1)
      V(2)=V(2)+DP(1)
      V(3)=V(3)+DP(1)
      DO 418 I = 1,6
      V(I)=-V(I)
418   CONTINUE
      NT=6
      RETURN
CRAY 1
      ENTRY OPAB1(L,M,N,GA,V,NT,D)
CDC 1
C      ENTRY OPAB1
      NT=3
      FNU(1)=OLAP(L,M,N,GA)
      FNU(2)=OLAP(L+1,M,N,GA)
      FNU(3)=OLAP(L,M+1,N,GA)
      FNU(4)=OLAP(L,M,N+1,GA)
      V(1)=FNU(2)+D(1)*FNU(1)
      V(2)=FNU(3)+D(2)*FNU(1)
      V(3)=FNU(4)+D(3)*FNU(1)
      DO 300 I=1,NT
300   V(I)=-V(I)
      RETURN
CRAY 1
      ENTRY OPAB2(L,M,N,GA,V,NT,D)
CDC 1
C      ENTRY OPAB2
      NT = 10
      FNU(1)=OLAP(L,M,N,GA)
      FNU(2)=OLAP(L+1,M,N,GA)
      FNU(3)=OLAP(L,M+1,N,GA)
      FNU(4)=OLAP(L,M,N+1,GA)
      FNU(5)=OLAP(L+2,M,N,GA)
      FNU(6)=OLAP(L,M+2,N,GA)
      FNU(7)=OLAP(L,M,N+2,GA)
      FNU(8)=OLAP(L+1,M+1,N,GA)
      FNU(9)=OLAP(L+1,M,N+1,GA)
      FNU(10)=OLAP(L,M+1,N+1,GA)
      V(1)=FNU(5)+2.*D(1)*FNU(2)+D(1)**2*FNU(1)
      V(2)=FNU(6)+2.*D(2)*FNU(3)+D(2)**2*FNU(1)
      V(3)=FNU(7)+2.*D(3)*FNU(4)+D(3)**2*FNU(1)
      V(4)=(FNU(8)+D(1)*FNU(3)+D(2)*FNU(2)+D(1)*D(2)*FNU(1))*1.5
      V(5)=(FNU(9)+D(1)*FNU(4)+ D(3)*FNU(2)+D(1)*D(3)*FNU(1))*1.5
      V(6)=(FNU(10)+D(2)*FNU(4)+D(3)*FNU(3)+D(2)*D(3)*FNU(1))*1.5
      V(7)=V(1)+V(2)+V(3)
      V(8) = V(1)
      V(9) = V(2)
      V(10) = V(3)
      V(1)=(3.*V(1)-V(7))/2.
      V(2)=(3.*V(2)-V(7))/2.
      V(3)=(3.*V(3)-V(7))/2.
      DO 301 I=1,NT
301   V(I)=-V(I)
      RETURN
CRAY 1
      ENTRY OPAB3(L,M,N,GA,V,NT,D)
CDC 1
C      ENTRY OPAB3
      NT=10
      FNU(1)=OLAP( L,M,N,GA)
      FNU(2)=OLAP(L+1,M,N,GA)
      FNU(3)=OLAP(L,M+1,N,GA)
      FNU(4)=OLAP(L,M,N+1,GA)
      FNU(5)=OLAP(L+2,M,N,GA)
      FNU(6)=OLAP(L,M+2,N,GA)
      FNU(7)=OLAP(L,M,N+2,GA)
      FNU(8)=OLAP(L+1,M+1,N,GA)
      FNU(9)=OLAP(L+1,M,N+1,GA)
      FNU(10)=OLAP(L,M+1,N+1,GA)
      FNU(11)=OLAP(L+3,M,N,GA)
      FNU(12)=OLAP(L,M+3,N,GA)
      FNU(13)=OLAP(L,M,N+3,GA)
      FNU(14)=OLAP(L+2,M+1,N,GA)
      FNU(15)=OLAP(L+2,M,N+1,GA)
      FNU(16)=OLAP(L+1,M+2,N,GA)
      FNU(17)=OLAP(L,M+2,N+1,GA)
      FNU(18)=OLAP(L+1,M,N+2,GA)
      FNU(19)=OLAP(L,M+1,N+2,GA)
      FNU(20)=OLAP(L+1,M+1,N+1,GA)
      DP(2)=D(1)
      DP(3)=D(2)
      DP(4)=D(3)
      DP(5)=D(1)**2
      DP(6)=D(2)**2
      DP(7)=D(3)**2
      DP(8)=DP(2)*DP(3)
      DP(9)=DP(2)*DP(4)
      DP(10)=DP(3)*DP(4)
      DP(11)=DP(5)*DP(2)
      DP(12)=DP(6)*DP(3)
      DP(13)=DP(7)*DP(4)
      DP(14)=DP(5)*DP(3)
      DP(15)=DP(5)*DP(4)
      DP(16)=DP(6)*DP(2)
      DP(17)=DP(6)*DP(4)
      DP(18)=DP(7)*DP(2)
      DP(19)=DP(7)*DP(3)
      DP(20)=DP(8)*DP(4)
      VS(1)=FNU(11)+3.*DP(2)*FNU(5)+3.*DP(5)*FNU(2)+DP(11)*FNU(1)
      VS(2)=FNU(12)+3.*DP(3)*FNU(6)+3.*DP(6)*FNU(3)+DP(12)*FNU(1)
      VS(3)=FNU(13)+3.*DP(4)*FNU(7)+3.*DP(7)*FNU(4)+DP(13)*FNU(1)
      VS(4)=FNU(14)+2.*DP(2)*FNU(8)+DP(3)*FNU(5)+DP(5)*FNU(3)+
     1 2.*DP(8)*FNU(2)+DP(14)*FNU(1)
      VS(5)=FNU(15)+2.*DP(2)*FNU(9)+DP(4)*FNU(5)+DP(5)*FNU(4)+
     1 2.*DP(9)*FNU(2)+DP(15)*FNU(1)
      VS(6)=FNU(16)+2.*DP(3)*FNU(8)+DP(2)*FNU(6)+DP(6)*FNU(2)+
     1 2.*DP(8)*FNU(3)+DP(16)*FNU(1)
      VS(7)=FNU(17)+2.*DP(3)*FNU(10)+DP(4)*FNU(6)+DP(6)*FNU(4)+
     1 2.0*DP(10)*FNU(3)+DP(17)*FNU(1)
      VS(8)=FNU(18)+2.*DP(4)*FNU(9)+DP(2)*FNU(7)+DP(7)*FNU(2)+
     1 2.*DP(9)*FNU(4)+DP(18)*FNU(1)
      VS(9)=FNU(19)+2.*DP(4)*FNU(10)+DP(3)*FNU(7)+DP(7)*FNU(3)+
     1 2.*DP(10)*FNU(4)+DP(19)*FNU(1)
      VS(10)=FNU(20)+DP(2)*FNU(10)+DP(3)*FNU(9)+DP(4)*FNU(8)+
     1 DP(10)*FNU(2)+DP(9)*FNU(3)+
     2 DP(8)*FNU(4)+DP(20)*FNU(1)
C
C....    NOW TRANSFORM FROM 3RD MOMENTS TO OCTOPOLE MOMENTS
C
      V(1) = VS(1) - 1.5*(VS(6) + VS(8))
      V(2) = VS(2) - 1.5*(VS(4) + VS(9))
      V(3) = VS(3) - 1.5*(VS(5) + VS(7))
      V(4) = 2.*VS(4) - 0.5*(VS(2) + VS(9))
      V(5) = 2.*VS(5) - 0.5*(VS(3) + VS(7))
      V(6) = 2.*VS(6) - 0.5*(VS(1) + VS(8))
      V(7) = 2.*VS(7) - 0.5*(VS(3) + VS(5))
      V(8) = 2.*VS(8) - 0.5*(VS(1) + VS(6))
      V(9) = 2.*VS(9) - 0.5*(VS(2) + VS(4))
      V(10) = VS(10)
      DO 303 I=1,NT
303   V(I)=-V(I)
      RETURN
CRAY 1
      ENTRY OPAB4(L,M,N,GA,V,NT,D)
CDC 1
C      ENTRY OPAB4
      NT=6
      FNU(1)=OLAP(L,M,N,GA)
      FNU(2)=OLAP(L+1,M,N,GA)
      FNU(3)=OLAP(L,M+1,N,GA)
      FNU(4)=OLAP(L,M,N+1,GA)
      FNU(5)=OLAP(L+2,M,N,GA)
      FNU(6)=OLAP(L,M+2,N,GA)
      FNU(7)=OLAP(L,M,N+2,GA)
      FNU(8)=OLAP(L+1,M+1,N,GA)
      FNU(9)=OLAP(L+1,M,N+1,GA)
      FNU(10)=OLAP(L,M+1,N+1,GA)
      VS(1)=FNU(5)+2.*D(1)*FNU(2)+D(1)**2*FNU(1)
      VS(2)=FNU(6)+2.*D(2)*FNU(3)+D(2)**2*FNU(1)
      VS(3)=FNU(7)+2.*D(3)*FNU(4)+D(3)**2*FNU(1)
      VS(4)=FNU(8)+D(1)*FNU(3)+D(2)*FNU(2)+D(1)*D(2)*FNU(1)
      VS(5)=FNU(9)+D(1)*FNU(4)+D(3)*FNU(2)+D(1)*D(3)*FNU(1)
      VS(6)=FNU(10)+D(2)*FNU(4)+D(3)*FNU(3)+D(2)*D(3)*FNU(1)
C
C....    NOW TRANSFORM FROM 2ND MOMENT TO DIAMAGNETIC SUSCEPTIBILITIES
C
      V(1) = 1.5*(VS(2) + VS(3))
      V(2) = 1.5*(VS(1) + VS(3))
      V(3) = 1.5*(VS(1) + VS(2))
      V(4) = -1.5*VS(4)
      V(5) = -1.5*VS(5)
      V(6) = -1.5*VS(6)
      DO 304 I=1,NT
304   V(I)=-V(I)
      RETURN
CRAY 1
      ENTRY OPAB5(L,M,N,GA,V,NT,D)
CDC 1
C      ENTRY OPAB5
      NT=6
      FNU(1)=OLAP(L,M,N,GA)
      FNU(2)=OLAP(L+1,M,N,GA)
      FNU(3)=OLAP(L,M+1,N,GA)
      FNU(4)=OLAP(L,M,N+1,GA)
      FNU(5)=OLAP(L+2,M,N,GA)
      FNU(6)=OLAP(L,M+2,N,GA)
      FNU(7)=OLAP(L,M,N+2,GA)
      FNU(8)=OLAP(L+1,M+1,N,GA)
      FNU(9)=OLAP(L+1,M,N+1,GA)
      FNU(10)=OLAP(L,M+1,N+1,GA)
      V(1)=FNU(5)+2.*D(1)*FNU(2)+D(1)**2*FNU(1)
      V(2)=FNU(6)+2.*D(2)*FNU(3)+D(2)**2*FNU(1)
      V(3)=FNU(7)+2.*D(3)*FNU(4)+D(3)**2*FNU(1)
      V(4)=FNU(8)+D(1)*FNU(3)+D(2)*FNU(2)+D(1)*D(2)*FNU(1)
      V(5)=FNU(9)+D(1)*FNU(4)+D(3)*FNU(2)+D(1)*D(3)*FNU(1)
      V(6)=FNU(10)+D(2)*FNU(4)+D(3)*FNU(3)+D(2)*D(3)*FNU(1)
c      DO 504 I=1,NT
c504   V(I)=-V(I)
      RETURN
CRAY 1
      ENTRY OPAB6(L,M,N,GA,V,NT,D)
CDC 1
C      ENTRY OPAB6
      NT=10
      FNU(1)=OLAP( L,M,N,GA)
      FNU(2)=OLAP(L+1,M,N,GA)
      FNU(3)=OLAP(L,M+1,N,GA)
      FNU(4)=OLAP(L,M,N+1,GA)
      FNU(5)=OLAP(L+2,M,N,GA)
      FNU(6)=OLAP(L,M+2,N,GA)
      FNU(7)=OLAP(L,M,N+2,GA)
      FNU(8)=OLAP(L+1,M+1,N,GA)
      FNU(9)=OLAP(L+1,M,N+1,GA)
      FNU(10)=OLAP(L,M+1,N+1,GA)
      FNU(11)=OLAP(L+3,M,N,GA)
      FNU(12)=OLAP(L,M+3,N,GA)
      FNU(13)=OLAP(L,M,N+3,GA)
      FNU(14)=OLAP(L+2,M+1,N,GA)
      FNU(15)=OLAP(L+2,M,N+1,GA)
      FNU(16)=OLAP(L+1,M+2,N,GA)
      FNU(17)=OLAP(L,M+2,N+1,GA)
      FNU(18)=OLAP(L+1,M,N+2,GA)
      FNU(19)=OLAP(L,M+1,N+2,GA)
      FNU(20)=OLAP(L+1,M+1,N+1,GA)
      DP(2)=D(1)
      DP(3)=D(2)
      DP(4)=D(3)
      DP(5)=D(1)**2
      DP(6)=D(2)**2
      DP(7)=D(3)**2
      DP(8)=DP(2)*DP(3)
      DP(9)=DP(2)*DP(4)
      DP(10)=DP(3)*DP(4)
      DP(11)=DP(5)*DP(2)
      DP(12)=DP(6)*DP(3)
      DP(13)=DP(7)*DP(4)
      DP(14)=DP(5)*DP(3)
      DP(15)=DP(5)*DP(4)
      DP(16)=DP(6)*DP(2)
      DP(17)=DP(6)*DP(4)
      DP(18)=DP(7)*DP(2)
      DP(19)=DP(7)*DP(3)
      DP(20)=DP(8)*DP(4)
      V(1)=FNU(11)+3.*DP(2)*FNU(5)+3.*DP(5)*FNU(2)+DP(11)*FNU(1)
      V(2)=FNU(12)+3.*DP(3)*FNU(6)+3.*DP(6)*FNU(3)+DP(12)*FNU(1)
      V(3)=FNU(13)+3.*DP(4)*FNU(7)+3.*DP(7)*FNU(4)+DP(13)*FNU(1)
      V(4)=FNU(14)+2.*DP(2)*FNU(8)+DP(3)*FNU(5)+DP(5)*FNU(3)+
     1 2.*DP(8)*FNU(2)+DP(14)*FNU(1)
      V(5)=FNU(15)+2.*DP(2)*FNU(9)+DP(4)*FNU(5)+DP(5)*FNU(4)+
     1 2.*DP(9)*FNU(2)+DP(15)*FNU(1)
      V(6)=FNU(16)+2.*DP(3)*FNU(8)+DP(2)*FNU(6)+DP(6)*FNU(2)+
     1 2.*DP(8)*FNU(3)+DP(16)*FNU(1)
      V(7)=FNU(17)+2.*DP(3)*FNU(10)+DP(4)*FNU(6)+DP(6)*FNU(4)+
     1 2.0*DP(10)*FNU(3)+DP(17)*FNU(1)
      V(8)=FNU(18)+2.*DP(4)*FNU(9)+DP(2)*FNU(7)+DP(7)*FNU(2)+
     1 2.*DP(9)*FNU(4)+DP(18)*FNU(1)
      V(9)=FNU(19)+2.*DP(4)*FNU(10)+DP(3)*FNU(7)+DP(7)*FNU(3)+
     1 2.*DP(10)*FNU(4)+DP(19)*FNU(1)
      V(10)=FNU(20)+DP(2)*FNU(10)+DP(3)*FNU(9)+DP(4)*FNU(8)+
     1 DP(10)*FNU(2)+DP(9)*FNU(3)+
     2 DP(8)*FNU(4)+DP(20)*FNU(1)
      DO 503 I=1,NT
503   V(I)=-V(I)
      RETURN
CRAY 1
      ENTRY OPAC1(L,M,N,GA,V,NT,D)
CDC 1
C      ENTRY OPAC1
      NT=3
      DO 105 I=1,3
      FNU(I)=-D(I)
      V(I)=0.
105   CONTINUE
      G=.5/GA
      LH=L/2
      MH=M/2
      NH=N/2
      PRE=2.*G*PI
      IF(2*NH-N) 9,3,9
3     IF(2*MH-M) 7,4,7
4     V(1)=1.
      IF(L) 39,41,39
39    V(1)=FNU(1)**L
41    V(1)=PRE*V(1)*EXP(-GA*FNU(1)**2)*G**(MH+NH)*DFTR(MH+1)*DFTR(NH+1)
      IF(2*LH-L) 50,5,50
5     V(2)=1.
      IF(M) 49,51,49
49    V(2)=FNU(2)**M
51    V(2)=PRE*V(2)*EXP(-GA*FNU(2)**2)*G**(LH+NH)*DFTR(LH+1)*DFTR(NH+1)
6     V(3)=1.
      IF(N) 59,61,59
59    V(3)=FNU(3)**N
61    V(3)=PRE*V(3)*EXP(-GA*FNU(3)**2)*G**(LH+MH)*DFTR(LH+1)*DFTR(MH+1)
      RETURN
7     IF(2*LH-L) 50,8,50
8     V(2)=1.
      IF(M) 79,81,79
79    V(2)=FNU(2)**M
81    V(2)=PRE*V(2)*EXP(-GA*FNU(2)**2)*G**(LH+NH)*DFTR(LH+1)*DFTR(NH+1)
      RETURN
9     IF(2*MH-M) 50,10,50
10    IF(2*LH-L) 50,6,50
50    RETURN
CRAY 1
      ENTRY OPAC2(L,M,N,GA,V,NT,D)
CDC 1
C      ENTRY OPAC2
      NT=3
      DO 201 I=1,3
      DD(I)=-D(I)
      DS(I)=D(I)**2
      V(I)=0.
201   CONTINUE
      G=.5/GA
      PRE = SQRPI*SQRT(G+G)
      LH=L/2
      MH=M/2
      NH=N/2
      IF(2*NH-N) 204,203,204
203   V(1)=1.
      IF(L) 228,229,228
228   V(1)=DD(1)**L
229   IF(M) 230,231,230
230   V(1)=V(1)*DD(2)**M
231   V(1)=PRE*G**NH*V(1)*EXP (-GA*(DS(1)+DS(2)))*DFTR(NH+1)
204   IF(2*MH-M) 206,205,206
205   V(2)=1.
      IF(L) 248,249,248
248   V(2)=DD(1)**L
249   IF(N) 250,251,250
250   V(2) = V(2)*DD(1)**N
251   V(2)=PRE*G**MH*V(2)*EXP (-GA*(DS(1)+DS(3)))*DFTR(MH+1)
206   IF(2*LH-L) 208,207,208
207   V(3)=1.
      IF(M) 268,269,268
268   V(3)=DD(2)**M
269   IF(N) 270,271,270
270   V(3)=V(3)*DD(3)**N
271   V(3)=PRE*G**LH*V(3)*EXP (-GA*(DS(2)+DS(3)))*DFTR(LH+1 )
208   RETURN
CRAY 1
      ENTRY OPAC3(L,M,N,GA,V,NT,D)
CDC 1
C      ENTRY OPAC3
      NT=1
      DDSQ=0.
      DO 111 I=1,3
      DD(I)=-D(I)
      DDSQ=DDSQ+D(I)**2
111   CONTINUE
      V(1) =1.
      IF(L) 112,113,112
112   V(1)=DD(1)**L
113   IF(M) 114,115,114
114   V(1)=V(1)*DD(2)**M
115   IF(N) 116,117,116
116   V(1)=V(1)*DD(3)**N
117   V(1) = V(1)*EXP(-GA*DDSQ)
      RETURN
CRAY 1
      ENTRY OPAD1(L,M,N,GA,V,NT,D)
CDC 1
C      ENTRY OPAD1
      NT=1
      V(1)=OLAP(L,M,N,GA)
      RETURN
      END
