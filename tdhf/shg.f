      SUBROUTINE SHG(F,FW,E,EW,UP,UWP,UM,UWM,NSIZ1,NSIZO)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*1 AXYZ 
      COMMON/CONST/ EBETA,EGAMMA 
      DIMENSION UP(NSIZ1,NSIZO,3),UM(NSIZ1,NSIZO,3),UWP(NSIZ1,NSIZO,3)
     X ,UWM(NSIZ1,NSIZO,3),F(NSIZ1,NSIZ1,3),FW(NSIZ1,NSIZ1,3),
     X  E(NSIZO,NSIZO,3),EW(NSIZO,NSIZO,3)
      DIMENSION ICOMP(6),JCOMP(6),IJCP(3,3),AXYZ(3)
      DIMENSION ENRG(3,6)
      DATA ZERO/0.D00/,TWO/2.D0/,HALF/5.D-01/,FOUR/4.D0/
      DATA ONE/1.D00/,OCTH/1.D-8/,THREE/3.D0/,SIX/6.D0/           
C   ........ (xx,yy,zz,xy,yz,zx) ...........
      DATA ICOMP/ 1,2,3,1,2,3/,
     X     JCOMP/ 1,2,3,2,3,1/,
     X     IJCP/ 1,4,6,4,2,5, 6,5,3 /,
     X     AXYZ/ 'x','y','z'/                 
      WRITE(6,*) ' >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> '
      WRITE(6,*) ' Beta:  Second Harmonic Generation '
      WRITE(6,*) '  by first order wavefunctions '
      WRITE(6,*) ' >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> '
      DO 10 I=1,3
      DO 10 JK=1,6
      J=ICOMP(JK)
      K=JCOMP(JK)
      UFU=ZERO
      DO 11 L=1,NSIZ1
      DO 11 M=1,NSIZ1
      DO 11 N=1,NSIZO
      UFU=UFU
     X +UWP(L,N,I)*F(L,M,J)*UP(M,N,K)+UM(L,N,K)*F(L,M,J)*UWM(M,N,I)
     X +UM(L,N,J)*F(L,M,K)*UWM(M,N,I)+UWP(L,N,I)*F(L,M,K)*UP(M,N,J)
     X +UM(L,N,K)*FW(M,L,I)*UP(M,N,J)+UM(L,N,J)*FW(M,L,I)*UP(M,N,K)
   11 CONTINUE
      UUE=ZERO
      DO 12 L=1,NSIZ1
      DO 12 M=1,NSIZO
      DO 12 N=1,NSIZO
      UUE=UUE
     X +UWP(L,M,I)*UP(L,N,K)*E(N,M,J)+UM(L,M,K)*UWM(L,N,I)*E(N,M,J)
     X +UM(L,M,J)*UWM(L,N,I)*E(N,M,K)+UWP(L,M,I)*UP(L,N,J)*E(N,M,K)
     X +UM(L,M,K)*UP(L,N,J)*EW(M,N,I)+UM(L,M,J)*UP(L,N,K)*EW(M,N,I)
   12 CONTINUE
      ENRG(I,JK)=  -(UFU - UUE)*TWO
      WRITE(6,100) AXYZ(I),AXYZ(J),AXYZ(K),ENRG(I,JK)
  100 FORMAT(1H ,3A1,' = ',F20.7)
   10 CONTINUE
      BXXX=3.D0*ENRG(1,1)+ENRG(1,2)+ENRG(2,4)+ENRG(2,4)
     X                   +ENRG(1,3)+ENRG(3,6)+ENRG(3,6)
      BXXX=BXXX/5.D0
      BYYY=3.D0*ENRG(2,2)+ENRG(2,1)+ENRG(1,4)+ENRG(1,4)
     X                   +ENRG(2,3)+ENRG(3,5)+ENRG(3,5)
      BYYY=BYYY/5.D0
      BZZZ=3.D0*ENRG(3,3)+ENRG(3,1)+ENRG(1,6)+ENRG(1,6)
     X                   +ENRG(3,2)+ENRG(2,5)+ENRG(2,5)
      BZZZ=BZZZ/5.D0
      BZPP=(TWO*(ENRG(3,1)+ENRG(3,2)+ENRG(3,3))
     X         -ENRG(1,6)-ENRG(2,5)-ENRG(3,3))/5.D0
      BXPP=(TWO*(ENRG(1,1)+ENRG(1,2)+ENRG(1,3))
     X         -ENRG(1,1)-ENRG(2,4)-ENRG(3,6))/5.D0
      BYPP=(TWO*(ENRG(2,1)+ENRG(2,2)+ENRG(2,3))
     X         -ENRG(1,4)-ENRG(2,2)-ENRG(3,5))/5.D0
      WRITE(6,*) '  Parallel x,y,z(a.u.) '
      WRITE(6,*) BXXX,BYYY,BZZZ
      BXESU=BXXX*EBETA
      BYESU=BYYY*EBETA
      BZESU=BZZZ*EBETA
      WRITE(6,*) ' Beta in 10**-32 e.s.u. '
      WRITE(6,*) BXESU,BYESU,BZESU
      WRITE(6,*) ' perpendicular (a.u.) '
      WRITE(6,*) BXPP,BYPP,BZPP
      RETURN
      END
