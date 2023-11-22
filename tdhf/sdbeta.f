      SUBROUTINE SDBETA(F,FS,E,ES,UP,UM,US,NSIZ1,NSIZO)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*1 AXYZ 
      COMMON/CONST/ EBETA,EGAMMA 
      DIMENSION UP(NSIZ1,NSIZO,3),UM(NSIZ1,NSIZO,3),US(NSIZ1,NSIZO,3)
     X ,F(NSIZ1,NSIZ1,3),FS(NSIZ1,NSIZ1,3),E(NSIZO,NSIZO,3)
     X ,ES(NSIZO,NSIZO,3)
      DIMENSION ICOMP(9),JCOMP(9),IJCP(3,3),AXYZ(3),IJSP(3,3)
      DIMENSION ENRG(3,9)
      DATA ZERO/0.D00/,TWO/2.D0/,HALF/5.D-01/,FOUR/4.D0/
      DATA ONE/1.D00/,OCTH/1.D-8/,THREE/3.D0/,SIX/6.D0/           
C   ........ (xx,yy,zz,xy,yz,zx,yx,zy,xz) ...........
      DATA ICOMP/ 1,2,3,1,2,3,2,3,1/,
     X     JCOMP/ 1,2,3,2,3,1,1,2,3/,
     X     AXYZ/ 'x','y','z'/                 
      DATA IJCP/ 1,4,6, 4,2,5, 6,5,3 /,IJSP/ 1,7,6, 4,2,8, 9,5,3 /
      WRITE(6,*) ' >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> '
      WRITE(6,*) ' Beta: Electro Optic Pockel Effect and '
      WRITE(6,*) ' Optical Rectification by first order wavefunctions '
      WRITE(6,*) ' >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> '
      DO 10 I=1,3
      DO 10 JK=1,9
      J=ICOMP(JK)
      K=JCOMP(JK)
      UFUPE=ZERO
      UFUKE=ZERO
      DO 11 L=1,NSIZ1
      DO 11 M=1,NSIZ1
      DO 11 N=1,NSIZO
      UFUPE=UFUPE 
     X +UP(L,N,I)*FS(L,M,J)*UP(M,N,K)+UM(L,N,K)*FS(L,M,J)*UM(M,N,I)
     X +US(L,N,J)*F(L,M,K)*UM(M,N,I)+UP(L,N,I)*F(L,M,K)*US(M,N,J)
     X +UM(L,N,K)*F(M,L,I)*US(M,N,J)+US(L,N,J)*F(M,L,I)*UP(M,N,K)
      UFUKE=UFUKE
     X +US(L,N,I)*F(L,M,J)*UM(M,N,K)+UP(L,N,K)*F(L,M,J)*US(M,N,I)
     X +UM(L,N,J)*F(M,L,K)*US(M,N,I)+US(L,N,I)*F(M,L,K)*UP(M,N,J)
     X +UP(L,N,K)*FS(L,M,I)*UP(M,N,J)+UM(L,N,J)*FS(L,M,I)*UM(M,N,K)
   11 CONTINUE
      UUEPE=ZERO
      UUEKE=ZERO
      DO 12 L=1,NSIZ1
      DO 12 M=1,NSIZO
      DO 12 N=1,NSIZO
      UUEPE=UUEPE
     X +UP(L,M,I)*UP(L,N,K)*ES(N,M,J)+UM(L,M,K)*UM(L,N,I)*ES(N,M,J)
     X +US(L,M,J)*UM(L,N,I)*E(N,M,K)+UP(L,M,I)*US(L,N,J)*E(N,M,K)
     X +UM(L,M,K)*US(L,N,J)*E(M,N,I)+US(L,M,J)*UP(L,N,K)*E(M,N,I)
      UUEKE=UUEKE
     X +US(L,M,I)*UM(L,N,K)*E(N,M,J)+UP(L,M,K)*US(L,N,I)*E(N,M,J)
     X +UM(L,M,J)*US(L,N,I)*E(M,N,K)+US(L,M,I)*UP(L,N,J)*E(M,N,K)
     X +UP(L,M,K)*UP(L,N,J)*ES(N,M,I)+UM(L,M,J)*UM(L,N,K)*ES(N,M,I)
   12 CONTINUE
      BETAPE =  -(UFUPE - UUEPE)*TWO
      BETAKE =  -(UFUKE - UUEKE)*TWO
      ENRG(I,JK) = BETAKE
      WRITE(6,100) AXYZ(I),AXYZ(J),AXYZ(K),BETAPE,BETAKE
  100 FORMAT(1H ,3A1,' = ',2F20.7)
   10 CONTINUE
      BXXX=3.D0*ENRG(1,1)+ENRG(1,2)+ENRG(2,4)+ENRG(2,7)
     X                   +ENRG(1,3)+ENRG(3,6)+ENRG(3,9)
      BXXX=BXXX/5.D0
      BYYY=3.D0*ENRG(2,2)+ENRG(2,1)+ENRG(1,4)+ENRG(1,7)
     X                   +ENRG(2,3)+ENRG(3,5)+ENRG(3,8)
      BYYY=BYYY/5.D0
      BZZZ=3.D0*ENRG(3,3)+ENRG(3,1)+ENRG(1,6)+ENRG(1,9)
     X                   +ENRG(3,2)+ENRG(2,5)+ENRG(2,8)
      BZZZ=BZZZ/5.D0
      BXPP=(TWO*(ENRG(1,1)+ENRG(2,7)+ENRG(3,6))
     X         -ENRG(1,1)-ENRG(2,4)-ENRG(3,9))/5.D0
      BYPP=(TWO*(ENRG(1,4)+ENRG(2,2)+ENRG(3,8))
     X         -ENRG(1,7)-ENRG(2,2)-ENRG(3,5))/5.D0
      BZPP=(TWO*(ENRG(1,9)+ENRG(2,5)+ENRG(3,3))
     X         -ENRG(1,6)-ENRG(2,8)-ENRG(3,3))/5.D0
      WRITE(6,*) ' Parallel x,y,z(a.u.) '
      WRITE(6,*)  BXXX,BYYY,BZZZ
      BXESU=BXXX*EBETA
      BYESU=BYYY*EBETA
      BZESU=BZZZ*EBETA
      WRITE(6,*) ' Beta in 10**-32 e.s.u. '
      WRITE(6,*) BXESU,BYESU,BZESU
      WRITE(6,*) ' Perpendicular(a.u.) '
      WRITE(6,*)  BXPP,BYPP,BZPP
      RETURN
      END
