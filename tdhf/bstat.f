      SUBROUTINE BSTAT(F,E,U,NSIZ1,NSIZO)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*1 AXYZ 
      COMMON/CONST/ EBETA,EGAMMA 
      DIMENSION U(NSIZ1,NSIZO,3),F(NSIZ1,NSIZ1,3),E(NSIZO,NSIZO,3)
      DIMENSION ICOMP(6),JCOMP(6),IJCP(3,3),AXYZ(3)
      DIMENSION ENRG(3,6)
      DATA ZERO/0.D00/,TWO/2.D0/,HALF/5.D-01/,FOUR/4.D0/
      DATA ONE/1.D00/,OCTH/1.D-8/,THREE/3.D0/,SIX/6.D0/           
C   ........ (xx,yy,zz,xy,yz,zx) ...........
      DATA ICOMP/ 1,2,3,1,2,3/, JCOMP/ 1,2,3,2,3,1/,IJCP/ 1,4,6,
     X 4,2,5, 6,5,3 /,AXYZ/ 'x','y','z'/                 
      WRITE(6,*) ' >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> '
      WRITE(6,*) ' Static Beta by first order wavefunctions '
      WRITE(6,*) ' >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> '
      DO 10 I=1,3
      DO 10 JK=1,6
      J=ICOMP(JK)
      K=JCOMP(JK)
      UFU=ZERO
      DO 11 L=1,NSIZ1
      DO 11 M=1,NSIZ1
      DO 11 N=1,NSIZO
      UFU=UFU 
     X  +U(L,N,I)*F(L,M,J)*U(M,N,K)+U(L,N,K)*F(L,M,J)*U(M,N,I)
     X  +U(L,N,J)*F(L,M,K)*U(M,N,I)+U(L,N,I)*F(L,M,K)*U(M,N,J)
     X  +U(L,N,K)*F(L,M,I)*U(M,N,J)+U(L,N,J)*F(L,M,I)*U(M,N,K)
   11 CONTINUE
      UUE=ZERO
      DO 12 L=1,NSIZ1
      DO 12 M=1,NSIZO
      DO 12 N=1,NSIZO
      UUE=UUE
     X  +U(L,M,I)*U(L,N,K)*E(N,M,J)+U(L,M,K)*U(L,N,I)*E(N,M,J)
     X  +U(L,M,J)*U(L,N,I)*E(N,M,K)+U(L,M,I)*U(L,N,J)*E(N,M,K)
     X  +U(L,M,K)*U(L,N,J)*E(N,M,I)+U(L,M,J)*U(L,N,K)*E(N,M,I)
   12 CONTINUE
      ENRG(I,JK) =  -(UFU - UUE)*TWO
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
      WRITE(6,*) ' Parallel (x,y,z)  '
      WRITE(6,*) BXXX,BYYY,BZZZ,' a.u. '
      WRITE(6,*) ' By Kleinman symmetry, perpendicular components
     X  are 1/3 of that of Parallel. '
      BXESU=BXXX*EBETA
      BYESU=BYYY*EBETA
      BZESU=BZZZ*EBETA
      WRITE(6,*) ' Beta in 10**-32 e.s.u. '
      WRITE(6,*) BXESU,BYESU,BZESU
      RETURN
      END
