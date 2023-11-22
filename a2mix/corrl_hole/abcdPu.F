      REAL*8 FUNCTION abcdPu(alpha,beta,gamma,delta,AL,BL,CL,DL,u,m)
      INTEGER AL(3),BL(3),CL(3),DL(3),m,ALtmp(3),BLtmp(3),CLtmp(3),DLtmp(3)
      INTEGER ALLoc(3),BLLoc(3),CLLoc(3),DLLoc(3)
      INTEGER xi,i,Xang,Yang,Zang
      REAL*8 alpha,beta,gamma,delta,ssssPu,zeta,eta,u
      EXTERNAL ssssPu

C  First check if there is anything to be calculated or can we eliminate 
C  integral by symmetry
      Xang=AL(1)+BL(1)+CL(1)+DL(1)
      Yang=AL(2)+BL(2)+CL(2)+DL(2)
      Zang=AL(3)+BL(3)+CL(3)+DL(3)

      IF(BTEST(Xang,0).OR.BTEST(Yang,0).OR.BTEST(Zang,0)) THEN
         abcdPu=0.0d0
         RETURN
      ENDIF

      IF(AL(1).LT.0.OR.AL(2).LT.0.OR.AL(3).LT.0.OR.BL(1).LT.0.OR.BL(2).LT.0
     &    .OR.BL(3).LT.0.OR.CL(1).LT.0.OR.CL(2).LT.0.OR.CL(3).LT.0
     &     .OR.DL(1).LT.0.OR.DL(2).LT.0.OR.DL(3).LT.0) THEN
         abcdPu=0.0d0
         RETURN
      ENDIF

      zeta=alpha+beta
      eta=gamma+delta
         
C  Need to have some local versions otherwise get overwritten in recursion
      DO i=1,3
         ALLoc(i)=AL(i)
         BLLoc(i)=BL(i)
         CLLoc(i)=CL(i)
         DLLoc(i)=DL(i)
         ALtmp(i)=AL(i)
         BLtmp(i)=BL(i)
         CLtmp(i)=CL(i)
         DLtmp(i)=DL(i)
      ENDDO

C  Start by figuring out which center is to worked on
      IF(DL(1).GT.0.OR.DL(2).GT.0.OR.DL(3).GT.0) THEN
         IF(DL(1).GT.0) THEN
            xi=1
            GOTO 10
         ELSEIF(DL(2).GT.0) THEN
            xi=2
            GOTO 10
         ELSE 
            xi=3
            GOTO 10
         ENDIF
10       DLLoc(xi)=DL(xi)-1
         DLtmp(xi)=DLLoc(xi)-1
         ALtmp(xi)=AL(xi)-1
         BLtmp(xi)=BL(xi)-1
         CLtmp(xi)=CL(xi)-1
         abcdPu=
     &       DLLoc(xi)*((1d0/(2d0*eta))*abcdPu(alpha,beta,gamma,delta,AL,BL,CL,DLtmp,u,m)-
     &       (1d0/(2d0*eta))**2*abcdPu(alpha,beta,gamma,delta,AL,BL,CL,DLtmp,u,m+1))+
     &       CL(xi)*((1d0/(2d0*eta))*abcdPu(alpha,beta,gamma,delta,AL,BL,CLtmp,DLLoc,u,m)-
     &       (1d0/(2d0*eta))**2*abcdPu(alpha,beta,gamma,delta,AL,BL,CLtmp,DLLoc,u,m+1))+
     &       AL(xi)*(1d0/(2d0*zeta))*(1d0/(2d0*eta))*abcdPu(alpha,beta,gamma,delta,ALtmp,BL,CL,DLLoc,u,m+1)+
     &       BL(xi)*(1d0/(2d0*zeta))*(1d0/(2d0*eta))*abcdPu(alpha,beta,gamma,delta,AL,BLtmp,CL,DLLoc,u,m+1)
         RETURN
      ELSEIF(CL(1).GT.0.OR.CL(2).GT.0.OR.CL(3).GT.0) THEN
         IF(CL(1).GT.0) THEN
            xi=1
            GOTO 11
         ELSEIF(CL(2).GT.0) THEN
            xi=2
            GOTO 11
         ELSE 
            xi=3
            GOTO 11
         ENDIF
11       CLLoc(xi)=CL(xi)-1
         CLtmp(xi)=CLLoc(xi)-1
         ALtmp(xi)=AL(xi)-1
         BLtmp(xi)=BL(xi)-1
         DLtmp(xi)=DL(xi)-1
         abcdPu=
     &       CLLoc(xi)*((1d0/(2d0*eta))*abcdPu(alpha,beta,gamma,delta,AL,BL,CLtmp,DL,u,m)-
     &       (1d0/(2d0*eta))**2*abcdPu(alpha,beta,gamma,delta,AL,BL,CLtmp,DL,u,m+1))+
     &       DL(xi)*((1d0/(2d0*eta))*abcdPu(alpha,beta,gamma,delta,AL,BL,CLLoc,DLtmp,u,m)-
     &       (1d0/(2d0*eta))**2*abcdPu(alpha,beta,gamma,delta,AL,BL,CLLoc,DLtmp,u,m+1))+
     &       AL(xi)*(1d0/(2d0*zeta))*(1d0/(2d0*eta))*abcdPu(alpha,beta,gamma,delta,ALtmp,BL,CLLoc,DL,u,m+1)+
     &       BL(xi)*(1d0/(2d0*zeta))*(1d0/(2d0*eta))*abcdPu(alpha,beta,gamma,delta,AL,BLtmp,CLLoc,DL,u,m+1)
         RETURN
      ELSEIF(BL(1).GT.0.OR.BL(2).GT.0.OR.BL(3).GT.0) THEN
         IF(BL(1).GT.0) THEN
            xi=1
            GOTO 12
         ELSEIF(BL(2).GT.0) THEN
            xi=2
            GOTO 12
         ELSE 
            xi=3
            GOTO 12
         ENDIF
12       BLLoc(xi)=BL(xi)-1
         BLtmp(xi)=BLLoc(xi)-1
         ALtmp(xi)=AL(xi)-1
         CLtmp(xi)=CL(xi)-1
         DLtmp(xi)=DL(xi)-1
         abcdPu=
     &       BLLoc(xi)*((1d0/(2d0*zeta))*abcdPu(alpha,beta,gamma,delta,AL,BLtmp,CL,DL,u,m)-
     &       (1d0/(2d0*zeta))**2*abcdPu(alpha,beta,gamma,delta,AL,BLtmp,CL,DL,u,m+1))+
     &       AL(xi)*((1d0/(2d0*zeta))*abcdPu(alpha,beta,gamma,delta,ALtmp,BLLoc,CL,DL,u,m)-
     &       (1d0/(2d0*zeta))**2*abcdPu(alpha,beta,gamma,delta,ALtmp,BLLoc,CL,DL,u,m+1))+
     &       CL(xi)*(1d0/(2d0*zeta))*(1d0/(2d0*eta))*abcdPu(alpha,beta,gamma,delta,AL,BLLoc,CLtmp,DL,u,m+1)+
     &       DL(xi)*(1d0/(2d0*zeta))*(1d0/(2d0*eta))*abcdPu(alpha,beta,gamma,delta,AL,BLLoc,CL,DLtmp,u,m+1)
         RETURN
      ELSEIF(AL(1).GT.0.OR.AL(2).GT.0.OR.AL(3).GT.0) THEN
         IF(AL(1).GT.0) THEN
            xi=1
            GOTO 13
         ELSEIF(AL(2).GT.0) THEN
            xi=2
            GOTO 13
         ELSE 
            xi=3
            GOTO 13
         ENDIF
13       ALLoc(xi)=AL(xi)-1
         ALtmp(xi)=ALLoc(xi)-1
         BLtmp(xi)=BL(xi)-1
         CLtmp(xi)=CL(xi)-1
         DLtmp(xi)=DL(xi)-1
         abcdPu=
     &       ALLoc(xi)*((1d0/(2d0*zeta))*abcdPu(alpha,beta,gamma,delta,ALtmp,BL,CL,DL,u,m)-
     &       (1d0/(2d0*zeta))**2*abcdPu(alpha,beta,gamma,delta,ALtmp,BL,CL,DL,u,m+1))+
     &       BL(xi)*((1d0/(2d0*zeta))*abcdPu(alpha,beta,gamma,delta,ALLoc,BLtmp,CL,DL,u,m)-
     &       (1d0/(2d0*zeta))**2*abcdPu(alpha,beta,gamma,delta,ALLoc,BLtmp,CL,DL,u,m+1))+
     &       CL(xi)*(1d0/(2d0*zeta))*(1d0/(2d0*eta))*abcdPu(alpha,beta,gamma,delta,ALLoc,BL,CLtmp,DL,u,m+1)+
     &       DL(xi)*(1d0/(2d0*zeta))*(1d0/(2d0*eta))*abcdPu(alpha,beta,gamma,delta,ALLoc,BL,CL,DLtmp,u,m+1)
         RETURN
      ELSE
         abcdPu=ssssPu(alpha,beta,gamma,delta,u,m)
         RETURN
      ENDIF
      END
