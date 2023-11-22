C************************************************************
      SUBROUTINE ker_corr_vwn(DENS3,DVCAA,DVCBB)
c NEEDS MODIFICATIONS!!!!!!!!!!
c ONLY OPEN SHELL!!!!!
c PROBLEMS!!!!!!


CJC
C     Calculation of V_c (VWN) for each point of the grid
C
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION CAP(4),CEP(4)
C
c   
c     CAP --> g pour alpha
c     CEP --> g pour epsilon paramagnetic
c    
      DATA CAP/-0.0337737,1.13107,13.0045,-0.0047584/
     &    ,CEP/ 0.0621814,3.72744,12.9352,-0.1049800/
      DATA THIRD,PI43/0.33333333333333,4.1887901/
      DATA ONE,TWO,FOUR,SIX/1.0,2.0,4.0,6.0/
C*************Definition of statement functions****************
c --------- g(x)
      GVWN(X,A,B,C,X0)=A*LOG(X*X/(X-X0)**2)-A/(X0*X0+B*X0+C)*
     & ((X0*X0+C)*LOG((X*X+B*X+C)/(X-X0)**2)+
     & (X0*X0-C)*ATAN(SQRT(FOUR*C-B*B)/(TWO*X+B))*TWO*B/
     & SQRT(FOUR*C-B*B))
c --------- dg /dx
      DGVWN(X,A,B,C,X0)=TWO*A*(ONE/X-(ONE+B/(X-X0))*X/(X*X+B*X+C))
c --------- d**2 g / d x**2
      XVWN(X,A,B,C,X0)=-(TWO*A*(ONE/(X*X)-(ONE/(X*X+B*X+C)*
     & ((B*X)/((X-X0)**2)+(ONE+(B/(X-X0)))*((C-X*X)/(X*X+B*X+C))))))
C*************** End statement functions **********************

      call callstack_push('ker_corr_vwn')
      RS=ONE/(PI43*DENS3)**THIRD
      RX=SQRT(RS)
c ==========
      DVC=1.0d0-(1/SIX)
      DVC=DVC*DGVWN(RX,CEP(1),CEP(2),CEP(3),CEP(4))/TWO
      DVC=DVC-(RX/SIX)*XVWN(RX,CEP(1),CEP(2),CEP(3),CEP(4))/TWO
      DVC=-(RX/SIX)*DVC
c     spin Dzeta = g_alpha/2
      Dzeta=GVWN(RX,CAP(1),CAP(2),CAP(3),CAP(4))/TWO
      DVCAA=(DVC+Dzeta)/DENS3
      DVCBB=DVCAA
      DVCAB=(DVC-Dzeta)/DENS3
      dvcaa=dvcaa+dvcab
      dvcbb=dvcbb+dvcab
      call callstack_pop
      RETURN
      END
