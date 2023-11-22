      SUBROUTINE GETPRODUCT(LI, MI, NI, LJ, MJ, NJ, MAXPRM, EXP1,
     &                      EXP2,  CENTER, CNTMU, CNTNU, PRDUCT)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL LIMINUS, LJMINUS
C
      DIMENSION  CENTER(3), CNTMU(3), CNTNU(3)
C
      COMMON /HIGHL/ LMNVAL(3,84), ANORM(84)
C
       LJMINUS = (LJ .LT. 0 .OR. MJ .LT. 0 .OR. NJ .LT. 0) 
       LIMINUS = (LI .LT. 0 .OR. MI .LT. 0 .OR. NI .LT. 0) 

       IF (LJMINUS .OR. LIMINUS) THEN
C     
       PRDUCT = 0.0D0
C
      ELSE
C
         DISTNCEI=((CENTER(1)-CNTMU(1))**2 + (CENTER(2)-CNTMU(2))**2 +
     &            (CENTER(3)-CNTMU(3))**2)
         DISTNCEJ=((CENTER(1)-CNTNU(1))**2 + (CENTER(2)-CNTNU(2))**2 +
     &            (CENTER(3)-CNTNU(3))**2)
      
         PREFCTI=((CENTER(1)-CNTMU(1))**LI)*((CENTER(2)-CNTMU(2))**MI)*
     &             ((CENTER(3)-CNTMU(3))**NI)
         PREFCTJ=((CENTER(1)-CNTNU(1))**LJ)*((CENTER(2)-CNTNU(2))**MJ)*
     &             ((CENTER(3)-CNTNU(3))**NJ)
C
         EXPFCTI = DEXP(-1.0D0*EXP1*DISTNCEI)
         EXPFCTJ = DEXP(-1.0D0*EXP2*DISTNCEJ)

         PRDUCT = PREFCTI*EXPFCTI*PREFCTJ*EXPFCTJ
      ENDIF
C
      RETURN
      END



























