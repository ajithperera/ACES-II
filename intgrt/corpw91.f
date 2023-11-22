c----------------------------------------------------------------------
c######################################################################
c----------------------------------------------------------------------
      SUBROUTINE CORpw91(RS,ZET,G,EC,ECRS,ECZET,T,UU,VV,WW,H,
     1                   DVCUP,DVCDN,LPOT)
C  pw91 CORRELATION, modified by K. Burke to put all arguments 
c  as variables in calling statement, rather than in common block
c  May, 1996.
C  lpot added by S. Ivanov. It tells if potential is needed - lpot>0
C  INPUT RS: SEITZ RADIUS
C  INPUT ZET: RELATIVE SPIN POLARIZATION
C  INPUT T: ABS(GRAD D)/(D*2.*KS*G)
C  INPUT UU: (GRAD D)*GRAD(ABS(GRAD D))/(D**2 * (2*KS*G)**3)
C  INPUT VV: (LAPLACIAN D)/(D * (2*KS*G)**2)
C  INPUT WW:  (GRAD D)*(GRAD ZET)/(D * (2*KS*G)**2
C  OUTPUT H: NONLOCAL PART OF CORRELATION ENERGY PER ELECTRON
C  OUTPUT DVCUP,DVCDN:  NONLOCAL PARTS OF CORRELATION POTENTIALS
      IMPLICIT REAL*8 (A-H,O-Z)
      parameter(xnu=15.75592D0,cc0=0.004235D0,cx=-0.001667212D0)
      parameter(alf=0.09D0)
      parameter(c1=0.002568D0,c2=0.023266D0,c3=7.389D-6,c4=8.723D0)
      parameter(c5=0.472D0,c6=7.389D-2,a4=100.D0)
      parameter(thrdm=-0.333333333333D0,thrd2=0.666666666667D0)
      BET = XNU*CC0
      DELT = 2.D0*ALF/BET
      G3 = G**3
      G4 = G3*G
      PON = -DELT*EC/(G3*BET)
      B = DELT/(DEXP(PON)-1.D0)
      B2 = B*B
      T2 = T*T
      T4 = T2*T2
      T6 = T4*T2
      RS2 = RS*RS
      RS3 = RS2*RS
      Q4 = 1.D0+B*T2
      Q5 = 1.D0+B*T2+B2*T4
      Q6 = C1+C2*RS+C3*RS2
      Q7 = 1.D0+C4*RS+C5*RS2+C6*RS3
      CC = -CX + Q6/Q7
      R0 = 0.663436444d0*rs
      R1 = A4*R0*G4
      COEFF = CC-CC0-3.D0*CX/7.D0
      R2 = XNU*COEFF*G3
      R3 = DEXP(-R1*T2)
      H0 = G3*(BET/DELT)*DLOG(1.D0+DELT*Q4*T2/Q5)
      H1 = R3*R2*T2
      H = H0+H1
      IF (LPOT.EQ.0) RETURN
C  LOCAL CORRELATION OPTION:
C     H = 0.0D0
C  ENERGY DONE. NOW THE POTENTIAL:
      CCRS = (C2+2.*C3*RS)/Q7 - Q6*(C4+2.*C5*RS+3.*C6*RS2)/Q7**2
      RSTHRD = RS/3.D0
      R4 = RSTHRD*CCRS/COEFF
      GZ = ((1.D0+ZET)**THRDM - (1.D0-ZET)**THRDM)/3.D0
      FAC = DELT/B+1.D0
      BG = -3.D0*B2*EC*FAC/(BET*G4)
      BEC = B2*FAC/(BET*G3)
      Q8 = Q5*Q5+DELT*Q4*Q5*T2
      Q9 = 1.D0+2.D0*B*T2
      H0B = -BET*G3*B*T6*(2.D0+B*T2)/Q8
      H0RS = -RSTHRD*H0B*BEC*ECRS
      FACT0 = 2.D0*DELT-6.D0*B
      FACT1 = Q5*Q9+Q4*Q9*Q9
      H0BT = 2.D0*BET*G3*T4*((Q4*Q5*FACT0-DELT*FACT1)/Q8)/Q8
      H0RST = RSTHRD*T2*H0BT*BEC*ECRS
      H0Z = 3.D0*GZ*H0/G + H0B*(BG*GZ+BEC*ECZET)
      H0T = 2.*BET*G3*Q9/Q8
      H0ZT = 3.D0*GZ*H0T/G+H0BT*(BG*GZ+BEC*ECZET)
      FACT2 = Q4*Q5+B*T2*(Q4*Q9+Q5)
      FACT3 = 2.D0*B*Q5*Q9+DELT*FACT2
      H0TT = 4.D0*BET*G3*T*(2.D0*B/Q8-(Q9*FACT3/Q8)/Q8)
      H1RS = R3*R2*T2*(-R4+R1*T2/3.D0)
      FACT4 = 2.D0-R1*T2
      H1RST = R3*R2*T2*(2.D0*R4*(1.D0-R1*T2)-THRD2*R1*T2*FACT4)
      H1Z = GZ*R3*R2*T2*(3.D0-4.D0*R1*T2)/G
      H1T = 2.D0*R3*R2*(1.D0-R1*T2)
      H1ZT = 2.D0*GZ*R3*R2*(3.D0-11.D0*R1*T2+4.D0*R1*R1*T4)/G
      H1TT = 4.D0*R3*R2*R1*T*(-2.D0+R1*T2)
      HRS = H0RS+H1RS
      HRST = H0RST+H1RST
      HT = H0T+H1T
      HTT = H0TT+H1TT
      HZ = H0Z+H1Z
      HZT = H0ZT+H1ZT
      COMM = H+HRS+HRST+T2*HT/6.D0+7.D0*T2*T*HTT/6.D0
      PREF = HZ-GZ*T2*HT/G
      FACT5 = GZ*(2.D0*HT+T*HTT)/G
      COMM = COMM-PREF*ZET-UU*HTT-VV*HT-WW*(HZT-FACT5)
      DVCUP = COMM + PREF
      DVCDN = COMM - PREF
C  LOCAL CORRELATION OPTION:
C     DVCUP = 0.0D0
C     DVCDN = 0.0D0
      RETURN
      END
