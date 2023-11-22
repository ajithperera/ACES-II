C                                                                       A  10720
      SUBROUTINE DCSHG(F1,F2,FS,FS2,FW,E1,E2,ES,ES2,EW,
     X UP1,UP2,UM1,UM2,US,USP2,USM2,UWP,UWM,
     X LCOMP,NCOMP,NSIZ1,NSIZ3,NSIZO)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*1 AXYZ 
      COMMON/INFOA/NBASIS,NUMSCF,NX,NMO2,NOC,NVT,NVO                    
      COMMON/CONST/ EBETA,EGAMMA 
      DIMENSION UP1(NSIZ1,NSIZO,3),UP2(NSIZ1,NSIZO,6),UM1(NSIZ1,NSIZO,3)
     X ,UM2(NSIZ1,NSIZO,6),UWP(NSIZ1,NSIZO,3),UWM(NSIZ1,NSIZO,3)
     X ,US(NSIZ1,NSIZO,3),USP2(NSIZ1,NSIZO,9),USM2(NSIZ1,NSIZO,9)
      DIMENSION F1(NSIZ1,NSIZ1,3),F2(NSIZ1,NSIZ1,6),E1(NSIZO,NSIZO,3)
     X ,E2(NSIZO,NSIZO,6),FW(NSIZ1,NSIZ1,3),EW(NSIZO,NSIZO,3)
     X ,FS(NSIZ1,NSIZ1,3),FS2(NSIZ1,NSIZ1,9),ES(NSIZO,NSIZO,3)
     X ,ES2(NSIZO,NSIZO,9)
      DIMENSION LCOMP(21,4),IJCP(3,3),IJSP(3,3),AXYZ(3)
      DATA ZERO/0.D00/,TWO/2.D0/,HALF/5.D-01/,FOUR/4.D0/
      DATA ONE/1.D00/,OCTH/1.D-8/,THREE/3.D0/,SIX/6.D0/           
      DATA IJCP/ 1,4,6, 4,2,5, 6,5,3 /,IJSP/ 1,7,6, 4,2,8, 9,5,3 /
      DATA AXYZ/'x','y','z'/
      WRITE(6,*) ' Now in DC-SHG' 
      WRITE(6,*) ' I J K L and GANMA  '
C     DO 100 MCOM=1,3
C     write(6,*) ' E1',MCOM
C     CALL OUTMXD(E1(1,1,MCOM),NSIZO,NOC,NOC)
C     write(6,*) ' EW',MCOM
C     CALL OUTMXD(EW(1,1,MCOM),NSIZO,NOC,NOC)
C 100 CONTINUE
      GPARA = ZERO
      GPERP = ZERO
      DO 50 LCOM=1,NCOMP
      U1U2F1=ZERO
      U1U2F2=ZERO
      U1U2E1=ZERO
      U1U2E2=ZERO
      U1U1F=ZERO
      U1U1E=ZERO
      L1 = LCOMP(LCOM,1)
      L2 = LCOMP(LCOM,2)
      L3 = LCOMP(LCOM,3)
      L4 = LCOMP(LCOM,4)
      DO 10 I=1,NUMSCF               
      DO 10 J=1,NUMSCF                   
      DO 10 K=1,NOC
      U1U2F1 = U1U2F1 + UWP(I,K,L1)*FS(I,J,L2)*UP2(J,K,IJCP(L3,L4))
     X + UM2(I,K,IJCP(L3,L4))*FS(I,J,L2)*UWM(J,K,L1)
     X + UWP(I,K,L1)*F1(I,J,L3)*USP2(J,K,IJSP(L2,L4))
     X + USM2(I,K,IJSP(L2,L4))*F1(I,J,L3)*UWM(J,K,L1)
     X + UWP(I,K,L1)*F1(I,J,L4)*USP2(J,K,IJSP(L2,L3))
     X + USM2(I,K,IJSP(L2,L3))*F1(I,J,L4)*UWM(J,K,L1)
C
      U1U2F2 = U1U2F2 +  US(I,K,L2)*FW(J,I,L1)*UP2(J,K,IJCP(L3,L4))
     X + UM2(I,K,IJCP(L3,L4))*FW(J,I,L1)*US(J,K,L2)
     X + UM1(I,K,L3)*FW(J,I,L1)*USP2(J,K,IJSP(L2,L4))
     X + USM2(I,K,IJSP(L2,L4))*FW(J,I,L1)*UP1(J,K,L3)
     X + UM1(I,K,L4)*FW(J,I,L1)*USP2(J,K,IJSP(L2,L3))
     X + USM2(I,K,IJSP(L2,L3))*FW(J,I,L1)*UP1(J,K,L4)
C
      U1U1F = U1U1F + UWP(I,K,L1)*F2(I,J,IJCP(L3,L4))*US(J,K,L2)
     X + US(I,K,L2)*F2(I,J,IJCP(L3,L4))*UWM(J,K,L1)
     X + UWP(I,K,L1)*FS2(I,J,IJSP(L2,L4))*UP1(J,K,L3)
     X + UM1(I,K,L3)*FS2(I,J,IJSP(L2,L4))*UWM(J,K,L1)
     X + UWP(I,K,L1)*FS2(I,J,IJSP(L2,L3))*UP1(J,K,L4)
     X + UM1(I,K,L4)*FS2(I,J,IJSP(L2,L3))*UWM(J,K,L1)
C
   10 CONTINUE
      DO 11 K=1,NOC
      DO 11 L=1,NOC
      DO 11 I=1,NUMSCF               
      U1U2E1 = U1U2E1 - UWP(I,K,L1)*UP2(I,L,IJCP(L3,L4))*ES(L,K,L2)
     X  - UM2(I,K,IJCP(L3,L4))*UWM(I,L,L1)*ES(L,K,L2)
     X  - UWP(I,K,L1)*USP2(I,L,IJSP(L2,L4))*E1(L,K,L3)
     X  - USM2(I,K,IJSP(L2,L4))*UWM(I,L,L1)*E1(L,K,L3)
     X  - UWP(I,K,L1)*USP2(I,L,IJSP(L2,L3))*E1(L,K,L4)
     X  - USM2(I,K,IJSP(L2,L3))*UWM(I,L,L1)*E1(L,K,L4)
C
      U1U2E2 = U1U2E2 - US(I,K,L2)*UP2(I,L,IJCP(L3,L4))*EW(K,L,L1)
     X  - UM2(I,K,IJCP(L3,L4))*US(I,L,L2)*EW(K,L,L1)
     X  - UM1(I,K,L3)*USP2(I,L,IJSP(L2,L4))*EW(K,L,L1)
     X  - USM2(I,K,IJSP(L2,L4))*UP1(I,L,L3)*EW(K,L,L1)
     X  - UM1(I,K,L4)*USP2(I,L,IJSP(L2,L3))*EW(K,L,L1)
     X  - USM2(I,K,IJSP(L2,L3))*UP1(I,L,L4)*EW(K,L,L1)
      U1U1E = U1U1E - UWP(I,K,L1)*US(I,L,L2)*E2(L,K,IJCP(L3,L4))
     X - US(I,K,L2)*UWM(I,L,L1)*E2(L,K,IJCP(L3,L4))
     X - UWP(I,K,L1)*UP1(I,L,L3)*ES2(L,K,IJSP(L2,L4))
     X - UM1(I,K,L3)*UWM(I,L,L1)*ES2(L,K,IJSP(L2,L4))
     X - UWP(I,K,L1)*UP1(I,L,L4)*ES2(L,K,IJSP(L2,L3))
     X - UM1(I,K,L4)*UWM(I,L,L1)*ES2(L,K,IJSP(L2,L3))
C
   11 CONTINUE
C  .. two;occupation
      GANMA = -(U1U2F1+U1U2F2+U1U1F+U1U2E1+U1U2E2+U1U1E)*TWO
C     WRITE(6,*) L1,',',L2,',',L3,',',L4,' gamma = ',GANMA
      WRITE(6,20) AXYZ(L1),AXYZ(L2),AXYZ(L3),AXYZ(L4),GANMA
   20 FORMAT(1H ,4A1,' = ',F20.5)
      GPARA = GPARA + GANMA
      IF(LCOM.EQ.1.OR.LCOM.EQ.2.OR.LCOM.EQ.3) GPARA=GPARA+GANMA+GANMA
      IF(LCOM.EQ.1.OR.LCOM.EQ.2.OR.LCOM.EQ.3) THEN
      GPERP=GPERP+GANMA
      ELSE
      IF(LCOM.EQ.6.OR.LCOM.EQ.7.OR.LCOM.EQ.12.OR.LCOM.EQ.13.OR.LCOM
     X .EQ.18.OR.LCOM.EQ.19)  THEN
      GPERP=GPERP+GANMA*TWO
      ELSE
      GPERP = GPERP - GANMA*HALF
      ENDIF
      ENDIF
C     WRITE(6,*) U1U2F1,U1U2F2,U1U1F,U1U2E1,U1U2E2,U1U1E
   50 CONTINUE
      GPARA= GPARA/15
      GPERP= GPERP/15
      WRITE(6,*) ' Parallel and Perpondecular (a.u.) '
      WRITE(6,*) GPARA,GPERP
      GPARAE = GPARA*EGAMMA
      GPERPE = GPERP*EGAMMA
      WRITE(6,*) '  in 10**-39 e.s.u. '
      WRITE(6,*) GPARAE,GPERPE
      RETURN
      END
