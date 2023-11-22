C      
      SUBROUTINE MKT2IAEAB(W, T, T1, Z, ISPIN, POP1, POP2, VRT1, VRT2,
     &                    NVRTSQ, DISSYW, DISSYT, NUMSYW, NUMSYT, 
     &                    IRREPWL, IRREPTL, IRREPWR, IRREPTR, IRREPX,
     &                    TMP, IOFFSET, IUHF, FACT)
C
C Calculate the following dontributions to the F(EA) intermediate
C depending the vslue of ISPIN. Also handles RHF spin adapted code.
C
C  Z(AE) = - SUM f SUM M,n Hbar(Mn,Ef)*T(Mn,Af) [ISPIN = 1]
C
C  Z(ae) = - SUM F SUM m,N Hbar(mN,eF)*T(mN,aF) [ISPIN = 2]
C
C RHF spin adapted expression.
C
C  Z(AE) = - SUM f SUM M,n {2*Hbar(Mn,Ef) - Hbar(Mn,Fe)}*T(Mn,Af)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      INTEGER DISSYT,DISSYW,DIRPRD,POP1,POP2,VRT1,VRT2
      DIMENSION W(NUMSYW,DISSYW),T(NUMSYT,DISSYT),Z(NVRTSQ)
      DIMENSION T1(DISSYW,NUMSYW)
      DIMENSION TMP(1)
      DIMENSION POP1(8),POP2(8),VRT1(8),VRT2(8),IWOFF(8)
C
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),DIRPRD(8,8)
C
C Common blocks used in the quadratic term
C 
      COMMON /QADINTI/ INTIMI, INTIAE, INTIME
      COMMON /QADINTG/ INGMNAA, INGMNBB, INGMNAB, INGABAA,
     &                 INGABBB, INGABAB, INGMCAA, INGMCBB,
     &                 INGMCABAB, INGMCBABA, INGMCBAAB,
     &                 INGMCABBA
      COMMON /APERTT1/ IAPRT1AA             
      COMMON /BPERTT1/ IBPRT1AA
      COMMON /APERTT2/ IAPRT2AA1, IAPRT2BB1, IAPRT2AB1, IAPRT2AB2, 
     &                 IAPRT2AB3, IAPRT2AB4, IAPRT2AA2, IAPRT2BB2,
     &                 IAPRT2AB5 
C
      DATA ONE,ONEM,TWO,HALF /1.0D0,-1.0D0,2.0D0,0.5D0/
C
      FACTM = - FACT
C
C Pick up first the relevent T2 and Hbar pieces     
C
      LISTW = 16
      LISTT = IAPRT2AB5  
C    
      CALL GETLST(T1, 1, NUMSYW, 2, IRREPWR, LISTW)
      CALL TRANSP(T1, W, NUMSYW, DISSYW)
C     
C Spin adapted code for RHF
C     
      IF(IUHF .EQ. 0) THEN
         CALL SPINAD1(IRREPWL, VRT1, NUMSYW, W, TMP, TMP(1 + NUMSYW))
      ENDIF
C     
      CALL GETLST (T1, 1, NUMSYT,1, IRREPTR, LISTT)
      CALL TRANSP (T1, T, NUMSYT, DISSYT)
C     
      IF (IUHF .EQ. 1 .AND. ISPIN .EQ. 1) THEN
         CALL SYMTR1(IRREPWL, VRT1, VRT2, NUMSYW, W, TMP, 
     &               TMP(1 + NUMSYW), TMP(1 + 2*NUMSYW))
         CALL SYMTR1(IRREPTL, VRT1, VRT2, NUMSYT, T, TMP, 
     &               TMP(1 + NUMSYT), TMP(1 + 2*NUMSYT))
      ENDIF
C     
      JOFF = 1
      KOFF = 1
      IWOFF(1) = 1
C
      DO 5000 IRREPE = 2, NIRREP
C
         IRREP = IRREPE - 1
         IRREPF = DIRPRD(IRREP, IRREPWL)
C
         IWOFF(IRREPE) = IWOFF(IRREPE - 1) + VRT1(IRREP)*VRT2(IRREPF)
C
 5000 CONTINUE
C
      DO 90 IRREPAO = 1, NIRREP
C     
C Determine irrepf and irrepa.
C     
         IRREPF = DIRPRD(IRREPAO, IRREPTL)
         IRREPE = DIRPRD(IRREPF, IRREPWL)
C     
C Get number of virtuals orbitals for irrepe
C     
         NVRTE = VRT1(IRREPE)
         NVRTF = VRT2(IRREPF)
         NVRTA = VRT1(IRREPAO)
C
         IOFF = IWOFF(IRREPE)
C
C If zero nothing to compute     
C     
         IF (NVRTF .GT. 0 .AND. NVRTE .GT. 0) THEN
C     
            CALL XGEMM('T', 'N', NVRTE, NVRTA, NUMSYW*NVRTF, FACTM,
     &                  W(1,IOFF), NVRTF*NUMSYW, T(1,JOFF),
     &                  NVRTF*NUMSYT, ONE, Z(KOFF), NVRTE)
         ENDIF
C     
C Update the offset for the intermediate Z
C     
         JOFF = JOFF + NVRTA*NVRTF
         KOFF = KOFF + NVRTE*NVRTA
C
 90   CONTINUE
C     
      RETURN
      END
