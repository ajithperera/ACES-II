C
      SUBROUTINE MKT2IAEAA(W, T, T1, Z, ISPIN, POP, VRT, NVRTSQ,
     &                     NVRT2SQ1, NVRT2SQ2, DISSYW, DISSYT, NUMSYW, 
     &                     NUMSYT, IRREPWL, IRREPTL, IRREPWR, IRREPTR,
     &                     IRREPX, TMP, IOFFSET, IUHF, FACT)
C 
C Calculates the following contributions to the F(AE) intermediate
C depending on the value of ISPIN.
C
C  Z(AE) = - SUM F SUM M<N Hbar(MN,EF)*T(MN,AF) [ISPIN = 1]
C
C  Z(ae) = - SUM f SUM m<n Hbar(mn,ef)*T(mn,af) [ISPIN = 2]

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DISSYT, DISSYW, DIRPRD, POP, VRT
      DIMENSION W(NUMSYW,NVRT2SQ1),T(NUMSYT,NVRT2SQ2)
      DIMENSION T1(DISSYT,NUMSYT),Z(NVRTSQ), IWOFF(8)
      DIMENSION TMP(1)
      DIMENSION POP(8),VRT(8)
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
      DATA ONEM,ONE,HALF,TWO /-1.0D0,1.D0,0.5D0,2.D0/
C
      IND(I,J)=((I-2)*(I-1))/2+J
      FACTM = - FACT
C
C Pick up first the relevent T2 and Hbar pieces
C
      LISTW = ISPIN + 13
      LISTT = ISPIN + (IAPRT2AA2 - 1)       
C
      CALL GETLST(T1, 1, NUMSYW, 2, IRREPWR, LISTW)
      CALL TRANSP(T1, W, NUMSYW, DISSYW)
C
      CALL GETLST(T1, 1, NUMSYT, 1, IRREPTR, LISTT)   
      CALL TRANSP(T1, T, NUMSYT, DISSYT)
C
      CALL SYMEXP(IRREPWL, VRT, NUMSYW, W)
      CALL SYMEXP(IRREPTL, VRT, NUMSYT, T)
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
         IWOFF(IRREPE) = IWOFF(IRREPE - 1) + VRT(IRREPF)*VRT(IRREP)
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
C Get number of virtual orbitals for irrepe
C     
         NVRTE = VRT(IRREPE)
         NVRTF = VRT(IRREPF)
         NVRTA = VRT(IRREPAO)
C
         IOFF = IWOFF(IRREPE)
C
C If zero nothing to compute
C     
         IF(NVRTF .GT. 0 .AND. NVRTE .GT. 0) THEN
C     
            CALL XGEMM('T', 'N', NVRTE, NVRTA, NVRTF*NUMSYW, FACTM,
     &                  W(1,IOFF), NVRTF*NUMSYW, T(1,JOFF),
     &                  NVRTF*NUMSYT, ONE, Z(KOFF), NVRTE)
         ENDIF
C
         JOFF = JOFF + NVRTA*NVRTF
         KOFF = KOFF + NVRTE*NVRTA
C     
 90   CONTINUE
C     
      RETURN
      END
