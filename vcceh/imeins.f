C
      SUBROUTINE IMEINS(ICORE, MAXCOR, IRREPX, IOFFSET, IUHF,ISPIN)
C
C This routine compute the contribution of I(ME) intermediate
C in the singles equation.
C
C  Z(i,a) = SUM ME T2(IM,AE) * I(ME) 
C 
C Spin integrated formulas are as follows and alpha and beta are 
C implicit in these equations.
C
C  Z(A,I) = I(E,M) * T2(IM,AE) + F(e,m) * T2(Im,Ae)
C
C  Z(a,i) = F(e,m) * T2(im,ae) + F(E,M) * T2(iM,aE)
C
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE,ONEM,ZILCH, SDOT
      DIMENSION ICORE(MAXCOR)
      CHARACTER*6 SPCASE(2)
C
      COMMON /INFO/ NOCA,NOCB,NVRTA,NVRTB
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),NTOT(18)
      COMMON /FLAGS/ IFLAGS(100)
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
      DATA ONE /1.0D0/ 
      DATA ONEM /-1.0D0/
      DATA ZILCH /0.0D0/
      DATA SPCASE /'AA =  ', 'BB =  '/
C
      LISTTA = (IAPRT2AA1 - 1) + ISPIN
      LISTTB = (IAPRT2AB1 - 1) + ISPIN
      IF(IUHF .EQ. 0) LISTTB = IAPRT2AB2
      LISTI = INTIME
C
C We only need to consider one DPD irrep, IRREPX the symmetry of
C the perturbation
C
      NLENT1 = IRPDPD(1, 8 +  ISPIN)
      NLENIA = IRPDPD(IRREPX, 8 + ISPIN)
      NLENIB = IRPDPD(IRREPX, 11 - ISPIN)
C
C I000 Holds Q(AI) target in eqaution above.
C I010 Holds I(ME) in equation above, 
C I020 Holds I(me) in equation above.
C I030 Holds individual IRREPX of blocked T2 list.
C
      I000 = 1
      I010 = I000 + NLENT1*IINTFP
C
      IF (IUHF .EQ. 0) THEN
         I020 = I010
      ELSE
         I020 = I010 + NLENIA*IINTFP
      ENDIF
C
      I030 = I020 + NLENIB*IINTFP
C
C Load the I(ME) from the disk
C
      CALL GETLST(ICORE(I010), 1, 1, 1, ISPIN, LISTI)
      IF(IUHF .NE. 0) CALL GETLST(ICORE(I020), 1, 1, 1, 3 - ISPIN,
     &                            LISTI)
C
      IRREPIEM  = IRREPX
      IRREPTAEM = IRREPIEM
      IRREPTAAI = DIRPRD(IRREPTAEM, IRREPX)
      IRREPTBEM = IRREPTAEM
      IRREPTBAI = DIRPRD(IRREPTBEM, IRREPX)
C
      NTDSZA = IRPDPD(IRREPTAEM, ISYTYP(1, LISTTA))
      NTDISA = IRPDPD(IRREPTAAI, ISYTYP(2, LISTTA))
      NTDSZB = IRPDPD(IRREPTBEM, ISYTYP(1, LISTTB))
      NTDISB = IRPDPD(IRREPTBAI, ISYTYP(2, LISTTB))
C
      I040 = I030 + NTDSZA*NTDISA*IINTFP
C
      IF (I040 .GT. MAXCOR) CALL INSMEM('IMEINS', I040, MAXCOR)
C
C Get appropriate T2(beta) amplitudes and take the product.
C
      CALL GETLST(ICORE(I030), 1, NTDISA, 1, IRREPTAAI, LISTTA)
C     
      CALL XGEMM('N','N', 1, NLENT1, NLENIA, ONEM, ICORE(I010),
     &            1, ICORE(I030), NTDSZA, ZILCH, ICORE(I000), 1)

      I040 = I030 + NTDSZB*NTDISB*IINTFP
C
      IF (I040 .GT. MAXCOR) CALL INSMEM('IMEINS', I040, MAXCOR)
C
      CALL GETLST(ICORE(I030), 1, NTDISB, 1, IRREPTBAI, LISTTB)
C
      CALL XGEMM('N', 'N', 1, NLENT1, NLENIB, ONE, ICORE(I020),
     &            1, ICORE(I030), NTDSZB, ONE, ICORE(I000), 1)
C         
C Update the Q(AI) singles list. Allocate memory for Q(AI)
C components already on the disk.
C
      IF (IFLAGS(1) .GE. 20) THEN
C
         IF (ISPIN .EQ. 1) THEN
            NSIZE = IRPDPD(1, 9)
         ELSE
            NSIZE = IRPDPD(1, 10)
         ENDIF
C
         CALL HEADER('Checksum @-IMEINS', 0, LUOUT)
C            
         WRITE(LUOUT, *) SPCASE(ISPIN), SDOT(NSIZE, ICORE(I000),
     &                    1, ICORE(I000), 1)
      ENDIF
C
      I020 = I010 + NLENT1*IINTFP
C
      CALL IZERO(ICORE(I010), NLENT1*IINTFP)
      CALL GETLST(ICORE(I010), 1, 1, 1, 2 + ISPIN, 90)
      CALL SAXPY(NLENT1, ONE , ICORE(I000), 1, ICORE(I010), 1) 
C
      CALL PUTLST(ICORE(I010), 1, 1, 1, 2 + ISPIN, 90)
C
      RETURN
      END
