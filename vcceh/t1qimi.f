C
      SUBROUTINE T1QIMI(ICORE, MAXCOR, IRREPX, IOFFSET, IUHF, ISPIN)
C
C This subroutine calculates the T1 contribution to the I(MI)
C intermediate. Spin orbital expression for the contribution 
C calculated in this routine is as  follows
C 
C  Z(MI) = SUM E,N T(E,N)*Hbar(MN,IE)
C
C Spin integrated formulas are as given below
C 
C UHF
C
C  Z(MI) = SUM E,N T(E,N)*Hbar(MN,[IE]) + SUM e,n T(e,n)*Hbar(Mn,Ie) [AA]
C
C  Z(mi) = SUM e,n T(e,n)*Hbar(mn,[ie]) + SUM E,N T(E,N)*Hbar(mN,iE) [BB]
C
C RHF 
C 
C  Z(MI) = SUM E,N T(E,N)*Hbar(MN,[IE]) + SUM e,n T(e,n)*Hbar(Mn,Ie) [AA]
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DOUBLE PRECISION ONE, ONEM, ZILCH
      INTEGER POP, VRT, DIRPRD, DISSYW1, DISSYW2, DISSYW
      DIMENSION ICORE(MAXCOR)
      CHARACTER*6 SPCASE(2)
C           
      COMMON /INFO/ NOCA,NOCB,NVRTA,NVRTB
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),NTOT(18)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
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
      DATA ONE /1.0D+00/
      DATA ONEM /-1.0D+0/
      DATA ZILCH /0.0D+0/ 
      DATA SPCASE /'AA =  ', 'BB =  '/
C     
      LISTT1 = IAPRT1AA
      LISTT2 = IAPRT1AA
      LISTW1 = 6 + ISPIN
      LISTW2 = 11 - ISPIN
      IF (IUHF .EQ. 0) LISTW2 = 10
C      
C Let's determine the size of the T1 vector and allocate memory for it.
C Also allocate memory for the target matrix Z(MI).
C
      NFOOTAR = IRPDPD(IRREPX, 20 + ISPIN)
      NFVOT1  = IRPDPD(IRREPX,  8 + ISPIN)
      NFVOT2  = IRPDPD(IRREPX, 11 - ISPIN)
C
C Allocate memory for the target and T1 amplitudes,      
C
      I000 = 1
      I010 = I000 + IINTFP*NFOOTAR
C      
      IF (IUHF .EQ. 0) THEN
         I020 = I010
      ELSE
         I020 = I010 + IINTFP*NFVOT1
      ENDIF
C
      I030 = I020 + IINTFP*NFVOT2
C
C The address I000 corresponds to the target I(MI) matrix.
C The address I010 and I020 corresponds to the perturb T1 amplitudes.
C The address I030 corresponds to the Hbar matrix elements.
C
C Read the perturb T1 amplitudes in to the memory.
C     
      CALL GETLST(ICORE(I010), 1, 1, 1, ISPIN, LISTT1)
C
      IF (IUHF .NE. 0) THEN
         CALL GETLST(ICORE(I020), 1, 1, 1, 3 - ISPIN, LISTT2)
      ENDIF
C
C Unfortunately, Hbar(IJ,KA), [I<J,KA] lists stored in disk (List 7-10)      
C does not have the correct structure for straight DPD multiplication.
C We have to resort these list such as Hbar(IJ,KA), [I<J,KA] for direct DPD
C multiplication. This is a memory intensive step and partial out-of core routines
C are available. Load the full LISTW1 into the memory. The W1 vector is  ordered
C as (M<N,IE) and (m<n,ie) depending on the value of ISPIN. Change the ordering
C to (NE,MI) [ISPIN = 1] and (ne,mi) [ISPIN = 2] respectively.
C
C Determine the total size of the resorted list
C     
      ISYTYPL =  15 + ISPIN
      ISYTYPR =  20 + ISPIN
C
      NSIZEOUT = IDSYMSZ(1, ISYTYPL, ISYTYPR)
C
C Allocate memory for the sorted list
C     
      I040 = I030 + NSIZEOUT*IINTFP
C
      IF (I040 .LT. MAXCOR) THEN
C
         MXCOR = (MAXCOR - I040 + 1)/IINTFP
C
         CALL QRESORT(ICORE(I040), ICORE(I030), MXCOR, POP(1,ISPIN), 
     &                POP(1,ISPIN), POP(1,ISPIN), VRT(1,ISPIN), 
     &                1, NSIZEOUT, '2413', LISTW1, .TRUE., 20+ISPIN,
     &                ISYTYPL, ISYTYPR)
C
      ELSE
         CALL INSMEM('T1QIMI', I040, MAXCOR)
      ENDIF
C
C We only need to consider the DPD of IRREPX.
C
      IRREPEN = IRREPX
      IRREPMI = IRREPEN
C
C Resorted list is in ICORE(I030). Get the address corresponds to 
C the IRREPX.
C
      ITOP = I030
C
      DO 100 IRREP = 1, (IRREPX - 1)
         DISSYW = IRPDPD(IRREP, ISYTYPL)
         NUMSYW = IRPDPD(IRREP, ISYTYPR)
         ITOP = ITOP + DISSYW*NUMSYW*IINTFP
 100  CONTINUE
C
      DISSYW1 = IRPDPD(IRREPEN, ISYTYPL)
      NUMSYW1 = IRPDPD(IRREPMI, ISYTYPR)
      IOFF    = ITOP
C
C Copy the elements correspond to irrepx to the memory locations
C begining from I030.
C
      I040 = I030 + DISSYW1*NUMSYW1*IINTFP
      I050 = I040 + DISSYW1*NUMSYW1*IINTFP
C
      IF (I050 .LT. MAXCOR) THEN
c YAU : old
c        CALL ICOPY(DISSYW1*NUMSYW1*IINTFP,ICORE(IOFF),1,ICORE(I030),1)
c YAU : new
         CALL DCOPY(DISSYW1*NUMSYW1,ICORE(IOFF),1,ICORE(I030),1)
c YAU : end
      ELSE
         CALL INSMEM('T1QIMI', I050, MAXCOR)
      ENDIF
C
C Change the ordering (NE,MI) to (EN,MI) [ISPIN = 1] and (ne,mi) to
C (en,mi) [ISPIN = 2] for those elements to be used in the following
C multiplication.
C
      CALL SYMTRA2(IRREPEN, POP(1, ISPIN), VRT(1, ISPIN), DISSYW1,
     &             NUMSYW1, ICORE(I030), ICORE(I040))
C
C Now do the first multiplication. The product is stored in ICORE(I000)
C and the T1 vector is stored in ICORE(I010). 
C
      CALL XGEMM('N', 'N', 1, NUMSYW1, DISSYW1, ONE, ICORE(I010),
     &            1, ICORE(I040), DISSYW1, ZLICH, ICORE(I000), 1)
C
C Now do the second term of the full contribution. The W2 vector
C is ordered (Mn,Ie) [ISPIN = 1] and (Nm,Ei) for [ISPIN = 2].
C Reorder them in such a way that (ne,MI) [ISPIN = 1] and (NE,mi)
C when [ISPIN = 2]. 
C  
      ISYTYPL = 15 + (3 - ISPIN)
      ISYTYPR = 20 + ISPIN
C
      NSIZEOUT = IDSYMSZ(1, ISYTYPL, ISYTYPR)
C
C Allocate memory for the sorted list 
C     
      I040 = I030 + NSIZEOUT*IINTFP
C
      IF (I040 .LT. MAXCOR) THEN
C     
         MXCOR = (MAXCOR - I040 + 1)/IINTFP
C
         IF (ISPIN .EQ. 1) THEN
C
            CALL QRESORT(ICORE(I040), ICORE(I030), MXCOR, POP(1,1), 
     &                   POP(1,2), POP(1,1), VRT(1,2), 1, NSIZEOUT, 
     &                  '2413', LISTW2, .FALSE., 14, ISYTYPL, ISYTYPR)
C
         ELSE IF (ISPIN .EQ. 2) THEN
C     
            CALL QRESORT(ICORE(I040), ICORE(I030), MXCOR, POP(1,1), 
     &                   POP(1,2), VRT(1,1), POP(1,2), 1, NSIZEOUT, 
     &                  '1324', LISTW2, .FALSE., 14, ISYTYPL, ISYTYPR)
         ENDIF
C      
      ELSE
         CALL INSMEM('T1QIMI', I040, MAXCOR)
      ENDIF
C
C Resorted list is in ICORE(I030). Get the address corresponds to 
C the IRREPX.
C
      ITOP = I030
C
      DO 200 IRREP = 1, (IRREPX - 1)
         DISSYW = IRPDPD(IRREP, ISYTYPL)
         NUMSYW = IRPDPD(IRREP, ISYTYPR)
         ITOP = ITOP + DISSYW*NUMSYW*IINTFP
 200  CONTINUE
C
      DISSYW2 = IRPDPD(IRREPEN, ISYTYPL)
      NUMSYW2 = IRPDPD(IRREPMI, ISYTYPR)
      IOFF    = ITOP
C
C Copy the elements correspond to irrepx to the memory locations
C begining from I030.
C
      I040 = I030 + DISSYW2*NUMSYW2*IINTFP
      I050 = I040 + DISSYW2*NUMSYW2*IINTFP
C
      IF (I050 .LT. MAXCOR) THEN
c YAU : old
c        CALL ICOPY(DISSYW2*NUMSYW2*IINTFP,ICORE(IOFF),1,ICORE(I030),1)
c YAU : new
         CALL DCOPY(DISSYW2*NUMSYW2,ICORE(IOFF),1,ICORE(I030),1)
c YAU : end
      ELSE
         CALL INSMEM('T1QIMI', I050, MAXCOR)
      ENDIF
C
C Change the ordering (ne,MI) to (en,MI) [ISPIN = 1] and (NE,mi) to
C (EN,mi) [ISPIN = 2] for those elements to be used in the following
C multiplication.
C
      CALL SYMTRA2(IRREPEN, POP(1, 3-ISPIN), VRT(1, 3-ISPIN), DISSYW2,
     &             NUMSYW2, ICORE(I030), ICORE(I040))
C     
C Now do the second multiplication. The product is stored in ICORE(I000)
C and the T1 vector is stored in ICORE(I020).
C
      CALL XGEMM('N', 'N', 1, NUMSYW2, DISSYW2, ONE, ICORE(I020),
     &            1, ICORE(I040), DISSYW2, ONE, ICORE(I000), 1)
C
C Write them to the disk. Note that this is the first contribution
C to the F(MI)/F(mi) interemediate. So we do not have to update.
C 
      CALL PUTLST(ICORE(I000), 1, 1, 1, ISPIN, INTIMI)
C
      IF (IFLAGS(1) .GE. 20) THEN

         CALL HEADER('Checksum @-T1QIMI', 0, LUOUT)
            
         WRITE(LUOUT, *) SPCASE(ISPIN), SDOT(NFOOTAR, ICORE(I000),
     &                   1, ICORE(I000), 1)
C         
      ENDIF

      RETURN
      END
