C
      SUBROUTINE T2QIMI(ICORE, MAXCOR, IRREPX, IOFFSET, IUHF, FACT)
C
C This routine computes,
C     
C      Z(MI) = 1/2 SUM N SUM E,F Hbar(MN,EF)*T(IN,EF)
C
C The multiplication is carried out using symmetry and the resulting
C product array is stored in a symmetry adapted way. Notice the
C difference from the CC code in which T(IN,EF) is in fact Tau(prime)(IN,EF).
C Spin integrated formulas are given as given below.
C
C UHF
C
C  Z(MI) =  SUM N SUM E<F {T(IN,EF)*Hbar(MN,EF)} [AAAA]
C
C  Z(mi) =  SUM n SUM e<f {T(in,ef)*Hbar(mn,ef)} [BBBB]
C
C  Z(MI) =  SUM n SUM E,f {T(In,Ef)*Hbar(Mn,Ef)} [ABAB]
C
C  Z(mi) =  SUM N SUM e,F {T(iN,eF)*Hbar(mN,eF)} [BABA]
C
C RHF; Spin adapted 
C
C  Z(MI) = SUM N SUM E,f {2*Hbar(Mn,Ef) - Hbar(Mn,Fe)}*T(In,Ef)
C
C Originaly coded by JG June/90 and modified by ajith 06/1994
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      INTEGER DIRPRD,DISSYW,DISSYT,POP,VRT
      DIMENSION ICORE(MAXCOR)
      CHARACTER*6 SPCASE(2)
C
      COMMON /SYM/ POP(8,2),VRT(8,2),NTAA,NTBB,NF1AA,NF1BB,NF2AA,NF2BB
      COMMON /INFO/ NOCCO(2),NVRTO(2)
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
      DATA SPCASE /'AA =  ', 'BB =  '/
C
      MXCOR = MAXCOR
C     
      DO 1000 ISPIN = 1, IUHF + 1 
C
         NOCCSQ = 0
         DO 50 IRREPM = 1, NIRREP
            IRREPI = DIRPRD(IRREPM, IRREPX) 
            NOCCSQ = NOCCSQ + POP(IRREPM, ISPIN)*POP(IRREPI, ISPIN)
 50      CONTINUE
C     
C AAAA and BBBB spin cases
C     
         I000 = MXCOR + 1 - IINTFP*NOCCSQ
C
C Load F(MI) contributions already in the disk in to the memory
C
         CALL GETLST(ICORE(I000), 1, 1, 1, ISPIN, INTIMI)
C
         IF (IUHF .EQ. 1) THEN
C
            LISTT = (IAPRT2AA2 - 1) + ISPIN
            LISTW = 13 + ISPIN
C     
C Loop over irreps of EF block
C     
            DO 100 IRREPEF = 1, NIRREP
C
               IRREPWEF = IRREPEF
               IRREPWMN = IRREPWEF
               IRREPTEF = IRREPWEF
               IRREPTIN = DIRPRD(IRREPTEF, IRREPX)
C    
C Retrieve Hbar integrals and T2 amplitudes
C     
               NOCC2SQ1 = 0
               NOCC2SQ2 = 0
               DO 110 IRREPM = 1, NIRREP
                  IRREPN  = DIRPRD(IRREPM, IRREPWMN)
                  NOCC2SQ1 = NOCC2SQ1 + POP(IRREPM, ISPIN)
     &                      *POP(IRREPN, ISPIN)
 110           CONTINUE 
C
               DO 120 IRREPN = 1, NIRREP
                  IRREPI  = DIRPRD(IRREPN, IRREPTIN)
                  NOCC2SQ2 = NOCC2SQ2 + POP(IRREPI, ISPIN)
     &                      *POP(IRREPN, ISPIN)
 120           CONTINUE 
C
               DISSYW = IRPDPD(IRREPWEF, ISYTYP(1, LISTW))
               NUMSYW = IRPDPD(IRREPWMN, ISYTYP(2, LISTW))  
               DISSYT = IRPDPD(IRREPTEF, ISYTYP(1, LISTT))
               NUMSYT = IRPDPD(IRREPTIN, ISYTYP(2, LISTT))
C     
               I001 = 1
               I002 = I001 + IINTFP*MAX(NOCC2SQ2*DISSYT,
     &                NOCC2SQ1*DISSYW)
               I003 = I002 + IINTFP*MAX(NOCC2SQ2*DISSYT,
     &                NOCC2SQ1*DISSYW)
C
               IF (MIN(NUMSYT, NUMSYW, DISSYT, DISSYW) .NE. 0) THEN
C
                  I004 = I003 + IINTFP*MAX(NUMSYT, NUMSYW)
C
                  IF(I004 .LT. MXCOR) THEN
C     
C In core version
C
                     CALL MKT2IMIAA(ICORE(I001), ICORE(I002), 
     &                              ICORE(I000), ISPIN, POP(1,ISPIN),
     &                              VRT(1,ISPIN), DISSYW, DISSYT,
     &                              NUMSYW, NUMSYT, NOCCSQ, NOCC2SQ1, 
     &                              NOCC2SQ2, IRREPWEF, IRREPWMN,
     &                              IRREPTEF, IRREPTIN, IRREPX,
     &                              ICORE(I003), IOFFSET, IUHF, FACT)
                  ELSE
                     CALL INSMEM('T2QIMIAA', I004, MXCOR)
                  ENDIF
C
               ELSE
               ENDIF                
C
 100        CONTINUE
C
         ENDIF
C     
C AB spin case
C     
         LISTT = IAPRT2AB5
         LISTW = 16
C     
C Loop over irreps 
C     
         DO 200 IRREPEF = 1, NIRREP
C
            IRREPWEF = IRREPEF
            IRREPWMN = IRREPWEF
            IRREPTEF = IRREPWEF
            IRREPTIN = DIRPRD(IRREPTEF, IRREPX)
C     
C Retrive integrals and T2 amplitudes
C     
            DISSYW = IRPDPD(IRREPWEF, ISYTYP(1, LISTW))
            NUMSYW = IRPDPD(IRREPWMN, ISYTYP(2, LISTW))
            DISSYT = IRPDPD(IRREPTEF, ISYTYP(1, LISTT))
            NUMSYT = IRPDPD(IRREPTIN, ISYTYP(2, LISTT))
C
            I001 = 1
            I002 = I001 + IINTFP*MAX(NUMSYW*DISSYW, NUMSYT*DISSYT)
            I003 = I002 + IINTFP*MAX(NUMSYW*DISSYW, NUMSYT*DISSYT)
C
            IF (MIN (NUMSYT, NUMSYW, DISSYT, DISSYW) .NE. 0) THEN
C
               I004 = I003 + 3*IINTFP*MAX(NUMSYT, NUMSYW, DISSYW,
     &                DISSYT)
               IF (I004 .LT. MXCOR) THEN
C     
C In core version
C     
                  CALL MKT2IMIAB(ICORE(I001), ICORE(I002), 
     &                           ICORE(I000), ISPIN, POP(1,ISPIN),
     &                           POP(1,3-ISPIN), VRT(1,ISPIN), 
     &                           VRT(1,3-ISPIN), DISSYW, DISSYT,
     &                           NUMSYW, NUMSYT, NOCCSQ, IRREPWEF,
     &                           IRREPTEF, IRREPWMN, IRREPTIN, IRREPX,
     &                           ICORE(I003), IOFFSET, IUHF, FACT)
               ELSE
                  CALL INSMEM('T2QIMIAA', I004, MXCOR)
               ENDIF
C
            ELSE
            ENDIF
C
 200     CONTINUE
C
C Update the F(MI) intermediate lists
C
         CALL PUTLST (ICORE(I000), 1, 1, 1, ISPIN, INTIMI)
C
         IF (IFLAGS(1) .GE. 20) THEN
C     
            NSIZE = NOCCSQ
C
            CALL HEADER('Checksum @-T2QIMI', 0, LUOUT)
            
            WRITE(LUOUT, *) SPCASE(ISPIN), SDOT(NSIZE, 
     &                      ICORE(I000), 1, ICORE(I000), 1)
         ENDIF
C
 1000 CONTINUE
C     
      RETURN
      END
