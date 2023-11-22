C
      SUBROUTINE T2QGMNIJ(ICORE, MAXCOR, IRREPX, IOFFSET, IUHF, FACT)
C
C This routine computes,
C
C       Z(MN,IJ) = 1/2 SUM E,F T(IJ, EF)*Hbar(MN, EF)
C
C The multiplication is carried out using symmetry and the resulting
C product array is stored in a symmetry adapted way. Notice the
C difference from the CC code in which T(IJ,EF) is in fact Tau(IJ,EF).
C
C Notice the factor half which is different from JCP, 94, 4334, 1991.
C The reason for this is that in the final expression (doubles)
C
C     Z(MN,IJ) = (1/4 SUM E,F T(IJ, EF)*Hbar(MN, EF)) and
C
C     Z(AB,EF) = (1/4 SUM M,N T(MN, AB)*Hbar(MN, EF))
C
C will be indentical. Because of that we need only one of them
C and now carry a factor half. 
C
C Note  ISPIN = 1 :  Hbar(MN,EF) (UHF only)
C       ISPIN = 2 :  Hbar(mn,ef) (UHF only)
C       ISPIN = 3 :  Hbar(Mn,Ef) (RHF and UHF)
C
C Originally coded by JG June/90 and modified by ajith 06/1994
C 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      INTEGER DIRPRD,DISSYW,DISSYT,POP1,POP2,VRT1,VRT2
      DIMENSION ICORE(MAXCOR)
C
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
      COMMON /SYM/ POP1(8),POP2(8),VRT1(8),VRT2(8),NTAA,NTBB,
     &             NF1AA,NF1BB,NF2AA,NF2BB
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
      MXCOR = MAXCOR
C
      IF (IUHF .EQ. 1) THEN
C
C AAAA and BBBB spin cases
C
         DO 1000 ISPIN = 1, 2   
C
            LISTT = ISPIN + (IAPRT2AA2 - 1)
            LISTG = ISPIN + (INGMNAA - 1)
            LISTW = ISPIN + 13
C
C Loop over Irreps
C
            DO 100 IRREPEF = 1, NIRREP
               IRREPWEF = IRREPEF
               IRREPWMN = IRREPWEF
               IRREPTEF = IRREPWEF
               IRREPTIJ = DIRPRD(IRREPTEF, IRREPX) 
               IRREPGMN = IRREPWMN
               IRREPGIJ = IRREPTIJ
C
C Retrive Hbar integrals and T2 amplitudes.
C
               DISSYW = IRPDPD(IRREPWEF, ISYTYP(1, LISTW))
               NUMSYW = IRPDPD(IRREPWMN, ISYTYP(2, LISTW))
               DISSYT = IRPDPD(IRREPTEF, ISYTYP(1, LISTT))
               NUMSYT = IRPDPD(IRREPTIJ, ISYTYP(2, LISTT))
C
               I001 = 1
               I002 = I001 + IINTFP*NUMSYW*DISSYW
               I003 = I002 + IINTFP*NUMSYT*DISSYT
               I004 = I003 + IINTFP*NUMSYT*NUMSYW
C
               IF (MIN(NUMSYT, NUMSYW, DISSYT, DISSYW) .NE. 0) THEN
                  I005 = I004 + IINTFP*MAX(DISSYT, DISSYW)
                   IF (I005 .LT. MXCOR) THEN
C  
C In core version        
C
                     CALL MKT2GMNIJ(ICORE(I001), ICORE(I002),
     &                              ICORE(I003), ISPIN, DISSYW, DISSYT,
     &                              NUMSYW, NUMSYT, IRREPWMN, IRREPTIJ, 
     &                              IRREPGIJ, ICORE(I004), IOFFSET, 
     &                              IRSTART, FACT)
                  ELSE
                     CALL INSMEM('T2QGMNIJ', I005, MXCOR)
                  ENDIF
               ELSE
C
C Put zero for this irrep. Note that in the usual CC program this
C is initialized by appropriate bare integrals.
C
                  CALL IZERO(ICORE(I003), IINTFP*NUMSYT*NUMSYW)
                  CALL PUTLST(ICORE(I003), 1, NUMSYT, 1, IRREPGIJ,
     &                        LISTG)
C
               ENDIF
C
 100        CONTINUE
 1000    CONTINUE
C
      ENDIF
C     
C AB spin case
C
      LISTT = IAPRT2AB5
      LISTG = INGMNAB
      LISTW = 16
C
C Loop over irreps.
C
      DO 110 IRREPEF = 1, NIRREP
         IRREPWEF = IRREPEF
         IRREPWMN = IRREPWEF
         IRREPTEF = IRREPWEF
         IRREPTIJ = DIRPRD(IRREPTEF, IRREPX) 
         IRREPGMN = IRREPWMN
         IRREPGIJ = IRREPTIJ
C 
C Retrive T2 amplitudes and Hbar integrals
C
         DISSYW = IRPDPD(IRREPWEF, ISYTYP(1, LISTW))
         NUMSYW = IRPDPD(IRREPWMN, ISYTYP(2, LISTW))
         DISSYT = IRPDPD(IRREPTEF, ISYTYP(1, LISTT))
         NUMSYT = IRPDPD(IRREPTIJ, ISYTYP(2, LISTT))
C
         I001 = 1
         I002 = I001 + IINTFP*NUMSYW*DISSYW
         I003 = I002 + IINTFP*NUMSYT*DISSYT
         I004 = I003 + IINTFP*NUMSYT*NUMSYW
C
         IF (MIN(NUMSYT, NUMSYW, DISSYT, DISSYW) .NE. 0) THEN
            I005 = I004 + IINTFP*MAX(DISSYT,DISSYW)
            IF (I005 .LT. MXCOR) THEN
C
C In core version
C
               CALL MKT2GMNIJ(ICORE(I001), ICORE(I002), ICORE(I003), 
     &                        3, DISSYW, DISSYT, NUMSYW, NUMSYT,
     &                        IRREPWMN, IRREPTIJ, IRREPGIJ,
     &                        ICORE(I004), IOFFSET, IRSTART, FACT)
            ELSE
                CALL INSMEM('T2QGMNIJ', I005, MXCOR)
            ENDIF
         ELSE
C
C Put zero for this irrep. Note that in the usual CC program this
C is initialized by appropriate bare integrals.
C
            CALL IZERO(ICORE(I003), IINTFP*NUMSYT*NUMSYW)
            CALL PUTLST(ICORE(I003), 1, NUMSYT, 1, IRREPGIJ, LISTG)
         ENDIF
C
 110  CONTINUE
C     
      RETURN
      END
