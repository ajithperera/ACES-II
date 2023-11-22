      SUBROUTINE DWMBEJNR(ICORE,MAXCOR,SPCASE,IUHF,C1,C2,C3,C4)
C
C THIS ROUTINE AND DEPENDENTS COMPUTE THE RING-TYPE W(MBEJ)
C  INTERMEDIATE FOR ALL SIX POSSIBLE SPIN CASES.  THE ALGORITHM
C  ASSUMES IN-CORE STORAGE OF SYMMETRY PACKED TARGET, T2 AND W
C  VECTORS AND USES THE DPD SYMMETRY APPROACH TO EVALUATE THE
C  CONTRACTIONS.
C
C SPIN ORBITAL EQUATION FOR THIS INTERMEDIATE (EVENTUAL MODIFICATION
C                     TO INCLUDE T1 IS NECESSARY FOR CCSD):
C
C         W(mbej)=<mb||ej>-(1/2)*SUM T(jn,fb) <mn||ef>
C                                n,f
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION C1,C2,C3,C4,X,TWO,ONE,ONEM,HALF,ZILCH,TWOM
      CHARACTER*4 SPCASE,TAUTYP1,TAUTYP2,SPNCASEW(2),SPNCASET(2),
     &            SPNTYPW, SPNTYP2, SPNTYPT
      DIMENSION ICORE(MAXCOR),IOFFT1(8,2)
      LOGICAL MBPT3,MBPT4,CC,TRPEND,SNGEND,GRAD,MBPTT,SING1,
     &                QCISD,UCC
      COMMON /SWITCH/ MBPT3,MBPT4,CC,TRPEND,SNGEND,GRAD,MBPTT,SING1,
     &                QCISD,UCC
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),NTOT(18)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT1AA,NT1BB,NF1AA,NF1BB,NF2AA,
     &             NF2BB
      COMMON /FLAGS/  IFLAGS(100)
      DATA ZILCH /0.0/
      DATA ONE /1.0/
      DATA HALF/0.5/
      DATA TWO /2.0/
      DATA TWOM /-2.0/
      DATA ONEM/-1.0/
      NNP1O2(I)=I*(I+1)/2
      X=1.0
      IF(C3.NE.0.0)CALL GETT1(ICORE,MAXCOR,MXCOR,IUHF,IOFFT1)
C
C********************************************************************
C
C DO SPIN CASE AAAA AND BBBB:
C
C       W(MBEJ)=<MB||EJ>-(1/2)*SUM [T(JN,FB)<MN||EF>-T(Jn,fB)<Mn|Ef>]
C
C       W(mbej)=<mb||ej>-(1/2)*SUM [T(jn,fb)<mn||ef>-T(jN,Fb)<mN|eF>]
C
C TERMS A AND B BOTH HANDLED BELOW.
C
      IF(SPCASE.EQ.'AAAA'.OR.SPCASE.EQ.'BBBB')THEN
         IF(SPCASE.EQ.'AAAA')THEN
            ISPIN  = 1
            SPNCASEW(1) = 'AABB'
            SPNCASET(1) = 'BBAA'
            TAUTYP1='AAAA'
         ELSEIF(SPCASE.EQ.'BBBB')THEN
            ISPIN  = 2
            SPNCASEW(2) = "BBAA"
            SPNCASET(2) = "AABB"
            TAUTYP1='BBBB'
         ENDIF
C
         REFLISTTA = 43 + ISPIN
         REFLISTWA = 13 + ISPIN
         REFLISTWB = 16
         REFLISTTB = 46
         TARLISTTA = 33 + ISPIN
         TARLISTWA = 18 + ISPIN
         TARLISTTB = 35 + ISPIN
         IF (IUHF .EQ. 0)  TARLISTTB = 37
         TARLISTWB = 19 - ISPIN
C
         ITARSIZWA = ISYMSZ(ISYTYP(1, TARLISTWA), 
     &                      ISYTYP(2, TARLISTWA))
         ITARSIZWB = ISYMSZ(ISYTYP(1, TARLISTWB), 
     &                      ISYTYP(2, TARLISTWB))
         ITARSIZTA = ISYMSZ(ISYTYP(1, TARLISTTA), 
     &                      ISYTYP(2, TARLISTTA))
         ITARSIZTB = ISYMSZ(ISYTYP(1, TARLISTTB), 
     &                      ISYTYP(2, TARLISTTB))
         IREFSIZWA = ISYMSZ(ISYTYP(1, REFLISTWA), 
     &                      ISYTYP(2, REFLISTWA))
         IREFSIZWB = ISYMSZ(ISYTYP(1, REFLISTWB), 
     &                      ISYTYP(2, REFLISTWB))
         IREFSIZTA = ISYMSZ(ISYTYP(1, REFLISTTA), 
     &                      ISYTYP(2, REFLISTTA))
         IREFSIZTB = ISYMSZ(ISYTYP(1, REFLISTTB),
     &                      ISYTYP(2, REFLISTTB))
         ISCRSIZW  = NOCCO(ISPIN)*NVRTO(ISPIN) + (NVRTO(ISPIN)*
     &               (NVRTO(ISPIN)-1))/2 + (NOCCO(ISPIN)*
     &               (NOCCO(ISPIN)-1))/2
         ISCRSIZT  = NOCCO(1)*NOCCO(2) + NVRTO(1)*NVRTO(2) + 
     &               NVRTO(1)*NOCCO(1) + NVRTO(2)*NOCCO(2)
C
         IBGN = 1
         I001 = IBGN + IINTFP*ITARSIZWA
         I002 = I001 + IINTFP*ITARSIZWB
         I003 = I002 + IINTFP*ITARSIZTA
         I004 = I003 + IINTFP*ITARSIZTB
         I005 = I004 + IINTFP*MAX(IREFSIZWA, ITARSIZWA)
         I006 = I005 + ISCRSIZW
         IF (I006 .GT. MAXCOR) CALL INSMEM('DWMBEJNR', I006, MAXCOR)
C
         CALL GETALL(ICORE(I004), IREFSIZWA, 1, REFLISTWA)
         CALL SST003(ICORE(I004), ICORE(IBGN), IREFSIZWA, ITARSIZWA, 
     &               ICORE(I005), SPCASE, "AIBJ")
C
         I005 = I004 + IINTFP*IREFSIZWB
         I006 = I005 + ISCRSIZT
         IF (I006 .GT. MAXCOR) CALL INSMEM('DWMBEJNR', I006, MAXCOR)
C
         CALL GETALL(ICORE(I004), IREFSIZWB, 1, REFLISTWB)
         CALL SST002(ICORE(I004), ICORE(I001), IREFSIZWB, ITARSIZWB, 
     &               ICORE(I005), SPNCASEW(ISPIN))
C
         I005 = I004 + IINTFP*MAX(IREFSIZTA, ITARSIZTA)
         I006 = I005 + ISCRSIZW
         IF (I006 .GT. MAXCOR) CALL INSMEM('DWMBEJNR', I006, MAXCOR)
C
         CALL GETALL(ICORE(I004), IREFSIZTA, 1, REFLISTTA)    
         CALL SST003(ICORE(I004), ICORE(I002), IREFSIZTA, ITARSIZTA, 
     &               ICORE(I005), SPCASE, "AJBI")
C
         I005 = I004 + IINTFP*IREFSIZTB
         I006 = I005 + ISCRSIZT
         IF (I006 .GT. MAXCOR) CALL INSMEM('DWMBEJNR', I006, MAXCOR)
C
         CALL GETALL(ICORE(I004), IREFSIZTB, 1, REFLISTTB)    
         IF (IUHF .NE. 0) THEN
            CALL SST002(ICORE(I004), ICORE(I003), IREFSIZTB,
     &                  ITARSIZTB, ICORE(I005), SPNCASET(ISPIN))
         ELSE
            CALL SST002(ICORE(I004), ICORE(I003), IREFSIZTB,
     &                  ITARSIZTB, ICORE(I005), "AABB")
         ENDIF
C
         MAXSIZ = 0
         LISTWB = 22 + ISPIN
         LISTQ  = 53 + ISPIN
C
         IOFFSTWA    = IBGN
         IOFFSTWB    = I001
         IOFFSTTA    = I002
         IOFFSTTB    = I003   
C
         DO 100 IRREP=1,NIRREP 
            
            DSSYWA = IRPDPD(IRREP, ISYTYP(1, TARLISTWA))
            NMSYWA = IRPDPD(IRREP, ISYTYP(2, TARLISTWA))
            DSSYTA = IRPDPD(IRREP, ISYTYP(1, TARLISTTA))
            NMSYTA = IRPDPD(IRREP, ISYTYP(2, TARLISTTA))
            DSSYWB = IRPDPD(IRREP, ISYTYP(1, TARLISTWB))
            NMSYWB = IRPDPD(IRREP, ISYTYP(2, TARLISTWB))
            DSSYTB = IRPDPD(IRREP, ISYTYP(1, TARLISTTB))
            NMSYTB = IRPDPD(IRREP, ISYTYP(2, TARLISTTB))
            DISSYQ = IRPDPD(IRREP, ISYTYP(1, LISTQ))
            NUMSYQ = IRPDPD(IRREP, ISYTYP(2, LISTQ))
C
            I000 = I004
            I010 = I000+IINTFP*DISSYQ*NUMSYQ
            I020 = I010+IINTFP*NMSYWA*DSSYWA
            I030 = I020+IINTFP*NMSYTA*DSSYTA
            I030 = MAX(I030,I010+DISSYQ*NUMSYQ*IINTFP)
            IF (I030 .GT. MAXCOR) CALL INSMEM('DWMBEJNR', I030, MAXCOR)
C
            CALL IZERO(ICORE(I000),DISSYQ*NUMSYQ*IINTFP)
C
            NSIZE=DISSYQ*NUMSYQ

            IF(C1.NE.0.0)THEN
C
               CALL GETLST(ICORE(I000), 1, NUMSYQ, 1, IRREP, LISTWB)
C
               IF(UCC) THEN
                  CALL SUMSYM2(ICORE(I000),ICORE(I010),NSIZE,1,
     &                         IRREP,LISTQ)
                  CALL IZERO(ICORE(I000),DISSYQ*NUMSYQ*IINTFP)
               ENDIF  
            ENDIF
C
            IF(C4.NE.0.0)THEN
C
               IF(C3.NE.0.0)THEN
                  CALL GETARGS(TARLISTTA,I1,I2)
                  CALL F2TAU(ICORE(IOFFSTTA),ICORE(IOFFT1(1,I1)),
     &                       ICORE(IOFFT1(1,I2)),DSSYTA,NMSYTA,
     &                       POP(1,I1),POP(1,I2),VRT(1,I1),VRT(1,I2),
     &                       IRREP,TWO,TAUTYP1)
               ENDIF
C
               CALL CNTRCTNR(ICORE(I000),ICORE(IOFFSTWA),
     &                       ICORE(IOFFSTTA),DISSYQ,NUMSYQ,
     &                       DSSYTA,NMSYTA,DSSYWA,NMSYWA,
     &                       'WxT',X,-C4,0)
C
               IF(UCC) THEN
                  CALL DT(TWO,ICORE(I000),DISSYQ,NUMSYQ,'ADD',
     &                    ICORE(I030))
                  CALL SUMSYM2(ICORE(I000),ICORE(IOFFSTWA),NSIZE,1,
     &                         IRREP,LISTQ)
                  CALL IZERO(ICORE(I000),DISSYQ*NUMSYQ*IINTFP)
               ENDIF  
C
               CALL CNTRCTNR(ICORE(I000),ICORE(IOFFSTWB),
     &                       ICORE(IOFFSTTB),DISSYQ,NUMSYQ,
     &                       DSSYTB,NMSYTB,DSSYWB,NMSYWB,
     &                       'WxT',X,C4,0)
C
               IF(UCC) THEN
                  CALL DT(TWO,ICORE(I000),DISSYQ,NUMSYQ,'ADD',
     &                    ICORE(I030))
               ENDIF  
            ENDIF
C
            CALL SUMSYM2(ICORE(I000),ICORE(I010),NSIZE,1,IRREP,
     &                   LISTQ)
C
            IOFFSTWA = IINTFP*NMSYWA*DSSYWA + IOFFSTWA
            IOFFSTWB = IINTFP*NMSYWB*DSSYWB + IOFFSTWB
            IOFFSTTA = IINTFP*NMSYTA*DSSYTA + IOFFSTTA
            IOFFSTTB = IINTFP*NMSYTB*DSSYTB + IOFFSTTB

 100     CONTINUE
C
C DO SPIN CASES ABAB AND BABA:
C
C       W(MbEj)=<Mb|Ej>-(1/2)*SUM [T(jn,fb)<Mn|Ef>-T(jN,bF)<MN||EF>]
C
C       W(mBeJ)=<mB|eJ>-(1/2)*SUM [T(JN,FB)<mN|eF>-T(Jn,Bf)<mn||ef>]
C    
      ELSEIF (SPCASE.EQ.'ABAB'.OR.SPCASE.EQ.'BABA')THEN

         IF(IUHF.NE.0)THEN
C
            IF(SPCASE.EQ.'ABAB')THEN
               LISTQ  = 56
               LISTQ1 = 57
C
               SPNTYPW   = 'AAAA'
               IF(IUHF.EQ.0)THEN
                  TARLISTTA = 34
                  REFLISTTA = 44 
                  SPNTYPT   = 'AAAA'
               ELSE
                  TARLISTTA = 35
                  REFLISTTA = 45
                  SPNTYPT   = 'BBBB'
               ENDIF
               TARLISTTB = 37
               TARLISTWA = 18
               TARLISTWB = 19
               REFLISTTB = 46
               REFLISTWA = 16
               REFLISTWB = 14
               ISPIN     = 1
               TAUTYP1   = 'BBBB'
               SPNTYP2   = 'AABB'
            ELSEIF(SPCASE.EQ.'BABA')THEN 
               LISTQ     = 57
               LISTQ1    = 56
               TARLISTTA = 34
               TARLISTWA = 17
               TARLISTWB = 20
               TARLISTTB = 36
               REFLISTTA = 44
               REFLISTWA = 16
               REFLISTWB = 15
               REFLISTTB = 46
               ISPIN     = 2
               TAUTYP1   = 'AAAA'
               SPNTYPT   = 'AAAA'
               SPNTYPW   = 'BBBB'
               SPNTYP2   = 'BBAA'
            ENDIF
C     
            ITARSIZWA = ISYMSZ(ISYTYP(1, TARLISTWA), ISYTYP(2, 
     &                  TARLISTWA))
            ITARSIZWB = ISYMSZ(ISYTYP(1, TARLISTWB), ISYTYP(2, 
     &                  TARLISTWB))
            ITARSIZTA = ISYMSZ(ISYTYP(1, TARLISTTA), ISYTYP(2, 
     &                  TARLISTTA))
            ITARSIZTB = ISYMSZ(ISYTYP(1, TARLISTTB), ISYTYP(2, 
     &                  TARLISTTB))
            IREFSIZWA = ISYMSZ(ISYTYP(1, REFLISTWA), ISYTYP(2, 
     &                  REFLISTWA))
            IREFSIZWB = ISYMSZ(ISYTYP(1, REFLISTWB), ISYTYP(2, 
     &                  REFLISTWB))
            IREFSIZTA = ISYMSZ(ISYTYP(1, REFLISTTA), ISYTYP(2,
     &                  REFLISTTA))
            IREFSIZTB = ISYMSZ(ISYTYP(1, REFLISTTB), ISYTYP(2, 
     &                  REFLISTTB))
            ISCRSIZW  = NOCCO(ISPIN)*NVRTO(ISPIN) + NVRTO(ISPIN)*
     &                  (NVRTO(ISPIN)-1)/2 +
     &                  NOCCO(ISPIN)*(NOCCO(ISPIN)-1)/2
            ISCRSIZT  = NOCCO(1)*NOCCO(2) + NVRTO(1)*NVRTO(2) + 
     &                  NVRTO(1)*NOCCO(1) + NVRTO(2)*NOCCO(2)
C
            IBGN = 1
            I001 = IBGN + IINTFP*ITARSIZWA
            I002 = I001 + IINTFP*ITARSIZWB
            I003 = I002 + IINTFP*ITARSIZTA
            I004 = I003 + IINTFP*ITARSIZTB
            I005 = I004 + IINTFP*MAX(IREFSIZWB, ITARSIZWB)
            I006 = I005 + ISCRSIZW
            IF (I006 .GT. MAXCOR) CALL INSMEM('DWMBEJNR', I006, MAXCOR)
C
            CALL GETALL(ICORE(I004), IREFSIZWB, 1, REFLISTWB)
            CALL SST003(ICORE(I004), ICORE(I001), IREFSIZWB,
     &                  ITARSIZWB, ICORE(I005), SPNTYPW, "AIBJ")
            
            I005 = I004 + IINTFP*IREFSIZWA
            I006 = I005 + ISCRSIZT
            IF (I006 .GT. MAXCOR) CALL INSMEM('DWMBEJNR', I006, MAXCOR)
C
            CALL GETALL(ICORE(I004), IREFSIZWA, 1, REFLISTWA)
            CALL SST002(ICORE(I004), ICORE(IBGN), IREFSIZWA, 
     &                  ITARSIZWA, ICORE(I005), SPNTYP2)

            I005 = I004 + IINTFP*MAX(IREFSIZTA, ITARSIZTA)
            I006 = I005 + ISCRSIZW
            IF (I006 .GT. MAXCOR) CALL INSMEM('DWMBEJNR', I006, MAXCOR)
C
            CALL GETALL(ICORE(I004), IREFSIZTA, 1, REFLISTTA)
            CALL SST003(ICORE(I004), ICORE(I002), IREFSIZTA, 
     &                  ITARSIZTA, ICORE(I005), SPNTYPT, "AJBI")
            
            I005 = I004 + IINTFP*IREFSIZTB
            I006 = I005 + ISCRSIZT
            IF (I006 .GT. MAXCOR) CALL INSMEM('DWMBEJNR', I006, MAXCOR)
C
            CALL GETALL(ICORE(I004), IREFSIZTB, 1, REFLISTTB)
            CALL SST002(ICORE(I004), ICORE(I003), IREFSIZTB, 
     &                  ITARSIZTB, ICORE(I005), SPNTYP2)
C     
            IOFFSTWA    = IBGN
            IOFFSTWB    = I001
            IOFFSTTA    = I002
            IOFFSTTB    = I003
            IOFFSTWQ    = IBGN
C
            DO 200 IRREP=1,NIRREP 
C
               DSSYWA = IRPDPD(IRREP, ISYTYP(1, TARLISTWA))
               NMSYWA = IRPDPD(IRREP, ISYTYP(2, TARLISTWA))
               DSSYTA = IRPDPD(IRREP, ISYTYP(1, TARLISTTA))
               NMSYTA = IRPDPD(IRREP, ISYTYP(2, TARLISTTA))
               DSSYWB = IRPDPD(IRREP, ISYTYP(1, TARLISTWB))
               NMSYWB = IRPDPD(IRREP, ISYTYP(2, TARLISTWB))
               DSSYTB = IRPDPD(IRREP, ISYTYP(1, TARLISTTB)) 
               NMSYTB = IRPDPD(IRREP, ISYTYP(2, TARLISTTB))
               DISSYQ = IRPDPD(IRREP, ISYTYP(1, LISTQ))
               NUMSYQ = IRPDPD(IRREP, ISYTYP(2, LISTQ))
C
               I000 = I004
               I010 = I000 + IINTFP*DISSYQ*NUMSYQ
               I020 = I010 + IINTFP*NMSYWA*DSSYWA
               IF (I020 .GT. MAXCOR) CALL INSMEM('DWMBEJNR', 
     &                               I020, MAXCOR)
C
               CALL IZERO(ICORE(I000),DISSYQ*NUMSYQ*IINTFP)
               NSIZE=DISSYQ*NUMSYQ
C
               IF(C1.NE.0.0)THEN
	          CALL DCOPY(NSIZE, ICORE(IOFFSTWQ), 1, 
     &                       ICORE(I000), 1)
                  IF(UCC) THEN
                     CALL SUMSYM2(ICORE(I000),ICORE(I010),NSIZE,1,
     &                            IRREP, LISTQ)
                     CALL IZERO(ICORE(I000),DISSYQ*NUMSYQ*IINTFP)
                  ENDIF  
               ENDIF
CSSS       call checksum("list 57:",icore(i000),nsize)
C
               IF(C4.NE.0.0)THEN
C
                  IF(C3.NE.0.0)THEN

                     CALL GETARGS(TARLISTTA,I1,I2)
                     CALL F2TAU(ICORE(IOFFSTTA),ICORE(IOFFT1(1,I1)),
     &                          ICORE(IOFFT1(1,I2)),DSSYTA,NMSYTA,
     &                          POP(1,I1),POP(1,I2),VRT(1,I1),
     &                          VRT(1,I2),IRREP,TWO,TAUTYP1)
                  ENDIF
C
                  CALL CNTRCTNR(ICORE(I000),ICORE(IOFFSTWA),
     &                          ICORE(IOFFSTTA),DISSYQ,NUMSYQ,DSSYTA,
     &                          NMSYTA,DSSYWA,NMSYWA,'WxT',X,C4,0) 
CSSS          call checksum("list 57 after 1st cntrctnr:",icore(i000),nsize)
C
                  IF(UCC) THEN
                     CALL SUMSYM2(ICORE(I000),ICORE(I010),NSIZE,
     &                            1,IRREP,LISTQ)
                     CALL DT(TWO,ICORE(I000),DISSYQ,NUMSYQ,'NSY',
     &                       ICORE(I010))
                     CALL SUMSYM2(ICORE(I010),ICORE(I000),NSIZE,1,
     &                            IRREP,LISTQ1)
                     CALL IZERO(ICORE(I000),DISSYQ*NUMSYQ*IINTFP)
                  ENDIF  
C
CSSS I020=I010+IINTFP*NMSYWB*DSSYWB
CSSS I030=I020+IINTFP*NMSYTB*DSSYTB
CSSS I030=MAX(I030,I010+DISSYQ*NUMSYQ*IINTFP)
CSSS CALL GETLST(ICORE(I020),1,NMSYTB,1,IRREP,LISTTB)
C
                  I020 = I010 + IINTFP*NMSYWB*DSSYWB
                  IF (I020 .GT. MAXCOR) CALL INSMEM('DWMBEJNR',
     &                                               I020, MAXCOR)
C
                  CALL CNTRCTNR(ICORE(I000),ICORE(IOFFSTWB),
     &                          ICORE(IOFFSTTB),DISSYQ,NUMSYQ,
     &                          DSSYTB,NMSYTB,DSSYWB,NMSYWB,
     &                          'WxT',X,-C4,0)
CSSS          call checksum("list 57 after 2nd cntrctnr:",icore(i000),nsize)
C
                  IF(UCC) THEN
                     CALL SUMSYM2(ICORE(I000),ICORE(I010),NSIZE,1,IRREP,
     &                            LISTQ)
                     CALL DT(TWO,ICORE(I000),DISSYQ,NUMSYQ,'NSY',
     &                       ICORE(I010))
                     CALL SUMSYM2(ICORE(I010),ICORE(I000),NSIZE,1,IRREP,
     &                            LISTQ1)
                  ENDIF  
               ENDIF
               IF(.NOT.UCC)
     &              CALL SUMSYM2(ICORE(I000),ICORE(I010),NSIZE,1,IRREP,
     &                           LISTQ)
C
	       IOFFSTWA = IINTFP*NMSYWA*DSSYWA + IOFFSTWA
	       IOFFSTWB = IINTFP*NMSYWB*DSSYWB + IOFFSTWB
	       IOFFSTTA = IINTFP*NMSYTA*DSSYTA + IOFFSTTA
               IOFFSTTB = IINTFP*NMSYTB*DSSYTB + IOFFSTTB
               IOFFSTWQ = IINTFP*DISSYQ*NUMSYQ + IOFFSTWQ

 200        CONTINUE
C
C SPIN ADAPTED CODE FOR RHF "ABAB" SPIN CASE
C
         ELSE
C
            LISTQ     = 56
            LSTINTB   = 25
CSSS            TARLSTINT = 18
            TARLISTTA = 37
            TARLISTTB = 39
            TARLISTWA = 18
            TARLISTWB = 21
            REFLISTWA = 16
            REFLISTTA = 46
            REFLISTTB = 46
            REFLISTWB = 16
            IASPIN    =  2
            TAUTYP1   = 'ABAB'

            ITARSIZWA = ISYMSZ(ISYTYP(1, TARLISTWA), 
     &                         ISYTYP(2, TARLISTWA))
            ITARSIZWB = ISYMSZ(ISYTYP(1, TARLISTWB), 
     &                         ISYTYP(2, TARLISTWB))
            ITARSIZTA = ISYMSZ(ISYTYP(1, TARLISTTA), 
     &                         ISYTYP(2, TARLISTTA))
            ITARSIZTB = ISYMSZ(ISYTYP(1, TARLISTTB), 
     &                         ISYTYP(2, TARLISTTB))
            IREFSIZWA = ISYMSZ(ISYTYP(1, REFLISTWA), 
     &                         ISYTYP(2, REFLISTWA))
            IREFSIZWB = ISYMSZ(ISYTYP(1, REFLISTWB), 
     &                         ISYTYP(2, REFLISTWB))
            IREFSIZTA = ISYMSZ(ISYTYP(1, REFLISTTA), 
     &                         ISYTYP(2, REFLISTTA))
            IREFSIZTB = ISYMSZ(ISYTYP(1, REFLISTTB), 
     &                         ISYTYP(2, REFLISTTB))
            ISCRSIZW  = NOCCO(ISPIN)*NVRTO(ISPIN) + NVRTO(ISPIN)*
     &                 (NVRTO(ISPIN)-1)/2 +
     &                  NOCCO(ISPIN)*(NOCCO(ISPIN)-1)/2
            ISCRSIZT  = NOCCO(1)*NOCCO(2) + NVRTO(1)*NVRTO(2) + 
     &                  NVRTO(1)*NOCCO(1) + NVRTO(2)*NOCCO(2)
C
          IBGN = 1
          I001 = IBGN + IINTFP*MAX(ITARSIZWA, ITARSIZWB)
          I002 = I001 + IINTFP*MAX(ITARSIZWB, ITARSIZWA)
          I003 = I002 + IINTFP*MAX(ITARSIZTA, ITARSIZTB)
          I004 = I003 + IINTFP*MAX(ITARSIZTB, ITARSIZTA)
          I005 = I004 + IINTFP*IREFSIZWA
          I006 = I005 + ISCRSIZT
          IF (I006 .GT. MAXCOR) CALL INSMEM('DWMBEJNR', I006, MAXCOR)
C
          CALL GETALL(ICORE(I004), IREFSIZWA, 1, REFLISTWA)
          CALL SST002(ICORE(I004), ICORE(IBGN), IREFSIZWA, ITARSIZWA, 
     &                ICORE(I005), "AABB")
          CALL SSTRNG(ICORE(IBGN), ICORE(I001), ITARSIZWA, ITARSIZWB, 
     &                ICORE(I005), "AABB")

          I004 = I003 + IINTFP*IREFSIZTA
          I006 = I005 + ISCRSIZT
          IF (I006 .GT. MAXCOR) CALL INSMEM('DWMBEJNR', I006, MAXCOR)
C
          CALL GETALL(ICORE(I004), IREFSIZTA, 1, REFLISTTA)
          CALL SST002(ICORE(I004), ICORE(I002), IREFSIZTA, ITARSIZTA, 
     &                ICORE(I005), "AABB")
          CALL SSTRNG(ICORE(I002), ICORE(I003), ITARSIZTA, ITARSIZTB,
     &                ICORE(I005), "AABB")
          
          IOFFSTWA    = IBGN
          IOFFSTWB    = I001
          IOFFSTTA    = I002
          IOFFSTTB    = I003
C
          DO 201 IRREP=1,NIRREP 
C
             DSSYWA = IRPDPD(IRREP, ISYTYP(1, TARLISTWA))
             NMSYWA = IRPDPD(IRREP, ISYTYP(2, TARLISTWA))
             DSSYTA = IRPDPD(IRREP, ISYTYP(1, TARLISTTA))
             NMSYTA = IRPDPD(IRREP, ISYTYP(2, TARLISTTA))
             DSSYWB = IRPDPD(IRREP, ISYTYP(1, TARLISTWB))
             NMSYWB = IRPDPD(IRREP, ISYTYP(2, TARLISTWB))
             DSSYTB = IRPDPD(IRREP, ISYTYP(1, TARLISTTB))
             NMSYTB = IRPDPD(IRREP, ISYTYP(2, TARLISTTB))
             DISSYQ = IRPDPD(IRREP, ISYTYP(1, LISTQ))
             NUMSYQ = IRPDPD(IRREP, ISYTYP(2, LISTQ))
C
             I000 = I004 
             I010 = I000 + IINTFP*DISSYQ*NUMSYQ
             I020 = I010 + IINTFP*NMSYWA*DSSYWA
             I030 = I020 + IINTFP*NMSYTA*DSSYTA
             I040 = I030 + IINTFP*NUMSYQ*DISSYQ
             IF (I040 .GT. MAXCOR) CALL INSMEM('DWMBEJNR', I040, 
     &                             MAXCOR)             
C
C SPIN ADAPT THE T AMPLITUDES
CSSS I020=I010+IINTFP*NMSYWA*DSSYWA
CSSS I030=I020+IINTFP*NMSYTA*DSSYTA
CSSS CALL GETLST(ICORE(I020),1,NMSYTA,1,IRREP,LISTTA)
CSSS CALL GETLST(ICORE(I010),1,NMSYTA,1,IRREP,LISTTB)
C
             CALL DCOPY(DSSYWA*NMSYWA, ICORE(IOFFSTWA), 1,
     &                  ICORE(I010), 1)
C
             IF(C4.NE.0.0)THEN
C
                IF(C3.NE.0.0)THEN
C
                   CALL GETARGS(TARLISTTA,I1,I2)
                   CALL F2TAU(ICORE(IOFFSTTB),ICORE(IOFFT1(1,I1)),
     &                        ICORE(IOFFT1(1,I2)),DSSYTA,NMSYTA,
     &                        POP(1,I1),POP(1,I2),VRT(1,I1),VRT(1,I2),
     &                        IRREP,TWO,TAUTYP1)
                ENDIF
C
                CALL SSCAL (NMSYTA*NMSYTA,TWO,ICORE(IOFFSTTA),1)
                CALL SAXPY (NMSYTA*NMSYTA,ONEM,ICORE(IOFFSTTB),1,
     &                      ICORE(IOFFSTTA),1)
C
                CALL CNTRCTNR(ICORE(IOFFSTWB),ICORE(IOFFSTWA),
     &                        ICORE(IOFFSTTA),DISSYQ,NUMSYQ,
     &                        DSSYTA,NMSYTA,DSSYWA,NMSYWA, 
     &                        'WxT',ZILCH,-C4,1)
C
	        CALL DCOPY(NUMSYQ*DISSYQ, ICORE(IOFFSTWB), 1,
     &                     ICORE(I000), 1) 
CSSS      call checksum("list 56 after cntrctnr:",icore(i000),numsyq*dissyq)
C
                IF(UCC) THEN
                   CALL DT(TWO,ICORE(I000),DISSYQ,NUMSYQ,'ADD',
     &                     ICORE(I030))
                ENDIF  
             ENDIF
C
             NSIZE=DISSYQ*NUMSYQ
C
C NOW INCLUDE THE SPIN-ADAPTED INTEGRALS
C
             IF(C1.NE.0.0)THEN
C
CSSS CALL GETLST(ICORE(I010),1,NUMSYQ,1,IRREP,LSTINT)
C
                CALL GETLST(ICORE(I020),1,NUMSYQ,1,IRREP,LSTINTB)
                CALL SSCAL (DISSYQ*DISSYQ,TWO,ICORE(I010),1)
                CALL SAXPY (DISSYQ*DISSYQ,ONEM,ICORE(I020),1,
     &                      ICORE(I010),1) 
C
                IF(C4.NE.0.0)THEN
                   CALL SAXPY (DISSYQ*DISSYQ,ONE,ICORE(I010),1,
     &                         ICORE(I000),1)
                ELSE
c YAU : old
c CALL ICOPY(DISSYQ*DISSYQ*IINTFP,ICORE(I010),1,ICORE(I000),1)
c YAU : new
                   CALL DCOPY(DISSYQ*DISSYQ,ICORE(I010),1,
     &                        ICORE(I000),1)
c YAU : end
                ENDIF
             ENDIF
CSSS      call checksum("list 56 after spina1:",icore(i000),numsyq*dissyq)
C
C WE ALSO NEED SPIN-ADAPTED T1RING STUFF FOR CCSD
C
c     IF(C3.NE.0.0)THEN
CJDW 3/20/98. The T1RING stuff must also be picked up when we have
C             non-HF SDQ-MBPT(4) or MBPT(4). Also, require C4.NE.0
C             so that we know we are calculating second-order MBEJ,
C             which implies that T1RING term has been calculated
C             and must be included.
C
             IF(C3.NE.0.0 .OR.
     &         ( IFLAGS(38).NE.0 .AND. 
     &         (IFLAGS(2).EQ.3.OR.IFLAGS(2).EQ.4) .AND.
     &         C4.NE.0)                                )THEN
C
                CALL GETLST(ICORE(I010), 1, NUMSYQ, 1, IRREP, 56)
                CALL GETLST(ICORE(I020), 1, NUMSYQ, 1, IRREP, 58)
                CALL SSCAL(NSIZE, TWO, ICORE(I010), 1)
                CALL SAXPY(NSIZE, ONEM, ICORE(I020), 1, ICORE(I010), 1)
                CALL SAXPY(NSIZE, ONE , ICORE(I010), 1, ICORE(I000), 1)
                CALL PUTLST(ICORE(I000), 1, NUMSYQ, 1, IRREP, 56)
             ELSE
                CALL SUMSYM2(ICORE(I000), ICORE(I010), NSIZE, 1, 
     &                       IRREP, LISTQ)
             ENDIF
CSSS      call checksum("list 56 after spina2:",icore(i000),numsyq*dissyq)
C  
               IOFFSTWA = IINTFP*NMSYWA*DSSYWA + IOFFSTWA
               IOFFSTWB = IINTFP*NMSYWB*DSSYWB + IOFFSTWB
               IOFFSTTA = IINTFP*NMSYTA*DSSYTA + IOFFSTTA
               IOFFSTTB = IINTFP*NMSYTB*DSSYTB + IOFFSTTB

 201         CONTINUE
          ENDIF
C
       ELSEIF(SPCASE.EQ.'ABBA'.OR.SPCASE.EQ.'BAAB')THEN
C
C DO SPIN CASES ABBA AND BAAB:
C
C       W(MbeJ) = -<Mb|Je>+(1/2)*SUM [T(Jn,Fb)<Mn|Fe>]
C
C       W(mBEj) = -<mB|jE>+(1/2)*SUM [T(jN,fB)<mN|fE>]
C
          IF(SPCASE.EQ.'ABBA')THEN
C
             LISTQ     = 58
             LSTINT    = 25
             TARLISTTA = 38
             TARLISTTB = 36
             TARLISTWA = 21
             TARLISTWB = 18
             REFLISTTA = 46
             REFLISTWA = 16
             TAUTYP1     = 'ABAB'
             SPNCASEW(1) = 'AABB'
             SPNCASET(1) = 'BBAA'
             ISPIN     = 1 
             IF (IUHF .EQ. 0) THEN 
                TARLISTTA = 39
                TARLISTTB = 37
                SPNCASET(1)  = 'AABB'
             ENDIF
C
          ELSEIF(SPCASE.EQ.'BAAB')THEN
C             
             LISTQ     = 59
             LSTINT    = 26
             TARLISTTA = 39
             TARLISTTB = 37
             TARLISTWA = 22
             TARLISTWB = 17
             REFLISTTA = 46             
             REFLISTWA = 16
             TAUTYP1      = 'BABA'
             SPNCASEW(2)  = 'BBAA'
             SPNCASET(2)  = 'AABB'
             ISPIN  = 2
          ENDIF
C
          ITARSIZWA = ISYMSZ(ISYTYP(1, TARLISTWA), 
     &                       ISYTYP(2, TARLISTWA))
          ITARSIZWB = ISYMSZ(ISYTYP(1, TARLISTWB), 
     &                       ISYTYP(2, TARLISTWB))
          ITARSIZTA = ISYMSZ(ISYTYP(1, TARLISTTA), 
     &                       ISYTYP(2, TARLISTTA))
          ITARSIZTB = ISYMSZ(ISYTYP(1, TARLISTTB), 
     &                       ISYTYP(2, TARLISTTB))
          IREFSIZWA = ISYMSZ(ISYTYP(1, TARLISTWA), 
     &                       ISYTYP(2, TARLISTWA))
          IREFSIZTA = ISYMSZ(ISYTYP(1, REFLISTTA), 
     &                       ISYTYP(2, REFLISTTA))
          ISCRSIZW  = NOCCO(ISPIN)*NVRTO(ISPIN) + NVRTO(ISPIN)*
     &                (NVRTO(ISPIN)-1)/2 +
     &                NOCCO(ISPIN)*(NOCCO(ISPIN)-1)/2
          ISCRSIZT  = NOCCO(1)*NOCCO(2) + NVRTO(1)*NVRTO(2) + 
     &                NVRTO(1)*NOCCO(1) + NVRTO(2)*NOCCO(2)
C
         IBGN = 1
         I001 = IBGN + IINTFP*MAX(ITARSIZWA, ITARSIZWB)
         I002 = I001 + IINTFP*MAX(ITARSIZWA, ITARSIZWB) 
         I003 = I002 + IINTFP*IREFSIZWA
         I004 = I003 + ISCRSIZT
         IF (I004 .GT. MAXCOR) CALL INSMEM('DWMBEJNR', I004, MAXCOR)
C
          CALL GETALL(ICORE(I002), IREFSIZWA, 1, REFLISTWA)
          CALL SST002(ICORE(I002), ICORE(IBGN), IREFSIZWA, ITARSIZWB, 
     &                ICORE(I003), SPNCASEW(ISPIN))
          CALL SSTRNG(ICORE(IBGN), ICORE(I001), ITARSIZWB, ITARSIZWA, 
     &                ICORE(I003), SPNCASEW(ISPIN))
          CALL DCOPY(ITARSIZWA, ICORE(I001), 1, ICORE(IBGN), 1)
          
          I002 = I001 + IINTFP*MAX(ITARSIZTA, ITARSIZTB) 
          I003 = I002 + IINTFP*MAX(ITARSIZTA, ITARSIZTB)
          I004 = I003 + ISCRSIZT
          IF (I004 .GT. MAXCOR) CALL INSMEM('DWMBEJNR', I004, MAXCOR)
C
          CALL GETALL(ICORE(I002), IREFSIZTA, 1, REFLISTTA)
          CALL SST002(ICORE(I002), ICORE(I001), IREFSIZTA, ITARSIZTB, 
     &                ICORE(I003), SPNCASET(ISPIN))
          CALL SSTRNG(ICORE(I001), ICORE(I002), ITARSIZTB, ITARSIZTA,
     &                ICORE(I003), SPNCASET(ISPIN))
          CALL DCOPY(ITARSIZTA, ICORE(I002), 1, ICORE(I001), 1)
C          
          IOFFSTWA    = IBGN
          IOFFSTTA    = I001
C
          DO 300 IRREP=1,NIRREP 
C
             DSSYWA = IRPDPD(IRREP, ISYTYP(1, TARLISTWA))
             NMSYWA = IRPDPD(IRREP, ISYTYP(2, TARLISTWA))
             DSSYTA = IRPDPD(IRREP, ISYTYP(1, TARLISTTA))
             NMSYTA = IRPDPD(IRREP, ISYTYP(2, TARLISTTA))
             DISSYQ = IRPDPD(IRREP, ISYTYP(1, LISTQ))
             NUMSYQ = IRPDPD(IRREP, ISYTYP(2, LISTQ))
C
             I000 = I002
             I010 = I000 + IINTFP*DISSYQ*NUMSYQ
             I020 = I010 + IINTFP*NMSYWA*DSSYWA
             I030 = I020 + IINTFP*DISSYQ*NUMSYQ
             IF (I030 .GT. MAXCOR) CALL INSMEM('DWMBEJNR', I030,
     &                             MAXCOR)
CSSS I030=I020+IINTFP*NMSYTA*DSSYTA
CSSS I030=MAX(I030,I010+DISSYQ*NUMSYQ*IINTFP)
C     
             CALL IZERO(ICORE(I000),DISSYQ*NUMSYQ*IINTFP)
             NSIZE=DISSYQ*NUMSYQ
C
C LIST ORDERING OF TARGET IS SAME AS INTEGRAL INVOLVED IN
C CONTRACTION, SO JUST USE THIS SET OF POINTERS TO FETCH 
C PACKED INTEGRAL LISTS.
C
             IF(C1.NE.0.0)THEN
                CALL GETLST(ICORE(I000),1,NUMSYQ,2,IRREP,LSTINT)
                IF(UCC) THEN
                   CALL SUMSYM2(ICORE(I000),ICORE(I010),NSIZE,1,
     &                          IRREP,LISTQ)
                   CALL IZERO(ICORE(I000),DISSYQ*NUMSYQ*IINTFP)
                ENDIF  
             ENDIF
C
             IF(C4.NE.0.0)THEN
C
CSSS CALL GETLST(ICORE(I020),1,NMSYTA,1,IRREP,LISTTA)
C
                IF(C3.NE.0.0)THEN
                   CALL GETARGS(TARLISTTA,I1,I2)
                   CALL F2TAU(ICORE(IOFFSTTA),ICORE(IOFFT1(1,I1)),
     &                        ICORE(IOFFT1(1,I2)),DSSYTA,NMSYTA,
     &                        POP(1,I1),POP(1,I2),VRT(1,I1),VRT(1,I2),
     &                        IRREP,TWO,TAUTYP1)
                ENDIF
C
                CALL CNTRCTNR(ICORE(I000),ICORE(IOFFSTWA),
     &                        ICORE(IOFFSTTA),DISSYQ,NUMSYQ,
     &                        DSSYTA,NMSYTA,DSSYWA,NMSYWA,
     &                        'WxT',X,C4,0)
                IF(UCC) THEN
C
CSSS CALL SUMSYM2(ICORE(I000),ICORE(I010),NSIZE,1,IRREP,LISTQ)
cSSS CALL DT(ONE,ICORE(I000),DISSYQ,NUMSYQ,'EXT',icore(i030))
C
                   CALL DT(TWO,ICORE(I000),DISSYQ,NUMSYQ,'ADD',
     &                     ICORE(I030))
                ENDIF  
             ENDIF
             CALL SUMSYM2(ICORE(I000),ICORE(I010),NSIZE,1,IRREP,LISTQ)

             IOFFSTWA = IINTFP*NMSYWA*DSSYWA + IOFFSTWA
             IOFFSTTA = IINTFP*NMSYTA*DSSYTA + IOFFSTTA
C
 300      CONTINUE
       ENDIF
C
       RETURN
       END
