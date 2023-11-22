      SUBROUTINE G4ALL(ICORE,MAXCOR,SPCASE,IUHF,bRedundant)
C
C  THIS ROUTINE AND DEPENDENTS COMPUTE THE RING-TYPE G(IA,JB)
C  INTERMEDIATE FOR ALL SIX POSSIBLE SPIN CASES.  THE ALGORITHM
C  ASSUMES IN-CORE STORAGE OF SYMMETRY PACKED TARGET AND T2 LISTS
C  VECTORS AND USES THE DPD SYMMETRY APPROACH TO EVALUATE THE
C  CONTRACTIONS.
C
C  MBPT(3) :
C
C  G(IA,JB) = - SUM M,E T[1](IM,BE) T[1](MJ,EA)
C
C  ROHF-MBPT(3) :
C
C  G(IA,JB) = - SUM M,E T[1](IM,BE) T[1](MJ,EA)
C
C             - T[1](I,B) T[1](J,A)]
C
C  MBPT(4) :
C
C  G(IA,JB) = - P(IA,JB) SUM M,E [2 T2(IM,BE) - T1(IM,BE) ] T1(MJ,EA)
C
C  CCD AND UCC :
C
C  G(IA,JB) = - 1/2 P(IA,JB) H(IB,JA)
C
C  QCISD : 
C
C  G(IA,JB) = - 1/2 P(IA,JB) H(IB,JA)
C
C             - 1/2 P(IA,JB) [T(I,B) L(J,A)]
C
C  CCSD :
C
C  G(IA,JB) = -1/2 P(IA,JB)   H(IB,JA) 
C 
C             -1/2 P(IA,JB) T(I,B) L(J,A)
C
C             +1/2 P(IA,JB) SUM M,E T(I,E) T(M,B) L(MJ,EA)
C
CEND
C
C  CODED AUGSUT/90  JG
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL ISAME,CIS,EOM,LTRP
      LOGICAL MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC
      LOGICAL CC,MBPT4
      LOGICAL DENS,GRAD,QRHF,NONHF,ROHF,SEMI,CANON,TRIP1,TRIP2,
     &        GABCD,RELAXED,TRULY_NONHF
      logical bRedundant
      INTEGER DIRPRD,DISSYG,DSSYT1A,DSSYT2A,DSSYT1B,DSSYT2B,
     &        DSSYT1,DSSYT2,DSSYT3A,DSSYT3B,DSSYT4A,DSSYT4B,
     &        DISSYG2,POP,VRT  
      CHARACTER*4 SPCASE
      DIMENSION ICORE(MAXCOR)
      COMMON/DERIV/DENS,GRAD,QRHF,NONHF,ROHF,SEMI,CANON,TRIP1,
     &              TRIP2,GABCD,RELAXED,TRULY_NONHF
      COMMON/METH/MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/FLAGS/IFLAGS(100)
      COMMON/SYMINF/NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
      COMMON/EXCITE/CIS,EOM
      COMMON/LTRIP/LTRP
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),NTOT(18)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NF1(2),NF2(2)
C
      DATA AZERO,HALF,ONE,ONEM,TWO /0.0D+0,0.5D0,1.0D+0,-1.D+0,2.D0/
C
      MXCOR=MAXCOR
C
      MBPT4=M4DQ.OR.M4SDQ.OR.M4SDTQ
      CC=CCD.OR.CCSD.OR.QCISD.OR.UCC
C
      IOFFLISTH=0
      IF(IFLAGS(3).EQ.2) IOFFLISTH=200
C
C DO SPIN CASE AAAA AND BBBB:
C
C  MBPT(3) :    G(IAJB) = - SUM M,E T[1](IM,BE) T[1](JM,AE) 
C                                   + T[1](Im,Be) T[1](Jm,Ae)
C
C  MBPT(4) :    G(IAJB) AS MBPT(3) WITH TWO DIFFERENT AMPLITUDES 
C
C  CCD  :    ....
C
C  QCISD : ... IN ADDITION - 1/2 P((IA),(JB) L(IB) T(JA)
C
C  STORAGE MODE :   BI,AJ    ACTUALLY  = T1(MI,BE) T1(MJ,AE)
C                                  STORAGE HERE   BI,ME AND ME, AJ
C                                  SIGN = ONEM 
C                            T1(Im,Be) T1(Jm,Ae)  STORAGE HERE BI, me
C                                   SIGN = ONEM
C  AAAA :     RHF AND UHF
C
C  BBBB :     UHF ONLY           
C
      IF(SPCASE.EQ.'AAAA'.OR.SPCASE.EQ.'BBBB')THEN
       IF(SPCASE.EQ.'AAAA')THEN
        ISPIN=1
       ELSEIF(SPCASE.EQ.'BBBB')THEN
        ISPIN=2
       ENDIF
       MAXSIZ=0
       LISTG =122+ISPIN
C
C CC OR CIS CODE
C
       IF(CC.OR.CIS) THEN
        LISTH=53+ISPIN+IOFFLISTH
        DO 90 IRREP=1,NIRREP
        DISSYG=IRPDPD(IRREP,ISYTYP(1,LISTG))
        NUMSYG=IRPDPD(IRREP,ISYTYP(2,LISTG))
        I000=1
        I010=I000+IINTFP*NUMSYG*DISSYG
        IF(.NOT.CIS)THEN
         CALL GETLST(ICORE(I000),1,NUMSYG,1,IRREP,LISTH)
        ELSE
         CALL ZERO(ICORE(I000),NUMSYG*DISSYG)
        ENDIF
C
C FOR QCISD AND CCSD ADD SINGLE CONTRIBUTION
C
        IF(CIS)THEN
         IONE=1
         CALL GETREC(20,'JOBARC','IOPTSYM ',IONE,IOPTSYM)
        ENDIF
        IF(QCISD.AND.IRREP.EQ.1.OR.CCSD.AND.IRREP.EQ.1
     &    .OR.CIS.AND.IRREP.EQ.IOPTSYM)THEN
         CALL T1ING4(ICORE(I000),DISSYG,NUMSYG,ICORE(I010),
     &               ONEM,ISPIN,ISPIN,90,190)
        ENDIF
        IF (EOM .OR. LTRP)THEN
         CALL GETLST(ICORE(I010),1,NUMSYG,1,IRREP,LISTG)
         CALL SAXPY (NUMSYG*DISSYG,ONE,ICORE(I010),1,ICORE(I000),1)
        ENDIF
C
C  SYMMETRIZE THE GAMMA INTERMEDIATE
C
        CALL SYMMET2(ICORE(I000),DISSYG)
C
C  SAVE GAMMA4 INTERMEDIATE ON DISK
C
        CALL PUTLST(ICORE(I000),1,NUMSYG,2,IRREP,LISTG)
C
90     CONTINUE
C
       ELSE
C
C MBPT CODE
C       
       IF(MBPT3) THEN
        ISAME=.TRUE.
        FACT=ONEM
        LIST1A=33+ISPIN
        LIST2A=33+ISPIN
        LIST1B=38-ISPIN
        LIST2B=38-ISPIN
       ELSE IF(MBPT4) THEN
        LIST1A=33+ISPIN
        LIST2A=133+ISPIN
        LIST1B=38-ISPIN
        LIST2B=138-ISPIN
        ISAME=.FALSE.
        FACT=ONEM
       ENDIF
       DO 100 IRREP=1,NIRREP 
        DSSYT1A=IRPDPD(IRREP,ISYTYP(1,LIST1A))
        DSSYT2A=IRPDPD(IRREP,ISYTYP(1,LIST2A))
        NMSYT1A=IRPDPD(IRREP,ISYTYP(2,LIST1A))
        NMSYT2A=IRPDPD(IRREP,ISYTYP(2,LIST2A))
        DSSYT1B=IRPDPD(IRREP,ISYTYP(1,LIST1B))
        DSSYT2B=IRPDPD(IRREP,ISYTYP(1,LIST2B))
        NMSYT1B=IRPDPD(IRREP,ISYTYP(2,LIST1B))
        NMSYT2B=IRPDPD(IRREP,ISYTYP(2,LIST2B))
        DISSYG=IRPDPD(IRREP,ISYTYP(1,LISTG))
        NUMSYG=IRPDPD(IRREP,ISYTYP(2,LISTG))
C
C FIRST TERM
C
        I000=1
        I010=I000+IINTFP*DISSYG*NUMSYG
        IF(ISAME) THEN
         I020=I010
        ELSE
         I020=I010+IINTFP*NMSYT1A*DSSYT1A
        ENDIF
        I030=I020+IINTFP*NMSYT2A*DSSYT2A
        IF(MXCOR.LT.I030) CALL INSMEM('G4ALL',I030,MXCOR)
        CALL IZERO(ICORE(I000),DISSYG*NUMSYG*IINTFP)
        IF (bRedundant) THEN
           CALL GETLST(ICORE(I010),1,NMSYT1A,1,IRREP,LIST1A)
           IF (.NOT.ISAME) CALL GETLST(ICORE(I020),1,NMSYT2A,2,
     &                                 IRREP,LIST2A)
        ELSE
           CALL GETLST_NR(ICORE(I010),ICORE(I030),MXCOR-I030+1,
     &                    LIST1A,IRREP)
           IF (.NOT.ISAME) CALL GETLST_NR(ICORE(I020),ICORE(I030),
     &                                    MXCOR-I030+1,LIST2A,IRREP)
        END IF
C
C  FOR MBPT(4) FORM APPROBIATE LINEAR COMBINATION OF AMPLITUDES
C
        IF(MBPT4) THEN
         CALL SSCAL(NMSYT2A*DSSYT2A,TWO,ICORE(I020),1)
         CALL SAXPY(NMSYT2A*DSSYT2A,ONE,ICORE(I010),1,ICORE(I020),1)
        ENDIF
C
        CALL XGEMM('N','N',DISSYG,NUMSYG,DSSYT2A,FACT,ICORE(I010),
     &             DSSYT1A,ICORE(I020),DSSYT2A,AZERO,ICORE(I000),
     &             DISSYG)
C
C SECOND TERM
C
        IF(ISAME) THEN
         I020=I010
        ELSE
         I020=I010+IINTFP*NMSYT1B*DSSYT1B
        ENDIF
        I030=I020+IINTFP*NMSYT2B*DSSYT2B
        IF (bRedundant) THEN
           CALL GETLST(ICORE(I010),1,NMSYT1B,1,IRREP,LIST1B)
           IF (.NOT.ISAME) CALL GETLST(ICORE(I020),1,NMSYT2B,2,
     &                                 IRREP,LIST2B)
        ELSE
           CALL GETLST_NR(ICORE(I010),ICORE(I030),MXCOR-I030+1,
     &                    LIST1B,IRREP)
           IF (.NOT.ISAME) CALL GETLST_NR(ICORE(I020),ICORE(I030),
     &                                    MXCOR-I030+1,LIST2B,IRREP)
        END IF
C
C  FOR MBPT(4) FORM APPROBIATE LINEAR COMBINATION OF AMPLITUDES
C
        IF(MBPT4) THEN
         CALL SSCAL(NMSYT2B*DSSYT2B,TWO,ICORE(I020),1)
         CALL SAXPY(NMSYT2B*DSSYT2B,ONE,ICORE(I010),1,ICORE(I020),1)
        ENDIF
        CALL XGEMM('N','T',DISSYG,NUMSYG,NMSYT2B,FACT,ICORE(I010),
     &             DSSYT1B,ICORE(I020),DSSYT2B,ONE,ICORE(I000),
     &             DISSYG)
C
C SINGLES CONTRIBUTION IN ROHF-MBPT(3)
C
        IF(ROHF.AND.MBPT3.AND.IRREP.EQ.1) THEN
         CALL T1ING4(ICORE(I000),DISSYG,NUMSYG,ICORE(I010),
     &               ONEM,ISPIN,ISPIN,90,90)
        ENDIF
        IF(EOM .OR. LTRP)THEN
         CALL GETLST(ICORE(I010),1,NUMSYG,1,IRREP,LISTG)
         CALL SAXPY (NUMSYG*DISSYG,ONE,ICORE(I010),1,ICORE(I000),1)
        ENDIF
C
C  SYMMETRIZE THE GAMMA INTERMEDIATE, NOT REQUIRED FOR MBPT3
C
        IF(.NOT.ISAME) CALL SYMMET2(ICORE(I000),DISSYG)
C
C  SAVE GAMMA4 INTERMEDIATE ON DISK
C
        CALL PUTLST(ICORE(I000),1,NUMSYG,2,IRREP,LISTG)
100    CONTINUE
C
       ENDIF
C
C DO SPIN CASES ABBA AND BAAB:
C
C MBPT(3) :  G(IajB) = - SUM M,E[ T1(IM,BE) T1(jM,aE) + T1(Im,Be) T1(jm,ae)]
C
C
C  FOR QCISD :  ... IN ADDITION -1/2 P((IA),(JB)) L(IB) T(JA)
C
C  THIS CORRESPONDS TO <Ia//jB> = <Ia/jB> = - <Ij//Ba> (LIST 18)
C
C STORAGE : BI,aj (ACTUALLY HERE WE ARE INTERESTED IN THE NEGATIVE OF
C                  THIS TERM (<Ia//jB> =  - <Ia//Bj> = - <Ba//Ij>)
C
C      ACTUALLY T1(IM,BE) T1(jM,aE) = - T1(MI,BE) T1(Mj,Ea)
C
C      STORAGE    BI, ME ; ME, aj    OVERALL SIGN :  ONEM
C
C      ACTUALLY  T1(Im,Be) T1(jm,ae) = - T1(Im,Be) T1(jm,ea)
C
C      STORAGE    Bi, me ;  me,  aj   OVERALL SIGN : ONEM
C
C  NOTE THE SECOND TERM IS FOR MBPT(3)  IN RHF  TRANSPOSITION OF THE
C  FIRST TERM AND THEREFORE NOT REQUIRED
C
C  ALSO AT LEAST FOR MBPT(3)  ONLY THE ABAB CASE IS REQUIRED, SINCE
C  WE HAVE SYMMETRY WITH RESPECT TO THE INDICES
C
C FOR ALL OTHER CASES : RHF CASE : CALCULATE ONLY THE 118 LIST AND
C                                  SYMMETRIZE
C                       UHF CASE : CALCULATE 118 LIST AND THEN
C                                  THE 117 LIST    
C    
      ELSEIF(SPCASE.EQ.'ABBA'.OR.SPCASE.EQ.'BAAB')THEN
       IF(SPCASE.EQ.'ABBA')THEN
        LISTG=118
        IF(MBPT3) THEN
         FACT=ONEM
         LIST1A=34
         LIST2A=37
         LIST1B=37
         LIST2B=34+IUHF
        ELSE IF(MBPT4) THEN
         FACT=ONEM
         LIST1A=34
         LIST2A=137
         LIST5A=37
         LIST1B=37
         LIST2B=134+IUHF
         LIST5B=34+IUHF
        ELSE IF(CC.OR.CIS) THEN
         LISTH=56+IOFFLISTH
        ENDIF
       ENDIF
       IF(SPCASE.EQ.'ABBA'.AND.IUHF.EQ.1.AND..NOT.MBPT3)THEN 
        LISTG2=117
        FACT=ONEM
        LIST3A=35
        LIST4A=136
        LIST6A=36
        LIST3B=36
        LIST4B=134
        LIST6B=34
       ENDIF
       IF(CC.OR.CIS) THEN
        DO 190 IRREP=1,NIRREP
        DISSYG=IRPDPD(IRREP,ISYTYP(1,LISTG))
        NUMSYG=IRPDPD(IRREP,ISYTYP(2,LISTG))
        I000=1
        I010=I000+IINTFP*NUMSYG*DISSYG
        IF(.NOT.CIS)THEN
         CALL GETLST(ICORE(I000),1,NUMSYG,1,IRREP,LISTH)
CSSS         call checksum("H4ABBA  :",ICORE(I000),DISSYG*NUMSYG)
        ELSE
         CALL ZERO(ICORE(I000),NUMSYG*DISSYG)
        ENDIF
C
C FOR QCISD AND CCSD ADD SINGLE CONTRIBUTION
C
        IF(CIS)THEN
         IONE=1
         CALL GETREC(20,'JOBARC','IOPTSYM ',IONE,IOPTSYM)
        ENDIF
        IF(QCISD.AND.IRREP.EQ.1.OR.CCSD.AND.IRREP.EQ.1
     &    .OR.CIS.AND.IRREP.EQ.IOPTSYM)THEN
         CALL T1ING4(ICORE(I000),DISSYG,NUMSYG,ICORE(I010),
     &               ONE,1,1+IUHF,90,190)
        ENDIF
        IF(IUHF.EQ.0) THEN
         IF(EOM .OR. LTRP)THEN
          CALL GETLST(ICORE(I010),1,NUMSYG,1,IRREP,LISTG)
          CALL SAXPY (NUMSYG*DISSYG,ONE,ICORE(I010),1,ICORE(I000),1)
         ENDIF
         CALL SYMMET2(ICORE(I000),DISSYG)
         CALL PUTLST(ICORE(I000),1,NUMSYG,2,IRREP,LISTG)
        ELSE
         CALL TRANSP(ICORE(I000),ICORE(I010),NUMSYG,DISSYG)
         IF(.NOT.CIS)THEN
          CALL GETLST(ICORE(I000),1,DISSYG,1,IRREP,LISTH+1)
         ELSE
          CALL ZERO(ICORE(I000),DISSYG*NUMSYG)
         ENDIF
         CALL SAXPY(NUMSYG*DISSYG,ONE,ICORE(I010),1,ICORE(I000),1)
C
C FOR QCISD AND CCSD ADD SINGLE CONTRIBUTION
C
        IF(CIS)THEN
         IONE=1
         CALL GETREC(20,'JOBARC','IOPTSYM ',IONE,IOPTSYM)
        ENDIF
         IF(QCISD.AND.IRREP.EQ.1.OR.CCSD.AND.IRREP.EQ.1
     &    .OR.CIS.AND.IRREP.EQ.IOPTSYM)THEN
          CALL T1ING4(ICORE(I000),NUMSYG,DISSYG,ICORE(I010),
     &                ONE,2,1,90,190)
         ENDIF
C
C  SCALE BY 1/2
C
         CALL SSCAL(NUMSYG*DISSYG,HALF,ICORE(I000),1)
C
C  SAVE 117 ON DISK
C
         IF(EOM .OR. LTRP)THEN
          CALL GETLST(ICORE(I010),1,DISSYG,1,IRREP,LISTG-1)
          CALL SAXPY (NUMSYG*DISSYG,ONE,ICORE(I010),1,ICORE(I000),1)
         ENDIF
         CALL PUTLST(ICORE(I000),1,DISSYG,2,IRREP,LISTG-1)
C
C  FORM 118 BY TRANSPOSITION AND SAVE ON DISK
C
         CALL TRANSP(ICORE(I000),ICORE(I010),DISSYG,NUMSYG)
         CALL PUTLST(ICORE(I010),1,NUMSYG,2,IRREP,LISTG)
        ENDIF
190    CONTINUE
C
       ELSE
C
       DO 200 IRREP=1,NIRREP 
        DSSYT1A=IRPDPD(IRREP,ISYTYP(1,LIST1A))
        DSSYT2A=IRPDPD(IRREP,ISYTYP(1,LIST2A))
        NMSYT1A=IRPDPD(IRREP,ISYTYP(2,LIST1A))
        NMSYT2A=IRPDPD(IRREP,ISYTYP(2,LIST2A))
        DSSYT1B=IRPDPD(IRREP,ISYTYP(1,LIST1B))
        DSSYT2B=IRPDPD(IRREP,ISYTYP(1,LIST2B))
        NMSYT1B=IRPDPD(IRREP,ISYTYP(2,LIST1B))
        NMSYT2B=IRPDPD(IRREP,ISYTYP(2,LIST2B))
        DISSYG=IRPDPD(IRREP,ISYTYP(1,LISTG))
        NUMSYG=IRPDPD(IRREP,ISYTYP(2,LISTG))
C
C  FIRST TERM
C
        I000=1
        I010=I000+IINTFP*DISSYG*NUMSYG
        I020=I010+IINTFP*MAX(NMSYT1A*DSSYT1A,NMSYT2A*DSSYT2A)
        I030=I020+IINTFP*NMSYT2A*DSSYT2A
        CALL IZERO(ICORE(I000),DISSYG*NUMSYG*IINTFP)
        IF (bRedundant) THEN
           CALL GETLST(ICORE(I020),1,NMSYT2A,1,IRREP,LIST2A)
        ELSE
           CALL GETLST_NR(ICORE(I020),ICORE(I030),MXCOR-I030+1,
     &                    LIST2A,IRREP)
        END IF
        IF(MBPT4) THEN
         IF (bRedundant) THEN
            CALL GETLST(ICORE(I010),1,NMSYT2A,1,IRREP,LIST5A)
         ELSE
            CALL GETLST_NR(ICORE(I010),ICORE(I030),MXCOR-I030+1,
     &                     LIST5A,IRREP)
         END IF
         CALL SSCAL(NMSYT2A*DSSYT2A,TWO,ICORE(I020),1)
         CALL SAXPY(NMSYT2A*DSSYT2A,ONE,ICORE(I010),1,ICORE(I020),1)
        ENDIF
        IF (bRedundant) THEN
           CALL GETLST(ICORE(I010),1,NMSYT1A,1,IRREP,LIST1A)
        ELSE
           CALL GETLST_NR(ICORE(I010),ICORE(I030),MXCOR-I030+1,
     &                    LIST1A,IRREP)
        END IF
        CALL XGEMM('N','N',DISSYG,NUMSYG,DSSYT2A,FACT,ICORE(I010),
     &             DSSYT1A,ICORE(I020),DSSYT2A,AZERO,ICORE(I000),
     &             DISSYG)
C
C IN MBPT3 WITH RHF-REFERENCE THE SECOND IS THE TRANSPOSE OF THE FIRST
C
        IF(IUHF.NE.0.OR..NOT.MBPT3) THEN
         I020=I010+IINTFP*MAX(NMSYT1B*DSSYT1B,NMSYT2B*DSSYT2B)
         I030=I020+IINTFP*NMSYT2B*DSSYT2B
         IF(MXCOR.LT.I030) CALL INSMEM('G4ALL',I030,MXCOR)
         IF (bRedundant) THEN
            CALL GETLST(ICORE(I020),1,NMSYT2B,1,IRREP,LIST2B)
         ELSE
            CALL GETLST_NR(ICORE(I020),ICORE(I030),MXCOR-I030+1,
     &                     LIST2B,IRREP)
         END IF
         IF(MBPT4) THEN
          IF (bRedundant) THEN
             CALL GETLST(ICORE(I010),1,NMSYT2B,1,IRREP,LIST5B)
          ELSE
             CALL GETLST_NR(ICORE(I010),ICORE(I030),MXCOR-I030+1,
     &                      LIST5B,IRREP)
          END IF
          CALL SSCAL(NMSYT2B*DSSYT2B,TWO,ICORE(I020),1)
          CALL SAXPY(NMSYT2B*DSSYT2B,ONE,ICORE(I010),1,ICORE(I020),1)
         ENDIF
         IF (bRedundant) THEN
            CALL GETLST(ICORE(I010),1,NMSYT1B,1,IRREP,LIST1B)
         ELSE
            CALL GETLST_NR(ICORE(I010),ICORE(I030),MXCOR-I030+1,
     &                     LIST1B,IRREP)
         END IF
         CALL XGEMM('N','N',DISSYG,NUMSYG,DSSYT2B,FACT,ICORE(I010),
     &              DSSYT1B,ICORE(I020),DSSYT2B,ONE,ICORE(I000),
     &              DISSYG)
        ELSE
         I020=I010+IINTFP*NUMSYG*DISSYG
         IF(MXCOR.LT.I020) CALL INSMEM('G4ALL',I020,MXCOR)
         CALL TRANSP(ICORE(I000),ICORE(I010),NUMSYG,DISSYG)
         CALL SAXPY(NUMSYG*DISSYG,ONE,ICORE(I010),1,ICORE(I000),1)
        ENDIF
C
C ADD SINGLES IN ROHF-MBPT(3)
C
        IF(ROHF.AND.MBPT3.AND.IRREP.EQ.1) THEN
         CALL T1ING4(ICORE(I000),DISSYG,NUMSYG,ICORE(I010),
     &               ONE,1,1+IUHF,90,90)
        ENDIF
C
C  FOR RHF :
C  SYMMETRIZE GAMMA 4 BUT NOT FOR MBPT3
C  NOTE THAT THIS REQUIRES IN UHF RUNS LIST 117
C  FOR RHF 117 AND 118 ARE THE TRANSPOSE OF EACH OTHER
C
        IF((EOM .OR. LTRP).AND. (IUHF.EQ.0))THEN
         CALL GETLST(ICORE(I010),1,NUMSYG,1,IRREP,LISTG)
         CALL SAXPY (NUMSYG*DISSYG,ONE,ICORE(I010),1,ICORE(I000),1)
        ENDIF
        IF(IUHF.EQ.0.AND..NOT.MBPT3) CALL SYMMET2(ICORE(I000),DISSYG)
        IF(IUHF.EQ.0.OR.MBPT3)  THEN
         CALL PUTLST(ICORE(I000),1,NUMSYG,2,IRREP,LISTG)
        ENDIF
        IF(IUHF.EQ.1) THEN
         IF(MBPT3) THEN
          CALL TRANSP(ICORE(I000),ICORE(I010),NUMSYG,DISSYG)
          CALL PUTLST(ICORE(I010),1,DISSYG,2,IRREP,LISTG-1)
         ELSE IF(.NOT.MBPT3) THEN
          CALL TRANSP(ICORE(I000),ICORE(I010),NUMSYG,DISSYG)
c YAU : old
c         CALL ICOPY(NUMSYG*DISSYG*IINTFP,ICORE(I010),1,ICORE(I000),1)
c YAU : new
          CALL DCOPY(NUMSYG*DISSYG,ICORE(I010),1,ICORE(I000),1)
c YAU : end
          DSSYT3A=IRPDPD(IRREP,ISYTYP(1,LIST3A))
          DSSYT4A=IRPDPD(IRREP,ISYTYP(1,LIST4A))
          NMSYT3A=IRPDPD(IRREP,ISYTYP(2,LIST3A))
          NMSYT4A=IRPDPD(IRREP,ISYTYP(2,LIST4A))
          DSSYT3B=IRPDPD(IRREP,ISYTYP(1,LIST3B))
          DSSYT4B=IRPDPD(IRREP,ISYTYP(1,LIST4B))
          NMSYT3B=IRPDPD(IRREP,ISYTYP(2,LIST3B))
          NMSYT4B=IRPDPD(IRREP,ISYTYP(2,LIST4B))
          DISSYG2=IRPDPD(IRREP,ISYTYP(1,LISTG2))
          NUMSYG2=IRPDPD(IRREP,ISYTYP(2,LISTG2))
C
C  FIRST TERM
C
          I020=I010+IINTFP*MAX(NMSYT3A*DSSYT3A,NMSYT4A*DSSYT4A)
          I030=I020+IINTFP*NMSYT4A*DSSYT4A
          IF (bRedundant) THEN
             CALL GETLST(ICORE(I020),1,NMSYT4A,1,IRREP,LIST4A)
          ELSE
             CALL GETLST_NR(ICORE(I020),ICORE(I030),MXCOR-I030+1,
     &                      LIST4A,IRREP)
          END IF
C
C FOR MBPT(4) FORM APPROBIATE LINEAR COMBINATION
C
          IF(MBPT4) THEN
           IF (bRedundant) THEN
              CALL GETLST(ICORE(I010),1,NMSYT4A,1,IRREP,LIST6A)
           ELSE
              CALL GETLST_NR(ICORE(I010),ICORE(I030),MXCOR-I030+1,
     &                       LIST6A,IRREP)
           END IF
           CALL SSCAL(NMSYT4A*DSSYT4A,TWO,ICORE(I020),1)
           CALL SAXPY(NMSYT4A*DSSYT4A,ONE,ICORE(I010),1,ICORE(I020),1)
          ENDIF
C
          IF (bRedundant) THEN
             CALL GETLST(ICORE(I010),1,NMSYT3A,1,IRREP,LIST3A)
          ELSE
             CALL GETLST_NR(ICORE(I010),ICORE(I030),MXCOR-I030+1,
     &                      LIST3A,IRREP)
          END IF
C
          CALL XGEMM('N','N',DISSYG2,NUMSYG2,DSSYT4A,FACT,ICORE(I010),
     &               DSSYT3A,ICORE(I020),DSSYT4A,ONE,ICORE(I000),
     &               DISSYG2)
          I020=I010+IINTFP*MAX(NMSYT3B*DSSYT3B,NMSYT4B*DSSYT4B)
          I030=I020+IINTFP*NMSYT4B*DSSYT4B
          IF(MXCOR.LT.I030) CALL INSMEM('G4ALL',I030,MXCOR)
          IF (bRedundant) THEN
             CALL GETLST(ICORE(I020),1,NMSYT4B,1,IRREP,LIST4B)
          ELSE
             CALL GETLST_NR(ICORE(I020),ICORE(I030),MXCOR-I030+1,
     &                      LIST4B,IRREP)
          END IF
C
C  FOR MBPT(4) FORM APPROBIATE LINEAR COMBINATION
C
          IF(MBPT4) THEN
           IF (bRedundant) THEN
              CALL GETLST(ICORE(I010),1,NMSYT4B,1,IRREP,LIST6B)
           ELSE
              CALL GETLST_NR(ICORE(I010),ICORE(I030),MXCOR-I030+1,
     &                       LIST6B,IRREP)
           END IF
           CALL SSCAL(NMSYT4B*DSSYT4B,TWO,ICORE(I020),1)
           CALL SAXPY(NMSYT4B*DSSYT4B,ONE,ICORE(I010),1,ICORE(I020),1)
          ENDIF
C
          IF (bRedundant) THEN
             CALL GETLST(ICORE(I010),1,NMSYT3B,1,IRREP,LIST3B)
          ELSE
             CALL GETLST_NR(ICORE(I010),ICORE(I030),MXCOR-I030+1,
     &                      LIST3B,IRREP)
          END IF
          CALL XGEMM('N','N',DISSYG2,NUMSYG2,DSSYT4B,FACT,ICORE(I010),
     &               DSSYT3B,ICORE(I020),DSSYT4B,ONE,ICORE(I000),
     &               DISSYG2)
C
C  SCALE BY 1/2
C
          CALL SSCAL(NUMSYG*DISSYG,HALF,ICORE(I000),1)
C
C  SAVE 117 ON DISK
C
          IF(EOM .OR. LTRP)THEN
           CALL GETLST(ICORE(I010),1,NUMSYG2,1,IRREP,LISTG2)
           CALL SAXPY (NUMSYG*DISSYG,ONE,ICORE(I010),1,ICORE(I000),1)
          ENDIF
          CALL PUTLST(ICORE(I000),1,NUMSYG2,2,IRREP,LISTG2)
C
C  FORM 118 BY TRANSPOSITION AND SAVE ON DISK
C
          CALL TRANSP(ICORE(I000),ICORE(I010),NUMSYG2,DISSYG2)
          CALL PUTLST(ICORE(I010),1,NUMSYG,2,IRREP,LISTG)
         ENDIF
        ENDIF
200    CONTINUE   
       ENDIF
C
      ELSEIF(SPCASE.EQ.'ABAB'.OR.SPCASE.EQ.'BABA')THEN
       IF(CIS)RETURN
C
C DO SPIN CASES ABAB AND BABA:
C
C       G(IaJb) =  -  SUM  M,E T1(Im,Eb) T1(Jm,Ea)
C
C  STORAGE  Bi Aj
       IF(SPCASE.EQ.'ABAB')THEN
        LISTG=125
        IF(MBPT3) THEN
         ISAME=.TRUE.
         FACT=ONEM
         LIST1=39
         LIST2=39
        ELSE IF(MBPT4) THEN
         LIST1=39
         LIST2=139
         FACT=ONEM
         ISAME=.FALSE.
        ELSE IF(CC) THEN
         LISTH=58+IOFFLISTH
        ENDIF
       ELSEIF(SPCASE.EQ.'BABA')THEN
        LISTG=126
        IF(MBPT3) THEN
         ISAME=.TRUE.
         FACT=ONEM
         LIST1=38
         LIST2=38
        ELSE IF(MBPT4) THEN
         FACT=ONEM
         ISAME=.FALSE. 
         LIST1=38
         LIST2=138
        ELSE IF(CC) THEN
         LISTH=59+IOFFLISTH
        ENDIF
       ENDIF
       IF(CC) THEN
       DO 290 IRREP=1,NIRREP
        DISSYG=IRPDPD(IRREP,ISYTYP(1,LISTG))
        NUMSYG=IRPDPD(IRREP,ISYTYP(2,LISTG))
        I000=1         
        I010=I000+IINTFP*NUMSYG*DISSYG
        CALL GETLST(ICORE(I000),1,NUMSYG,1,IRREP,LISTH)
        IF(EOM .OR. LTRP)THEN
         CALL GETLST(ICORE(I010),1,NUMSYG,1,IRREP,LISTG)
         CALL SAXPY (NUMSYG*DISSYG,ONE,ICORE(I010),1,ICORE(I000),1)
        ENDIF
C
C SYMMETRIZE THE GAMMA 4 INTERMEDIATE
C
        CALL SYMMET2(ICORE(I000),DISSYG)
C
C SAVE GAMMA 4 INTERMEDIATE ON DISK
C
        CALL PUTLST(ICORE(I000),1,NUMSYG,2,IRREP,LISTG)
CSSS        call checksum("H4ABAB  :",ICORE(I000),DISSYG*NUMSYG)
C
290    CONTINUE
C
       ELSE
C 
       DO 300 IRREP=1,NIRREP 
        DSSYT1=IRPDPD(IRREP,ISYTYP(1,LIST1))
        DSSYT2=IRPDPD(IRREP,ISYTYP(1,LIST2))
        NMSYT1=IRPDPD(IRREP,ISYTYP(2,LIST1))
        NMSYT2=IRPDPD(IRREP,ISYTYP(2,LIST2))
        DISSYG=IRPDPD(IRREP,ISYTYP(1,LISTG))
        NUMSYG=IRPDPD(IRREP,ISYTYP(2,LISTG))
        I000=1
        I010=I000+IINTFP*DISSYG*NUMSYG
        IF(ISAME) THEN
         I020=I010
        ELSE
         I020=I010+IINTFP*NMSYT1*DSSYT1
        ENDIF
        I030=I020+IINTFP*NMSYT2*DSSYT2
        IF(MXCOR.LT.I030) CALL INSMEM('G4ALL',I030,MXCOR)
        CALL IZERO(ICORE(I000),DISSYG*NUMSYG*IINTFP)
        IF (bRedundant) THEN
           CALL GETLST(ICORE(I010),1,NMSYT1,1,IRREP,LIST1)
           IF (.NOT.ISAME) CALL GETLST(ICORE(I020),1,NMSYT2,1,
     &                                 IRREP,LIST2)
        ELSE
           CALL GETLST_NR(ICORE(I010),ICORE(I030),MXCOR-I030+1,
     &                    LIST1,IRREP)
           IF (.NOT.ISAME) CALL GETLST_NR(ICORE(I020),ICORE(I030),
     &                                    MXCOR-I030+1,LIST2,IRREP)
        END IF
C
C FOR MBPT(4) FORM APPROBIATE LINEAR COMBINATION
C
        IF(MBPT4) THEN
         CALL SSCAL(NMSYT2*DSSYT2,TWO,ICORE(I020),1)
         CALL SAXPY(NMSYT2*DSSYT2,ONE,ICORE(I010),1,ICORE(I020),1)
        ENDIF
C
        CALL XGEMM('N','T',DISSYG,NUMSYG,NMSYT1,FACT,ICORE(I010),
     &             DSSYT1,ICORE(I020),DSSYT2,AZERO,ICORE(I000),DISSYG)
        IF(EOM .OR. LTRP)THEN
         CALL GETLST(ICORE(I010),1,NUMSYG,1,IRREP,LISTG)
         CALL SAXPY (NUMSYG*DISSYG,ONE,ICORE(I010),1,ICORE(I000),1)
        ENDIF
C
C SYMMETRIZE THE GAMMA 4 INTERMEDIATE BUT ONLY IF ISAME = FALSE
C
        IF(.NOT.ISAME) CALL SYMMET2(ICORE(I000),DISSYG)
C
C SAVE GAMMA 4 INTERMEDIATE ON DISK
C
        CALL PUTLST(ICORE(I000),1,NUMSYG,2,IRREP,LISTG)
CSSS        call checksum("H4ABAB  :",ICORE(I000),DISSYG*NUMSYG)
300    CONTINUE
       ENDIF
      ENDIF
      RETURN
      END
