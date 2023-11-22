      SUBROUTINE S2PHPHS2(ICORE,MAXCOR,IUHF, ISIDE, ISPIN)
C
C   THIS SUBROUTINE CALCULATES THE CONTRIBUTION OF S2 
C   TO S2 INVOLVING THE PHPH BLOCK OF EFFECTIVE INTEGRALS.
C
C CONTRACTIONS :
C
C    Z(IJ,AP) = -Sum Hbar(MA,JE) * S2(MI,EP) (ABAB and BABA terms) 
C                E,M
C    Z(IJ,AP) = -Sum Hbar(MA,IE) * S2(JM,EP)
C                E,M
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE,ONEM,ZILCH
C
      DIMENSION ICORE(MAXCOR)
      DIMENSION LISTW0(2,2), NUMSZS(8)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SINFO/ NS(8), SIRREP
      COMMON/SLISTS/LS1IN, LS1OUT, LS2IN(2,2), LS2OUT(2,2)
C
      DATA ONE  /1.0/
      DATA ZILCH/0.0/
      DATA ONEM/-1.0/
C
      Write(6,"(a)") "-----------------S2PHPHS2--------------------"
      LISTW0(1,1) = 54
      LISTW0(2,1) = 58
      LISTW0(1,2) = 59
      LISTW0(2,2) = 55
C
C  THE FIRST INDEX IN LISTW0 LABELS THE SPIN OF K AND I, THE SECOND INDEX
C  LABELS THE SPIN OF B AND C AND EQUALS ISPIN. THE FIRST INDEX DETERMINES 
C  ALSO IMIXSPIN IN S.      
C
C  Z(JI,AP) = Hbar(MA,JE)*S2(MI,EP)
C  
      DO 5 XIRREP = 1, NIRREP
         IRREPE = DIRPRD(XIRREP,SIRREP)
         NUMSZS(XIRREP) = VRT(IRREPE,ISPIN) * NS(SIRREP)
 5    CONTINUE      
C
      IF (IUHF.EQ.0) THEN
         IMIXSPIN = 1
      ELSE
         IMIXSPIN = 3 - ISPIN
      ENDIF
      ASPIN = IMIXSPIN
      MSPIN = ISPIN
      LISTS2IN = LS2IN(ISPIN,IMIXSPIN + 1 - IUHF)
      LISTW = LISTW0(IMIXSPIN+1-IUHF,ISPIN)

C UHF (LISTW=58) Second term: Hbar(Ma,Ie)*S2(jM,eP)
C RHF (LISTW,58 First and second terms)

       CALL GETLEN_HHP(LENS,POP(1,IMIXSPIN),POP(1,ISPIN), 
     $                 VRT(1,IMIXSPIN),NS)

       NMAX = MAX(POP(1,ISPIN)*POP(1,IMIXSPIN),
     $            VRT(1,IMIXSPIN)*NS(SIRREP))
       I000=1
       I010=I000+LENS*IINTFP
       I020=I010+LENS*IINTFP
       IEND=I020+NMAX

       IF (IEND .GT. MAXCOR) CALL INSMEM("S2PHPHS2",IEND,MAXCOR)

       IF ((IUHF.NE.0) .AND. ISPIN.EQ.IMIXSPIN) THEN       
          CALL GETEXP2_HHP(ICORE(I010),LENS,NUMSZS, 
     $                     IRPDPD(1,20+ISPIN),LISTS2IN, 
     $                     POP(1,ISPIN), 
     $                     IRPDPD(1,ISYTYP(1,LISTS2IN)))
       ELSE
          CALL GETALLS2_HHP(ICORE(I010),LENS,VRT(1,IMIXSPIN),
     $                      NS(1),1,LISTS2IN)
       ENDIF
C
C REORDER THE S COEFFICIENT ACCORDING TO S2IN(MI,EP) -> S(IM,EP) -> S(EM,IP)
C UHF = S2IN(Ma,Ie)*S2(jM,ep) No symtra (only sstgen). 
     
       IF (IUHF .EQ. 0) THEN
          SOFF = I010
          DO XIRREP = 1, NIRREP
             IRREPE = DIRPRD(XIRREP,SIRREP)
             DISSYS = IRPDPD(XIRREP,14)
             NUMDSS = VRT(IRREPE,IMIXSPIN)*NS(SIRREP)

             MAXS  = MAX(DISSYS,NUMDSS)
             ITMP1 = I020
             ITMP2 = ITMP1 + MAXS*IINTFP
             ITMP3 = ITMP2 + MAXS*IINTFP
             IEND  = ITMP3 + MAXS*IINTFP
             IF (IEND .GT. MAXCOR) CALL INSMEM("S2PHPHS2",IEND,MAXCOR)
             CALL SYMTR3(XIRREP,POP(1,ISPIN),POP(1,IMIXSPIN),DISSYS,
     &                   NUMDSS,ICORE(SOFF),ICORE(ITMP1),ICORE(ITMP2),
     &                   ICORE(ITMP3))
             SOFF = SOFF + DISSYS*NUMDSS
          ENDDO
       ENDIF 
C
C REORDER THE S COEFFICIENT ACCORDING TO S2(IM,EP) -> S2(EM,IP)
C
       CALL SSTGEN(ICORE(I010),ICORE(I000),LENS,POP(1,IMIXSPIN), 
     $             POP(1,IMIXSPIN),VRT(1,IMIXSPIN),NS,ICORE(I020),1,
     $             '3214')

#ifdef _DEBUG_LVL0
      call checksum("->s2    :",icore(i000),lens)
#endif
       CALL ZERO(ICORE(I010),LENS)
C
C LOOP OVER IRREPS OF RHS OF Hbar(MA,JE) stored as (Em,Aj) (RHF) 
C AND (Em,Aj),(eM,aJ) for UHF. 
C
       ICOUNT1 = I000
       ICOUNT2 = I010
       DO 100 XIRREP=1,NIRREP
          IRREPE  = DIRPRD(XIRREP,SIRREP)
          IRREPJ  = IRREPE 
          IRREPW  = XIRREP

          DISSYW=IRPDPD(IRREPW,ISYTYP(1,LISTW)) 
          NUMDSW=IRPDPD(IRREPW,ISYTYP(2,LISTW)) 
          DISSYS=DISSYW
          NUMDSS=POP(IRREPJ,IMIXSPIN) * NS(SIRREP)

          I030=I020+IINTFP*DISSYW*NUMDSW

          IF (I030.GT.MAXCOR) CALL INSMEM("S2PHPHS2",I030,MAXCOR)
C
          CALL GETLST(ICORE(I020),1,NUMDSW,1,IRREPW,LISTW)
#ifdef _DEBUG_LVL0
       write(6,"(a,i4)") "Listw:", Listw
       call checksum("->w     :",icore(i020),dissyw*numdsw)
#endif
C
C PERFORM MATRIX MULTIPLICATION

C  Z(Aj,IP) =  -SUM W(Em,Aj) * S2(Em,IP)  [ISIDE = 1]
C               ME
C  Z(Aj,IP) =  -SUM W(Aj,Em) * S2(Em,IP)  [ISIDE = 2]
C               ME

          NROW=NUMDSW
          NCOL=NUMDSS
          NSUM=DISSYW

#ifdef _DEBUG_LVLM
             write(6,*) ' W-integrals'
             call output(icore(i020), 1, dissyw, 1, numdsw, dissyw,
     $          numdsw, 1)
             write(6,*) ' s-coefficients'
             call output(icore(icount1), 1, dissys, 1, numdss, dissys,
     $          numdss, 1)
#endif 
          IF (ISIDE.EQ.1) THEN
            CALL XGEMM('T','N',NROW,NCOL,NSUM,ONEM,ICORE(I020),NSUM,
     &               ICORE(ICOUNT1),NSUM,ONE,ICORE(ICOUNT2),NROW)
          ELSE
            CALL XGEMM('N','N',NROW,NCOL,NSUM,ONEM,ICORE(I020),NROW,
     &               ICORE(ICOUNT1),NSUM,ONE,ICORE(ICOUNT2),NROW)
          ENDIF

#ifdef _DEBUG_LVLM
             write(6,*) ' output s-coefficients'
             call output(icore(icount2), 1, dissys, 1, numdss,
     $          dissys, numdss, 1)
#endif 
          ICOUNT1 = ICOUNT1 + DISSYS * NUMDSS * IINTFP
          ICOUNT2 = ICOUNT2 + DISSYS * NUMDSS * IINTFP
100    CONTINUE

#ifdef _DEBUG_LVL0
       call checksum("s2->s2  :",icore(i010),lens)
#endif
C
C  REORDER S COEFFICIENTS BACK TO ORIGINAL FORM S2EX[AJ,IP] -> S2EX[IJ,AP]
C
       CALL SSTGEN(ICORE(I010),ICORE(I000),LENS,VRT(1,IMIXSPIN),
     $             POP(1,ISPIN),POP(1,IMIXSPIN),NS(1),ICORE(I020),
     $             1,'3214')

C  S2EX(IJ,AP) -> S2EX(JI,AP), For UHF (S2EX(JI,AP) and no symtr3 is needed)

       IF (IUHF .EQ. 0) THEN
          SOFF  = I000
          ITMP1 = I020
          DO XIRREP = 1, NIRREP
             IRREPA = DIRPRD(XIRREP,SIRREP)
             DISSYS = IRPDPD(XIRREP,14)
             NUMDSS = VRT(IRREPA,IMIXSPIN)*NS(SIRREP)

             MAXS  = MAX(DISSYS,NUMDSS)
             ITMP2 = ITMP1 + MAXS*IINTFP
             ITMP3 = ITMP2 + MAXS*IINTFP
             IEND  = ITMP3 + MAXS*IINTFP
             IF (IEND .GT. MAXCOR) CALL INSMEM("S2PHPHS2",IEND,MAXCOR)
              CALL SYMTR3(XIRREP,POP(1,IMIXSPIN),POP(1,ISPIN),DISSYS,
     &                    NUMDSS,ICORE(SOFF),ICORE(ITMP1),ICORE(ITMP2),
     &                    ICORE(ITMP3))
               SOFF = SOFF + DISSYS*NUMDSS
          ENDDO
       ENDIF 
C
C ADD THE S INCREMENTS TO LISTS2EX
C
       LISTS2EX = LS2OUT(ISPIN, IMIXSPIN + 1 - IUHF)

       IF ((IUHF.NE.0) .AND. (IMIXSPIN.EQ.ISPIN)) THEN
          CALL ASSYMALL_HHP(ICORE(I000),LENS,NUMSZS, 
     $                      IRPDPD(1,20+ISPIN),
     $                      POP(1,ISPIN),ICORE(I020), 
     $                      MAXCOR-I020+1)          
          CALL GETEXP2_HHP(ICORE(I010),LENS,NUMSZS,
     $                     IRPDPD(1,20+ISPIN),LISTS2EX,
     $                     POP(1,ISPIN),IRPDPD(1,ISYTYP(1,LISTS2EX)))
       ELSE
          CALL GETALLS2_HHP(ICORE(I010),LENS,VRT(1,IMIXSPIN),NS(1),
     $                      1,LISTS2EX)
       ENDIF       
#ifdef _DEBUG_LVL0
      call checksum("->1hs2  :",icore(i010),lens)
      call checksum("s2->s2f :",icore(i000),lens)
#endif
C
       CALL SAXPY(LENS,ONE,ICORE(I000),1,ICORE(I010),1)
C
#ifdef _DEBUG_LVL0
      call checksum("1s2->fs2:",icore(i010),lens)
#endif
       IF ((IUHF.NE.0).AND. (IMIXSPIN.EQ.ISPIN)) THEN
          CALL PUTSQZ_HHP(ICORE(I010),LENS,NUMSZS, 
     $                    IRPDPD(1,20+ISPIN),
     $                    LISTS2EX,POP(1,ISPIN),
     $                    IRPDPD(1,ISYTYP(1,LISTS2EX)),
     $                    ICORE(I020),MAXCOR-I020+1)
       ELSE
          CALL PUTALLS2_HHP(ICORE(I010),LENS,VRT(1,IMIXSPIN),
     $                      NS(1),1,LISTS2EX)
       ENDIF
C
      IF (IUHF.EQ.0) THEN
C
C DO THE FIRST TERM FOR THE RHF CASE
C
C   Z(IJ,AP) =  -Hbar(MA,IE) * S2(JM,EP)
C
       CALL GETLEN_HHP(LENS,POP(1,1),POP(1,1),VRT(1,1),NS(1))
       I000 = 1
       I010 = I000 + LENS * IINTFP
       I020 = I010 + LENS * IINTFP
       IEND = I020
       IF (IEND .GT. MAXCOR) CALL INSMEM("S2PHPHS2",IEND,MAXCOR)
       LISTS2IN = LS2IN(1,2)
       LISTS2EX = LS2OUT(1,2)
       LISTW    = LISTW0(2,1)
       CALL GETALLS2_HHP(ICORE(I010),LENS,VRT(1,1),NS(1),1,
     &                   LISTS2IN)

C S2IN(JM,EP) -> S2IN(EM,JP)

       CALL SSTGEN(ICORE(I010),ICORE(I000),LENS,POP(1,1),POP(1,1),
     $             VRT(1,1),NS,ICORE(I020),1,'3214')

#ifdef _DEBUG_LVL0
      call checksum("->2s2   :",icore(i000),lens)
#endif
       CALL ZERO(ICORE(I010),LENS)
C
C LOOP OVER IRREPS OF RHS OF S VECTOR
C
       ICOUNT1 = I000
       ICOUNT2 = I010
       DO 200 XIRREP=1,NIRREP
        IRREPJ=DIRPRD(XIRREP,SIRREP)
        IRREPW=XIRREP
        DISSYW=IRPDPD(IRREPW,ISYTYP(1,LISTW)) 
        NUMDSW=IRPDPD(IRREPW,ISYTYP(2,LISTW)) 
        NUMDSS=POP(IRREPJ,ISPIN)*NS(SIRREP)
        I030=I020+IINTFP*DISSYW*NUMDSW
        IF(I030.GT.MAXCOR)THEN
C
C OUT-OF-CORE ALGORITHM
C
         WRITE(6,*)' out-of-core AB not coded '
         call errex
C
        ENDIF
C
C READ W INTO W(Em,Aj)
C     
         CALL GETLST(ICORE(I020),1,NUMDSW,1,IRREPW,LISTW)
#ifdef _DEBUG_LVLM
         call checksum("->2s2   :",icore(i000),lens)
         call checksum("->w2    :",icore(i020),dissyw*numdsw)
#endif
C
C PERFORM MATRIX MULTIPLICATION
C
C Z(Ai,JP) = - SUM W(Em,Ai) * S(Em,Jp)   [ISIDE = 1]
C              Em 
C Z(Ai,JP) = - SUM W(Ai,Em) * S(Em,JP)   [ISIDE = 2]
C              Em
          NROW=NUMDSW
          NCOL=NUMDSS
          NSUM=DISSYW
          IF (ISIDE.EQ.1) THEN
            CALL XGEMM('T','N',NROW,NCOL,NSUM,ONEM,ICORE(I020),NSUM,
     &               ICORE(ICOUNT1),NSUM,ONE,ICORE(ICOUNT2),NROW)
          ELSE
            CALL XGEMM('N','N',NROW,NCOL,NSUM,ONEM,ICORE(I020),NROW,
     &               ICORE(ICOUNT1),NSUM,ONE,ICORE(ICOUNT2),NROW)
          ENDIF
          ICOUNT1 = ICOUNT1 + DISSYW * NUMDSS * IINTFP
          ICOUNT2 = ICOUNT2 + NUMDSW * NUMDSS * IINTFP
 200   CONTINUE
 
#ifdef _DEBUG_LVL0
       call checksum("2s2->s2 :",icore(i010),lens)
#endif 
C
C  REORDER S COEFFICIENTS BACK TO ORIGINAL FORM S[Ai,JP] -> S[Ji,AP]
C
       CALL SSTGEN(ICORE(I010),ICORE(I000),LENS,VRT(1,1),POP(1,1),
     &                   POP(1,1),NS(1),ICORE(I020),1,'3214')

C ADD THE S INCREMENTS TO LISTS2EX
C
       CALL GETALLS2_HHP(ICORE(I010),LENS,VRT(1,1),NS(1),1,LISTS2EX)

#ifdef _DEBUG_LVL0
      call checksum("->2hs2  :",icore(i010),lens)
      call checksum("2s2->s2 :",icore(i000),lens)
#endif
       CALL DAXPY(LENS,ONE,ICORE(I000),1,ICORE(I010),1)
#ifdef _DEBUG_LVL0
      call checksum("2s2->fs2:",icore(i010),lens)
#endif
       CALL PUTALLS2_HHP(ICORE(I010),LENS,VRT(1,1),NS(1),1,LISTS2EX)
C
       ENDIF
C#endif 
C       
       RETURN
       END
