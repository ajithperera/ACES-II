      SUBROUTINE S2HPPHS2(ICORE,MAXCOR,IUHF, ISIDE, ISPIN)   
C
C THIS SUBROUTINE CALCULATES THE CONTRIBUTION OF S2 TO
C S2 INVOLVING THE HPPH BLOCK OF EFFECTIVE INTEGRALS W (RINGS).
C
C CONTRACTIONS :
C
C     Z(JI,AP) = SUM <MA|EI> (2 S(MJ,EP) - S(JM,EP)) [ RHF CASE, ISIDE=1]
C                M,E
C     Z(JI,AP) = SUM <IE|AM> (2 S(JM,EP) - S(MJ,EP)) [ RHF CASE, ISIDE=2]
C                M,E
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE,ONEM,ZILCH
      logical prints, print, skip
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
      DATA ONE  /1.0D0/
      DATA ZILCH/0.0D0/
      DATA ONEM/-1.0D0/
C
      Write(6,"(a)") "-----------------S2HPPHS2--------------------"
      LISTW0(1,1) = 54
      LISTW0(1,2) = 56
      LISTW0(2,1) = 57
      LISTW0(2,2) = 55
C
C LISTS 54 AND 55 CONTAIN HBAR(MB,IE) AND HBAR(mb,ie) STORED AS (EM,BI) or 
C (em,bi). 
C
      DO 5 XIRREP = 1, NIRREP
         IRREPE = DIRPRD(XIRREP,SIRREP)
         NUMSZS(XIRREP) = VRT(IRREPE,ISPIN) * NS(SIRREP)
 5    CONTINUE      
C
      DO 25 MIXSPIN = 1, 1+IUHF
         LISTS2IN = LS2IN(ISPIN,MIXSPIN + 1 - IUHF)
         CALL GETLEN_HHP(LENS,POP(1,MIXSPIN),POP(1,ISPIN),
     $                   VRT(1,MIXSPIN),NS)
       I000=1
       I010=I000+LENS*IINTFP
       IEND=I010
       IF (IEND .GT. MAXCOR) CALL INSMEM("S2PHPHS2",IEND,MAXCOR)
       IF ((IUHF.NE.0) .AND. (MIXSPIN.EQ.ISPIN)) THEN
          CALL GETEXP2_HHP(ICORE(I000),LENS,NUMSZS, 
     $                     IRPDPD(1,20+ISPIN),
     $                     LISTS2IN,POP(1,MIXSPIN),
     $                     IRPDPD(1,ISYTYP(1,LISTS2IN)))
       ELSE
          CALL GETALLS2_HHP(ICORE(I000),LENS,VRT(1,MIXSPIN),NS(1),1,
     $                      LISTS2IN)
       ENDIF

C S2IN(MJ,EP) -> S2IN(JM,EP)
      SOFF = I000
      DO XIRREP = 1,NIRREP
         IRREPE = DIRPRD(XIRREP,SIRREP)
         IF (ISPIN .EQ. MIXSPIN) THEN
            DISSYS = IRPDPD(XIRREP,20+ISPIN)
         ELSE
            DISSYS = IRPDPD(XIRREP,ISYTYP(1,LISTS2IN))
         ENDIF
         NUMDSS = VRT(IRREPE,MIXSPIN) * NS(SIRREP)

         MAXS  = MAX(DISSYS,NUMDSS)
         ITMP1 = I010
         ITMP2 = ITMP1 + MAXS*IINTFP
         ITMP3 = ITMP2 + MAXS*IINTFP
         IEND  = ITMP3 + MAXS*IINTFP
         IF (IEND .GT. MAXCOR) CALL INSMEM("S2PHPHS2",IEND,MAXCOR)
         CALL SYMTR3(XIRREP,POP(1,MIXSPIN),POP(1,ISPIN),DISSYS,NUMDSS,
     &                ICORE(SOFF),ICORE(ITMP1),ICORE(ITMP2),
     &                ICORE(ITMP3))
         SOFF = SOFF + DISSYS*NUMDSS
       ENDDO
#ifdef _DEBUG_LVL0          
       write(6,"(a,2i4)") 'soff-i000,lens:',soff-i000,lens
#endif 
C
C SPINADAPT THE S COEFFICIENTS BY IRREP S2(JM,Ep) (i.e. S2(Jm,Ep))
C
       IF (IUHF.EQ.0) THEN
          ICOUNT = I000
          DO 10 XIRREP = 1,NIRREP
             IRREPE = DIRPRD(XIRREP,SIRREP)
             DISSYS=IRPDPD(XIRREP,ISYTYP(1,LISTS2IN))
             NUMDSS= VRT(IRREPE,1) * NS(SIRREP)
             I020 = I010 + NUMDSS * IINTFP
             I030 = I020 + NUMDSS * IINTFP
             CALL SPINAD3(XIRREP,POP(1,1),DISSYS,NUMDSS,
     $                    ICORE(ICOUNT),ICORE(I010),ICORE(I020))
             ICOUNT = ICOUNT + NUMDSS * DISSYS * IINTFP
 10       CONTINUE
       ENDIF

C S2IN(JM,EP) -> S2IN(EM,JP)
       I020 = I010 + LENS*IINTFP
       IEND = I020  
       IF (IEND .GT. MAXCOR) CALL INSMEM("S2PHPHS2",IEND,MAXCOR)
       CALL SCOPY(LENS,ICORE(I000),1,ICORE(I010),1)
       CALL SSTGEN(ICORE(I010),ICORE(I000),LENS,POP(1,MIXSPIN),
     $             POP(1,ISPIN),VRT(1,MIXSPIN),NS,ICORE(I020),1, 
     $            '3214')
C
C ICORE(I000) CONTAINS ALL [SPIN ADDAPTED] S COEFFICIENTS [2S(MJ,EP) -
C S(JM,EP)], ORDERED AS EM;JP
C
C NOW CARRY OUT THE SUMMATION WHICH USES THE [ SPIN ADAPTED ]
C TRANSPOSED S-ELEMENTS.
C
       DO 50 SOUTSPIN = 1, 1 + IUHF
          CALL GETLEN_HHP(LENS,POP(1,SOUTSPIN),POP(1,ISPIN),
     $                    VRT(1,SOUTSPIN),NS(1))

       I020 = I010 + LENS *IINTFP
       CALL ZERO(ICORE(I010),LENS)

       IF (ISIDE.EQ.1) THEN
          LISTW = LISTW0(MIXSPIN, SOUTSPIN + 1 - IUHF)
C LIST-56, HBAR(MA|EI) -> HBAR(EM,ai) (RHF)
C LIST-54(EM,AI), LIST-55(em,ai),LIST-57 (em,AI)(UHF)
       ELSE
          LISTW = LISTW0(SOUTSPIN, MIXSPIN + 1 - IUHF)
C LIST-57, HBAR(IE|AM) -> HBAR(ai,EM) (RHF)
C LIST-54(AI,EM), LIST-55(ai,em),LIST-57 (ai,EM)(UHF)
       ENDIF
C
C LOOP OVER RHS IRREPS OF HBAR(MA,EI) stored as (EM,AI) and at this
C point S2 is stored as S2(EM,JP)
C
       ICOUNT1 = I000
       ICOUNT2 = I010
#ifdef _DEBUG_LVL0
      Write(6,"(a,3I4)") "ispin.mixspin,soutspin:",ispin,mixspin,
     &                     soutspin
      write(6,"(a,2i4)") 'listw,lens:',listw,lens
      call checksum("->s2    :",icore(i000),lens)
c      ioff = i000
c      do xirrep = 1, nirrep
c         irrepe = dirprd(xirrep,sirrep)
c         if (ispin .eq. mixspin) then
c            dissys = irpdpd(xirrep,20+ispin)
c         else
c            dissys = irpdpd(xirrep,isytyp(1,lists2in))
c         endif 
c         numdss = vrt(irrepe,soutspin)*Ns(sirrep)
c         call output(icore(ioff),1,dissys,1,numdss,dissys,
c     &               numdss,1)
c         ioff = ioff + dissys*numdss 
c      enddo
#endif

       DO 100 XIRREP=1,NIRREP 
          IRREPI  = DIRPRD(XIRREP,SIRREP)
          IRREPW  = XIRREP

          DISSYW=IRPDPD(IRREPW,ISYTYP(1,LISTW)) 
          NUMDSW=IRPDPD(IRREPW,ISYTYP(2,LISTW)) 
          IF (ISPIN .EQ. SOUTSPIN) THEN
             DISSYS=IRPDPD(IRREPW,8+ISPIN)
          ELSE
             DISSYS=IRPDPD(IRREPW,8+SOUTSPIN)
          ENDIF 
          NUMDSS=POP(IRREPI,ISPIN)*NS(SIRREP)
          I030=I020+IINTFP*DISSYW*NUMDSW

          IF(I030.GT.MAXCOR)THEN
             WRITE(6,*)' out-of-core AB not coded, s2hpphs2'
             CALL INSMEM("@-S2HPPHS2",I030,MAXCOR)
             CALL ERREX
          ENDIF
C
C READ Hbar INTO Hbar(MA,EI) stored as Hbar(EM,ai) (ISIDE=1,RHF,UHF)
C READ Hbar INTO Hbar(IE,AM) stored as Hbar(ai,EM) (ISIDE=2,RHF,UHF)
C     
         CALL GETLST(ICORE(I020),1,NUMDSW,1,IRREPW,LISTW)
#ifdef _DEBUG_LVL0
      call checksum("->w     :",icore(i020),dissyw*numdsw)
#endif
         IF ((LISTW.EQ.54) .OR. (LISTW.EQ.55)) THEN
C STORED AS HBAR(EM,BI) or HBAR(em,bj) (ELEMENTS ARE <MA||IE> AND ma||ie> 
            CALL SSCAL(NUMDSW*DISSYW,ONEM,ICORE(I020),1)
         ENDIF
C
C PERFORM MATRIX MULTIPLICATION
C
C Z(AI,Jp) = SUM W(EM,AI) S(EM,Jp) [ISIDE = 1] 
C                ME
C Z(AI,Jp) = SUM W(AI,EM) S(EM,Jp) [ISIDE = 2] 
C                ME

         IF (ISIDE.EQ.1) THEN
            NROW = NUMDSW
         ELSE
            NROW=DISSYW
         ENDIF

         NCOL=NUMDSS
         NSUM=DISSYS

#ifdef _DEBUG_LVLM
             write(6,*)' sinspin', sinspin
             write(6,*)' soutspin', soutspin
             write(6,*)'  input w-coefficients, irrepw =',irrepw
             call output(icore(i020), 1, dissyw, 1, numdsw,
     $          dissyw, numdsw, 1)
             write(6,*)'  input s-coefficients'
             call output(icore(icount1), 1, dissys, 1, numdss,
     $          dissys, numdss,1)
#endif 
          IF (ISIDE.EQ.1) THEN
            CALL XGEMM('T','N',NROW,NCOL,NSUM,ONE,ICORE(I020),NSUM,
     &               ICORE(ICOUNT1),NSUM,ONE,ICORE(ICOUNT2),NROW)
          ELSE
            CALL XGEMM('N','N',NROW,NCOL,NSUM,ONE,ICORE(I020),NROW,
     &               ICORE(ICOUNT1),NSUM,ONE,ICORE(ICOUNT2),NROW)
          ENDIF
#ifdef _DEBUG_LVLM
          write(6,*)' output s-coefficients'
          call output(icore(icount2),1,nrow,1,numdss,nrow,numdss,1)
#endif 
          ICOUNT1 = ICOUNT1 + DISSYS * NUMDSS * IINTFP
          ICOUNT2 = ICOUNT2 + NROW   * NUMDSS * IINTFP
100    CONTINUE

#ifdef _DEBUG_LVL0
       call checksum("s2->s2  :",icore(i010),lens)
#endif
C
C REORDER S COEFFICIENTS BACK TO ORIGINAL FORM S[AI,JP] -> S[JI,AP]
C
       I030 = I020 + LENS *IINTFP
       IEND = I030
       IF (IEND .GT. MAXCOR) CALL INSMEM("S2PHPHS2",IEND,MAXCOR)
       CALL SSTGEN(ICORE(I010),ICORE(I020),LENS,VRT(1,SOUTSPIN),
     $      POP(1,ISPIN),POP(1,SOUTSPIN),NS(1),ICORE(I030),1,
     $      '3214')

C S2EX(JI,AP) -> S2EX(IJ,AP)

      SOFF = I020
      DO XIRREP = 1,NIRREP
         IRREPA = DIRPRD(XIRREP,SIRREP)
         DISSYS = IRPDPD(XIRREP,20+ISPIN)
         NUMDSS = VRT(IRREPA,MIXSPIN) * NS(SIRREP)

         MAXS  = MAX(DISSYS,NUMDSS)
         ITMP1 = I030
         ITMP2 = ITMP1 + MAXS*IINTFP
         ITMP3 = ITMP2 + MAXS*IINTFP
         IEND  = ITMP3 + MAXS*IINTFP
         IF (IEND .GT. MAXCOR) CALL INSMEM("S2PHPHS2",IEND,MAXCOR)
         CALL SYMTR3(XIRREP,POP(1,MIXSPIN),POP(1,ISPIN),DISSYS,
     &               NUMDSS,ICORE(SOFF),ICORE(ITMP1),ICORE(ITMP2),
     &               ICORE(ITMP3))
         SOFF = SOFF + DISSYS*NUMDSS
      ENDDO 
C
C ADD THE S INCREMENTS TO LISTS2EX
C
       LISTS2EX = LS2OUT(ISPIN,SOUTSPIN+1-IUHF)
       IF ((IUHF.NE.0) .AND. (SOUTSPIN.EQ.ISPIN)) THEN
          CALL ASSYMALL_HHP(ICORE(I020), LENS, NUMSZS, 
     $                      IRPDPD(1,20+ISPIN),POP(1,ISPIN),
     $                      ICORE(I030), MAXCOR-I030+1)
          CALL GETEXP2_HHP(ICORE(I010),LENS,NUMSZS,IRPDPD(1,20+ISPIN),
     $                     LISTS2EX,POP(1,ISPIN),
     $                     IRPDPD(1,ISYTYP(1, LISTS2EX)))
       ELSE
          CALL GETALLS2_HHP(ICORE(I010),LENS,VRT(1,SOUTSPIN),NS(1),
     $                      1,LISTS2EX)
       ENDIF       
#ifdef _DEBUG_LVL0
      call checksum("->hs2   :",icore(i010),lens)
      call checksum("s2->s2  :",icore(i020),lens)
#endif
       CALL DAXPY (LENS,ONE,ICORE(I020),1,ICORE(I010),1)
#ifdef _DEBUG_LVL0
      call checksum("s2->fs2 :",icore(i010),lens)
#endif
       IF ((IUHF.NE.0).AND. (SOUTSPIN.EQ.ISPIN)) THEN
          CALL PUTSQZ_HHP(ICORE(I010),LENS,NUMSZS,IRPDPD(1,20+ISPIN),
     $                    LISTS2EX,POP(1,ISPIN),
     $                    IRPDPD(1,ISYTYP(1,LISTS2EX)),
     $                    ICORE(I020),MAXCOR-I020+1)
       ELSE
          CALL PUTALLS2_HHP(ICORE(I010),LENS,VRT(1,SOUTSPIN),NS(1),1,
     $       LISTS2EX)
       ENDIF
C
 50    CONTINUE
 25   CONTINUE
C
      RETURN
      END
