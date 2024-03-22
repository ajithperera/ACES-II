      SUBROUTINE S2HPS1(ICORE,MAXCOR,IUHF, ISIDE, ISPIN)
C
C   THIS SUBROUTINE CALCULATES THE CONTRIBUTION OF S2 TO
C   S1 INVOLVING THE HP BLOCK OF THE EFFECTIVE ONE-PARTICLE 
C   INTEGRALS.
C
C CONTRACTION :
C
C     Z(I,P) = + SUM [Hbar(m,e) (2 S(Mi,ep) - S(iM,ep)) [ RHF CASE ]
C                m,e
C     ONLY FOR ISIDE = 1
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE,ONEM,ZILCH
      LOGICAL SKIP
C
      DIMENSION ICORE(MAXCOR), NUMSZS(8)
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
      IF (ISIDE.EQ.2) THEN
         RETURN
      ENDIF
CSSS      Write(6,"(a)") "-----------------S2HPS1-------------------"
C
      DO 5 XIRREP = 1, NIRREP
         IRREPE = DIRPRD(XIRREP,SIRREP)
         NUMSZS(XIRREP) = VRT(IRREPE,ISPIN) * NS(SIRREP)
 5    CONTINUE      
C
C MIXSPIN indicates the spin of EM and ISPIN indicates the spin
C of I.
C
      DO 25 MIXSPIN = 1, 1+IUHF
         LISTS2IN = LS2IN(ISPIN, MIXSPIN + 1 - IUHF)
         CALL GETLEN_HHP(LENS, POP(1,MIXSPIN), POP(1,ISPIN),
     $                   VRT(1,MIXSPIN), NS)

       I000=1
       I010=I000+LENS*IINTFP
       IF ((IUHF.NE.0) .AND. (MIXSPIN.EQ.ISPIN)) THEN
          CALL GETEXP2_HHP(ICORE(I000),LENS,NUMSZS,
     +                     IRPDPD(1,20+ISPIN),
     +                     LISTS2IN,POP(1,ISPIN), 
     +                     IRPDPD(1,ISYTYP(1,LISTS2IN)))
       ELSE
          CALL GETALLS2_HHP(ICORE(I000),LENS,VRT(1,MIXSPIN),NS(1),1,
     $                      LISTS2IN)
       ENDIF
#ifdef _DEBUG_LVL0
      call checksum("->s2    :",icore(i000),Lens)
#endif 
C
C SPINADAPT THE S COEFFICIENTS BY IRREP
C
       IF (IUHF.EQ.0) THEN
          ICOUNT = I000
          DO 10 XIRREP = 1,NIRREP
             IRREPE=DIRPRD(XIRREP,SIRREP)
             DISSYS=IRPDPD(XIRREP,ISYTYP(1,LISTS2IN))
             NUMDSS=VRT(IRREPE,1) * NS(SIRREP)
             I020 = I010 + NUMDSS * IINTFP
             IEND = I020 + NUMDSS * IINTFP
             IF (IEND .GT. MAXCOR) CALL INSMEM("S2HPS1",IEND,MAXCOR)
             CALL SPINAD3(XIRREP,POP(1,1),DISSYS,NUMDSS,
     $                   ICORE(ICOUNT), ICORE(I010), ICORE(I020))
             ICOUNT = ICOUNT + NUMDSS * DISSYS * IINTFP
 10       CONTINUE
       ENDIF
C
C REORDER THE S COEFFICIENT ACCORDING TO S(MI,EP) -> S(ME,IP)
C
       I020 = I010 + LENS*IINTFP
       IEND = I020
       IF(IEND.GT.MAXCOR)THEN
C
C OUT-OF-CORE ALGORITHM
C
          CALL INSMEM("@-S2HPS1",IEND,MAXCOR)
          CALL ERREX 
        ENDIF
C
C S2(MI,EP) --> S2(ME,IP)
C
       CALL SCOPY(LENS,ICORE(I000),1,ICORE(I010),1)

#ifdef _DEBUG_LVLM
       call checksum("->s2    :",icore(i010),Lens)
#endif 
       CALL SSTGEN(ICORE(I010),ICORE(I000),LENS,POP(1,MIXSPIN),
     $      POP(1,ISPIN),VRT(1,MIXSPIN),NS,ICORE(I020),1,'1324')

C ICORE(I000) CONTAINS ALL [SPIN ADAPTED] S COEFFICIENTS [2S(Me,iP) -
C S(eM,iP)], ORDERED AS Me,iP
C
      I020 = I010 + NT(MIXSPIN) * IINTFP
      IEND = I020 + NT(MIXSPIN) * IINTFP
      IF (IEND .GT. MAXCOR) CALL INSMEM("@-S2HPS1",IEND,MAXCOR) 

       CALL GETLST(ICORE(I020),1,1,1,MIXSPIN,93)
C
C ONLY THE FIRST IRREP OF S VECTOR IS NEEDED (IRREP(ME) = 1)
C
       DISSYS= IRPDPD(1,15+MIXSPIN)
       IF (IUHF .EQ. 0) DISSYS = IRPDPD(1,18)
       NUMDSS= POP(SIRREP,ISPIN) * NS(SIRREP)

       IOFF = I020
       JOFF = I010
       DO IRREP_R = 1, NIRREP
          IRREP_L = DIRPRD(IRREP_R,1)
          NROW = VRT(IRREP_L,MIXSPIN)
          NCOL = POP(IRREP_R,MIXSPIN)
          CALL TRANSP(ICORE(IOFF),ICORE(JOFF),NCOL,NROW)
          IOFF = IOFF + NROW*NCOL
          JOFF = JOFF + NROW*NCOL
       ENDDO

       CALL GETLST(ICORE(I020), 1, 1, 1, ISPIN, LS1OUT)

#ifdef _DEBUG_LVL0
      call checksum("->hs1   :",icore(i020),numdss) 
C      call output(icore(i020),1,1,1,numdss,1,numdss,1)
#endif 
C
C PERFORM MATRIX MULTIPLICATION
C
C  Z(1,IP) = +SUM U(ME,1) * S(ME,IP) [ISIDE = 1 OR 2]
C             ME
C  
#ifdef _DEBUG_LVL0
      call checksum("fme     :",icore(i010),dissys)
C      ioff = i010
C      do irrep_r = 1, Nirrep
C         irrep_l = Dirprd(Irrep_r,1)
C         nrow = pop(irrep_l,mixspin)
C         ncol = vrt(irrep_r,mixspin)
C         call output(icore(ioff),1,nrow,1,ncol,nrow,ncol,1)
C         ioff = ioff + Nrow*Ncol
C      enddo 
#endif 
          NROW=1
          NCOL=NUMDSS
          NSUM=DISSYS
          CALL XGEMM('T','N',NROW,NCOL,NSUM,ONE,ICORE(I010),NSUM,
     &                 ICORE(I000),NSUM,ONE,ICORE(I020),NROW)
C
           CALL PUTLST(ICORE(I020),1, 1, 1, ISPIN, LS1OUT)

#ifdef _DEBUG_LVL0
      call checksum("s2>s1   :",icore(i020),numdss) 
c      call output(icore(i020),1,numdss,1,1,numdss,1,1)
#endif 
C
 25   CONTINUE
C
      RETURN
      END
