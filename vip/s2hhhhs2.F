      SUBROUTINE S2HHHHS2(ICORE,MAXCOR,IUHF,ISIDE,ISPIN)
C
C  THIS SUBROUTINE CALCULATES THE CONTRIBUTION OF S2 TO
C  S2 INVOLVING THE HHHH BLOCK OF EFFECTIVE INTEGRALS W.
C
C CONTRACTION :
C
C     Z(JI,AP) =  SUM <MN|JI> S(NM,AP)
C                 N,M
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE,ONEM,ZILCH,HALF
C
      DIMENSION ICORE(MAXCOR)
      DIMENSION LISTW0(2,2)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SINFO/ NS(8), SIRREP
      COMMON/SLISTS/LS1IN, LS1OUT, LS2IN(2,2), LS2OUT(2,2)
C
      DATA HALF  /0.5/
      DATA ONE  /1.0D0/
      DATA ZILCH/0.0/
      DATA ONEM/-1.0/
C
      Write(6,"(a)") "-----------------S2HHHHS2-------------------"
      LISTW0(1,1) = 51
      LISTW0(1,2) = 53
      LISTW0(2,1) = 53
      LISTW0(2,2) = 52
C
C  IF ISPIN = 1  AND MIXSPIN =2 : SPECIAL CASE
C
C LOOP OVER IRREPS OF RHS OF S VECTOR
C
      DO 50 MIXSPIN = 1, 1+IUHF
         LISTW = LISTW0(ISPIN, MIXSPIN + 1 - IUHF)

         LISTS2EX = LS2OUT(ISPIN, MIXSPIN + 1 - IUHF)
         LISTS2IN = LS2IN(ISPIN,  MIXSPIN + 1 - IUHF)

         DO 100 IRREPA=1,NIRREP
            XIRREP =DIRPRD(IRREPA,SIRREP)
            IRREPW =XIRREP

            DISSYW=IRPDPD(IRREPW,ISYTYP(1,LISTW)) 
            NUMDSW=IRPDPD(IRREPW,ISYTYP(2,LISTW)) 

            DISSYS=IRPDPD(XIRREP,ISYTYP(1,LISTS2IN))
            DISSYSEX=IRPDPD(XIRREP,ISYTYP(1, LISTS2EX))
            NUMDSS=VRT(IRREPA,MIXSPIN)*NS(SIRREP)
            MAXS = MAX(NUMDSS, DISSYS)

            I000=1
            I010=I000
            I020=I010+IINTFP*DISSYS*NUMDSS
            I030=I020+IINTFP*DISSYSEX*NUMDSS
            I040 = I030 + MAXS * IINTFP
            I050 = I040 + MAXS * IINTFP
            I060 = I050 + MAXS * IINTFP
            ITOP = I060

            I070=ITOP+IINTFP*DISSYW*NUMDSW
            CALL DZERO(ICORE(I020),DISSYSEX*NUMDSS)
C
            IF(I070.GT.MAXCOR) CALL INSMEM("@-S2HHHHS2",I070,MAXCOR)

            CALL GETLST(ICORE(ITOP),1,NUMDSW,1,IRREPW,LISTW)

#ifdef _DEBUG_LVLM
            call checksum("W(ijkl) :",icore(itop),numdsw*dissyw)
            write(6,*)' W-integrals hhhh-block'
            write(6,*)' irrepw   = ', irrepw
            write(6,*)' imixspin =' , mixspin
            call output(icore(iTOP), 1, dissyw, 1, numdsw,
     $                  dissyw, numdsw,1)
#endif 
            CALL GETLST(ICORE(I010),1,NUMDSS,1,XIRREP,LISTS2IN)

#ifdef _DEBUG_LVL0
            call checksum("->s2    :",icore(i010),numdss*dissys)
            call checksum("->w     :",icore(itop),numdsw*dissyw)
#endif 
C
C IF ISPIN = 1 AND MIXSPIN = 2 THEN SPECIAL MULTIPLICATION:
C Hbar(Mn,Ij):S2(nM,aP)
            IF ((ISPIN.EQ.1.AND.MIXSPIN.EQ.2) .AND. IUHF.NE.0) THEN
               CALL SYMTR3(IRREPW,POP(1,ISPIN),POP(1,MIXSPIN),
     $                     DISSYS,NUMDSS,ICORE(I010),
     $                     ICORE(I030),ICORE(I040),ICORE(I050))
           ENDIF 
C PERFORM MATRIX MULTIPLICATION HBAR(NM,IJ)*S2IN(NM,AP)->S2EX(IJ,AP)
C
            NROW=NUMDSW 
            NCOL=NUMDSS
            NSUM=DISSYW
            IF (NROW*NCOL*NSUM .NE. 0) THEN

            IF (ISIDE.EQ.1) THEN
               CALL XGEMM('T','N',NROW,NCOL,NSUM,ONE,ICORE(ITOP),NROW,
     &                     ICORE(I010),NSUM,ONE,ICORE(I020),NROW)
            ELSE
C
C ISIDE = 2, LEFT HAND SIDE
C
               CALL XGEMM('N','N',NROW,NCOL,NSUM,ONE,ICORE(ITOP),NROW,
     &                    ICORE(I010),NSUM,ONE,ICORE(I020),NROW)       
             ENDIF
          ENDIF 
C
C Z(Ij,aP) --> Z(jI,aP)

            IF (ISPIN.EQ.1.AND.MIXSPIN.EQ.2 .AND. IUHF .NE.0) THEN
                CALL SYMTR3(IRREPW,POP(1,ISPIN),POP(1,MIXSPIN),
     $                      DISSYS,NUMDSS,ICORE(I020),
     $                      ICORE(I030),ICORE(I040),ICORE(I050))
            ENDIF 

            CALL GETLST(ICORE(I010),1,NUMDSS,1,XIRREP,LISTS2EX)
#ifdef _DEBUG_LVL0
          call checksum("s2->s2  :",icore(i020),numdss*dissys)
          call checksum("->hs2   :",icore(i010),numdss*dissys)
#endif 
           CALL SAXPY(NUMDSS*DISSYS,ONE,ICORE(I020),1, 
     $                ICORE(I010),1)
           CALL PUTLST(ICORE(I010),1,NUMDSS,1,XIRREP,LISTS2EX)

#ifdef _DEBUG_LVL0
           call checksum("s2->s2f :",icore(i010),numdss*dissys)
#endif 
C
 100   CONTINUE
 50    CONTINUE
C
       RETURN
       END
