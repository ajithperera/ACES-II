      SUBROUTINE S2HPHHS1(ICORE,MAXCOR,IUHF,ISIDE, ISPIN)
C
C   THIS SUBROUTINE CALCULATES THE CONTRIBUTION OF S2 TO
C   S1 INVOLVING THE HPHH BLOCK OF EFFECTIVE INTEGRALS W.
C
C CONTRACTION :
C
C     Z(I,P) = -SUM  [2<Mn|Ei>-<Mn|iE>] S(Nm,EP) [SPIN ADAPTED RHF]
C               Nm,E
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE,ONEM,ZILCH, FACTOR
C
      DIMENSION ICORE(MAXCOR),NI(8),NE(8)
      DIMENSION LISTW0(2,2)
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
C Hbar(NM,EI), lists 7-10;stored as NM,IE

CSSS      Write(6,"(a)") "-----------------S2HPHHS1-----------------"
      LISTW0(1,1) = 7
      LISTW0(1,2) = 10
      LISTW0(2,1) = 9
      LISTW0(2,2) = 8 
C
C  LIST 9 CONTAINS INTEGRALS IN THE ORDER (MN,EI) AND NEEDS SPECIAL
C  TREATMENT
C
       LENTAR= NS(SIRREP) * POP(SIRREP,ISPIN)

       I000=1
       I010=I000+LENTAR*IINTFP

       CALL ZERO(ICORE(I000),LENTAR)
C
C LOOP OVER IRREPS OF RHS OF S VECTOR
C
       DO 50 MIXSPIN =1, 1 + IUHF 
          LISTS2IN = LS2IN(ISPIN,MIXSPIN + 1 - IUHF)
          LISTS1EX = LS1OUT
          LISTW = LISTW0(ISPIN,MIXSPIN + 1 - IUHF) + 100 * (ISIDE-1)

          DO 100 IRREPE=1,NIRREP
             IRREPMN=DIRPRD(IRREPE,SIRREP)  
             XIRREP =IRREPMN
             IRREPEI=IRREPMN
             IRREPW =IRREPEI

             DISSYW=IRPDPD(IRREPMN,ISYTYP(1,LISTW)) 
             NUMDSW=VRT(IRREPE,MIXSPIN) * POP(SIRREP,ISPIN)

             DISSYS=IRPDPD(IRREPMN,ISYTYP(1,LISTS2IN))
             NUMDSS=VRT(IRREPE,MIXSPIN) * NS(SIRREP)

             MAXW=MAX(NUMDSW,DISSYW)
             I020=I010+IINTFP*DISSYW*NUMDSW
             IEND=I020+IINTFP*MAX(DISSYS*NUMDSS,3*MAXW)
             
             IF (IEND.GT.MAXCOR)THEN
                WRITE(6,*)' out-of-core AB not coded, s2hphhs1 '
                CALL INSMEM("@-S2HHPHS1",IEND,MAXCOR)
             ENDIF
C
C READ Hbar(MN,IE) as (MN,IE) amd transpose ket indices to (MN,EI)
C
             ITMP1=I020
             ITMP2=ITMP1+IINTFP*MAXW
             ITMP3=ITMP2+IINTFP*MAXW
             IEND =ITMP3+IINTFP*MAXW
             IF (IEND .GT. MAXCOR) CALL INSMEM("S2HPHHS1",IEND,MAXCOR)
C
        IF (MOD(LISTW,100). NE. 9) THEN

           NFIRST = 1
           DO 110 IRREPE0 = 1, IRREPE - 1
              IRREPI0 = DIRPRD(IRREPE0,IRREPEI)
              NFIRST = NFIRST + VRT(IRREPE0,MIXSPIN)*POP(IRREPI0,ISPIN)
 110       CONTINUE

           CALL GETLST(ICORE(I010),NFIRST,NUMDSW,1,IRREPW,LISTW)

           IF ( IUHF .EQ. 0) THEN
              CALL SPINAD3(IRREPMN,POP(1,1),DISSYW,NUMDSW,
     &                     ICORE(I010),ICORE(ITMP1),ICORE(ITMP2))
           ENDIF
C
C TRANSPOSE KET INDICES TO Hbar(MN,IE) -> Hbar(NM,EI)
C
          CALL IZERO(NE,8)
          CALL IZERO(NI,8)
          NI(SIRREP) = POP(SIRREP,ISPIN)
          NE(IRREPE) = VRT(IRREPE,MIXSPIN)
          CALL SYMTR1(IRREPEI,NI,NE,DISSYW,ICORE(I010),
     &                ICORE(ITMP1),ICORE(ITMP2),ICORE(ITMP3))

C TRANSPOSE BRA INDICES TO Hbar(MN,EI) -> Hbar(NM,EI)

         IF ((IUHF.EQ.0) .OR. (MIXSPIN .NE. ISPIN) ) THEN
            FACTOR = ONEM
            CALL SYMTR3(IRREPMN,POP(1,ISPIN),POP(1,MIXSPIN), DISSYW,
     $      NUMDSW,ICORE(I010),ICORE(ITMP1),ICORE(ITMP2), ICORE(ITMP3))
         ELSE
            FACTOR = ONE
         ENDIF
        ELSE
C
C Read list 9 AND 109, Hbar(Nm,Ei) is stored as (Nm,Ei)
C
           FACTOR = ONEM
           NFIRST = 1
           IRREPI = SIRREP
           DO 120 IRREPI0 = 1, IRREPI - 1
              IRREPE0 = DIRPRD(IRREPI0,IRREPEI)
              NFIRST = NFIRST + VRT(IRREPE0,MIXSPIN)*POP(IRREPI0,ISPIN)
 120       CONTINUE
           CALL GETLST(ICORE(I010),NFIRST,NUMDSW,1,IRREPW,LISTW)
          ENDIF
C
C READ S2 VECTOR; stored as S(Nm,eP) 
C
         CALL GETLST(ICORE(I020),1,NUMDSS,1,XIRREP,LISTS2IN)
C
C PERFORM MATRIX MULTIPLICATION
C
C                              +          
C         Z(I,P) = SUM W(Nm,Ei) * S(NmE,P)
C                  Efm
C
          NROW=POP(SIRREP,ISPIN)
          NCOL=NS(SIRREP)
          NSUM=DISSYW*VRT(IRREPE,MIXSPIN)

#ifdef _DEBUG_LVL0
      call checksum("s2-in   :",icore(i020),nsum)
      call checksum("W       :",icore(i010),NUMDSW*DISSYW)
#endif
#ifdef _DEBUG_LVLM
         if ((ncol*nsum*nrow).gt.0) then
            write(6,*)' s2phhhs1'
            write(6,"(a,3i3)")' iside, ispin, mspin ', iside, ispin, 
     $      mixspin
            write(6,*)' input s before multiplication, irrep', xirrep
            call output(icore(i020), 1, nsum, 1, ncol,
     $         nsum, ncol,1)
            write(6,*)' listw,input integrals ', listw
            call output(icore(i010), 1,nsum, 1, nrow,nsum, nrow, 1)
         endif
#endif 
          CALL XGEMM('T','N',NROW,NCOL,NSUM,FACTOR,ICORE(I010),NSUM,
     &               ICORE(I020),NSUM,ONE,ICORE(I000),NROW)

100    CONTINUE
#ifdef _DEBUG_LVL0
      call checksum("s2->s1  :",icore(i000),lentar)
c      write(6,"(a,2I4)") "Ispin,Mspin :", Ispin,Mixspin
c      call output(icore(i000),1,lentar,1,1,lentar,1,1)
#endif 
 50   CONTINUE
C
       CALL GETLST(ICORE(I010),1,1,1,ISPIN,LISTS1EX)
#ifdef _DEBUG_LVL0
      ndim= POP(SIRREP,ISPIN)*NS(SIRREP)*IINTFP
      call checksum("s1-in   :",icore(i010),Ndim)
C      call output(icore(i010),1,lentar,1,1,lentar,1,1)
#endif 
       CALL SAXPY (LENTAR,ONE,ICORE(I000),1,ICORE(I010),1)
       CALL PUTLST(ICORE(I010),1,1,1,ISPIN,LISTS1EX)

#ifdef _DEBUG_LVL0
      call checksum("s2->s1f :",icore(i010),lentar)
c      call output(icore(i010),1,lentar,1,1,lentar,1,1)
#endif
C
      RETURN
      END
