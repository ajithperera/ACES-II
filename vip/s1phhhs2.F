      SUBROUTINE S1PHHHS2(ICORE,MAXCOR,IUHF,ISIDE, ISPIN)
C
C   THIS SUBROUTINE CALCULATES THE CONTRIBUTION OF S1 TO
C   S2 INVOLVING THE PHHH BLOCK OF EFFECTIVE INTEGRALS W.
C
C CONTRACTION :
C
C     Z(JI,AP) = -SUM  <AM||JI> s1(M,Ip)
C                 M
C     Z(JI,AP) = -SUM  (2.0<Am|Ji>-<Am|iJ>) s1(M,Ip)
C                 M
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE,ONEM,ZILCH, FACTOR
      logical print
C
      DIMENSION ICORE(MAXCOR),NA(8),NM(8)
      DIMENSION LISTW0(2,2)
C
C  LISTW0 IS PARAMETRIZED AS (IVRTSPIN=ISPIN, IMIXSPIN)
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
CSSS      Write(6,"(a)") "-----------------S1PHHHS2-----------------"
C
C Hbar(KA,IJ), lists 107-110, stored as I<J,KA, Hbar(Ak,Ij) Ij,Ak (109),
C Hbar(Ka,Ij) Ij,K,a (110)
C     
      LISTW0(1,1) = 07 
      LISTW0(1,2) = 10
      LISTW0(2,1) = 09
      LISTW0(2,2) = 08
C
      LISTS1IN = LS1IN
C
C READ IN S1(I,P) COEFFICIENTS
C
      I000 = 1
      I010 = I000 + POP(SIRREP,ISPIN) * NS(SIRREP) * IINTFP
C
C GET S1 COEFFIENTS
C
      CALL GETLST(ICORE(I000),1,1,1,ISPIN,LISTS1IN)
C
      DO 50 IMIXSPIN = 1, 1+IUHF
         MSPIN = ISPIN
         ASPIN = IMIXSPIN
         LISTS2EX = LS2OUT(ISPIN, IMIXSPIN+1-IUHF)
         LISTW = LISTW0(ISPIN, IMIXSPIN+1-IUHF) + (2 - ISIDE) * 100
C
C LOOP OVER IRREPS OF RHS OF W VECTOR (JI,AM) (S1 is (M,Ip))
C
         DO 100 IRREPA=1,NIRREP

            XIRREP=DIRPRD(IRREPA,SIRREP)
            IRREPW=XIRREP

            DISSYW=IRPDPD(IRREPW,ISYTYP(1,LISTW)) 
            NUMDSW=VRT(IRREPA,ASPIN)*POP(SIRREP,MSPIN)
            DISSYS=IRPDPD(XIRREP,ISYTYP(1,LISTS2EX))
            NUMDSS=VRT(IRREPA,ASPIN) * NS(SIRREP)

            MAXW=MAX(NUMDSW,DISSYW)
            I020=I010+IINTFP*DISSYW*NUMDSW
            I030=I020+IINTFP*MAX(DISSYS*NUMDSS,3*MAXW)
            IEND=I030
            IF(IEND.GT.MAXCOR)THEN
C
C OUT-OF-CORE ALGORITHM

               Write(6,"(2a)") "@-S1PHHHS2 out-of-core algorithm is ",
     +                         "not coded" 
               CALL INSMEM("S1PHHHS2",IEND,MAXCOR)
            ENDIF
C
C DO IN-CORE ALGORITHM.  READ W  
C
        IF (MOD(LISTW,100).EQ.9) THEN
C
C SPECIAL TREATMENT OF LISTS 109 <Ji|Am>. ONLY THE INTEGRALS CIRREP=SIRREP, 
C MIRREP ARE NEEDED THEREFORE CALCULATE OFFSET.
C
           NFIRST = 1
           DO 110 IRREPM0 = 1, (SIRREP -1)
              IRREPA0 = DIRPRD(IRREPM0,XIRREP)
              NFIRST = NFIRST + VRT(IRREPA0,ASPIN) *
     +                          POP(IRREPM0,MSPIN)
 110          CONTINUE
           CALL GETLST(ICORE(I010),NFIRST,NUMDSW,1,IRREPW,LISTW)
           FACTOR = ONEM
C     
         ELSE
C
C READ <JI|AM>=<IJ|MA> stored as <IJ,MA> -> (JI,AM)
C AND TRANSPOSE BOTH BRA AND KET INDICES TO YIELD W(JI,AM)
C ONLY THE INTEGRALS CIRREP=SIRREP, MIRREP ARE NEEDED
C THEREFORE CALCULATE OFFSET.
C
           NFIRST = 1
           DO 115 IRREPA0 = 1, (IRREPA -1)
              IRREPM0 = DIRPRD(IRREPA0, XIRREP)
              NFIRST  = NFIRST + VRT(IRREPA0,ASPIN) *
     +                           POP(IRREPM0,MSPIN)
 115       CONTINUE
          CALL GETLST(ICORE(I010),NFIRST,NUMDSW,1,IRREPW,LISTW)

          ITMP1=I020
          ITMP2=ITMP1+IINTFP*MAXW
          ITMP3=ITMP2+IINTFP*MAXW
          IEND=ITMP3+IINTFP*MAXW
          IF (IEND .GT. MAXCOR) CALL INSMEM("S1PHHHS2",IEND,MAXCOR)
C
C W(IJ,MA) -> W(IJ,AM)
C
          CALL IZERO(NA,8)
          CALL IZERO(NM,8)
          NA(IRREPA) =VRT(IRREPA,ASPIN)
          NM(SIRREP) =POP(SIRREP,MSPIN)
          CALL SYMTR1(IRREPW,NM,NA,DISSYW,ICORE(I010),
     &                ICORE(ITMP1),ICORE(ITMP2),ICORE(ITMP3))
C 
C W(IJ,AM) -> W(JI,AM).IF ISPIN = IMIXSPIN THIS AMOUNTS TO A MINUS SIGN
C
          IF ((IUHF.EQ.0) .OR. (ISPIN.NE.IMIXSPIN)) THEN
             FACTOR = ONEM
             CALL SYMTR3(IRREPW,POP(1,MSPIN),POP(1,ASPIN),DISSYW,
     $          NUMDSW,ICORE(I010),ICORE(ITMP1),ICORE(ITMP2),
     $          ICORE(ITMP3))
          ELSE
             FACTOR = ONE
          ENDIF
       ENDIF
C
C READ S2 VECTOR TO BE UPDATED
C
         CALL GETLST(ICORE(I020),1,NUMDSS,1,XIRREP,LISTS2EX)
C
C PERFORM MATRIX MULTIPLICATION
C
C         
C Z(JI,AM) = SUM W(JIA,M) * S(M,P) 
C            M
C  
          NROW=DISSYW*VRT(IRREPA,ASPIN)
          NCOL=NS(SIRREP)
          NSUM=POP(SIRREP,ISPIN)

#ifdef _DEBUG_LVL0
         print = .false.
         if ((print) .and.( (ncol*nsum*nrow).gt.0)) then
            write(6,*)' output s before multiplication, irrep', xirrep
            write(6,*)' iside, ispin, mspin ', iside, ispin, mspin
            write(6,*)' listw, lists2ex', listw, lists2ex
            call output(icore(i020), 1, nrow, 1, ncol,
     $         nrow, ncol,1)
            write(6,*)' input integrals'
            call output(icore(i010), 1,nrow, 1, nsum, nrow, nsum, 1)
         endif
#endif 
#ifdef _DEBUG_LVL0
      call checksum("w       :",icore(i010),dissyw*numdsw)
      call checksum("s1-in   :",icore(i000),Ncol*Nsum)
      call checksum("s2-in   :",icore(i020),Nrow)
#endif 
          CALL XGEMM('N','N',NROW,NCOL,NSUM,FACTOR,ICORE(I010),NROW,
     &               ICORE(I000),NSUM,ONE,ICORE(I020),NROW)
          CALL PUTLST(ICORE(I020),1,NUMDSS,1,XIRREP,LISTS2EX)
#ifdef _DEBUG_LVL0
      call checksum("s1->s2  :",icore(i020),Nrow*Ncol)
c      call output(icore(i020),1,nrow,1,ncol,nrow,ncol,1)
#endif 
 100   CONTINUE
       If (IUHF.EQ. 0 .AND. ISPIN .EQ. MSPIN) THEN
CSSS           CALL S2AB_2S2AA(ICORE(I000),MAXCOR,ISPIN)
       ENDIF 
 50   CONTINUE
C
      RETURN
      END
