      SUBROUTINE S2PHPPS1(ICORE,MAXCOR,IUHF,ISIDE, ISPIN)
C
C   THIS SUBROUTINE CALCULATES THE CONTRIBUTION OF S2 TO
C   S1 INVOLVING THE PHPP BLOCK OF EFFECTIVE INTEGRALS W.
C
C CONTRACTION :
C
C     Z(A,P) = SUM  [2<Am|Ef>-<Am|Fe>] S(fE,mP) [SPIN ADAPTED RHF]
C              mef
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE,ONEM,ZILCH, FACTOR
      logical print
C
      DIMENSION ICORE(MAXCOR),NA(8),NM(8)
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
C     30: AbCi [V -> W(Ci, Ab)], ordered according to Ab;Ci
C     130: AbCi, Wintegrals, ordered according to Ab; Ci
C
      LISTW0(1,1) = 27
      LISTW0(1,2) = 30
      LISTW0(2,1) = 29
      LISTW0(2,2) = 28
C
C  LIST 29 CONTAINS INTEGRALS IN THE ORDER AB, IC AND NEEDS SPECIAL
C  TREATMENT
C
       LENTAR= NS(SIRREP) * VRT(SIRREP,ISPIN)
       I000=1
       I010=I000+LENTAR*IINTFP
       CALL ZERO(ICORE(I000),LENTAR)
C
C LOOP OVER IRREPS OF RHS OF S VECTOR
C
       DO 50 IMIXSPIN = 1, 1 + IUHF
          MSPIN = IMIXSPIN
          LISTS2IN = LS2IN(ISPIN, IMIXSPIN + 1 - IUHF)
          LISTS1EX = LS1OUT
          LISTW = LISTW0(ISPIN, IMIXSPIN + 1 - IUHF) + 100 * (ISIDE-1)
       DO 100 MIRREP=1,NIRREP
        XIRREP=DIRPRD(MIRREP,SIRREP)
        IRREPW=XIRREP
        DISSYW=IRPDPD(IRREPW,ISYTYP(1,LISTW)) 
        NUMDSW=VRT(SIRREP,ISPIN) * POP(MIRREP,MSPIN)
        DISSYS=IRPDPD(XIRREP,ISYTYP(1,LISTS2IN))
        NUMDSS=POP(MIRREP,MSPIN) * NS(SIRREP)
        MAXW=MAX(NUMDSW,DISSYW)
        I020=I010+IINTFP*DISSYW*NUMDSW
        I030=I020+IINTFP*MAX(DISSYS*NUMDSS,3*MAXW)
        IF(I030.GT.MAXCOR)THEN
C
C OUT-OF-CORE ALGORITHM
C
         WRITE(6,*)' out-of-core AB not coded, s2phpps1 '
         WRITE(6,1000) I030, MAXCOR
 1000    FORMAT('  REQUIRED MEMORY :', I12, '  AVAILABLE : ', I12)
         call errex
C
        ENDIF
C
C
C DO IN-CORE ALGORITHM
C
C
C READ W INTO 
C     
C        W(efam) =  2 W(Ef,Am) - W(FeAm)    [ RHF CASE ]
C
C        ONLY THE AIRREP=SIRREP, AND MIRREP INTEGRALS ARE NEEDED.
C        THEREFORE CALULATE OFFSET
C
          ITMP1=I020
          ITMP2=ITMP1+IINTFP*MAXW
          ITMP3=ITMP2+IINTFP*MAXW
          ITMP4=ITMP3+IINTFP*MAXW
C
        IF (MOD(LISTW,100). NE.29) THEN
           NFIRST = 1
           DO 110 IRREPM = 1, MIRREP - 1
              AIRREP = DIRPRD(IRREPM, XIRREP)
              NFIRST = NFIRST + VRT(AIRREP,ISPIN)*POP(IRREPM,MSPIN)
 110          CONTINUE
         CALL GETLST(ICORE(I010),NFIRST,NUMDSW,1,IRREPW,LISTW)
C
          IF ( IUHF .EQ. 0) THEN
             CALL SPINAD3(IRREPW,VRT(1,1),DISSYW,NUMDSW,
     &                  ICORE(I010),ICORE(ITMP1),ICORE(ITMP2))
          ENDIF
C
C TRANSPOSE KET INDICES TO W(Efm,A) 
C WE ONLY HAVE 1 PARTICULAR SYMMETRY TYPE OF THE A AND M INDICES, AND
C THIS SHOULD BE PASSED TO SYMTR1
C
          CALL IZERO(NA, 8)
          CALL IZERO(NM, 8)
              NA(SIRREP) =VRT(SIRREP,ISPIN)
              NM(MIRREP)=POP(MIRREP,MSPIN)
          CALL SYMTR1(IRREPW,NA,NM,DISSYW,ICORE(I010),
     &                ICORE(ITMP1),ICORE(ITMP2),ICORE(ITMP3))
C
C TRANSPOSE BRA INDICES TO W(fEm,A), TO GET ORDERING AS IN S2EX
C
          IF ((IUHF.EQ.0) .OR. (IMIXSPIN .NE. ISPIN) ) THEN
            FACTOR = ONE
            CALL SYMTR3(IRREPW,VRT(1,ISPIN), VRT(1,MSPIN), DISSYW,
     $      NUMDSW,ICORE(I010),ICORE(ITMP1),ICORE(ITMP2), ICORE(ITMP3))
         ELSE
            FACTOR = ONEM
         ENDIF
       ELSE
C
C  TREAT 29 AND 129 AS SPECIAL CASE
C
             FACTOR = ONE
             NFIRST = 1
             AIRREP = SIRREP
             DO 120 IRREPA = 1, AIRREP - 1
                IRREPM = DIRPRD(IRREPA, XIRREP)
                NFIRST = NFIRST + VRT(IRREPA,ISPIN)*POP(IRREPM,MSPIN)
 120         CONTINUE
             CALL GETLST(ICORE(I010),NFIRST,NUMDSW,1,IRREPW,LISTW)
          ENDIF
C
C READ S2 VECTOR INTO 
C
C           S(fEm,P) 
C
         CALL GETLST(ICORE(I020),1,NUMDSS,1,XIRREP,LISTS2IN)
C
C PERFORM MATRIX MULTIPLICATION
C
C                              +          
C         Z(A,P) = SUM W(fEm,A) * S(fEm,P)
C                  Efm
          NROW=VRT(SIRREP,ISPIN)
          NCOL=NS(SIRREP)
          NSUM=DISSYW*POP(MIRREP,MSPIN)
          write(6,"(a,3i4)") "nrow,ncol.nsum: ",nrow,ncol,nsum
          call checksum("s2-in   :",icore(i020),nsum)

         print = .false.
         if ((print) .and.( (ncol*nsum*nrow).gt.0)) then
            write(6,*)' s2phpps1'
            write(6,*)' input s before multiplication, irrep', xirrep
            write(6,*)' iside, ispin, mspin ', iside, ispin, mspin
            write(6,*)' listw ', listw
            call output(icore(i020), 1, nsum, 1, ncol,
     $         nsum, ncol,1)
            write(6,*)' input integrals'
            call output(icore(i010), 1,nsum, 1, nrow,nsum, nrow, 1)
         endif
          CALL XGEMM('T','N',NROW,NCOL,NSUM, FACTOR, ICORE(I010),NSUM,
     &               ICORE(I020),NSUM,ONE,ICORE(I000),NROW)
           call checksum("s2-s1   :",icore(i000),lentar)
100    CONTINUE
          print*, "ispin :", ispin, mspin
          call output(icore(i000),1,lentar,1,1,lentar,1,1)
 50   CONTINUE
C
       CALL GETLST(ICORE(I010),1,1,1,ISPIN,LISTS1EX)
       print = .true.
       if (print) then
          write(6,*)' output s1ex before leaving routine'
          nr = vrt(sirrep, ispin)
          nc = ns(sirrep)
          call output(icore(i010), 1, nr, 1, nc, nr, nc, 1)
          write(6,*)' correction calculated in s2phpps1'
          call output(icore(i000), 1, nr, 1, nc, nr, nc, 1)
       endif
       CALL SAXPY (LENTAR,ONE,ICORE(I000),1,ICORE(I010),1)
       CALL PUTLST(ICORE(I010),1,1,1,ISPIN,LISTS1EX)
       if (print) then
          write(6,*)' @-s2phpps1 exit'
          call output(icore(i010), 1, nr, 1, nc, nr, nc, 1)
       endif 
C
      RETURN
      END
