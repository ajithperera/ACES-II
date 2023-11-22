      SUBROUTINE BUILTDIR(EAMAT, NEA, ICORE, DMAXCOR, IUHF,
     $   ISIDE, ISPIN, NEAST)
C
C  THE EA-EOM MATRIX IS CONSTRUCTED IN A DIRECT FASHION 
C  BY APPLYING EADIR TO THE IDENTITY MATRIX IN EA-EOM SPACE, 
C  EITHER FROM THE LEFT (ISIDE = 1) OR FROM THE RIGHT (ISIDE =2)
C
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE, ONEM,ZILCH, EAMAT, ICORE
      logical prints, printea, diag, lefthand, excicore,
     $   singonly, DROPCORE, ADDST, skip,skip1
C
C    ICORE IS USED AS DOUBLE PRECISION ARRAY
C
      DIMENSION ICORE(DMAXCOR), EAMAT(NEA, NEA), NDUMS(8)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SINFO/ NS(8), SIRREP
      COMMON/SLISTS/LS1IN, LS1OUT, LS2IN(2,2), LS2OUT(2,2)
      COMMON/EACALC/LEFTHAND, EXCICORE, SINGONLY, DROPCORE
      COMMON/COREINFO/IREPCORE, SPINCORE, IORBCORE, IORBOCC
      COMMON/STINFO/ITOTALS, STMODE, LISTST, STCALC, NSIZEST
C
      DATA ONE, ZILCH, ONEM /1.0D0, 0.0D0,  -1.0D0/
C
      MAXCOR = IINTFP* DMAXCOR
C
C
C  PUT LS1OUT AND LS2OUT TO ZERO
C
      CALL SETZERO(IUHF, ISPIN, ICORE, MAXCOR, LS1OUT,
     $   LS2OUT)
      CALL SETZERO(IUHF, 3-ISPIN, ICORE, MAXCOR, LS1OUT,
     $   LS2OUT)
C
C  PUT LS1IN AND LS2IN TO ZERO
C
      CALL SETZERO(IUHF, ISPIN, ICORE, MAXCOR, LS1IN,
     $   LS2IN)
      CALL SETZERO(IUHF, 3-ISPIN, ICORE, MAXCOR, LS1IN,
     $   LS2IN)
C
      CALL IZERO(NS,8)
      NS(SIRREP) = NEA
      ADDST = .FALSE.
      IF (EXCICORE) THEN
         ADDST = (STMODE.EQ.1 .AND. (ISPIN .EQ. SPINCORE))
      ENDIF
C
C PUT THE ELEMENTS IN LISTS1IN IN IDENTITIY FORM      
C
      PRINTS = .FALSE.
      LISTS1IN = LS1IN
      I000 = 1
      NS1 = VRT(SIRREP,ISPIN) * NEA
      I010 = I000 + NS1
      CALL ZERO(ICORE(I000), NS1)
      ICOUNT = I000
      DO 10 I = 1, VRT(SIRREP,ISPIN)
         ICORE(ICOUNT) = ONE
         ICOUNT = ICOUNT + VRT(SIRREP,ISPIN) + 1
 10   CONTINUE
      IF (PRINTS) THEN
         WRITE(6,*)' S1 PART OF IDENTITY MATIX'
         CALL OUTPUT(ICORE(I000), 1, VRT(SIRREP, ISPIN), 1, NEA,
     $      VRT(SIRREP, ISPIN), NEA, 1)
      ENDIF
      CALL PUTLST(ICORE(I000), 1,1,1, ISPIN, LISTS1IN)
C
      IF (ADDST) THEN
C
C  ADD SPIN-FLIP COMPONENTS TO IDENTITY MATRIX
C
         DSPIN = 3 - ISPIN
         NS1 = VRT(SIRREP,DSPIN) * NEA
         I010 = I000 + NS1
         CALL ZERO(ICORE(I000), NS1)
         ICOUNT = I000 + VRT(SIRREP,DSPIN) * (NEA-NEAST)
         DO 11 I = 1, VRT(SIRREP,DSPIN)
         ICORE(ICOUNT) = ONE
         ICOUNT = ICOUNT + VRT(SIRREP,DSPIN) + 1
 11   CONTINUE
      IF (PRINTS) THEN
         WRITE(6,*)' S1 PART OF IDENTITY MATIX'
         CALL OUTPUT(ICORE(I000), 1, VRT(SIRREP, DSPIN), 1, NEA,
     $      VRT(SIRREP, DSPIN), NEA, 1)
      ENDIF
      CALL PUTLST(ICORE(I000), 1,1,1, DSPIN, LISTS1IN)
C
      ENDIF
C
C  NOW PUT THE ELEMENTS IN LISTS2IN IN IDENTITY FORM.
C  THE REVERSE OF THE SCHEME BELOW, TO CONSTRUCT EAMAT FROM THE
C  INFORMATION ON LISTS2EX. THIS SCHEME IS 'NECESSARY' DUE TO THE
C  COMPLICATED STRUCTURE OF THE S2-LISTS.
C
      CALL ZERO(EAMAT, NEA*NEA)
      DO 2 I = 1 , NEA
         EAMAT(I,I) = ONE
 2    CONTINUE
      skip1 = .false.
      if (skip1) then
         write(6,*) ' original higher excitation operators are skipped '
      else
      NC = VRT(SIRREP, ISPIN)
      CALL IZERO(NDUMS, 8)
      NDUMS(SIRREP) = 1
      IF (IUHF . NE. 0) THEN
         CALL GETLEN2(LEN2, IRPDPD(1, ISPIN), POP(1,ISPIN), NDUMS)
      ENDIF
      DO 20 MSPIN = 1, 1 + IUHF
         IF (IUHF.EQ.0) THEN
            ICOL = NC + 1
            ITYP = 19
         ELSE
            IF (ISPIN .EQ.  MSPIN) THEN
               ICOL = NC + 1
               ITYP = ISPIN
            ELSE
               ICOL = NC + LEN2 + 1
               ITYP = 15
            ENDIF
         ENDIF
         LISTS2IN = LS2IN(ISPIN, MSPIN + 1 - IUHF)
         DO 30 XIRREP = 1, NIRREP
         IIRREP = DIRPRD(XIRREP, SIRREP)
         NCOL = POP(IIRREP,MSPIN) * IRPDPD(XIRREP, ITYP)
         IF (NCOL.GT.0) THEN
         NUMDSS = POP(IIRREP,MSPIN) * NEA
         I010 = I000 + NCOL * NEA
         I020 = I010 + NCOL * NEA
         IF (I020.GT.DMAXCOR) THEN
            WRITE(6,*)'  OUT OF CORE ALGORITM NEEDED'
            CALL ERREX
         ENDIF
         CALL ZERO(ICORE(I000), NEA* NCOL)
         CALL SAXPY(NEA*NCOL, ONE, EAMAT(1, ICOL),1, ICORE(I000), 1)
         CALL TRANSP(ICORE(I000), ICORE(I010), NCOL, NEA)
         CALL PUTLST(ICORE(I010),1, NUMDSS, 1, XIRREP, LISTS2IN)
         ENDIF
         ICOL = ICOL + NCOL
 30   CONTINUE
 20   CONTINUE
      endif
C************
      IF (ADDST) THEN
         skip = .false.
         if (skip) then
            write(6,*) ' higher spin flip operators are put to zero'
         else
         DSPIN = 3 -ISPIN
      NC = VRT(SIRREP, DSPIN) + (NEA - NEAST)
      CALL IZERO(NDUMS, 8)
      NDUMS(SIRREP) = 1
      IF (IUHF . NE. 0) THEN
         CALL GETLEN2(LEN2ST, IRPDPD(1, DSPIN), POP(1,DSPIN), NDUMS)
      ENDIF
      DO 21 MSPIN = 1, 1 + IUHF
         IF (IUHF.EQ.0) THEN
            ICOL = NC + 1
            ITYP = 19
         ELSE
            IF (DSPIN .EQ.  MSPIN) THEN
               ICOL = NC + 1
               ITYP = DSPIN
            ELSE
               ICOL = NC + LEN2ST + 1
               ITYP = 15
            ENDIF
         ENDIF
         LISTS2IN = LS2IN(DSPIN, MSPIN + 1 - IUHF)
         DO 31 XIRREP = 1, NIRREP
         IIRREP = DIRPRD(XIRREP, SIRREP)
         NCOL = POP(IIRREP,MSPIN) * IRPDPD(XIRREP, ITYP)
         IF (NCOL.GT.0) THEN
         NUMDSS = POP(IIRREP,MSPIN) * NEA
         I010 = I000 + NCOL * NEA
         I020 = I010 + NCOL * NEA
         IF (I020.GT.DMAXCOR) THEN
            WRITE(6,*)'  OUT OF CORE ALGORITM NEEDED'
            CALL ERREX
         ENDIF
         CALL ZERO(ICORE(I000), NEA* NCOL)
         CALL SAXPY(NEA*NCOL, ONE, EAMAT(1, ICOL),1, ICORE(I000), 1)
         CALL TRANSP(ICORE(I000), ICORE(I010), NCOL, NEA)
         CALL PUTLST(ICORE(I010),1, NUMDSS, 1, XIRREP, LISTS2IN)
         ENDIF
         ICOL = ICOL + NCOL
 31   CONTINUE
 21   CONTINUE
      endif
      ENDIF
C
C   NOW PERFORM MULTIPLICATION -> EAMAT IMPLICIT ON LISTS(I)EX
C
      MAXCOR = DMAXCOR * IINTFP
C
      CALL PROJECTS(ISPIN, IUHF, ICORE, DMAXCOR, LS1IN, LS2IN,
     $   STMODE.NE.0, DROPCORE)
C      IF (ADDST) THEN
C         CALL PROJECTS(DSPIN, IUHF, ICORE, DMAXCOR, LS1IN, LS2IN,
C     $      ADDST,DROPCORE)
C      ENDIF
C
C        CALL SINGTRIP(ICORE, MAXCOR, IUHF, ISIDE, ISPIN)
C      IF (addst) then
C         CALL SINGTRIP(ICORE, MAXCOR, IUHF, ISIDE, DSPIN)
C       endif
C
      CALL EADIR(ICORE, MAXCOR, IUHF, ISIDE, ISPIN)
C
      CALL PROJECTS(ISPIN, IUHF, ICORE, DMAXCOR, LS1OUT, LS2OUT,
     $   STMODE.NE.0,DROPCORE)
C      IF (ADDST) THEN
C         CALL PROJECTS(DSPIN, IUHF, ICORE, DMAXCOR, LS1OUT, LS2OUT,
C     $      ADDST, DROPCORE)
C      ENDIF
C
C  NOW CONSTRUCT EAMAT FROM THE INFORMATION ON LISTS(I)EX  
C
      NC = VRT(SIRREP, ISPIN)
      I010 = I000 + NC * NEA
      I020 = I010 + NC * NEA
      IF (I020.GT.DMAXCOR) THEN
         WRITE(6,*)'  OUT OF CORE ALGORITM NEEDED'
         CALL ERREX
      ENDIF
      CALL ZERO(EAMAT, NEA*NEA)
      CALL GETLST(ICORE(I000), 1,1,1, ISPIN, LS1OUT)
      if (prints) then
         write(6,*)' s1 part of EA matix'
         call output(icore(i000), 1, vrt(sirrep, ISPIN), 1, nea,
     $      vrt(sirrep, ISPIN), nea, 1)
      endif
C
C      A[C, ABI] = S1[C, ABI]   [ISIDE = 1]
C      A[ABI, C] = S1[C, ABI]   [ISIDE = 2]
C
C      IT IS MOST CONVENIENT TO CALCULATE A(TRANSPOSE) FOR ISIDE = 1
C      THEN ALL ROWS ARE FULL SIZE (NEA).
C
      CALL TRANSP(ICORE(I000), ICORE(I010), NEA, NC)
      CALL SAXPY(NEA*NC, ONE, ICORE(I010),1, EAMAT(1,1),1)
C  *******
C
      if (skip1) then
         write(6,*) ' higher original 2ph operators are skipped'
      else
      DO 40 MSPIN = 1, 1 + IUHF
         IF (IUHF.EQ.0) THEN
            ICOL = NC + 1
            ITYP = 19
         ELSE
            IF (ISPIN .EQ.  MSPIN) THEN
               ICOL = NC + 1
               ITYP = ISPIN
            ELSE
               ICOL = NC + LEN2 + 1
               ITYP = 15
            ENDIF
         ENDIF
         LISTS2EX = LS2OUT(ISPIN, MSPIN + 1 - IUHF)
         DO 50 XIRREP = 1, NIRREP
         IIRREP = DIRPRD(XIRREP, SIRREP)
         NCOL = POP(IIRREP,MSPIN) * IRPDPD(XIRREP, ITYP)
         IF (NCOL.GT.0) THEN
         NUMDSS = POP(IIRREP,MSPIN) * NEA
         I010 = I000 + NCOL * NEA
         I020 = I010 + NCOL * NEA
         IF (I020.GT.DMAXCOR) THEN
            WRITE(6,*)'  OUT OF CORE ALGORITM NEEDED'
            CALL ERREX
         ENDIF
         CALL GETLST(ICORE(I000), 1, NUMDSS, 1, XIRREP, LISTS2EX)
         CALL TRANSP(ICORE(I000), ICORE(I010), NEA, NCOL)
         CALL SAXPY(NEA*NCOL, ONE, ICORE(I010), 1, EAMAT(1, ICOL),1)
         ENDIF
         ICOL = ICOL + NCOL
 50   CONTINUE
 40   CONTINUE
      endif
C*******
      IF(ADDST) THEN
C
C  EXTEND THE MATRIX WITH SPIN-FLIP COMPONENTS
C
      NC = VRT(SIRREP, DSPIN)
      I010 = I000 + NC * NEA
      I020 = I010 + NC * NEA
      IF (I020.GT.DMAXCOR) THEN
         WRITE(6,*)'  OUT OF CORE ALGORITM NEEDED'
         CALL ERREX
      ENDIF
      CALL GETLST(ICORE(I000), 1,1,1, DSPIN, LS1OUT)
      if (prints) then
         write(6,*)' s1 part of EA matix', DSPIN
         call output(icore(i000), 1, vrt(sirrep, DSPIN), 1, nea,
     $      vrt(sirrep, DSPIN), nea, 1)
      endif
C
C      A[C, ABI] = S1[C, ABI]   [ISIDE = 1]
C      A[ABI, C] = S1[C, ABI]   [ISIDE = 2]
C
C      IT IS MOST CONVENIENT TO CALCULATE A(TRANSPOSE) FOR ISIDE = 1
C      THEN ALL ROWS ARE FULL SIZE (NEA).
C
      ICOL = NEA - NEAST + 1
      CALL TRANSP(ICORE(I000), ICORE(I010), NEA, NC)
      CALL SAXPY(NEA*NC, ONE, ICORE(I010),1, EAMAT(1,ICOL),1)
C  *****
      if (skip) then
         write(6,*) ' higher spin flip operators are not included'
      else
      NC = VRT(SIRREP, DSPIN) + NEA - NEAST
      DO 41 MSPIN = 1, 1 + IUHF
         IF (DSPIN .EQ.  MSPIN) THEN
            ICOL = NC + 1
            ITYP = DSPIN
         ELSE
            ICOL = NC + LEN2ST + 1
            ITYP = 15
         ENDIF
         LISTS2EX = LS2OUT(DSPIN, MSPIN + 1 - IUHF)
         DO 51 XIRREP = 1, NIRREP
         IIRREP = DIRPRD(XIRREP, SIRREP)
         NCOL = POP(IIRREP,MSPIN) * IRPDPD(XIRREP, ITYP)
         IF (NCOL.GT.0) THEN
         NUMDSS = POP(IIRREP,MSPIN) * NEA
         I010 = I000 + NCOL * NEA
         I020 = I010 + NCOL * NEA
         IF (I020.GT.DMAXCOR) THEN
            WRITE(6,*)'  OUT OF CORE ALGORITM NEEDED'
            CALL ERREX
         ENDIF
         CALL GETLST(ICORE(I000), 1, NUMDSS, 1, XIRREP, LISTS2EX)
         CALL TRANSP(ICORE(I000), ICORE(I010), NEA, NCOL)
         CALL SAXPY(NEA*NCOL, ONE, ICORE(I010), 1, EAMAT(1, ICOL),1)
         ENDIF
         ICOL = ICOL + NCOL
 51   CONTINUE
 41   CONTINUE
      endif
C******
      ENDIF
C
C    NOW TRANSPOSE EAMAT IF ISIDE = 1
C
      IF (ISIDE.EQ.1) THEN
         CALL SCOPY(NEA*NEA,EAMAT, 1, ICORE(I000), 1)
         CALL TRANSP(ICORE(I000), EAMAT, NEA, NEA)
      ENDIF
      printea = .false.
      if (printea) then
         write(6,*)' directly constructed ea matrix'
         call output(eamat, 1, nea, 1, nea, nea, nea, 1)
      endif
C
C  CALCULATE EIGENVALUES AND EIGENVECTORS
C
      diag =.true.
      if (diag) then
      EIGVEC = 0
      I000 = 1
      I010 = I000 + NEA*NEA
      CALL SCOPY(NEA*NEA, EAMAT(1,1),1, ICORE(I000), 1)
      CALL DIAGEA(ICORE(I000), NEA, ICORE(I010), MAXCOR, EIGVEC)
      endif
C
      RETURN
      END
