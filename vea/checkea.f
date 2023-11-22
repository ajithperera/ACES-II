C
C ********************************************************
C  THE FOLLOWING PROCEDURES CHECK THE CODE
C  BY COMPARING THE EA-EOM MATRICES CONSTRUCTED FROM 
C  LEFT AND RIGHT HAND MULTIPLICATIONS RESPECTIVELY.
C ********************************************************
C
      SUBROUTINE CHECKEA(ICORE, MAXCOR, IUHF)
C
C THIS ROUTINE CHECKS THE DIFFERENT PROCEDURES TO CONSTRUCT
C THE FULL EA-EOM MATRICES.
C
      IMPLICIT INTEGER (A-Z)
       logical lefthand, excicore, singonly, DROPCORE, ADDST
C
      DIMENSION ICORE(MAXCOR)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SINFO/ NS(8), SIRREP
      COMMON/STINFO/ITOTALS, STMODE, LISTST, STCALC, NSIZEST
      COMMON/EACALC/LEFTHAND, EXCICORE, SINGONLY, DROPCORE
      COMMON/EAINFO/NUMROOT(8,3)
      COMMON/COREINFO/IREPCORE, SPINCORE, IORBCORE, IORBOCC
C
      STCALC = 1
       DO 100 IRREP = 1, NIRREP
          SIRREP = IRREP
         DO 50 ISPIN = 1, 1 + IUHF
C
            ADDST = .FALSE.
            IF (EXCICORE) THEN
               ADDST = ((STMODE .EQ. 1) .AND. (ISPIN .EQ.SPINCORE))
            ENDIF
            ISTART = 1
            IEND = 1
            IF ((EXCICORE) .AND. (ISPIN .EQ. SPINCORE) .AND.
     $         (STMODE .EQ. 2)) THEN
               ISTART = 0
            ENDIF
            DO 25 IS = ISTART, IEND
                  ITOTALS = IS
            IF (NUMROOT(IRREP,ISPIN+1-IS) .GT. 0) THEN                  
C
C FIRST FIND THE DIMENSION OF THE MATRICES
C
         CALL IZERO(NS, 8)
         NS(IRREP) = 1
         IF (ADDST) THEN
            WRITE(6,*) ' DIMENSIONS IN EAMAT: B, BB, AB, A, AA, BA '
         ELSE
            WRITE(6,*) ' DIMENSIONS IN EAMAT: VRT, AA, AB'
         ENDIF
         IF (ADDST) THEN
            DSPIN = 3 - ISPIN
            MSPIN = 3 - DSPIN
            LENA = VRT(SIRREP, DSPIN)
            CALL GETLEN2(LENAA, IRPDPD(1, DSPIN), POP(1,DSPIN), NS(1))
            CALL GETLEN2(LENAB, IRPDPD(1,13), POP(1,MSPIN), NS(1))
            MSPIN = 3 - ISPIN
            LENB = VRT(SIRREP, ISPIN)
            CALL GETLEN2(LENBB, IRPDPD(1, ISPIN), POP(1,ISPIN), NS(1))
            CALL GETLEN2(LENBA, IRPDPD(1,13), POP(1,MSPIN), NS(1))
            NEA = LENA + LENAA + LENAB + LENB + LENBB + LENBA
            NEAST = LENA + LENAA + LENAB
            WRITE(6,*) LENB, LENBB, LENBA, LENA, LENAA, LENAB
         ELSE
            IF (IUHF.EQ.0) THEN
               CALL GETLEN(LENAB, VRT(1,1), VRT(1,1), POP(1,1), NS(1))
               LENAA = 0
            ELSE
               MSPIN = 3 - ISPIN
            CALL GETLEN2(LENAA, IRPDPD(1, ISPIN), POP(1,ISPIN), NS(1))
               CALL GETLEN2(LENAB, IRPDPD(1,13), POP(1,MSPIN), NS(1))
            ENDIF
            WRITE(6,*) VRT(IRREP, ISPIN), LENAA, LENAB
            NEA = VRT(IRREP, ISPIN) + LENAA + LENAB
            NEAST = 0
         ENDIF
         CALL GETNEA(NSIZEST,IRREP,3-ISPIN,IUHF)
         NS(SIRREP) = NEA
         I000 = 1
         I010 = I000 
         I020 = I010 + NEA*NEA*IINTFP
         I030 = I020 + NEA*NEA*IINTFP
         NEWCOR = MAXCOR - I030 + 1
         IF (NEWCOR.LT.0) THEN
            WRITE(6,*) 'CHECKEA REQUIRES TO MUCH MEMORY'
            CALL ERREX
         ENDIF
         write(6,1001) nea
 1001    format(' dimension of current EA matrix : ', i10)
C
         IF (NEA.EQ.0) GO TO 100
C
C  CONSTRUCT NEW LIST POSITIONS        
C
         CALL NEWLIST(SIRREP, NEA, IUHF, ISPIN,1)
C
C         CONSTRUCT MATRIX IN STRAIGHTFORWARD WAY
C
         write(6,*) '**********************'
         write(6,*) ' irrep        :  ', sirrep
         write(6,*) ' spin         :  ', ispin
         if (excicore) then
            write(6,*) ' multiplicity :  ', 2*itotals + 1
         endif
         write(6,*) '**********************'
C
C  CONSTRUCT THE MATRIX BY RIGHT MULTIPLYING THE IDENTITY
C
         ISIDE = 1
         DNEWCOR = NEWCOR/IINTFP
         CALL BUILTDIR(ICORE(I010), NEA,ICORE(I030),
     $     DNEWCOR, IUHF, ISIDE, ISPIN, NEAST)
C
C  ALTERNATIVELY CONSTRUCT THE MATRIX BY LEFT MULTIPLICATION
C
         if (.not. lefthand) then
            write(6,*) ' left hand matrix not calculated'
         else
         ISIDE = 2
         WRITE(6,*) ' LEFT HAND PROBLEM'
         CALL BUILTDIR(ICORE(I020), NEA,ICORE(I030),
     $     DNEWCOR, IUHF, ISIDE, ISPIN, neast)
            WRITE(6,*) ' EQUALITY OF LEFT AND RIGHT MULTIPLICATION'
            CALL COMPAREM(ICORE(I010), ICORE(I020), NEA,
     $         IUHF, ISPIN)
         endif
C
      ENDIF
 25   CONTINUE
 50   CONTINUE
 100  CONTINUE
      RETURN
      END
