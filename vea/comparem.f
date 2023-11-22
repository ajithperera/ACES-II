      SUBROUTINE COMPAREM(EAMDIR, EAMAT, NEA, IUHF, ISPIN)
C
C COMPARES THE MATRICES EAMDIR AND EAMAT IN THE VARIOUS BLOCKS
C THE LARGEST DISCREPANCY (IF NOT ZERO) IS PRINTED
C
      IMPLICIT INTEGER (A - Z)
      DOUBLE PRECISION EAMAT, EAMDIR, X, ZILCH, DEL, Y
      DIMENSION EAMAT(NEA, NEA), EAMDIR(NEA,NEA), BLOKOFF(8), NSZ(8)
      LOGICAL LEFTHAND, EXCICORE, SINGONLY, DROPCORE, ADDST,
     $   ALLZERO, print
C
      COMMON/EACALC/LEFTHAND, EXCICORE, SINGONLY, DROPCORE
      COMMON/COREINFO/IREPCORE, SPINCORE, IORBCORE, IORBOCC
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SINFO/ NS(8), SIRREP
      COMMON/STINFO/ITOTALS, STMODE, LISTST, STCALC, NSIZEST
C
      CALL IZERO(NSZ,8)
      NSZ(SIRREP) = 1
      DATA ZILCH /0.0D0/
      ADDST = .FALSE.
      IF (EXCICORE) THEN
         ADDST = ((STMODE.EQ.1) .AND. (ISPIN .EQ.SPINCORE))
      ENDIF
C
      IF (ADDST) THEN
         DSPIN = 3 - ISPIN
         MSPIN = 3 - DSPIN
         LENA = VRT(SIRREP, DSPIN)
         CALL GETLEN2(LENAA, IRPDPD(1, DSPIN), POP(1,DSPIN), NSZ(1))
         CALL GETLEN2(LENAB, IRPDPD(1,13), POP(1,MSPIN), NSZ(1))
         MSPIN = 3 - ISPIN
         LENB = VRT(SIRREP, ISPIN)
         CALL GETLEN2(LENBB, IRPDPD(1, ISPIN), POP(1,ISPIN), NSZ(1))
         CALL GETLEN2(LENBA, IRPDPD(1,13), POP(1,MSPIN), NSZ(1))
         BLOKOFF(1) = 1
         BLOKOFF(2) = BLOKOFF(1) + LENB
         BLOKOFF(3) = BLOKOFF(2) + LENBB
         BLOKOFF(4) = BLOKOFF(3) + LENBA
         BLOKOFF(5) = BLOKOFF(4) + LENA
         BLOKOFF(6) = BLOKOFF(5) + LENAA
         BLOKOFF(7) = BLOKOFF(6) + LENAB
         NBLOCK = 6
      ELSE
         IF (IUHF.EQ.0) THEN
            CALL GETLEN(LENAB, VRT(1,1), VRT(1,1), POP(1,1), NSZ(1))
            LENA = VRT(SIRREP, ISPIN)
            BLOKOFF(1) = 1
            BLOKOFF(2) = BLOKOFF(1) + LENA
            BLOKOFF(3) = BLOKOFF(2) + LENAB
            NBLOCK = 2
         ELSE
            MSPIN = 3 - ISPIN
            CALL GETLEN2(LENAA, IRPDPD(1, ISPIN), POP(1,ISPIN),
     $         NSZ(1))
            CALL GETLEN2(LENAB, IRPDPD(1,13), POP(1,MSPIN), NSZ(1))
            LENA = VRT(SIRREP, ISPIN)
            BLOKOFF(1) = 1
            BLOKOFF(2) = BLOKOFF(1) + LENA
            BLOKOFF(3) = BLOKOFF(2) + LENAA
            BLOKOFF(4) = BLOKOFF(3) + LENAB
            NBLOCK = 3
         ENDIF
      ENDIF
      print = .true.
      if (print) then
         do 33 iblock = 1, nblock + 1
            write(6,*) iblock, blokoff(iblock)
 33      continue
      endif
C
C CHECK ALL POSSIBLE BLOCKS
C
      ALLZERO = .TRUE.
      DO 10 IBLOCK = 1, NBLOCK
         DO 20 JBLOCK = 1, NBLOCK            
            X = ZILCH
            Y = ZILCH
            DO 30 I = BLOKOFF(IBLOCK), BLOKOFF(IBLOCK+1) - 1
               DO 40 J = BLOKOFF(JBLOCK), BLOKOFF(JBLOCK+1) - 1
               DEL = ABS(EAMAT(I,J) - EAMDIR(I,J))
               X = MAX(DEL, X)
               Y = MAX((ABS(EAMAT(I,J))+ABS(EAMDIR(I,J))), Y)
 40            CONTINUE
 30         CONTINUE
            IF (X .GT. 1.E-7) THEN
               WRITE(6,'(3I5,3X,3I5,F14.7)') IBLOCK, BLOKOFF(IBLOCK),
     $            BLOKOFF(IBLOCK+1) -1,
     $            JBLOCK, BLOKOFF(JBLOCK), BLOKOFF(JBLOCK+1) -1, X
               ALLZERO = .FALSE.
            ENDIF
            IF (Y .GT. 1.E-7) THEN
               print = .false.
               if (print) then
                  write(6,*) 'right hand block ',iblock, jblock
                  ilow = blokoff(iblock)
                  ihigh= blokoff(iblock+1) - 1
                  jlow = blokoff(jblock)
                  jhigh= blokoff(jblock+1) - 1
                  call output(eamdir, ilow, ihigh, jlow, jhigh,
     $               nea, nea, 1)
                  write(6,*) 'left hand block ', iblock, jblock
                  ilow = blokoff(iblock)
                  ihigh= blokoff(iblock+1) - 1
                  jlow = blokoff(jblock)
                  jhigh= blokoff(jblock+1) - 1
                  call output(eamat, ilow, ihigh, jlow, jhigh,
     $               nea, nea, 1)
               endif
            ENDIF
 20      CONTINUE
 10   CONTINUE
      IF (ALLZERO) THEN
      WRITE(6,*) ' ALL BLOCKS ARE EQUAL TO ZERO'
      ENDIF
C
      RETURN
      END
