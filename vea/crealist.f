
      SUBROUTINE CREALIST(IUHF,ITYP)
C
C  CREATES LISTS USED IN EA-EOM CALCULATIONS.
C
      IMPLICIT INTEGER (A-Z)
      DIMENSION NS(8), MDIM(2)
      LOGICAL PRINT, LEFTHAND, EXCICORE, SINGONLY, DROPCORE
C     
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMPOP/ IRPDPD(8,22), ISYTYP(2,500), ID(18)
      COMMON/EXTRAP/MAXEXP,NREDUCE,NTOL,NSIZEC
      COMMON/SLISTS/LS1IN, LS1OUT, LS2IN(2,2), LS2OUT(2,2)
      COMMON/LISTDIP/DIPOO(2), DIPOV(2), DIPVV(2)
      COMMON/LISTDENS/DENSOO(2), DENSOV(2), DENSVV(2)
      COMMON/LISTDAV/LISTC, LISTHC, LISTH0
      COMMON/LISTF/LISTFJA, LISTFAJ, LISTFAB, LISTFIJ,
     $   LISTF0, NDIMFLST(2,4)
      COMMON/LISTU/LISTUIJ, LISTUAB, LISTUAI
      COMMON/EACALC/LEFTHAND,EXCICORE, SINGONLY, DROPCORE
      COMMON/STINFO/ITOTALS, STMODE, LISTST, STCALC, NSIZEST
C
C  IF ITYP = 1 or ITYP = 3 THEN
C  THE MAXIMUM SIZES OF THE LISTS ARE DETERMINED BY PUTTING NS TO THE 
C  MAXIMUM DIMENSION OF EAMAT, AS THIS IS DONE FIRST.
C
C  IF ITYP = 2
C  THE MAXIMUM DIMENSION OF THE LISTS ARE DETERMINED BY PUTTING NS TO 1
C
C  IF ITYP = 2 OR ITYP = 3
C  IN ADDITION THE LISTDAV, LISTDIP AND LISTDENS ARE INITIALIZED.
C
C  IF EXCICORE: LISTFJA, LISTFAJ, LISTFIJ AND LISTFAB ARE INITIALIZED
C  (USED TO ACCOUNT FOR PROPER SINGLET-TRIPLET SPLITTING)
C  IN THE SAME WAY LISTUIJ, LISTUAB AND LISTUAI ARE CREATED
C
C  CREATE S1 LISTS (MAXIMUM VALUES OF THE POINTERS)
C 
      MAXDIM = 0
      DO 5 ISPIN = 1, 1+IUHF
      DO 10 IRREP = 1, NIRREP
         CALL IZERO(NS, 8)
         NS(IRREP) = 1
         IF ((ITYP .EQ. 1) .OR. (ITYP .EQ. 3)) THEN
            LEN = 0
            DO 12 KSPIN = 1, 1 + IUHF
               CALL GETNEA(NEA, IRREP, KSPIN, IUHF)
               LEN = MAX(LEN, NEA)
 12         CONTINUE
         ELSEIF (ITYP.EQ.2) THEN
            LEN = 1
         ENDIF
         SIZES1 = VRT(IRREP, ISPIN) * LEN
         MAXDIM = MAX(SIZES1, MAXDIM)
 10   CONTINUE
      MDIM(ISPIN) = MAXDIM
 5    CONTINUE
c      INEWFIL = 1
      call aces_io_remove(54,'DERGAM')
      INEWFIL = 0
C
C  THE VERY FIRST TIME LIST 4?? IS ADDRESSED, INEWFIL = 1 OTHERWHISE 0
C
      CALL UPDMOI(1, MDIM(1), 1, LS1IN, INEWFIL, 0)
      INEWFIL = 0
      IF (IUHF.NE.0) THEN
         CALL UPDMOI(1, MDIM(2), 2, LS1IN, INEWFIL, 0)
      ENDIF
      DO 7 ISPIN = 1, 1+IUHF
         CALL UPDMOI(1, MDIM(ISPIN), ISPIN, LS1OUT, INEWFIL, 0)
 7    CONTINUE
C
C  NOW CREATE S2 LISTS
C
      DO 15 ISPIN = 1, 1 + IUHF
C
      DO 25 MSPIN = 1, 1 + IUHF
         IF (IUHF .EQ. 0 ) THEN
            DISTYP = 19
         ELSE
            IF (ISPIN .EQ.  MSPIN) THEN
               DISTYP = ISPIN
            ELSE
               DISTYP = 15
            ENDIF
         ENDIF
C
C  DISTYP DENOTES AB LABELS, RIGHT HAND LABELS DO NOT EXIST.
C
      ISYTYP(1, LS2IN(ISPIN, MSPIN)) = DISTYP
      ISYTYP(1, LS2OUT(ISPIN, MSPIN)) = DISTYP
C
C  FIRST DETERMINE MAXIMUM NUMBER OF S-COEFFICIENTS AND CORRESPONDING IRREP
C  FOR THE SPIN VARIANT UNDER CONSIDERATION.
C
      MAXIRREP = 0
      MAXDIM = 0
C
      DO 20 IRREP = 1, NIRREP
C
C  DETERMINE MAX NEA OVER THE TWO SPIN VARIANTS (SING-TRIP SPLITTING)
C
         IF ((ITYP .EQ. 1) .OR. (ITYP .EQ. 3)) THEN
            CALL IZERO(NS, 8)
            NS(IRREP) = 1
            LEN = 0
            DO 22 KSPIN = 1, 1 + IUHF
               CALL GETNEA(NEA, IRREP, KSPIN, IUHF)
               LEN = MAX(LEN, NEA)
 22         CONTINUE
         ELSEIF (ITYP .EQ. 2) THEN
            NEA = 1
         ENDIF
C
C  NOW DETERMINE NUMBER OF S-COEFFICIENTS
C
      CALL IZERO(NS, 8)
      NS(IRREP) = NEA
      NUMS = 0
      DO 30 XIRREP = 1, NIRREP
         MIRREP = DIRPRD(IRREP, XIRREP)
         NUMDSS = POP(MIRREP, MSPIN)*NEA
         DISSYS = IRPDPD(XIRREP, DISTYP)
         NUMS = NUMS + NUMDSS * DISSYS
 30   CONTINUE
      IF (NUMS .GT. MAXDIM) THEN
         MAXDIM = NUMS
         MAXIRREP = IRREP
      ENDIF
 20   CONTINUE
C
C NOW CREATE LISTPOINTERS CORRESPONDING TO THE MAXIMUM DIMENSION
C
      IRREP = MAXIRREP
C
C  DETERMINE NEA
C
      IF ((ITYP .EQ. 1) .OR. (ITYP .EQ. 3)) THEN
         CALL IZERO(NS, 8)
         NS(IRREP) = 1
            LEN = 0
            DO 32 KSPIN = 1, 1 + IUHF
               CALL GETNEA(NEA, IRREP, KSPIN, IUHF)
               LEN = MAX(LEN, NEA)
 32         CONTINUE
         ELSEIF (ITYP.EQ.2) THEN
            NEA = 1
      ENDIF
      DO 35 XIRREP = 1, NIRREP
         MIRREP = DIRPRD(MAXIRREP, XIRREP)
         NUMDSS = POP(MIRREP, MSPIN)*NEA
         DISSYS = IRPDPD(XIRREP, DISTYP)
         CALL UPDMOI(NUMDSS, DISSYS, XIRREP,
     $      LS2IN(ISPIN, MSPIN+1-IUHF), 0, 0)
 35   CONTINUE
      DO 40 XIRREP = 1, NIRREP
         MIRREP = DIRPRD(MAXIRREP, XIRREP)
         NUMDSS = POP(MIRREP, MSPIN)*NEA
         DISSYS = IRPDPD(XIRREP, DISTYP)
         CALL UPDMOI(NUMDSS, DISSYS, XIRREP,
     $      LS2OUT(ISPIN, MSPIN+1-IUHF), 0, 0)
 40   CONTINUE
 25   CONTINUE
 15   CONTINUE
C
      IF ((ITYP .EQ. 2) .OR. (ITYP .EQ. 3)) THEN
C
C   INITIALIZE LISTDAV
C
C  FIRST DETERMINE MAXIMUM DIMENSION
C
      MAXDIM = 0
C
      DO 45 ISPIN = 1, 1 + IUHF
C
      DO 50 IRREP = 1, NIRREP
C
C  DETERMINE NEA
C
         CALL GETNEA(NEA, IRREP, ISPIN, IUHF)
         IF (NEA .GT. MAXDIM) THEN
            MAXDIM = NEA
         ENDIF
 50   CONTINUE
 45   CONTINUE
C
      IF (LEFTHAND) THEN
         NUMLST = 2
      ELSE
         NUMLST = 1
      ENDIF
C
      DO 60 I=1, NUMLST
         CALL UPDMOI(MAXEXP, MAXDIM, I, LISTC, 0, 0)
         CALL UPDMOI(MAXEXP, MAXDIM, I, LISTHC, 0, 0)
 60   CONTINUE
         CALL UPDMOI(1+NREDUCE+8, MAXDIM, 1, LISTH0, 0, 0)
         CALL UPDMOI(1,MAXDIM,1,LISTST,0,0)
C
C  DETERMINE POINTERS TO DIPLIST AND DENSLIST
C
      DO 70 ISPIN = 1, 1+IUHF
         ISYTYP(1, DIPOO(ISPIN)) = 20 + ISPIN
         ISYTYP(1, DIPOV(ISPIN)) = 15 + ISPIN
         ISYTYP(1, DIPVV(ISPIN)) = 18 + ISPIN         
         ISYTYP(1, DENSOO(ISPIN)) = 20 + ISPIN
         ISYTYP(1, DENSOV(ISPIN)) = 15 + ISPIN
         ISYTYP(1, DENSVV(ISPIN)) = 18 + ISPIN
 70   CONTINUE
C
      DO 80 ISPIN = 1, 1+IUHF
      DO 90 XIRREP = 1, NIRREP
         NUMDSS = 3
         DISSYS = IRPDPD(XIRREP, ISYTYP(1,DIPOO(ISPIN)))
         CALL UPDMOI(NUMDSS, DISSYS, XIRREP,
     $      DIPOO(ISPIN), 0, 0)
 90   CONTINUE
 80   CONTINUE
C
      DO 81 ISPIN = 1, 1+IUHF
      DO 91 XIRREP = 1, NIRREP
         NUMDSS = 3
         DISSYS = IRPDPD(XIRREP, ISYTYP(1,DIPOV(ISPIN)))
         CALL UPDMOI(NUMDSS, DISSYS, XIRREP,
     $      DIPOV(ISPIN), 0, 0)
 91   CONTINUE
 81   CONTINUE
C
      DO 82 ISPIN = 1, 1+IUHF
      DO 92 XIRREP = 1, NIRREP
         NUMDSS = 3
         DISSYS = IRPDPD(XIRREP, ISYTYP(1,DIPVV(ISPIN)))
         CALL UPDMOI(NUMDSS, DISSYS, XIRREP,
     $      DIPVV(ISPIN), 0, 0)
 92   CONTINUE
 82   CONTINUE
C
      DO 83 ISPIN = 1, 1+IUHF
      DO 93 XIRREP = 1, NIRREP
         NUMDSS = 2
         DISSYS = IRPDPD(XIRREP, ISYTYP(1,DENSOO(ISPIN)))
         CALL UPDMOI(NUMDSS, DISSYS, XIRREP,
     $      DENSOO(ISPIN), 0, 0)
 93   CONTINUE
 83   CONTINUE
C
      DO 84 ISPIN = 1, 1+IUHF
      DO 94 XIRREP = 1, NIRREP
         NUMDSS = 2
         DISSYS = IRPDPD(XIRREP, ISYTYP(1,DENSOV(ISPIN)))
         CALL UPDMOI(NUMDSS, DISSYS, XIRREP,
     $      DENSOV(ISPIN), 0, 0)
 94   CONTINUE
 84   CONTINUE
C
      DO 85 ISPIN = 1, 1+IUHF
      DO 95 XIRREP = 1, NIRREP
         NUMDSS = 2
         DISSYS = IRPDPD(XIRREP, ISYTYP(1,DENSVV(ISPIN)))
         CALL UPDMOI(NUMDSS, DISSYS, XIRREP,
     $      DENSVV(ISPIN), 0, 0)
 95   CONTINUE
 85   CONTINUE
C
      ENDIF
c
      IF (EXCICORE) THEN
C
C      WE HAVE TO USE A DIFFERENT STRUCTURE HERE THAN USUAL IN ACES2
C      BECAUSE DIFFERENT COLUMS ON THE LIST CORREPOND TO DIFFERENT
C      SPIN VARIANTS.
C
         DO 410 I = 1, 2
            NDIMFLST(I,LISTFAJ-LISTF0) = IRPDPD(1, 10+I)
            NUMDSS = 1
            DISSYS = NDIMFLST(I,LISTFAJ-LISTF0)
            CALL UPDMOI(NUMDSS, DISSYS, I,LISTFAJ, 0, 0)
 410     CONTINUE
c
         DO 420 I = 1, 2
            NDIMFLST(I,LISTFJA-LISTF0) = IRPDPD(1, 10+3-I)
            NUMDSS = 1
            DISSYS = NDIMFLST(I,LISTFJA-LISTF0)
            CALL UPDMOI(NUMDSS, DISSYS, I,LISTFJA, 0, 0)
 420     CONTINUE
C
         DO 430 I = 1, 2
            NDIMFLST(I,LISTFAB-LISTF0) = IRPDPD(1, 15)
            NUMDSS = 1
            DISSYS = NDIMFLST(I,LISTFAB-LISTF0)
            CALL UPDMOI(NUMDSS, DISSYS, I,LISTFAB, 0, 0)
 430     CONTINUE
C
         DO 435 I = 1, 2
            NDIMFLST(I,LISTFIJ-LISTF0) = IRPDPD(1, 14)
            NUMDSS = 1
            DISSYS = NDIMFLST(I,LISTFIJ-LISTF0)
            CALL UPDMOI(NUMDSS, DISSYS, I,LISTFIJ, 0, 0)
 435     CONTINUE
C
C  CONSTRUCT U LISTS (SAME AS LISTS 91 - 93 FOR FOCK-MATRICES)      
C
      DO 440 ISPIN = 1, 2
         DISSYS = IRPDPD(1, 20 + ISPIN)
         CALL UPDMOI(1, DISSYS, ISPIN, LISTUIJ,0,0)
 440  CONTINUE
C
      DO 450 ISPIN = 1, 2
         DISSYS = IRPDPD(1, 18 + ISPIN)
         CALL UPDMOI(1, DISSYS, ISPIN, LISTUAB,0,0)
 450  CONTINUE
C
      DO 460 ISPIN = 1, 2
         DISSYS = IRPDPD(1, 8 + ISPIN)
         CALL UPDMOI(1, DISSYS, ISPIN, LISTUAI,0,0)
 460  CONTINUE
C
      ENDIF

c      CALL LISTINFO(IUHF,ITYP)

      RETURN
      END