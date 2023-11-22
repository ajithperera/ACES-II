      SUBROUTINE PRINTS(ICORE, MAXCOR, IUHF, ISPIN, LS1, LS2)
C
C  ALL S-COEFFICIENTS CORREPONDING TO ISPIN ARE PRINTED OUT
C
      IMPLICIT INTEGER(A-Z)
      DIMENSION ICORE(MAXCOR), LS2(2,2)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SINFO/ NS(8), SIRREP
C
      WRITE(6,*) ' S1 COEFFICIENTS ', ISPIN, LS1
      CALL GETLST(ICORE, 1, 1, 1, ISPIN, LS1)
      CALL PRINTT1(ICORE, VRT(1, ISPIN), NS)
C
      NS1 = NS(SIRREP)
      DO 50 KSPIN = 1, 1+IUHF
         IF (IUHF. NE.0) THEN
            MSPIN = KSPIN
         ELSE
            MSPIN = 3-ISPIN
         ENDIF
         LISTS2EX = LS2(ISPIN, MSPIN)
         WRITE(6,*) ' S2 COEFFICIENTS ', ISPIN, MSPIN, LISTS2EX
         DO 100 XIRREP = 1 , NIRREP
            DISSYS = IRPDPD(XIRREP, ISYTYP(1, LISTS2EX))
            MIRREP = DIRPRD(SIRREP, XIRREP)
            NUMDSS = POP(MIRREP, MSPIN) * NS1
            IF ((NUMDSS * DISSYS) .GT. 0) THEN
               CALL GETLST(ICORE, 1, NUMDSS, 1, XIRREP, LISTS2EX)
               WRITE(6,*) ' IRREP IS ', XIRREP
               CALL OUTPUT(ICORE,1,DISSYS,1,NUMDSS,DISSYS,NUMDSS,1)
            ENDIF
 100     CONTINUE
 50   CONTINUE
C
      RETURN
      END