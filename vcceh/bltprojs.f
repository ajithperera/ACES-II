C
      SUBROUTINE BLTPROJS(IUHF, SCR, DMAXCOR, IRREPX, IRREP, INDEX)
C
      IMPLICIT INTEGER (A - Z)
      DOUBLE PRECISION SCR(DMAXCOR), FACTOR
      CHARACTER*4 PHTYPE
      DIMENSION LS2OUT(2, 2),INDEX(2)
      LOGICAL MATCHS, MATCHD
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
C
      LS1OUT = 490
C
      FACTOR = 1.0D+00
C
C First consider the single excitation coefficients.
C
      DO 10 SSPIN = 1, 1 + IUHF
C
C         IF (SSPIN .EQ. 1)THEN
C            Write(6,*)"       Alpha Spin    "
C            Write(6, 310)
C            Write(6, *)
C         ELSE
C            Write(6, *)"      Beta  Spin    "
C            Write(6, 310)
C            Write(6, *)
C         ENDIF
C
         CALL GETLST(SCR, 1, 1, 1, SSPIN, LS1OUT)
         ICOUNT = 1
C
         DO 20 IIRREP = 1, NIRREP
C
            AIRREP = DIRPRD(IIRREP, IRREPX)
C
            DO 30 I = 1, POP(IIRREP, SSPIN)
C
               DO 40 A = 1, VRT(AIRREP, SSPIN)
C
                  IF (MATCHS(I, IIRREP, SSPIN, INDEX, IRREP)) THEN
                     SCR(ICOUNT) = SCR(ICOUNT) + FACTOR
                  ENDIF
C     
C                  IF (SSPIN .EQ. 1) THEN
C              Write(6, 300) I, A, IRREP, IIRREP, Index(1), Scr(ICOUNT)
C                  ELSE
C              Write(6, 300) I, A, IRREP, IIRREP, Index(2), Scr(ICOUNT)
C                  ENDIF
C
                  ICOUNT = ICOUNT + 1
C
 40            CONTINUE
 30         CONTINUE
 20      CONTINUE
C
         CALL PUTLST(SCR, 1, 1, 1, SSPIN, LS1OUT)
C
 10   CONTINUE
C
 300  Format(5X, I2, 5X, I2, 4X, I1, 4X, I1, 8X, I1, 4X, F4.2) 
 310  Format(6X,'I',5X,'A',4X,'IRREP',2X,'IIRREP',2X,'Indx',3X,'Scr')
C
      RETURN
      END
