C
      SUBROUTINE BLTPROJD(IUHF, SCR, DMAXCOR, IRREPX, IRREP, INDEX,
     &                    ICASE, IOFFSET)
C
C This routine bulild the double excitation projection vectors.
C                        
C LIST 444:    C(IJ,AB )       A<B ; I<J     AA AA
C      445:    C(ij,ab )       a<b ; i<j     BB BB
C      446:    C(Ij,Ab )       A,b ; I,j     AB AB
C
      IMPLICIT INTEGER (A - Z)
      DOUBLE PRECISION SCR(DMAXCOR), FACTOR
      CHARACTER*4 PHTYPE
      DIMENSION LS2OUT(2, 2),INDEX(2), IOFFSET(8, 2)
      LOGICAL MATCHS, MATCHD
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
C
      LS2OUT(1,1) = 444
      LS2OUT(1,2) = 446
      LS2OUT(2,1) = 446
      LS2OUT(2,2) = 445
C
      FACTOR = 1.0D+00
C
C Consider double excitation coefficients
C
      IF (ICASE .EQ. 1) THEN
         ASPIN = 1
         BSPIN = 2
C         Write(6,*) "      A-B spin case"
C         Write(6,360)
      ELSEIF (ICASE .EQ. 2) THEN
         ASPIN = 1
         BSPIN = 1
C         Write(6,*) "      A-A spin case"
C         Write(6,360)
      ELSEIF (ICASE .EQ. 3) THEN
C         Write(6,*) "      B-B spin case"
C         Write(6,360)
         ASPIN = 2
         BSPIN = 2
      ENDIF
C
      LISTS2EX = LS2OUT(ASPIN, BSPIN)
C
      ICOUNT = 1
C
      DO 60 RIRREP = 1, NIRREP
C
         LIRREP = DIRPRD(RIRREP, IRREPX)
C
         IF (RIRREP .EQ. 1) THEN
            IOFFSET(1, ICASE) = 0
         ELSE
            IOFFSET(RIRREP, ICASE) = ICOUNT - 1 
         ENDIF
C
         IF (ICASE .EQ. 1) THEN
            DISSYS = IRPDPD(LIRREP, ISYTYP(1, LISTS2EX))
            NUMDSS = IRPDPD(RIRREP, ISYTYP(2, LISTS2EX))
            DISSYA = DISSYS
            NUMDSA = NUMDSS
         ELSE
            DISSYA = IRPDPD(LIRREP, ISYTYP(1, LISTS2EX))
            NUMDSA = IRPDPD(RIRREP, ISYTYP(2, LISTS2EX))
            DISSYS = IRPDPD(LIRREP, 18 + ASPIN)
            NUMDSS = IRPDPD(RIRREP, 20 + ASPIN)
         ENDIF
C
         IF (DISSYA*NUMDSA .GT. 0) THEN
            I000 = 1
            I010 = I000 + DISSYS*NUMDSS*IINTFP
C
            DO 70 JIRREP = 1, NIRREP
               IIRREP = DIRPRD(JIRREP, RIRREP)
               DO 80 J= 1, POP(JIRREP, BSPIN)
                  DO 90 I = 1, POP(IIRREP, ASPIN)
                     DO 100 BIRREP = 1, NIRREP
                        AIRREP = DIRPRD(LIRREP, BIRREP)
                        DO 110 B = 1, VRT(BIRREP, BSPIN)
                           DO 120 A = 1, VRT(AIRREP, ASPIN)
C     
                              IF (MATCHD(I, J, A, B, IIRREP, 
     &                                   JIRREP, ASPIN, BSPIN, 
     &                                   INDEX, IRREP)) THEN
                                 SCR(ICOUNT) = SCR(ICOUNT) + FACTOR
                              ENDIF
C
C      Write(6,350) I, J, A, B, IIRREP, JIRREP, AIRREP, BIRREP, IRREP,
C     &             INDEX(ASPIN), INDEX(BSPIN), ICOUNT, SCR(ICOUNT)
C
                              ICOUNT = ICOUNT + 1
C
 120                       CONTINUE
 110                    CONTINUE
 100                 CONTINUE
 90               CONTINUE
 80            CONTINUE
 70         CONTINUE
C
         ENDIF
C
 60   CONTINUE
C
C 350  FORMAT(2X,I2,1X,I2,1X,I2,1X,I2,3X,I1,6X,I1,7X,I1,7X,I1,5X,I1,
C     &       7X,I1,7X,I1,5X,I2,5X,F6.2)
C 360  FORMAT(3X,'I',2X,'J',2X,'A',2X,'B',2X,'I-Irr',2X,'J-Irr',2x,
C     &       'A-Irr',2X,'B-Irr',2X,'Irrep',2X,'Indx1',2X,'Indx2',2X,
C     &       'Icount',2X, 'Scr')
C
      RETURN
      END
