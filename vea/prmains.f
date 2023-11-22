      SUBROUTINE PRMAINS(IUHF, ISPIN, SCR, DMAXCOR,TYPE)
C
C  S-COEFFICIENTS GREATER THAN THRESH ARE PRINTED AND 
C  CHARACTERIZED THROUGH THEIR ORBITAL INDICES
C
C  THE LABELS OF A COEFFICIENT S(AbI) ARE PRINTED AS
C  A [AIRREP]  ; B [BIRREP]  ; I [IIRREP]  ;     COEFFICIENT   
C
C  THE COEFFICIENTS ARE ASSUMED TO BE GIVEN ON LS2OUT AND LS1OUT
C
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION SCR(DMAXCOR), THRESH
      LOGICAL ORDERED, FIRST
      CHARACTER*14 STRINGA(2), STRINGB(2), STRINGI(2)
      CHARACTER*2 TYPE
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SINFO/ NS(8), SIRREP
      COMMON/SLISTS/LS1IN, LS1OUT, LS2IN(2,2), LS2OUT(2,2)
C
      STRINGA(1)='A  [A_IRREP]  '
      STRINGB(1)='B  [B_IRREP]; '
      STRINGI(1)='I  [I_IRREP]  '
      STRINGA(2)='a  [a_IRREP]  '
      STRINGB(2)='b  [b_IRREP]; '
      STRINGI(2)='i  [i_IRREP]  '
C
      THRESH = 0.1D0
C
      IF (NS(SIRREP).GT. 1) THEN
         WRITE(6,*) ' PRINT MAIN EIGENVECTOR COMPONENTS :'
         WRITE(6,*) ' NS > 1 NOT IMPLEMENTED'
         RETURN
      ENDIF
C
C  FIRST CONSIDER SINGLE EA COEFFICIENTS
C
      FIRST =.TRUE.
      CALL GETLST(SCR, 1, 1, 1, ISPIN, LS1OUT)
      DO 1 A = 1, VRT(SIRREP, ISPIN)
         IF (ABS(SCR(A)).GT.THRESH) THEN
            IF (FIRST) THEN
               WRITE(6,1300)
	       IF (TYPE .EQ. "EA") WRITE(6,*)'   SINGLE EA COEFFIENTS '
               IF (TYPE .EQ. "IP") WRITE(6,*)'   SINGLE IP COEFFIENTS '
               WRITE(6,999) STRINGA(ISPIN)
 999           FORMAT(7X,A14)
               WRITE(6,*)
            ENDIF
            FIRST =.FALSE.
            WRITE(6, 1001) A, SIRREP, SCR(A)
         ENDIF
 1    CONTINUE
 1001 FORMAT(4X, I4,4X,' [',I1,'] ', ';', 18X, F12.6)
C
      DO 100 MSPIN = 1, 1+IUHF
         FIRST = .TRUE.
         LISTS2EX = LS2OUT(ISPIN, MSPIN + 1 - IUHF)
         DO 200 XIRREP = 1, NIRREP
            MIRREP = DIRPRD(SIRREP, XIRREP)
            NUMDSS = POP(MIRREP, MSPIN) * NS(SIRREP)
            IF (ISPIN.EQ.MSPIN) THEN
               DISSYS = IRPDPD(XIRREP, 18 + ISPIN)
               DISSYA = IRPDPD(XIRREP, ISPIN)
            ELSE
               DISSYS = IRPDPD(XIRREP, 13)
            ENDIF
            I000 = 1
            I010 = I000 + DISSYS*NUMDSS *IINTFP
            IF ((IUHF.NE.0) .AND. (ISPIN.EQ.MSPIN)) THEN
               CALL GETLST(SCR(I010), 1, NUMDSS, 1, XIRREP,
     $            LISTS2EX)
               CALL SYMEXP2(XIRREP, VRT(1, ISPIN), DISSYS,
     $            DISSYA, NUMDSS, SCR(I000), SCR(I010))
            ELSE
               CALL GETLST(SCR(I000), 1, NUMDSS, 1, XIRREP,
     $            LISTS2EX)
            ENDIF
            ICOUNT = I000
            DO 10 S= 1, NS(SIRREP)
               DO 20 M = 1, POP(MIRREP, MSPIN)
                  DO 30 BIRREP = 1, NIRREP
                     AIRREP = DIRPRD(XIRREP, BIRREP)
                     ASYMCASE = 1
                 IF ((IUHF. NE.0) .AND. (MSPIN.EQ.ISPIN)) THEN
                        IF (AIRREP .EQ. BIRREP) ASYMCASE = 2
                        IF (AIRREP .GT. BIRREP) ASYMCASE = 3
                     ENDIF
                     DO 40 B = 1, VRT(BIRREP, ISPIN)
                        DO 50 A = 1, VRT(AIRREP, MSPIN)
                           IF(ABS(SCR(ICOUNT)).GT. THRESH) THEN
                              IF (ORDERED(ASYMCASE, A, B)) THEN
                                 IF (FIRST) THEN
                                    WRITE(6,1300)
                           IF (ISPIN .EQ. MSPIN+1-IUHF) THEN
               WRITE(6,*)'   AA SPIN COEFFICIENTS '
            ELSE
               WRITE(6,*)'   AB SPIN COEFFICIENTS '
            ENDIF
            WRITE(6,998) STRINGA(MSPIN+1-IUHF), STRINGB(ISPIN),
     $         STRINGI(MSPIN+1-IUHF)
 998        FORMAT(7X,A14,2X,A14,2X,A14)
            WRITE(6,*)
                                 ENDIF
                                 WRITE(6,1000) A, AIRREP,
     $                              B, BIRREP, M, MIRREP, SCR(ICOUNT)
                                 FIRST = .FALSE.
                              ENDIF
                           ENDIF
                              ICOUNT = ICOUNT + 1
 50                        CONTINUE
 40                     CONTINUE
 30                  CONTINUE
 20               CONTINUE
 10            CONTINUE
               IF (ICOUNT .NE. I000 + NUMDSS*DISSYS) THEN
                  WRITE(6,*)' SOMETHING WRONG IN PRMAINS',
     $               ICOUNT, I000+NUMDSS*DISSYS
               ENDIF
C
 200        CONTINUE
 100     CONTINUE
C
 1300       FORMAT(/)
 1000       FORMAT(4X,I4,4X, ' [',I1,'],  ',1X, I4,4X, ' [',I1,']   ;',
     $         I4,4X, ' [',I1,']    ', F12.6)
C
            RETURN
            END
