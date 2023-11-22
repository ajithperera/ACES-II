      SUBROUTINE LOCATE_HBARC(IUHF, SCR, DMAXCOR, IRREPX)
C
C LIST 444:    C(IJ,AB)       A<B ; I<J     AA AA
C      445:    C(ij,ab)       a<b ; i<j     BB BB
C      446:    C(Ij,Ab)       A,b ; I,j     AB AB
C      490:    C(i,a)
C
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION SCR(DMAXCOR), THRESH
      LOGICAL ORDERED, FIRST
      CHARACTER*12 STRINGI(2), STRINGJ(2), STRINGA(2),STRINGB(2)
      DIMENSION LS2OUT(2,2)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
C
      LS1OUT      = 490
      LS2OUT(1,1) = 444
      LS2OUT(1,2) = 446
      LS2OUT(2,1) = 446
      LS2OUT(2,2) = 445
C
      STRINGI(2)='i  [i_SYM]  '
      STRINGJ(2)='j  [j_SYM]  '
      STRINGA(2)='a  [a_SYM]  '
      STRINGB(2)='b  [b_SYM]  '
      STRINGI(1)='I  [I_SYM]  '
      STRINGJ(1)='J  [J_SYM]  '
      STRINGA(1)='A  [A_SYM]  '
      STRINGB(1)='B  [B_SYM]  '
C
C  FIRST CONSIDER SINGLE EXCITATION COEFFICIENTS
C
      DO SSPIN = 1, 1+IUHF
C
        WRITE(6,1300)
        IF (SSPIN .EQ. 1) THEN
          WRITE(6,*)'   SINGLE EXCITATION COEFFICIENTS AA'
        ELSE
          WRITE(6,*)'   SINGLE EXCITATION COEFFICIENTS BB'
        ENDIF
        FIRST =.TRUE.
        CALL GETLST(SCR, 1, 1, 1, SSPIN, LS1OUT)
        ICOUNT = 1
C
        DO IIRREP = 1, NIRREP
          AIRREP=DIRPRD(IIRREP,IRREPX)
C
          DO  I = 1, POP(IIRREP, SSPIN)
            DO  A = 1, VRT(AIRREP,SSPIN)
C
              IF (ABS(SCR(ICOUNT)).GT.THRESH) THEN
                IF (FIRST) THEN
                  WRITE(6,999) STRINGI(SSPIN), STRINGA(SSPIN)
  999             FORMAT(6X,A12,3X,A12)
                  WRITE(6,*)
                ENDIF
                FIRST =.FALSE.
                WRITE(6, 1001) I, IIRREP, A, AIRREP, SCR(ICOUNT)
              ENDIF
              ICOUNT = ICOUNT + 1
            ENDDO
          ENDDO
C
 1001   FORMAT(3X,I4,3X, ' [',I1,'],  ',1X, I4,3X, ' [',I1,']   ;',
     $       10X, F12.6)
        
C
        ENDDO
      ENDDO
C
      DO ICASE = 1, 1 + 2 * IUHF
         IF (ICASE .EQ. 1) THEN
            ASPIN = 1
            BSPIN = 2
         ELSEIF (ICASE .EQ. 2) THEN
            ASPIN = 1
            BSPIN = 1
         ELSEIF (ICASE .EQ. 3) THEN
            ASPIN = 2
           BSPIN = 2
         ENDIF
C
         FIRST = .TRUE.
         LISTS2EX = LS2OUT(ASPIN, BSPIN)
C
        DO RIRREP = 1, NIRREP
           LIRREP = DIRPRD(RIRREP,IRREPX)
C
           IF (ICASE .EQ. 1) THEN
              DISSYS = IRPDPD(LIRREP, ISYTYP(1,LISTS2EX))
              NUMDSS = IRPDPD(RIRREP, ISYTYP(2,LISTS2EX))
              DISSYA = DISSYS
             NUMDSA = NUMDSS
           ELSE
             DISSYA = IRPDPD(LIRREP, ISYTYP(1,LISTS2EX))
             NUMDSA = IRPDPD(RIRREP, ISYTYP(2,LISTS2EX))
             DISSYS = IRPDPD(LIRREP, 18 + ASPIN)
             NUMDSS = IRPDPD(RIRREP, 20 + ASPIN)
           ENDIF
C
           IF (DISSYA * NUMDSA .GT. 0) THEN
              I000 = 1
              I010 = I000 + DISSYS*NUMDSS *IINTFP
C
              IF (ICASE .EQ. 1) THEN
                 CALL GETLST(SCR(I000), 1, NUMDSS, 1, RIRREP,
     &                       LISTS2EX)
              ELSE
                 CALL GETLST(SCR(I000), 1, NUMDSA, 1, RIRREP,
     $                       LISTS2EX)            
                 CALL SYMEXP(RIRREP, POP(1,ASPIN),DISSYA,SCR(I000))
                 CALL SYMEXP2(LIRREP,VRT(1,ASPIN),DISSYS, DISSYA,
     &                        NUMDSS, SCR(I000), SCR(I000))
              ENDIF
C
              ICOUNT = I000
C
              DO JIRREP = 1, NIRREP
                 IIRREP = DIRPRD(JIRREP, RIRREP)
C
              DO J= 1, POP(JIRREP,BSPIN)
                 DO I = 1, POP(IIRREP, ASPIN)
C
                    DO BIRREP = 1, NIRREP
                       AIRREP = DIRPRD(LIRREP, BIRREP)
C
                       DO B = 1, VRT(BIRREP, BSPIN)
                          DO A = 1, VRT(AIRREP, ASPIN)
  
                          IF (ABS(SCR(ICOUNT)).GT. THRESH) THEN
                             IF (ORDERED(I,J,A,B,IIRREP,JIRREP,
     &                           AIRREP,BIRREP,ICASE)) THEN
C
                             IF (FIRST) THEN
                                 WRITE(6,1300)
                                 IF (ICASE .EQ. 2) THEN
                     WRITE(6,*)'   AA DOUBLE EXCITATION COEFFICIENTS '
                                 ELSEIF(ICASE .EQ. 3) THEN
                     WRITE(6,*)'   BB DOUBLE EXCITATION COEFFICIENTS '
                                 ELSEIF(ICASE .EQ. 1) THEN
                     WRITE(6,*)'   AB DOUBLE EXCITATION COEFFICIENTS '
                                 ENDIF
C
                                 WRITE(6,998) STRINGI(ASPIN),
     &                                        STRINGJ(BSPIN),
     &                                        STRINGA(ASPIN), 
     &                                        STRINGB(BSPIN)
  998                     FORMAT(6X,A12,3X,A12,3X,A12,3X,A12)
                          WRITE(6,*)
C
                             ENDIF
C
                             WRITE(6,1000) I, IIRREP, J, JIRREP, A, 
     &                                     AIRREP, B, BIRREP,
     &                                     SCR(ICOUNT)
                             FIRST = .FALSE.
                             ENDIF
                          ENDIF
C
                          ICOUNT = ICOUNT + 1
C
                          ENDDO 
                       ENDDO
C
                    ENDDO
                ENDDO
C
              ENDDO
            ENDDO
C
            IF (ICOUNT .NE. I000 + NUMDSS*DISSYS) THEN
                WRITE(6,*)' SOMETHING WRONG IN PRMAINX', ICASE,
     $          ICOUNT, I000+NUMDSS*DISSYS
            ENDIF
C
          ENDIF
C
        ENDDO
      ENDDO
C
      WRITE(6,1300)
C
 1300 FORMAT(/)
 1000 FORMAT(3X,I4,3X, ' [',I1,'],  ',1X, I4,3X, ' [',I1,']   ;',
     $   I4,3X, ' [',I1,']    ', I4,3X, ' [',I1,']  ', F12.6)
C
      RETURN
      END
                






