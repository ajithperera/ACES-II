      SUBROUTINE PRMAINX(IUHF, SCR, DMAXCOR, IRREPX, IROOT)
C
C  COEFFICIENTS GREATER THAN THRESH ARE PRINTED AND 
C  CHARACTERIZED THROUGH THEIR ORBITAL INDICES
C
C  THE LABELS OF A COEFFICIENT C(IjAb) ARE PRINTED AS
C  I [IIRREP]  A [AIRREP] ; J [JIRREP]  B [BIRREP]:   COEFFICIENT 
C
C  THE COEFFICIENTS ARE ASSUMED TO BE GIVEN ON LS2OUT AND LS1OUT
C
C
C                           all; each
C LIST 444:    C(IJ,AB )       A<B ; I<J     AA AA
C      445:    C(ij,ab )       a<b ; i<j     BB BB
C      446:    C(Ij,Ab )       A,b ; I,j     AB A
C
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION SCR(DMAXCOR),THRESH,EIGVAL,OSCSTR,EIGVAL_T
      LOGICAL ORDERED, FIRST

      Character*1 Nature(100,8)
      CHARACTER*12 STRINGI(2), STRINGJ(2), STRINGA(2),STRINGB(2)

      DIMENSION LS2OUT(2,2)
      
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/ROOTS/EIGVAL(100,8),EIGVAL_T(100,8),OSCSTR(100,8),
     &             BGN(100,8),BGN_IRP(100,8),END(100,8),
     &             END_IRP(100,8),NATURE
C
      LS1OUT = 490
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
c
      THRESH = 0.05D0
C
C  FIRST CONSIDER SINGLE EXCITATION COEFFICIENTS
C
      DO 500 SSPIN = 1, 1+IUHF
        WRITE(6,1300)
        IF (SSPIN .EQ. 1) THEN
          WRITE(6,*)'   SINGLE EXCITATION COEFFICIENTS AA'
        ELSE
          WRITE(6,*)'   SINGLE EXCITATION COEFFICIENTS BB'
        ENDIF
        FIRST =.TRUE.
        CALL GETLST(SCR, 1, 1, 1, SSPIN, LS1OUT)
        ICOUNT = 1
        DO IIRREP = 1, NIRREP
          AIRREP=DIRPRD(IIRREP,IRREPX)
          DO 1 I = 1, POP(IIRREP, SSPIN)
            DO 2 A = 1, VRT(AIRREP,SSPIN)
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
    2       CONTINUE
    1     CONTINUE
        ENDDO
 1001   FORMAT(3X,I4,3X, ' [',I1,'],  ',1X, I4,3X, ' [',I1,']   ;',
     $       10X, F12.6)
C
  500 CONTINUE
      If (IUHF .NE.0) THEN
         CALL ASSIGN_STATES_UHF(SCR,NATURE,BGN,END,BGN_IRP,END_IRP,
     &                          IUHF,IRREPX,IROOT)
      Else
         CALL ASSIGN_STATES_RHF(SCR,NATURE,BGN,END,BGN_IRP,END_IRP,
     &                          IUHF,IRREPX,IROOT)
      Endif 
C
C  PRINT OUT DOUBLE EXCITATION COEFFICIENTS
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
      FIRST = .TRUE.
      LISTS2EX = LS2OUT(ASPIN, BSPIN)
      DO 10 RIRREP = 1, NIRREP
        LIRREP = DIRPRD(RIRREP,IRREPX)
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
        IF (DISSYA * NUMDSA .GT. 0) THEN
          I000 = 1
          I010 = I000 + DISSYS*NUMDSS *IINTFP
          IF (ICASE .EQ. 1) THEN
            CALL GETLST(SCR(I000), 1, NUMDSS, 1, RIRREP,
     $       LISTS2EX)
          ELSE
            CALL GETLST(SCR(I000), 1, NUMDSA, 1, RIRREP,
     $       LISTS2EX)            
            CALL SYMEXP(RIRREP, POP(1,ASPIN),DISSYA,SCR(I000))
            CALL SYMEXP2(LIRREP,VRT(1,ASPIN),DISSYS, DISSYA,
     &         NUMDSS, SCR(I000), SCR(I000))
          ENDIF
          ICOUNT = I000
          DO 20 JIRREP = 1, NIRREP
            IIRREP = DIRPRD(JIRREP, RIRREP)
            DO 30 J= 1, POP(JIRREP,BSPIN)
              DO 40 I = 1, POP(IIRREP, ASPIN)
                DO 50 BIRREP = 1, NIRREP
                  AIRREP = DIRPRD(LIRREP, BIRREP)
                  DO 60 B = 1, VRT(BIRREP, BSPIN)
                    DO 70 A = 1, VRT(AIRREP, ASPIN)
                      IF(ABS(SCR(ICOUNT)).GT. THRESH) THEN
                  IF (ORDERED(I,J,A,B,IIRREP,JIRREP,AIRREP,BIRREP,
     &                     ICASE)) THEN
                        IF (FIRST) THEN
                          WRITE(6,1300)
                          IF (ICASE .EQ. 2) THEN
                     WRITE(6,*)'   AA DOUBLE EXCITATION COEFFICIENTS '
                          ELSEIF(ICASE .EQ. 3) THEN
                     WRITE(6,*)'   BB DOUBLE EXCITATION COEFFICIENTS '
                          ELSEIF(ICASE .EQ. 1) THEN
                     WRITE(6,*)'   AB DOUBLE EXCITATION COEFFICIENTS '
                          ENDIF
                          WRITE(6,998) STRINGI(ASPIN),STRINGJ(BSPIN),
     &                       STRINGA(ASPIN), STRINGB(BSPIN)
  998                     FORMAT(6X,A12,3X,A12,3X,A12,3X,A12)
                          WRITE(6,*)
                        ENDIF
                        WRITE(6,1000) I, IIRREP,
     $                     J, JIRREP, A, AIRREP, B, BIRREP,SCR(ICOUNT)
                        FIRST = .FALSE.
                      ENDIF
                    ENDIF
                    ICOUNT = ICOUNT + 1
   70             CONTINUE
   60           CONTINUE
   50         CONTINUE
   40       CONTINUE
   30     CONTINUE
   20   CONTINUE
        IF (ICOUNT .NE. I000 + NUMDSS*DISSYS) THEN
          WRITE(6,*)' SOMETHING WRONG IN PRMAINX', ICASE,
     $       ICOUNT, I000+NUMDSS*DISSYS
        ENDIF
      ENDIF
   10 CONTINUE
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
                






