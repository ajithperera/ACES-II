      SUBROUTINE GET_TOPAMPS(SCR,DMAXCOR,IMAP_A,IMAP_B,IUHF,NBAS, 
     &                       IRREPX,IUNIT,R1,R2)
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
      DOUBLE PRECISION R1,R2 
      DIMENSION IMAP_A(NBAS),IMAP_B(NBAS),NBFIRR(8)
      DIMENSION OCA_OFF(8),OCB_OFF(8),VTA_OFF(8),VTB_OFF(8)
      LOGICAL ORDERED, FIRST

      Character*1 Nature(100,8),BLANK
      CHARACTER*12 STRINGI(2), STRINGJ(2), STRINGA(2),STRINGB(2)
      CHARACTER*4 SPCASE(5)

      DIMENSION LS2OUT(2,2)
      
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/ROOTS/EIGVAL(100,8),EIGVAL_T(100,8),OSCSTR(100,8),
     &             BGN(100,8),BGN_IRP(100,8),END(100,8),
     &             END_IRP(100,8),NATURE
#include "info.com" 

      DATA SPCASE/'AA  ','BB  ','ABAB','AAAA','BBBB'/ 
      DATA BLANK/''/ 
C
      write(IUNIT,*)
      write(IUNIT,*)' Principal elements of right eigenvector: '
      write(iunit,*)' Orbitals given in HF energy ordering     '
      write(IUNIT,*)
      write(IUNIT,*)'------------------------------------------------'
      write(IUNIT,*)' CASE     I    J    A    B  C(ab,ij)  Spin Case '
      write(IUNIT,*)'------------------------------------------------'
c
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

      OCA_OFF(1) = 0
      OCB_OFF(1) = 0
      VTA_OFF(1) = NOCCO(1) 
      VTB_OFF(1) = NOCCO(2) 

      Call Getrec(20,"JOBARC",'NUMBASIR',Nirrep, Nbfirr)
      DO I = 1, (NIRREP-1) 
         OCA_OFF(I+1) = OCA_OFF(I) + POP(I,1)
         OCB_OFF(I+1) = OCB_OFF(I) + POP(I,2)
         VTA_OFF(I+1) = VTA_OFF(I) + VRT(I,1)
         VTB_OFF(I+1) = VTB_OFF(I) + VRT(I,2)
      ENDDO 
      R1 = 0.0D0
      R2 = 0.0D0
C
C  FIRST CONSIDER SINGLE EXCITATION COEFFICIENTS
C
      DO 500 SSPIN = 1, 1+IUHF

CSSS        WRITE(6,1300)
CSSS        IF (SSPIN .EQ. 1) THEN
CSSS          WRITE(6,*)'   SINGLE EXCITATION COEFFICIENTS AA'
CSSS        ELSE
CSSS          WRITE(6,*)'   SINGLE EXCITATION COEFFICIENTS BB'
CSSS        ENDIF

        CALL GETLST(SCR, 1, 1, 1, SSPIN, LS1OUT)
        ICOUNT = 1
        DO IIRREP = 1, NIRREP
          AIRREP=DIRPRD(IIRREP,IRREPX)
          DO 1 I = 1, POP(IIRREP, SSPIN)
            DO 2 A = 1, VRT(AIRREP,SSPIN)
             R1 = R1 + SCR(ICOUNT)*SCR(ICOUNT)
             IF (ABS(SCR(ICOUNT)).GT.THRESH) THEN
                 IF (SSPIN .EQ. 1) THEN
                    II = IMAP_A(I+OCA_OFF(IIRREP))  
                    AA = IMAP_A(A+VTA_OFF(AIRREP))
                 ELSE
                    II = IMAP_B(I+OCB_OFF(IIRREP)) 
                    AA = IMAP_B(A+VTB_OFF(AIRREP))
                 ENDIF 
                 WRITE(IUNIT,'(I4, 3x, 4I5,F15.10,1X,A4)') SSPIN,II,0,
     &                 AA,0,SCR(ICOUNT) ,SPCASE(SSPIN) 
              ENDIF
              ICOUNT = ICOUNT + 1
    2       CONTINUE
    1     CONTINUE
        ENDDO

  500 CONTINUE

      WRITE(IUNIT,"(A)") BLANK
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

          IF (I010 .GE. DMAXCOR) CALL INSMEM("@-GET_TOPAMPS",I010,
     &                                        DMAXCOR)
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
                      R2 = R2 + SCR(ICOUNT)*SCR(ICOUNT)
                      IF(ABS(SCR(ICOUNT)).GT. THRESH) THEN
                      IF (ORDERED(I,J,A,B,IIRREP,JIRREP,AIRREP,BIRREP,
     &                    ICASE)) THEN

                          IF (ICASE .EQ. 2) THEN

                             II = IMAP_A(I+OCA_OFF(IIRREP)) 
                             JJ = IMAP_A(J+OCA_OFF(JIRREP)) 
                             AA = IMAP_A(A+VTA_OFF(AIRREP))
                             BB = IMAP_A(B+VTA_OFF(BIRREP))
                             WRITE(IUNIT,'(I4, 3x, 4I5,F15.10,1X,A4)') 
     &                              1,II,JJ,AA,BB,SCR(ICOUNT) ,SPCASE(4) 

                          ELSEIF(ICASE .EQ. 3) THEN

                             II = IMAP_B(I+OCB_OFF(IIRREP))
                             JJ = IMAP_B(J+OCB_OFF(JIRREP))
                             AA = IMAP_B(A+VTB_OFF(AIRREP))
                             BB = IMAP_B(B+VTB_OFF(BIRREP)) 
                             WRITE(IUNIT,'(I4, 3x, 4I5,F15.10,1X,A4)') 
     &                              2,II,JJ,AA,BB,SCR(ICOUNT) ,SPCASE(5) 

                          ELSEIF(ICASE .EQ. 1) THEN

                             II = IMAP_A(I+OCA_OFF(IIRREP)) 
                             JJ = IMAP_B(J+OCB_OFF(JIRREP))
                             AA = IMAP_A(A+VTA_OFF(AIRREP))
                             BB = IMAP_B(B+VTB_OFF(BIRREP))
                              WRITE(IUNIT,'(I4, 3x, 4I5,F15.10,1X,A4)') 
     &                              3,II,JJ,AA,BB,SCR(ICOUNT) ,SPCASE(3) 
                          ENDIF 
                        ENDIF
                      ENDIF
                    ICOUNT = ICOUNT + 1
   70             CONTINUE
   60           CONTINUE
   50         CONTINUE
   40       CONTINUE
   30     CONTINUE
   20   CONTINUE
        IF (ICOUNT .NE. 1 + NUMDSS*DISSYS) THEN
          WRITE(6,*)' SOMETHING WRONG IN PRMAINX', ICASE,
     $       ICOUNT, 1+NUMDSS*DISSYS
        ENDIF
      ENDIF
   10 CONTINUE
      ENDDO
      
      write(IUNIT,*)'------------------------------------------------'
C
      RETURN
      END
                






