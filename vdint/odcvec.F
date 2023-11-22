      SUBROUTINE ODCVEC(NUCA,NUCB,NCONTA,NCONTB,NSETA,NSETB,
     *                  NPRIMA,NPRIMB,NRCOA,NRCOB,JSTRA,JSTRB,
     *                  NUCAB,TPRIAB,
     *                  SIGNAX,SIGNAY,SIGNAZ,
     *                  SIGNBX,SIGNBY,SIGNBZ,
     *                  THRESH,ITYPE,IPRINT)
C
CEND
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (LUCMD = 5, LUPRI = 6)
#include "baslims.par"
      PARAMETER (ZERO = 0.00 D00, HALF = 0.50 D00, ONE = 1.00 D00)
      PARAMETER (PI     = 3.14159 26535 89793 D00,
     *           R2PI52 = 5.91496 71727 95612 D00)
      LOGICAL TPRIAB
      COMMON /PRIMIT/ PRIEXP(MXPRIM), PRICCF(MXPRIM,MXCONT),
     *                PRICRX(MXPRIM), PRICRY(MXPRIM), PRICRZ(MXPRIM)
      COMMON /ODCCOM/ COORPX(MXAOSQ), COORQX(MXAOSQ), COORPY(MXAOSQ),
     *                COORQY(MXAOSQ), COORPZ(MXAOSQ), COORQZ(MXAOSQ),
     *                EXP12(MXAOSQ),  EXP34(MXAOSQ),
     *                FAC12(MXAOSQ),  FAC34(MXAOSQ)
      COMMON /TWOVEC/ DIFPAX(MXAOSQ), DIFPAY(MXAOSQ), DIFPAZ(MXAOSQ),
     *                DIFPBX(MXAOSQ), DIFPBY(MXAOSQ), DIFPBZ(MXAOSQ),
     *                TEXP1(MXAOSQ),  TEXP2(MXAOSQ),  HEXPPI(MXAOSQ)
      COMMON /CCFCOM/ CONT1 (MXCONT*MXAOVC), CONT2 (MXCONT*MXAOVC),
     *                CONT3 (MXCONT*MXAOVC), CONT4 (MXCONT*MXAOVC),
     *                CONTT1(MXCONT*MXAOVC), CONTT2(MXCONT*MXAOVC),
     *                CONTT3(MXCONT*MXAOVC), CONTT4(MXCONT*MXAOVC),
     *                NUC1X (MXAOVC),        NUC2X (MXAOVC),
     *                NUC3X (MXAOVC),        NUC4X (MXAOVC),
     *                NRC1X (MXAOVC),        NRC2X (MXAOVC),
     *                NRC3X (MXAOVC),        NRC4X (MXAOVC)
      COMMON /MAXOLD/ OLDMAX
      DIMENSION NPRIMA(MXAOVC), NPRIMB(MXAOVC),
     *          NRCOA(MXAOVC),  NRCOB(MXAOVC),
     *          JSTRA(MXAOVC),  JSTRB(MXAOVC)
C
C  Set up overlap distribution vectors 
C
      IF (ITYPE .EQ. 12) THEN
         IADD = 0
         SRFMAX = ZERO
C        THRSH = THRESH
         THRSH = ZERO
         SIGN = ONE
      ELSE
         IADD = MXAOSQ
C        THRSH = THRESH/OLDMAX
         THRSH = ZERO
         SIGN = - ONE
      END IF
      IPRIM = 0
      DO 100 ISETB = 1, NSETB
         ISTRB = JSTRB(ISETB) + 1
         IENDB = JSTRB(ISETB) + NPRIMB(ISETB)
         DO 110 IPRIMB = ISTRB, IENDB
            EXPB  = PRIEXP(IPRIMB)
            CRXB  = SIGNBX*PRICRX(IPRIMB)
            CRYB  = SIGNBY*PRICRY(IPRIMB)
            CRZB  = SIGNBZ*PRICRZ(IPRIMB)
            IF (TPRIAB) THEN
               JENDA = ISETB
            ELSE
               JENDA = NSETA
            END IF
            DO 200 ISETA = 1, JENDA
               ISTRA = JSTRA(ISETA) + 1
               IENDA = JSTRA(ISETA) + NPRIMA(ISETA)
               IF (TPRIAB .AND. ISTRA .EQ. ISTRB) THEN
                  IF (IENDA .GT. IPRIMB) IENDA = IPRIMB
               END IF
               DO 210 IPRIMA = ISTRA, IENDA
                  EXPA  = PRIEXP(IPRIMA)
                  CRXA  = SIGNAX*PRICRX(IPRIMA)
                  CRYA  = SIGNAY*PRICRY(IPRIMA)
                  CRZA  = SIGNAZ*PRICRZ(IPRIMA)
C
                  EXPP  = EXPA + EXPB
                  EXPPI = ONE/EXPP
                  DIFABX = CRXA - CRXB
                  DIFABY = CRYA - CRYB
                  DIFABZ = CRZA - CRZB
                  DISAB2 = DIFABX*DIFABX + DIFABY*DIFABY + DIFABZ*DIFABZ
                  FAC12I = R2PI52*EXPPI*DEXP(- EXPA*EXPB*DISAB2*EXPPI)
                  ABSFAC = ABS(FAC12I)
C                 IF (ABSFAC .GT. THRSH) THEN
                     IPRIM = IPRIM + 1
                     SRFMAX = DMAX1(SRFMAX,ABSFAC)
                     EXPAPI = EXPA*EXPPI
                     EXPBPI = EXPB*EXPPI
                     CORPX = EXPAPI*CRXA + EXPBPI*CRXB
                     CORPY = EXPAPI*CRYA + EXPBPI*CRYB
                     CORPZ = EXPAPI*CRZA + EXPBPI*CRZB
                     COORPX(IPRIM + IADD) = CORPX
                     COORPY(IPRIM + IADD) = CORPY
                     COORPZ(IPRIM + IADD) = CORPZ
                     EXP12 (IPRIM + IADD) = EXPP
                     FAC12 (IPRIM + IADD) = FAC12I
                     DIFPAX(IPRIM) = CORPX - CRXA
                     DIFPAY(IPRIM) = CORPY - CRYA
                     DIFPAZ(IPRIM) = CORPZ - CRZA
                     DIFPBX(IPRIM) = CORPX - CRXB
                     DIFPBY(IPRIM) = CORPY - CRYB
                     DIFPBZ(IPRIM) = CORPZ - CRZB
                     TEXP1 (IPRIM) = EXPA + EXPA
                     TEXP2 (IPRIM) = EXPB + EXPB
                     HEXPPI(IPRIM) = SIGN*HALF*EXPPI
C                 END IF
  210          CONTINUE
  200       CONTINUE
  110    CONTINUE
  100 CONTINUE
      NUCAB = IPRIM
      IF (ITYPE .EQ. 12) OLDMAX = SRFMAX
C
C     *******************************************************
C     ***** Set up matrices of contraction coefficients *****
C     *******************************************************
C
      IF (ITYPE .EQ. 12) THEN
         IOFFM = 0
      ELSE
         IOFFM = 2*MXCONT*MXAOVC
      END IF
C
C     Orbital A
C
      IADR1  = 1
      IADR20 = 1
      DO 300 I = 1, NSETA
         NPRMA = NPRIMA(I)
         NCNTA = NRCOA(I)
         IADR21 = IADR20
         DO 310 J = JSTRA(I) + 1, JSTRA(I) + NPRMA
            IADR2 = IADR21
            DO 320 K = 1, NCNTA
               CCFA = PRICCF(J,K)
               CONT1 (IOFFM + IADR1)  = CCFA
               CONTT1(IOFFM + IADR2) = CCFA
               IADR1 = IADR1 + 1
               IADR2 = IADR2 + NPRMA
  320       CONTINUE
            IADR21 = IADR21 + 1
  310    CONTINUE
         IADR20 = IADR20 + NPRMA*NCNTA
  300 CONTINUE
C
C     Orbital B
C
      IADR1  = 1
      IADR20 = 1
      DO 400 I = 1, NSETB
         NPRMB = NPRIMB(I)
         NCNTB = NRCOB(I)
         IADR21 = IADR20
         DO 410 J = JSTRB(I) + 1, JSTRB(I) + NPRMB
            IADR2 = IADR21
            DO 420 K = 1, NCNTB
               CCFB = PRICCF(J,K)
               CONT2 (IOFFM + IADR1)  = CCFB
               CONTT2(IOFFM + IADR2) = CCFB
               IADR1 = IADR1 + 1
               IADR2 = IADR2 + NPRMB
  420       CONTINUE
            IADR21 = IADR21 + 1
  410    CONTINUE
         IADR20 = IADR20 + NPRMB*NCNTB
  400 CONTINUE
C
C     *******************************************************
C     ***** Number of primitive and contracted orbitals *****
C     *******************************************************
C
      IF (ITYPE .EQ. 12) THEN
         IOFFV = 0
      ELSE
         IOFFV = 2*MXAOVC
      END IF
C
C     Orbital A
C
      NUCA   = 0
      NCONTA = 0
      DO 500 I = 1, NSETA
         NPRMA  = NPRIMA(I)
         NCNTA  = NRCOA(I)
         NUCA   = NUCA + NPRMA
         NCONTA = NCONTA + NCNTA
         NUC1X(IOFFV + I) = NPRMA
         NRC1X(IOFFV + I) = NCNTA
  500 CONTINUE
C
C     Orbital B
C
      NUCB   = 0
      NCONTB = 0
      DO 510 I = 1, NSETB
         NPRMB  = NPRIMB(I)
         NCNTB  = NRCOB(I)
         NUCB   = NUCB + NPRMB
         NCONTB = NCONTB + NCNTB
         NUC2X(IOFFV + I) = NPRMB
         NRC2X(IOFFV + I) = NCNTB
  510 CONTINUE
C
      IF (IPRINT .LT. 10) RETURN
C
C     *************************
C     ***** Print Section *****
C     *************************
C
      IELECT = 1
      IF (ITYPE .EQ. 34) IELECT = 2
      WRITE (LUPRI, 1000) IELECT
      WRITE (LUPRI, '(A,2I5)')   ' NUC   ', NUCA, NUCB
      WRITE (LUPRI, '(A,2I5)')   ' NCONT ', NCONTA, NCONTB
      WRITE (LUPRI, '(A,2I5)')   ' NSET  ', NSETA, NSETB
      WRITE (LUPRI, '(A,I5)')    ' NUCAB ', NUCAB
      WRITE (LUPRI, '(A,F4.1)' ) ' SIGN  ', SIGN
      WRITE (LUPRI, '(A,(5F12.6))') ' TEXP1 ', (TEXP1(I), I = 1, NUCAB)
      WRITE (LUPRI, '(A,(5F12.6))') ' TEXP2 ', (TEXP2(I), I = 1, NUCAB)
      WRITE (LUPRI, '(A,(12I5))') ' NPRIMA ', (NPRIMA(I), I = 1, NSETA)
      WRITE (LUPRI, '(A,(12I5))') ' NPRIMB ', (NPRIMB(I), I = 1, NSETB)
      WRITE (LUPRI, '(A,(12I5))') ' NRCOA  ', (NRCOA (I), I = 1, NSETA)
      WRITE (LUPRI, '(A,(12I5))') ' NRCOB  ', (NRCOB (I), I = 1, NSETB)
      WRITE (LUPRI, '(A,(12I5))') ' JSTRA  ', (JSTRA (I), I = 1, NSETA)
      WRITE (LUPRI, '(A,(12I5))') ' JSTRB  ', (JSTRB (I), I = 1, NSETB)
      IADD = 1
      DO 600 I = 1, NSETA
         WRITE (LUPRI, '(/A,I2,A)')
     *             ' Contraction matrices for set ',I,'.'
         NPRMA = NUC1X(IOFFV + I)
         NRCA  = NRC1X(IOFFV + I)
         WRITE (LUPRI, '(A,2I5)')
     *    ' Number of primitive and contracted functions: ',NPRMA,NRCA
         CALL PRIREC(CONT1 (IOFFM + IADD),NRCA,NPRMA,' CONT1  ')
         CALL PRIREC(CONTT1(IOFFM + IADD),NPRMA,NRCA,' CONTT1 ')
         IADD = IADD + NPRMA*NRCA
  600 CONTINUE
      IADD = 1
      DO 610 I = 1, NSETB
         WRITE (LUPRI, '(/A,I2,A)')
     *             ' Contraction matrices for set ',I,'.'
         NPRMB = NUC2X(IOFFV + I)
         NRCB  = NRC2X(IOFFV + I)
         WRITE (LUPRI, '(A,2I5)')
     *    ' Number of primitive and contracted functions: ',NPRMB,NRCB
         CALL PRIREC(CONT2 (IOFFM + IADD),NRCB,NPRMB,' CONT2  ')
         CALL PRIREC(CONTT2(IOFFM + IADD),NPRMB,NRCB,' CONTT2 ')
         IADD = IADD + NPRMB*NRCB
 610  CONTINUE
 1000 FORMAT(//,' <<<<<<<<<< SUBROUTINE ODCVEC >>>>>>>>>>',
     *       //,' Overlap distribution arrays for electron ',I1,/)
      RETURN
      END
