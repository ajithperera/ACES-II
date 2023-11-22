      SUBROUTINE DRSYM2(SO,AO,HKABCD,ISYMR,ISYMT,ISYMTS,MULATM,MULE,
     &                  SHABAB,NOABCD,LWRKAO,IPRINT)
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL TRANS
      LOGICAL TWOCEN, THRCEN, FOUCEN, DERONE, DERTWO
      LOGICAL         AEQB,   CEQD,   DIAGAB, DIAGCD, DIACAC,
     &                ONECEN, PQSYM,  DTEST,
     &                TPRIAB, TPRICD, TCONAB, TCONCD
      LOGICAL SHABAB
      INTEGER AND,OR,XOR
      DOUBLE PRECISION
     &        SIGN1X, SIGN1Y, SIGN1Z, SIGN2X, SIGN2Y, SIGN2Z,
     &        SIGN3X, SIGN3Y, SIGN3Z, SIGN4X, SIGN4Y, SIGN4Z
C
      PARAMETER (LUCMD = 5, LUPRI = 6)
#include "baslims.par"
#include "mxcent.par"
C
      DIMENSION AO(LWRKAO), SO(*)
C
      COMMON /DIRPRT/ SIGNDR(4,3), TRANS, NATOMS, IATOMS(4), ISOPDR(4)
      COMMON /EXPCOM/ SIGN1X, SIGN1Y, SIGN1Z, SIGN2X, SIGN2Y, SIGN2Z,
     &                SIGN3X, SIGN3Y, SIGN3Z, SIGN4X, SIGN4Y, SIGN4Z,
     &                NCENT1, NCENT2, NCENT3, NCENT4,
     &                ISO1,   ISO2,   ISO3,   ISO4,
     &                DERONE, DERTWO, TWOCEN, THRCEN, FOUCEN,
     &                NINTYP, NCCINT
      COMMON /INTINF/ THRESH,
     &                NHKTA,  NHKTB,  NHKTC,  NHKTD,
     &                MAXAB,  MAXCD,  JMAX0,
     &                KHKTA,  KHKTB,  KHKTC,  KHKTD,
     &                KHKTAB, KHKTCD, KHABCD,
     &                MHKTA,  MHKTB,  MHKTC,  MHKTD,
     &                MULA,   MULB,   MULC,   MULD,
     &                NORBA,  NORBB,  NORBC,  NORBD, NORBAB, NORBCD,
     &                NUCA,   NUCB,   NUCC,   NUCD,  NUCAB,  NUCCD,
     &                NSETA,  NSETB,  NSETC,  NSETD,
     &                ISTEPA, ISTEPB, ISTEPC, ISTEPD,
     &                NSTRA,  NSTRB,  NSTRC,  NSTRD,
     &                AEQB,   CEQD,
     &                DIAGAB, IAB0X,  IAB0Y,  IAB0Z,
     &                DIAGCD, ICD0X,  ICD0Y,  ICD0Z,
     &                DIACAC, ONECEN, PQSYM,  IPQ0X, IPQ0Y, IPQ0Z,
     &                TPRIAB, TPRICD, TCONAB, TCONCD,
     &                MAXDER, DTEST
      COMMON /SYMMET/ FMULT(0:7), PT(0:7),
     &                MAXLOP, MAXLOT, MULT(0:7), ISYTYP(3),
     &                ITYPE(8,36), NPARSU(8), NPAR(8), NAOS(8),
     &                NPARNU(8,8), IPTSYM(MXCORB,0:7),
     &                IPTCNT(3*MXCENT,0:7), NCRREP(0:7),
     &                IPTCOR(MXCENT*3), NAXREP(0:7), IPTAX(3),
     &                IPTXYZ(3,0:7)
C
      IBTAND(I,J) = AND(I,J)
      IBTOR(I,J)  = OR(I,J)
      IBTSHL(I,J) = ISHFT(I,J)
      IBTSHR(I,J) = ISHFT(I,-J)
      IBTXOR(I,J) = XOR(I,J)
C
      IF (IPRINT .GT. 10) THEN
       CALL HEADER('Subroutine DRSYM2',-1)
       WRITE (LUPRI,'(A,3I5)') ' ISYMR, ISYMT, ISYMTS ',
     &                             ISYMR, ISYMT, ISYMTS
      END IF
      IF (NATOMS .GT. 1) THEN
       NWRKAO = 3*(NATOMS + 1)*NCCINT
       IF (NWRKAO .GT. LWRKAO) THEN
        WRITE (LUPRI,'(/,A,2(/,A,I10))')
     &       ' Work space exceeded in DRSYM2.',
     &       ' NWRKAO:', NWRKAO, ' LWRKAO:', LWRKAO
        WRITE (LUPRI,'(A)') ' Increase dimension of WORK1.'
        CALL ERREX 
       END IF
      END IF
      ISOSTR = 1
      DO 100 ICOOR = 1, 3
       ISYTYE = ISYTYP(ICOOR)
       DO 200 IREPE = 0, MAXLOP
        IF (IBTAND(MULE,IBTXOR(IREPE,ISYTYE)) .EQ. 0) THEN
         IF (NATOMS .EQ. 1) THEN
          FACSYM = PT(IBTAND(ISOPDR(1),IREPE))
     &             *SIGNDR(1,ICOOR)*HKABCD
          IAOSTR = (ICOOR - 1)*NCCINT
         ELSE
          FACSYM = HKABCD
          IAOSTR = 3*NATOMS*NCCINT
          DO 300 ICENT = 1, NATOMS
           FACTOR = PT(IBTAND(ISOPDR(ICENT),IREPE))
     &              *SIGNDR(ICENT,ICOOR)
           IOFF   = (3*ICENT + ICOOR - 4)*NCCINT
           IF (ICENT .EQ. 1) THEN
Cparallelization note: no synchronization problems
Cparallelization note: no synchronization problems
           DO 400 INT = 1, NCCINT
            AO(IAOSTR + INT) = FACTOR*AO(IOFF + INT)
  400      CONTINUE
          ELSE
Cparallelization note: no synchronization problems
Cparallelization note: no synchronization problems
          DO 450 INT = 1, NCCINT
           AO(IAOSTR + INT) = AO(IAOSTR + INT)
     &                      + FACTOR*AO(IOFF + INT)
  450     CONTINUE
         END IF
  300   CONTINUE
       END IF
       IF (IPRINT .GT. 20) THEN
        WRITE (LUPRI,'(A,3I5)') ' ICOOR, IREPE, ISYTYE ',
     &                            ICOOR, IREPE, ISYTYE
        WRITE (LUPRI,'(A,F12.6)') ' FACSYM ', FACSYM
       END IF
       CALL SYM2(SO(ISOSTR),AO(IAOSTR+1),ISYMR,ISYMT,ISYMTS,
     &           SHABAB,NOABCD,FACSYM,IREPE,NSOINT,IPRINT)
       ISOSTR = ISOSTR + NSOINT
      END IF
  200 CONTINUE
  100 CONTINUE
      RETURN
      END
