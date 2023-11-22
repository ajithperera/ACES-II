
      SUBROUTINE EHT_PRJSTO3G_GUESS(DENS, CZAO, CSO, ZAOSO,
     &                              SCR1, SCRN, DCORE, IUHF,
     &                              NBAS, NBASX, LDIM1, MAXDCOR)

      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      LOGICAL MBS
C
#include "machsp.com"
#include "symm2.com"
#include "flags.com"
C
       DIMENSION DCORE(MAXDCOR),DENS((IUHF+1)*LDIM1),SCR1(NBAS*NBAS),
     &           SCRN(NBAS*NBAS), CZAO(NBASX*NBAS,(IUHF+1)),
     &           ZAOSO(NBAS*NBASX)

C-----------------------------------------------------------------------
C     JDW, Summer 2004. This is a fudge to enable two guess options.
C
C     IFLAGS(45) = 6 The program will read a file called MOFILE to
C                    get MOs in a minimal basis set and then generate
C                    an initial guess for another basis set. This has
C                    only been coded for certain basis sets. See the
C                    routine GETDENS3 and the SPLBAS file. MOFILE must
C                    have been created previously. MOFILE is identical
C                    to AOBASMOS. IFLAGS(45) = 6 is synonymous with the
C                    2-step procedure.
C
C     IFLAGS(45) = 7 The program will generate an EHT guess and write
C                    a file called EHTMOFILE. If this is NOT a minimal
C                    basis set calculation, GETDENS3 will then try to
C                    create an initial guess for the actual basis set
C                    being used. This is the one-step procedure. Alter-
C                    natively, if this option is run with a minimal
C                    basis set, GETDENS3 will be bypassed and an STO-3G
C                    initial density generated from the EHT coeff-
C                    icients. An STO-3G scf calculation will then be
C                    performed, created AOBASMOS, which can be renamed
C                    to MOFILE for the second job of the two-step proc-
C                    edure.
C-----------------------------------------------------------------------
C     Projecting an STO-3G set of MOs into current basis set.
C     Limitations for ROHF right now --- only creates occupied MOs and
C     and alpha and beta density matrices.
C-----------------------------------------------------------------------
C
      IF(IFLAGS(45).EQ.6)THEN
       call getdens3(dens,czao,cso,zaoso,scr1,scrn,
     &               nbasx,nbas,natoms,
     &               nirrep,ldim1,itriof,ireps,nbfirr,iuhf,2)
      ENDIF
C
C-----------------------------------------------------------------------
C     Extended Huckel Guess. Only available for H-Cl and STO-3G.
C     Limitations for ROHF right now --- only creates occupied MOs and
C     and alpha and beta density matrices.
C-----------------------------------------------------------------------
C
      IF(IFLAGS(45).EQ.7)THEN
       CALL EHTDIM(NBASEHT,NBASXEHT)
       CALL CALC_S(DCORE,DCORE(1+NBASXEHT*NBASXEHT),
     &                   DCORE(1+NBASXEHT*NBASXEHT+3*NATOMS),
     &                   DCORE(1+NBASXEHT*NBASXEHT+3*NATOMS),
     &             NBASEHT,NBASXEHT,NATOMS)
C
C     CALC_S puts overlap integrals in DCORE. Length: NBASXEHT*NBASXEHT.
C
cccc This has potential      CALL SYMTRANS
       CALL EHTGSS(DENS,DCORE(1+NBASX*NBASX),DCORE,
     &             DCORE(1+2*NBASX*NBASX),DCORE(1+3*NBASX*NBASX),SCR1,
     &             EVAL,NBASEHT,NBASXEHT,LDIM1,LDIM2,IUHF)
C
       MBS = (NBASX.EQ.NBASXEHT)
C
C     Read information from EHTMOFILE to create initial guess in
C     current basis set.
C
C     EHTGSS prepares EHTMOFILE. MOFILE has to be manually prepared:
C     it is the same as AOBASMOS for a minimal basis set calculation.
C
       IF(.NOT.MBS)THEN
        call getdens3(dens,czao,cso,zaoso,scr1,scrn,
     &                nbasx,nbas,natoms,
     &                nirrep,ldim1,itriof,ireps,nbfirr,iuhf,1)
       ENDIF
      ENDIF

      Return
      End

