#include "flags.h"
      SUBROUTINE INIT_IPEOM(ICORE, MAXCOR, IUHF, ICALC)
C
      IMPLICIT INTEGER (A-Z)
      DIMENSION ICORE(MAXCOR)
      CHARACTER TITLE*80
      CHARACTER*4 BOOL(2)

      LOGICAL  EMINFOL,EVECFOL,LEFTHAND,SINGONLY,DROPCORE
      LOGICAL  CCSD,HBAR_ABCD,HBAR_ABCI,ABCD_TYPE 
      LOGICAL  TRANABCD,TRANABCI

#include "flags.com"
#include "flags2.com"
#include "machsp.com"
#include "sym.com"
#include "syminf.com"
#include "sympop.com"
C
      COMMON/LISTDAV/LISTC,LISTHC,LISTH0
      COMMON/SLISTS/LS1IN,LS1OUT,LS2IN(2,2),LS2OUT(2,2)
      COMMON/LISTDIP/DIPOO(2),DIPOV(2),DIPVV(2)
      COMMON/LISTDENS/DENSOO(2),DENSOV(2),DENSVV(2)
      COMMON/EXTRAP/MAXEXP,NREDUCE,NTOL,NSIZEC
      COMMON/CNVRGE/EMINFOL,EVECFOL
      COMMON/IPCALC/LEFTHAND,SINGONLY,DROPCORE
      COMMON/IPINFO/NUMROOT(8,3)
      COMMON/IPCORE/COREIP 
      COMMON/LISTF/LISTFJA,LISTFAJ,LISTFAB,LISTFIJ,LISTF0
     &             NDIMFLST(2,4)
       COMMON/STINFO/ITOTALS,STMODE,LISTST,STCALC,NSIZEST
C
      DATA BOOL /'FALSE', ' TRUE'/
      DATA IONE,ONE,ONEM /1,1.0D0,-1.0D0/

      LISTC   = 470
      LISTHC  = 471
      LISTH0  = 472
      LISTF0  = 472
      LISTST  = 469
      LISTFAJ = 473
      LISTFJA = 474
      LISTFAB = 475
      LISTFIJ = 476
      LISTUIJ = 477
      LISTUAB = 478
      LISTUAI = 479

      LS1IN   = 420
      LS1OUT  = 430
C
      DIPOO(1) = 441
      DIPOO(2) = 442
      DIPOV(1) = 443
      DIPOV(2) = 444
      DIPVV(1) = 445
      DIPVV(2) = 446
C
      DENSOO(1) = 451
      DENSOO(2) = 452
      DENSOV(1) = 453
      DENSOV(2) = 454
      DENSVV(1) = 455
      DENSVV(2) = 456
C
      LS2IN(1,1) = 460
      LS2IN(1,2) = 461
      LS2IN(2,1) = 462
      LS2IN(2,2) = 463
C
      LS2OUT(1,1) = 480
      LS2OUT(1,2) = 481
      LS2OUT(2,1) = 482
      LS2OUT(2,2) = 483
C
      MAXEXP  =  30
      NREDUCE = 12
      NTOL    = IFLAGS(h_IFLAGS_estate_tol) 
      SIDE    = IFLAGS2(h_IFLAGS2_eom_solve)
      STMODE  = 0
C
      EMINFOL  = .TRUE.
      EVECFOL  = .FALSE.
      IF (SIDE .EQ. 1) LEFTHAND = .TRUE.
      SINGONLY = .FALSE.
      DROPCORE = .FALSE.
      EXCICORE = .FALSE.
C
      DO ISPIN = 1, 3
        DO IRREP = 1, NIRREP
           NUMROOT(IRREP,ISPIN) = 0
        ENDDO
      ENDDO 
C
      WRITE(6,1300)
      WRITE(6,"(a)")' **********************************************'
      WRITE(6,"(a)")' *                                            *'
      WRITE(6,"(a)")' *   ELECTRON DETACHMENT EOM-CC CALCULATION   *'
      WRITE(6,"(a)")' *                                            *'
      WRITE(6,"(a)")' *            CODED BY A. Perera              *'
      WRITE(6,"(a)")' *               FALL 2016                    *'
      WRITE(6,"(a)")' **********************************************'
      WRITE(6,1300)
C
C
C NUMROOT INFO IS OBTAINED FROM JOBARC
C
      CALL GETREC(20,'JOBARC','IPSYM_A ',NIRREP,NUMROOT(1,1))
      IF (IUHF .NE. 0) THEN
         CALL GETREC(20,'JOBARC','IPSYM_B ',NIRREP,NUMROOT(1,2))
      ENDIF         
C
      IF (SINGONLY) THEN

        DO IRREP = 1, NIRREP
          DO ISPIN = 1, 1+IUHF
              NUMROOT(IRREP,ISPIN) = MIN(POP(IRREP,ISPIN),
     +        NUMROOT(IRREP,ISPIN))
          ENDDO
        ENDDO
      ENDIF
C
C IF ALL NUMROOT ARE ZERO, ESTIMATE DESIRED ROOTS FROM SCF EIGENVALUES
C
      IROOT = 0
      DO IRREP = 1, NIRREP
        DO ISPIN = 1, 1 + IUHF
          IROOT = IROOT + NUMROOT(IRREP,ISPIN)
        ENDDO
      ENDDO

      IF (IROOT .GT. 0) THEN
         CALL GETREC(-1,'JOBARC','IP_IRREP',IONE,NIRREP2)         

         IF (NIRREP2 .NE. NIRREP) THEN
             WRITE(6,"(a,a,2I2)") ' NUMBER OF IRREPS FROM IP_SYM',
     +                  ' DOES NOT MATCH NIRREP', NIRREP2, NIRREP
          WRITE(6,"(a)") ' NUMBER OF ROOTS IS ESTIMATED AUTOMATICALLY'
          CALL IZERO(NUMROOT, 3*8)
        ENDIF
      ENDIF
C
      CALL INITNUMR_IP(ICORE, MAXCOR/IINTFP, IUHF)
C
      WRITE(6,*)
      WRITE(6,"(a)") ' NUMBER OF ROOTS PER IRREP, SPIN = ALPHA'
      WRITE(6,1100) (NUMROOT(i,1), i=1,nirrep)

      IF (IUHF.NE.0) THEN
        IF (STMODE .EQ. 2) THEN
          WRITE(6,"(2a)") ' NUMBER OF ROOTS PER IRREP ',
     $       'SPIN = BETA, TRIPLET'
          write(6,1100) (numroot(i,2), i=1,nirrep)
          WRITE(6,"(2a)") ' NUMBER OF ROOTS PER IRREP ',
     $       'SPIN = BETA, SINGLET'
          write(6,1100) (numroot(i,3), i=1,nirrep)
        ELSE
          WRITE(6,"(2a)") ' NUMBER OF ROOTS PER IRREP ',
     $       'SPIN = BETA'
          write(6,1100) (numroot(i,2), i=1,nirrep)
        ENDIF
      ENDIF

 1100 FORMAT(8x, 8I5)
      WRITE(6,1300)
 1300 FORMAT(/)
C
C NEW 400 LISTS ARE CREATED => DELETE POSSIBLE EXISTING DERGAM
C
      ICALC = 2
      CALL ACES_IO_REMOVE(54,'DERGAM')
      CALL CREATE_IPLISTS(IUHF,ICALC)
C
C INCLUDE THE DIAGONAL PART IN LISTS 91 AND 92
C
      CALL MODF_VIP(ICORE, MAXCOR, IUHF, 1)
C
      IF (IFLAGS(93) .EQ. 2) THEN
        WRITE(6,"(a)") ' AO-TYPE ABCD INTEGRALS NOT SUPPORTED'
        CALL ERREX
      ENDIF
C
      IF (ISYTYP(1, 233) .EQ. 5) THEN
        WRITE(6,"(a)")
        WRITE(6,"(a)")' COMPRESSED ABCD INTEGRALS'
        WRITE(6,"(a)")' NOT SUPPORTED IN CURRENT VERSION OF IP-EOMCC'
        CALL ERREX
      ENDIF
c
      CCSD      = .FALSE.
      HBAR_ABCD = .FALSE.
      HBAR_ABCI = .FALSE.
      ABCD_TYPE = .FALSE.
C
      CCSD       = (iFlags2(h_IFLAGS2_ref).eq.1)
      HBAR_ABCD  = (iFlags2(h_IFLAGS2_hbarabcd).eq.1)
      HBAR_ABCI  = (iFlags2(h_IFLAGS2_hbarabci).eq.1)
      ABCD_TYPE  = (iFlags(h_IFLAGS_abcdtype).eq.2)
C
      TRANABCD = (ISYTYP(1,233).NE.5 .AND. ABCD_TYPE .AND. HBAR_ABCD)
C
      TRANABCI = (CCSD .AND. HBAR_ABCI .AND. (ISYTYP(1,233) .NE. 5)
     +            .AND. (.NOT. ABCD_TYPE))
C
      IF (CCSD .AND. HBAR_ABCI) THEN
        write(6,"(a,a)") ' ABCD contribution to HBARABCI',
     +                   ' is not included !'
      ENDIF
C
      CALL MODHBAR(ICORE, MAXCOR, IUHF, TRANABCD, TRANABCI)
C
C ABCD INTEGRALS NEED TO BE MODIFIED BEFORE PHPH INTEGRALS !!
C
C  INCLUDE DENOMINATORS IN PHPH INTEGRALS
C
CSSS      CALL MODWPHPH_IPEOM(ICORE, MAXCOR, IUHF, 1.0D0, 91, 92)
C
      RETURN
      END
