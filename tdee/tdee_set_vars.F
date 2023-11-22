#include "flags.h"
      SUBROUTINE TDEE_SET_VARS(ICORE,IUHF)

      IMPLICIT INTEGER (A-Z)
      CHARACTER*8 STRING

      LOGICAL CIS,EOMCC,FULDIAG,INCORE,READGUES,ESPROP,NONSTD
      LOGICAL EMINFOL,EVECFOL,VPROP,LAMTHERE,GABCD,CISD
      LOGICAL MBPT2,CC,SS,SD,DS,DD,CCD,RPA,DRPA,CCD_DUMMY,MBPT2_DUMMY
      LOGICAL LCCD,LCCSD,CC2,CCSD
      LOGICAL RESTART_PRESENT

      CHARACTER*11 STRING1(6)
      CHARACTER*6  STRING2(4)

      DIMENSION ICORE(*)
      COMMON /FLAGS/  IFLAGS(100)
      COMMON /FLAGS2/ IFLAGS2(500)
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      COMMON /EOM_METH/ CIS,EOMCC,CISD,FULDIAG,INCORE,READGUES
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /METH/MBPT2_DUMMY,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD_DUMMY,
     +             QCISD,CCSD,UCC
      COMMON /REFTYPE/ MBPT2,CC,CCD,RPA,DRPA,LCCD,LCCSD,CC2
      COMMON /ABCD/ GABCD
      COMMON /CNVRGE/ EMINFOL,EVECFOL
      COMMON /LISTPROJ/ LISTH0, ICOLPR1, ICOLPR2
      COMMON /PROJECT/ IPROJECT, IPATTERN, NCALC, ICALC
      COMMON /DRVHBAR/ SS, SD, DS, DD

      DATA STRING1 /"NO", "SEARCH_ONLY","PRJCT_ALL","PRJCT_NOISE",
     &              "CORE","NTO"/
      DATA STRING2 /"LOWEST", "CORE", "LUMO", "HOMO"/
C
C PROCESS INPUT FOR EXCITATION ENERGY CALCULATION
C
      EOMCC=IFLAGS(87).EQ.3 .OR. IFLAGS(87) .GE. 7 .OR.
     &      IFLAGS(87).EQ.2 .OR. IFLAGS(87) .EQ. 4

      FULLDIAG = .FALSE.
      INCORE   = .FALSE.
      READGUESS= .FALSE.

      CIS   = .FALSE.
      CISD  = .FALSE.
      MBPT3 = .FALSE.
      M4DQ  = .FALSE.
      M4SDQ = .FALSE.
      QCISD = .FALSE.
      UCC   = .FALSE.
     
C  determine reference state (IFLAGS2(117) refer to eomref key-word)
C
      IF (EOMCC) THEN

        CCSD  = IFLAGS2(117) .EQ. 1
        MBPT2 = IFLAGS2(117) .EQ. 2
        CCD   = IFLAGS2(117) .EQ. 3
        RPA   = IFLAGS2(117) .EQ. 4
        DRPA  = IFLAGS2(117) .EQ. 5
        LCCD  = IFLAGS2(117) .EQ. 6
        LCCSD = IFLAGS2(117) .EQ. 7
        CC2   = IFLAGS2(117) .EQ. 9

        CCD_DUMMY   = CCD
        MBPT2_DUMMY = MBPT2
C
C Here CC meant Lambda is not zero. RPA and DRPA is assumed to
C be variational, Ajith Perera, 04, 2014.

        CC = .NOT. (MBPT2 .OR. RPA .OR. DRPA .OR. LCCD .OR. LCCSD)

        IF (CC .AND. .NOT.
     &     (IFLAGS2(117) .EQ. 1 .OR. IFLAGS2(117) .EQ.9)) THEN
          WRITE(6,*) ' SOMETHING STRANGE IN SETMET '
        ENDIF

      ENDIF
      IF (LCCD .OR. LCCSD) THEN
          WRITE(6,*) 
          WRITE(6,"(2a)") " LCCD or LCCSD is not currently supported",
     &                    " with TD-EOM-CC. Instead use EOM-CC."
          WRITE(6,*) 
          CALL ERREX 
      ENDIF 
C
      WRITE(6,*)
      IF (CCSD) THEN
      WRITE(6,103)
  103 FORMAT(T3,'EOM-CCSD excitation energies will be evaluated.')
      ELSE IF (CC2) THEN
      WRITE(6,104)
  104 FORMAT(T3,'EOM-CC2 excitation energies will be evaluated.')
      ELSE IF (LCCD) THEN
      WRITE(6,105)
  105 FORMAT(T3,'EOM-LCCD excitation energies will be evaluated.')
      ELSE IF (LCCSD) THEN
      WRITE(6,106)
  106 FORMAT(T3,'EOM-LCCSD excitation energies will be evaluated.')
      ELSE IF (RPA) THEN
      WRITE(6,107)
  107 FORMAT(T3,'EOM-RPA excitation energies will be evaluated.')
      ELSE IF (DRPA) THEN
      WRITE(6,108)
  108 FORMAT(T3,'EOM-DRPA excitation energies will be evaluated.')
      ELSE IF (CCD) THEN
      WRITE(6,109)
  109 FORMAT(T3,'EOM-CCD excitation energies will be evaluated.')
      ELSE IF (MBPT2) THEN
      WRITE(6,110)
  110 FORMAT(T3,'EOM-MBPT(2) excitation energies will be evaluated.')
       ENDIF
C
      GABCD=IFLAGS(100).EQ.0
C
C  LOGIC TO DETERMINE COMMON BLOCK /DRVHBAR/
C
      SS = .TRUE.
      SD = .TRUE.
      DS = .TRUE.
C
C EOM triples are introduced. DD is false only for partitioned
C and EOM-BWPT2 methods, 03/2014. Ajith Perera

      IF (IFLAGS(87) .EQ. 7 .OR. IFLAGS(87) .EQ. 8) THEN
        DD = .FALSE.
      ELSE
        DD = .TRUE.
      ENDIF

      IRREPX = 1
      CALL GETAOINF(IUHF,IRREPX)
C
C CHECK HBARLISTS
C
      IF (DD) THEN
        IF (IFLAGS2(122) .EQ. 2) THEN
          WRITE(6,*) ' HBARABCD SHOULD CORRESPOND TO BARE INTEGRALS'
          CALL ERREX
        ENDIF
        IF (CC .AND. IFLAGS2(123) .EQ. 2) THEN
          WRITE(6,*) ' T1*ABCD SHOULD NOT BE INCLUDED IN HBARABCI'
          CALL ERREX
        ENDIF
      ENDIF
C
C   LOGIC TO DETERMINE COMMON BLOCK /PROJECT/
C
      IPROJECT = Iflags2(h_IFLAGS2_eom_prjct)
      LISTH0 = 472
      ICOLPR1 = 5
      ICOLPR2 = 6
      IPATTERN =  Iflags2(h_IFLAGS2_ee_search)
C
      IF (IPATTERN .NE. 1) THEN
        IF (IPROJECT .EQ. 3) IPROJECT = 3
      ENDIF
C
C Projection must be ensured when IPROJECT=4 (NTO)
C
      IF (IPATTERN .EQ. 0 .AND. .NOT. (IPROJECT .EQ. 5)) IPROJECT = 0
C
      IF (Iflags2(h_IFLAGS2_core_window) .EQ. 0) THEN
         CALL TDEE_DETNCALC(ICORE,IUHF)
      ELSE
         NCALC = 1
      ENDIF
C
C  write out some stuff, related to projection
C
      IF (IPATTERN .NE. 0 .OR. IPROJECT .NE. 0) THEN
      write(6,*) ' Projection used in eom-cc calculation: ',
     &             String1(iproject+1)
      write(6,*) ' Excitation pattern: ',
     &             String2(ipattern+1)
      write(6,*) ' Number of different eom-cc calculations', ncalc
      WRITE(6,*)
      ENDIF

      EVECFOL = .FALSE.
      EMINFOL = .NOT. EVECFOL

      IF (EVECFOL) THEN
        WRITE(6,"(a,a)") "  Root following procedure: Time propagation",
     +               " or Lanczos"
      ELSE
        WRITE(6,"(a,a)") "  Root following procedure: Time propagation",
     +               " or Lanczos"
      ENDIF
C
      RETURN
      END
