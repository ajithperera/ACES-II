#include "flags.h"
      SUBROUTINE PSI4DBG_DRDTGSS(SCR, MAXCOR, IUHF, NAME, NLIST, 
     +                           READT)
C Driver for the rdtgss routine.
C
C SG 7/22/98
C
      IMPLICIT NONE
C
      INTEGER MAXCOR, IUHF, NLIST
      DOUBLE PRECISION SCR(MAXCOR)
      LOGICAL READT
      CHARACTER *8 NAME
C
      INTEGER POP, VRT, NT, NFMI, NFEA
      COMMON /SYM/ POP(8,2), VRT(8,2), NT(2), NFMI(2), NFEA(2)
      INTEGER IRPDPD, ISYTYP, ID
      COMMON /SYMPOP/ IRPDPD(8,22), ISYTYP(2,500), ID(18)
      INTEGER NSTART, NIRREP, IRREPS, DIRPRD
      COMMON /SYMINF/ NSTART, NIRREP, IRREPS(255,2), DIRPRD(8,8)
C
#include "flags2.com"
#include "flags.com"
C
      INTEGER ISYMSZ
C
      INTEGER I, OCC(8)
      INTEGER NAMLEN, I000, I010, I020, I030, I040, I050
      LOGICAL EXISTS
      CHARACTER *80 FULNAM
      CHARACTER *4  CPGRP, TPGRP
C
      CALL GFNAME(NAME, FULNAM, NAMLEN)
      INQUIRE(FILE=FULNAM(1:NAMLEN), EXIST=EXISTS)
      CALL GETCREC(01,'JOBARC','COMPPTGP',4,CPGRP)
      CALL GETCREC(-1,'JOBARC','PREVPTGP',4,TPGRP)
C
      READT = .FALSE.
c   o Make sure the occupation has not changed from the previous point.
      CALL GETREC(-1,'JOBARC','TGSSOCCA',NIRREP,OCC)
      DO I = 1, NIRREP
         IF (OCC(I).NE.POP(I,1).AND.OCC(I).NE.0) EXISTS=.FALSE.
      END DO
      IF (IUHF.NE.0) THEN
         CALL GETREC(-1,'JOBARC','TGSSOCCB',NIRREP,OCC)
         DO I = 1, NIRREP
            IF (OCC(I).NE.POP(I,2).AND.OCC(I).NE.0) EXISTS=.FALSE.
         END DO
      END IF
c   o If the PREVPTGP record does not exist, then getcrec will return
c     the null string. Let's assume that the TGUESS file has the same
c     symmetry as the current calculation.
      IF (EXISTS.AND.(CPGRP.EQ.TPGRP.OR.ICHAR(TPGRP(1:1)).EQ.0)) THEN
        I000 = 1
        I010 = I000 + NT(1)
        IF (IUHF .NE. 0) THEN
          I020 = I010 + NT(2)
          I030 = I020 + ISYMSZ(ISYTYP(1,45+NLIST), ISYTYP(2,45+NLIST))
        ELSE
          I030 = I010
        ENDIF
        I040 = I030 + ISYMSZ(ISYTYP(1,44+NLIST), ISYTYP(2,44+NLIST))
        I050 = I040 + ISYMSZ(ISYTYP(1,46+NLIST), ISYTYP(2,46+NLIST))
        IF (I050 .GT. MAXCOR) THEN
          WRITE(6,1000)
 1000     FORMAT(T3, '@DRTGSS-I, not enough memory to read T',
     &       ' as initial guess')
C
        ELSE
          CALL RDTGSS(SCR(I000), I050-I000, NAME, EXISTS)
CSSS          call checksum("rdtgss", scr(i000), I050-i000)

          If (IFLAGS(h_IFLAGS_calc) .EQ.8) Then
              IF (IUHF .EQ.0) THEN
                 Write(6,*) "Zeroing out T1"
                 CALL DZERO(SCR(I000),NT(1)) 
              ELSE
                 CALL DZERO(SCR(I000),NT(1)+NT(2)) 
              ENDIF
          ENDIF 

          CALL PUTLST(SCR(I000), 1, 1, 1, 1, 90+NLIST)

          IF (IUHF.NE.0) THEN
            CALL PUTLST(SCR(I010), 1, 1, 1, 2, 90+NLIST)          
            CALL PUTALL(SCR(I020),
     &         ISYMSZ(ISYTYP(1,45+NLIST), ISYTYP(2,45+NLIST)),
     &         1, 45+NLIST)
          ENDIF
          CALL PUTALL(SCR(I030),
     &       ISYMSZ(ISYTYP(1,44+NLIST), ISYTYP(2,44+NLIST)),
     &       1, 44+NLIST)
          CALL PUTALL(SCR(I040),
     &       ISYMSZ(ISYTYP(1,46+NLIST), ISYTYP(2,46+NLIST)),
     &       1, 46+NLIST)
          READT = .TRUE.
      Write(6,"(2x,a)"),' Starting T amplitudes are read from TGUESS'
        ENDIF
      ENDIF
C
      RETURN
      END
