#include "flags.h"
      SUBROUTINE DRDRGSS(SCR,MAXCOR,IUHF,NONSTD)
C Driver for the rdtgss routine.
C
C SG 7/22/98
C
      IMPLICIT NONE
C
      INTEGER MAXCOR, IUHF, NLIST, IEND
      DOUBLE PRECISION SCR(MAXCOR)
      DOUBLE PRECISION S
      LOGICAL NONSTD
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
      INTEGER NAMLEN, I000, I010, I020, I030, I040, I050
C
      I000 = 1
      I010 = I000 + NT(1)
      IF (IUHF .NE. 0) THEN
          I020 = I010 + NT(2)
          I030 = I020 + ISYMSZ(ISYTYP(1,445), ISYTYP(2,445))
      ELSE
          I030 = I010
      ENDIF
      I040 = I030 + ISYMSZ(ISYTYP(1,44), ISYTYP(2,44))
      IEND = I040 + ISYMSZ(ISYTYP(1,46), ISYTYP(2,46))

      IF (IEND .GT. MAXCOR) THEN
          WRITE(6,1000)
 1000     FORMAT(T3, '@-DRDRGSS, not enough memory to read R',
     &       ' as initial guess')
C         CALL INSMEM("DRDRGSS",IEND,MAXCOR)
      ELSE
          CALL RDRGSS(SCR(I000),IEND-I000,"RGUESS  ",.TRUE.)

#ifdef _DEBUG_LVL0
        call checksum("DRDRGSS", scr(i000),Iend-i000,S)
#endif 
          CALL PUTLST(SCR(I000), 1, 1, 1, 1, 490)

          IF (IUHF.NE.0) THEN
            CALL PUTLST(SCR(I010), 1, 1, 1, 2, 490)          
            CALL PUTALL(SCR(I020),ISYMSZ(ISYTYP(1,445),ISYTYP(2,445)),
     &                  1,445)
          ENDIF
          IF (.NOT. NONSTD) THEN
          CALL PUTALL(SCR(I030),ISYMSZ(ISYTYP(1,444),ISYTYP(2,444)),
     &                1,444)
          ENDIF 
          CALL PUTALL(SCR(I040),ISYMSZ(ISYTYP(1,446),ISYTYP(2,446)),
     &                1,446)
#ifdef _DEBUG_LVL0
          Write(6,"(2x,2a)") ' Starting R amplitudes are read from',
     &                       ' RGUESS'
#endif 
      ENDIF
C
      RETURN
      END
