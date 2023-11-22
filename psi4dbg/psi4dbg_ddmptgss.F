
#include "flags.h"
#include "disttype.h"

      SUBROUTINE PSI4DBG_DDMPTGSS(SCR, MAXCOR, IUHF, NLIST, NAME)
C
C Driver for DMPTGSS.  This routine figures out if it is safe to dump a T
C vector to use as a later guess.  If so, it loads T into core and calls DMPTGSS.
C
C SG 7/22/98
C
      IMPLICIT NONE
C
      INTEGER MAXCOR, IUHF, NLIST, LEngth
      DOUBLE PRECISION SCR(MAXCOR)
      CHARACTER *8 NAME
C
      INTEGER IFLAGS
      COMMON /FLAGS/ IFLAGS(100)
#include "flags2.com"
      INTEGER POP, VRT, NT, NFMI, NFEA
      COMMON /SYM/ POP(8,2), VRT(8,2), NT(2), NFMI(2), NFEA(2)
      INTEGER IRPDPD, ISYTYP, ID
      COMMON /SYMPOP/ IRPDPD(8,22), ISYTYP(2,500), ID(18)
C
      INTEGER ISYMSZ
      CHARACTER *4 FPGRP, CPGRP
C
      INTEGER I000, I010, I020, I030, I040, I050
C
C Dump the T vector if this is a finite difference calculation, geometry
C optimization, or Brueckner calculation. (And only if we have the memory.)
C
      I000 = 1
      I010 = I000 + NT(1)
      I020 = I010 + IUHF*NT(2)
      I030 = I020 + IUHF*ISYMSZ(h_DistType_altb,h_DistType_iltj)
      I040 = I030 + ISYMSZ(h_DistType_AltB,h_DistType_IltJ)
      I050 = I040 + ISYMSZ(h_DistType_Ab,h_DistType_Ij)
      IF (I050 .GT. MAXCOR) THEN
         WRITE(6,1000)
 1000    FORMAT(T3, '@DDMPTGSS-I, not enough memory to dump T',
     &              ' as initial guess')
      ELSE
         CALL GETCREC(0,'JOBARC','PREVPTGP',I000,CPGRP)
         IF (I000.LT.0) THEN
            CALL GETCREC(20,'JOBARC','COMPPTGP',4,CPGRP)
            CALL PUTCREC(20,'JOBARC','PREVPTGP',4,CPGRP)
         END IF
         CALL PUTREC(1,'JOBARC','TGSSOCCA',8,POP(1,1))
         IF (IUHF.NE.0) THEN
            CALL PUTREC(1,'JOBARC','TGSSOCCB',8,POP(1,2))
         END IF
         I000 = 1
         CALL GETLST(SCR(I000), 1, 1, 1, 1, 90+NLIST)
         IF (IUHF .NE. 0) THEN
            CALL GETLST(SCR(I010), 1, 1, 1, 2, 90+NLIST)
            CALL GETALL(SCR(I020),
     &                  ISYMSZ(ISYTYP(1,45+NLIST), ISYTYP(2,45+NLIST)),
     &                  1, 45)
         ENDIF
         CALL GETALL(SCR(I030),
     &               ISYMSZ(ISYTYP(1,44+NLIST), ISYTYP(2,44+NLIST)),
     &               1, 44)
         CALL GETALL(SCR(I040),
     &               ISYMSZ(ISYTYP(1,46+NLIST), ISYTYP(2,46+NLIST)),
     &               1, 46)
         CALL DMPTGSS(SCR(I000), I050-I000, NAME)
      ENDIF
C
      RETURN
      END
