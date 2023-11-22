
#include "flags.h"
#include "disttype.h"

      SUBROUTINE DDMPRGSS(SCR, MAXCOR, IUHF, NLIST, NAME)
C
C Driver for DMPRGSS. This routine figures out if it is safe to dump a R
C vector to use as a later guess.  If so, it loads R into core and calls 
c DMPRGSS.

      IMPLICIT NONE
C
      INTEGER MAXCOR, IUHF, LEngth, NLIST
      DOUBLE PRECISION SCR(MAXCOR), S
      CHARACTER *8 NAME
C
      INTEGER IFLAGS
      COMMON /FLAGS/ IFLAGS(100)
#include "flags2.com"
#include "syminf.com"
      INTEGER POP, VRT, NT, NFMI, NFEA, i
      COMMON /SYM/ POP(8,2), VRT(8,2), NT(2), NFMI(2), NFEA(2)
      INTEGER IRPDPD, ISYTYP, ID

      INTEGER IRREPR, IRREPL, DISSZT, NUMDST, DISSZF, NUMDSF
      INTEGER IOFFF, IOFFT, NSIZE, I060, IRREPX, ITMP1, ITMP2, IDIS

      COMMON /SYMPOP/ IRPDPD(8,22), ISYTYP(2,500), ID(18)
C
      INTEGER ISYMSZ
      CHARACTER *4 FPGRP, CPGRP
C
      INTEGER I000, I010, I020, I030, I040, I050, I070
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
      I060 = I050 + ISYMSZ(h_DistType_Ab,h_DistType_Ij)

      IF (I060 .GT. MAXCOR) THEN
         WRITE(6,1000)
 1000    FORMAT(T3, '@DDMPRGSS-I, not enough memory to dump T',
     &              ' as initial guess')
      ELSE

         I000 = 1
         CALL GETLST(SCR(I000), 1, 1, 1, 1, 90+NLIST)
         IF (IUHF .NE. 0) THEN
            CALL GETLST(SCR(I010), 1, 1, 1, 2, 90+NLIST)
            CALL GETALL(SCR(I020),
     &                  ISYMSZ(ISYTYP(1,45+NLIST), ISYTYP(2,45+NLIST)),
     &                  1, 45+NLIST)
         ENDIF

         IF (IUHF .NE. 0) THEN

         CALL GETALL(SCR(I030),
     &               ISYMSZ(ISYTYP(1,44+NLIST), ISYTYP(2,44+NLIST)),
     &               1, 44+NLIST)

         ELSE
 
C
C Lets generate R2(A<B,I<J) from R(Ab;Ij) for RHF
c
         CALL GETALL(SCR(I040), 
     &               ISYMSZ(ISYTYP(1,46+NLIST), ISYTYP(2,46+NLIST)),
     &               1, 46+NLIST)

#ifdef _DEBUG_LVLM
        NSIZE = ISYMSZ(ISYTYP(1,46), ISYTYP(2,46))
        CALL CHECKSUM("RAbIj",SCR(I040),NSIZE,S)
#endif 
         IRREPX = 1
         IOFFF  = I040
	 IOFFT  = I050

         DO IRREPR = 1, NIRREP
            IRREPL = DIRPRD(IRREPR,IRREPX)

            DISSZF=IRPDPD(IRREPL,13)
            NUMDSF=IRPDPD(IRREPR,14)


            CALL ASSYM2(IRREPR, POP(1,1), DISSZF, SCR(IOFFF))

            NUMDST=IRPDPD(IRREPR,3)

#ifdef _DEBUG_LVLM
            CALL CHECKSUM("R(AB;I<J)",SCR(IOFFF),NUMDST*DISSZF,S)
#endif 
            IOFFF = IOFFF + NUMDSF*DISSZF

         ENDDO 

         IOFFF  = I040
         DO IRREPR = 1, NIRREP
            IRREPL = DIRPRD(IRREPR,IRREPX)

            DISSZT=IRPDPD(IRREPL,1)
            NUMDST=IRPDPD(IRREPR,3)

            DISSZF=IRPDPD(IRREPL,13)
            NUMDSF=IRPDPD(IRREPR,14)

            ITMP1 = IOFFF 
            ITMP2 = IOFFT

            DO IDIS = 1, NUMDST

               CALL SQSYM(IRREPL,VRT(1,1),DISSZT,DISSZF,1,SCR(ITMP2),
     &                    SCR(ITMP1))
#ifdef _DEBUG_LVLM
               CALL CHECKSUM("R(A<B;I<J)",SCR(ITMP2),DISSZT,S)
#endif 

               ITMP1 = ITMP1 + DISSZF 
               ITMP2 = ITMP2 + DISSZT 
            ENDDO

            IOFFF = IOFFF + NUMDSF*DISSZF 

        ENDDO 
c
C NSIZE = R2(A<B; I<J) length
C
        NSIZE = ISYMSZ(ISYTYP(1,44), ISYTYP(2,44))

        CALL DCOPY(NSIZE, SCR(I050), 1, SCR(I030), 1)

        ENDIF 
C 
         CALL GETALL(SCR(I040),
     &               ISYMSZ(ISYTYP(1,46+NLIST), ISYTYP(2,46+NLIST)),
     &               1, 46+NLIST)
         CALL DMPRGSS(SCR(I000), I050-I000, NAME)
      ENDIF
C
      RETURN
      END
