#include "flags.h"
       SUBROUTINE DRIVE_CORE_TDA_WINDOW(SCR,MAXCOR,IRREPX,IUHF)
C
C The core EEs are obtained by setting the R(ij,ab) = 0 when i and
C j do  not belong to the core region. The present version works
C state specifc fashion.
c
C                           all; each
C LIST 444:    C(IJ,AB )       A<B ; I<J     AA AA
C      445:    C(ij,ab )       a<b ; i<j     BB BB
C      446:    C(Ij,Ab )       A,b ; I,j     AB A
C 
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION SCR(MAXCOR)
C
#include "maxbasfn.par"
C
      LOGICAL PROJECT_SINGLES
      CHARACTER*12 STRINGI(2), STRINGJ(2), STRINGA(2),STRINGB(2)
      DIMENSION LS2OUT(2,2), IMAP_A(MAXBASFN),IMAP_B(MAXBASFN) 
      DIMENSION OCA_OFF(8),OCB_OFF(8)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/EXTRAP/MAXEXP,NREDUCE,NTOL,NSIZEC
      COMMON/ROOTS/EIGVAL(100,8),EIGVAL_T(100,8),OSCSTR(100,8),
     &             BGN(100,8),BGN_IRP(100,8),END(100,8),
     &             END_IRP(100,8),NATURE
      COMMON/PROJECT/IPROJECT, IPATTERN, NCALC, ICALC, IWINDOW(8)
#include "info.com"
#include "flags.com"
#include "flags2.com"
C
      PROJECT_SINGLES = .TRUE.

#ifdef _DEBUG_LVLM
      Write(6,"(a,a)") "  Using the core projection for", 
     &                 " core-window option for TDA"
      Write(6,"(a,8i2)") "  The core window: ", 
     &                     (iwindow(i),i=1,NirreP)
#endif 
C
      LS1OUT      = 490
     
      NBAS   = NOCCO(1) + NVRTO(1)
      CALL GETREC(20,"JOBARC","ORBMAP_A",NBAS,IMAP_A)
      CALL GETREC(20,"JOBARC","ORBMAP_A",NBAS,IMAP_B)
      IF (IUHF.NE.0) CALL GETREC(20,"JOBARC","ORBMAP_B",NBAS,
     &                           IMAP_B)
      
      OCA_OFF(1) = 0
      OCB_OFF(1) = 0

      DO I = 1, (NIRREP-1)
         OCA_OFF(I+1) = OCA_OFF(I) + POP(I,1)
         OCB_OFF(I+1) = OCB_OFF(I) + POP(I,2)
      ENDDO

      CALL ZERO(SCR,NSIZEC)
c
C  FIRST CONSIDER SINGLE EXCITATION COEFFICIENTS
C
      IF (PROJECT_SINGLES) THEN

      dO 500 SSPIN = 1, 1+IUHF
        CALL GETLST(SCR, 1, 1, 1, SSPIN, LS1OUT)
        ICOUNT = 1
        DO IIRREP = 1, NIRREP
          AIRREP=DIRPRD(IIRREP,IRREPX)
          CORE_WINDOW = IWINDOW(IIRREP)
          DO 1 I = 1, POP(IIRREP, SSPIN)
            DO 2 A = 1, VRT(AIRREP,SSPIN)

                 IF (SSPIN .EQ. 1) THEN
                    II = IMAP_A(I+OCA_OFF(IIRREP))
                 ELSE
                    II = IMAP_B(I+OCB_OFF(IIRREP))
                 ENDIF
                 IF (I .LE. CORE_WINDOW) THEN
                    SCR(ICOUNT) = 1.0D0
                 ELSE 
                    SCR(ICOUNT) = 0.0D0
                 ENDIF 

#ifdef _DEBUG_LVLM
      Write(6,"(1x,a,5i3,F15.8)") "PH: ",SSPIN,IIRREP,AIRREP,I,A,
     +                             SCR(ICOUNT) 
#endif 
              ICOUNT = ICOUNT + 1
    2       CONTINUE
    1     CONTINUE
        ENDDO

        CALL PUTLST(SCR, 1, 1, 1, SSPIN, LS1OUT)
  500 CONTINUE
C
      ENDIF 
C
      RETURN
      END

