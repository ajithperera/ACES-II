
#include "flags.h"

      SUBROUTINE SETMET(pCCD,pCCDS,pCCDTS,pCCDTSD)
      IMPLICIT INTEGER(A-Z)
      CHARACTER*7 NAME
      CHARACTER*11 NAME1
      CHARACTER*4 NAME2
      CHARACTER*5 NAME3
      CHARACTER*6 NAME4
      CHARACTER*12 NAME5
      CHARACTER*8 NAME6,NAME8
      LOGICAL ROHF
      LOGICAL MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC
      LOGICAL DENS,GRAD,IHBAR,TRIPNI1,TRIPNI2,TRIPIT,T3STOR,
     &        IHBAR2,IHBAR3,CC2
      LOGICAL pCCD,pCCDS,pCCDTS,pCCDTSD
      COMMON/METH/MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC,CC2
      COMMON/TRIPLES/TRIPNI1,TRIPNI2,TRIPIT,T3STOR
      COMMON/DERIV/DENS,GRAD
      COMMON/REF/ROHF
      COMMON/HBAR/IHBAR,IHBAR2,IHBAR3
      COMMON/FLAGS/IFLAGS(100)
      COMMON/FLAGS2/IFLAGS2(500)
C
      EQUIVALENCE(METHOD,IFLAGS(2))
C
3000  FORMAT('  The Lambda equations are solved for ',A7,'.')
3001  FORMAT('  The Lambda equations are solved for ',A11,'.')
3002  FORMAT('  The Lambda equations are solved for ',A3,'.')
3003  FORMAT('  The Lambda equations are solved for ',A4,'.')
3004  FORMAT('  The Lambda equations are solved for ',A5,'.')
3025  FORMAT('  The Lambda equations are solved for ',A8,'.')
3026  FORMAT('  The Lambda equations are solved for ',A9,'.')
3005  FORMAT('  Hbar is calculated for a CCSD wave function.')
3015  FORMAT('  Hbar is calculated for a CC2 wave function.')
3020  FORMAT('  Hbar is calculated for a CCD wave function.')
3021  FORMAT('  Hbar is calculated for a pCCD wave function.')
3006  FORMAT('  The Lambda equations are solved for ',A12,'.')
3007  FORMAT('  The Lambda equations are solved for ',A8)
3008  FORMAT('  Hbar and lambda solutions will be evaluated.')
3009  FORMAT('  Hbar is calculated through second order in ',
     &       'correlation.')
3018  FORMAT('  Hbar and lambda is calculated through CCD in ',
     &       'correlation.')
3017  FORMAT('  Hbar is calculated through RPA in ',
     &       'correlation (RPA as reference).')
3016  FORMAT('  Hbar is calculated through DRPA in ',
     &       'correlation (dRPA as reference)')
3050  FORMAT('  Hbar is calculated through LCCD in ',
     &       'correlation (LCCD as reference)')
3051  FORMAT('  Hbar is calculated through LCCSD in ',
     &       'correlation (LCCSD as reference)')
C
CJDW  6/16/95. Add PROP=J_FC, PROP=J_SD, PROP=J_SO to IHBAR2 logic.
CJDW 10/25/95. Add PROP=JSC_ALL (PROP=13).
C IFLAGS(43):HBAR 
C IFLAGS(18):PROPS (Response properties)
C IFLASG(54):VIB (excited state optimizations and frequencies)
C IFLAGS(91):excited state properties as an expectation
C IFLAGS2(117):EOMREF 
C IFLGAS2(151):Resonance Raman 
      
      IHBAR  = .FALSE. 
      IHBAR2 = .FALSE.

      IHBAR=IFLAGS(43).NE.0 
      IHBAR2=IHBAR.AND.(IFLAGS(91).NE.0.OR.IFLAGS(18).GT.100
     &                                 .OR.IFLAGS(18).EQ. 8
     &                                 .OR.IFLAGS(18).EQ. 9
     &                                 .OR.IFLAGS(18).EQ.10
     &                                 .OR.IFLAGS(18).EQ.13
     &                                 .OR.IFLAGS(18).EQ.11)
     &       .AND.IFLAGS2(h_IFLAGS2_geom_opt).EQ.0.AND.IFLAGS(54).EQ.0
     &       .AND.IFLAGS2(3).EQ.0
CMN
C At this stage we do not solve lambda for RPA and DRPA. So
C There are no transition moments. If I am convinced that RPA
C is variational then we do not need it at all. 
C Ajith Perera, 07/2013.
C
      IF(IFLAGS2(117) .EQ. 2 .OR. IFLAGS2(117) .EQ. 4 .OR.
     &   IFLAGS2(117) .EQ. 5) 
     &   IHBAR2 = .FALSE.
C 
C Also, do not solve lamada equations of LCCD and LCCSD based
C EOM calculations (for that matter for gradients). We can
C also use the IFLAGS2(117) 6 and 7 to identify LCCD and 
C LCCSD.
C
      IF (IFLAGS(2) .EQ. 5 .OR. IFLAGS(2) .EQ. 6) IHBAR2 = .FALSE.
 
      IF (IFLAGS2(151) .EQ. 1) IHBAR2= .TRUE.
C
CMN END
C The following is simply to avoid doing the same thing twice. 
C
      IF(IHBAR2)IHBAR=.FALSE.

C      IHBAR = .FALSE.
C      IHBAR2 = .TRUE.
C
      MBPT2=.FALSE.
      MBPT3=.FALSE.
      M4DQ=.FALSE.
      M4SDQ=.FALSE.
      M4SDTQ=.FALSE.
      CCD=.FALSE.
      QCISD=.FALSE.
      CCSD=.FALSE.
      TRIPNI1=.FALSE.
      TRIPNI2=.FALSE.
      TRIPIT=.FALSE.
      ROHF=.FALSE.
      CC2 =.FALSE.
C
C Extensions to CCD, RPA and dRPA added (also see above). 
C Ajith Perera, 07/2013 (I am not sure why Marcel set CCSD
C is true for EOMREF=MBPT(2). It may be simply to get 
C lamabd code to make W(MB,EJ) intermediates. Also, added 
C LCCD and LCCSD (IFLAGS2(117,EOMREF)=6 and 7).

      IF(IHBAR) THEN

        IF(IFLAGS2(117) .EQ. 2) THEN
           WRITE(6,3009)
        ELSE IF (IFLAGS2(117) .EQ. 4) THEN 
           WRITE(6,3017)
        ELSE IF (IFLAGS2(117) .EQ. 5) THEN 
           WRITE(6,3016)
        ELSE IF (IFLAGS2(117) .EQ. 6) THEN 
           WRITE(6,3050)
           CCD = .TRUE. 
        ELSE IF (IFLAGS2(117) .EQ. 7) THEN 
           WRITE(6,3051)
           CCSD = .TRUE.
        ELSE 
           IF (METHOD .EQ. 47) THEN
              write(6,3015)
              CC2 = .TRUE.
           ELSE IF (METHOD .EQ. 8) THEN
               IF (pCCD.or.pCCDS.or.pCCDTS.or.pCCDTSD) THEN
                   write(6,3021)
               ELSE
                   write(6,3020)
               ENDIF 
               CCD=.TRUE.
           ELSE 
              write(6,3005)
              CCSD = .TRUE. 
           ENDIF 
           CCSD = .TRUE.
        ENDIF

      ELSEIF(IHBAR2)THEN

C METHOD = 10 (CCSD), = 5 (LCCD), = 7 (LCCSD), 4= 46 (CC2). 
C METHOD = 8  (CCD) 

       IF (METHOD .EQ.  10) THEN 
          write(6,3008)
          CCSD=.TRUE.
       ELSE IF (METHOD .EQ.5 .OR. METHOD.EQ.7 .OR. 
     &          METHOD .EQ. 47) THEN
          write(6,3008)
          CCSD=.TRUE.
          CC2 =.TRUE.
       ELSE IF (METHOD .EQ. 8) THEN
          write(6,3018)
          CCD=.TRUE.
       ENDIF
      
      ELSE 

      ROHF=IFLAGS(11).EQ.2
C
C CHECK WHETHER LAMBDA REALLY HAS TO BE CALLED
C
      IF ( METHOD.LE.1                 .OR.
     &    (METHOD.EQ.2.AND.(.NOT.ROHF)).OR.
     &     METHOD.EQ.7                 .OR.
     &     METHOD.EQ.9                
     &   ) THEN
         call aces_fin
         stop
      END IF

      IF(METHOD.EQ.2) THEN
       MBPT3=.TRUE.
       NAME='MBPT(3)'
       write(6,3000)NAME
      ELSE IF(METHOD.EQ.3) THEN
       M4SDQ=.TRUE.
       NAME1='SDQ-MBPT(4)'
       write(6,3001) NAME1   
      ELSE IF(METHOD.EQ.4) THEN
       M4SDTQ=.TRUE.
       NAME='MBPT(4)'
       TRIPNI1=.TRUE.
       write(6,3000)NAME
      ELSE IF(METHOD.EQ.8) THEN
       IF (pCCD) THEN
          write(6,3002) "pCCD"
          CCD=.TRUE.
       ELSE
          CCD=.TRUE.
          NAME2='CCD'
          write(6,3002)NAME2
       ENDIF
      ELSE IF(METHOD.EQ.10) THEN
      IF (pCCDS) THEN
          CCSD=.TRUE.
          WRITE(6,3000) "pCCD-(S)"
       ELSE IF (pCCDTS) THEN
          CCSD=.TRUE.
          WRITE(6,3000) "pCCD-(TS)"
       ELSE
          CCSD=.TRUE.
          NAME3='CCSD'
          write(6,3000)NAME3
       ENDIF
      ELSE IF(METHOD.EQ.11) THEN
       CCSD=.TRUE.
       NAME5='CCSD+T(CCSD)'
       TRIPNI1=.TRUE.
       write(6,3006) NAME5
      ELSE IF(METHOD.EQ.13) THEN
       CCSD=.TRUE.
       NAME8=' CCSDT-1.'
       IF(.NOT.IHBAR)THEN
        TRIPIT = .TRUE.
        WRITE(6,3007) NAME8
       ENDIF
      ELSE IF(METHOD.EQ.14) THEN
       CCSD=.TRUE.
       NAME8='CCSDT-1b.'
       IF(.NOT.IHBAR)THEN
        TRIPIT = .TRUE.
        WRITE(6,3007) NAME8
       ENDIF
      ELSE IF(METHOD.EQ.15) THEN
       CCSD=.TRUE.
       NAME8=' CCSDT-2.'
       IF(.NOT.IHBAR)THEN
        TRIPIT = .TRUE.
        WRITE(6,3007) NAME8
       ENDIF
      ELSE IF(METHOD.EQ.16) THEN
       CCSD=.TRUE.
       NAME8=' CCSDT-3.'
       IF(.NOT.IHBAR)THEN
        TRIPIT = .TRUE.
        WRITE(6,3007) NAME8
       ENDIF
      ELSE IF(METHOD.EQ.17) THEN
       CCSD=.TRUE.
       NAME8=' CCSDT-4.'
       IF(.NOT.IHBAR)THEN
        TRIPIT = .TRUE.
        WRITE(6,3007) NAME8
       ENDIF
      ELSE IF(METHOD.EQ.18) THEN
       CCSD=.TRUE.
       NAME8='  CCSDT.'
       IF(.NOT.IHBAR)THEN
        TRIPIT = .TRUE.
        WRITE(6,3007) NAME8
       ENDIF
      ELSE IF(METHOD.EQ.21) THEN
       QCISD=.TRUE.
       NAME6='QCISD(T).'
       TRIPNI1=.TRUE.
       TRIPNI2=.TRUE.
       write(6,3007) NAME6
      ELSE IF(METHOD.EQ.22) THEN
       CCSD=.TRUE.
       NAME='CCSD(T)'
       TRIPNI1=.TRUE.
       TRIPNI2=.TRUE.
       write(6,3000) NAME
      ELSE IF(METHOD.EQ.23) THEN
       QCISD=.TRUE.
       NAME4='QCISD'
       write(6,3004)NAME4
      ELSE IF(METHOD.EQ.33) THEN
       CCSD=.TRUE.
       NAME8='   CC3.'
       IF(.NOT.IHBAR)THEN
        TRIPIT = .TRUE.
        WRITE(6,3007) NAME8
       ENDIF
      ELSE IF (METHOD.EQ.42) THEN
        CCSD=.TRUE.
        NAME8='aCCSD(T).'
        METHOD=10
        WRITE(6,3007) NAME8
      ELSE IF (METHOD.EQ.47) THEN
        CCSD=.TRUE.
        CC2=.TRUE.
        NAME8='CC2.'
        WRITE(6,3007) NAME8
      ELSE IF (METHOD .EQ. 25) THEN
        CISD=.TRUE.
        CCSD=.TRUE.
        NAME8='CISD.'
        WRITE(6,3007) NAME8
      ENDIF
C
      IHBAR3=IFLAGS(3).EQ.2.AND.CCSD
      ENDIF

      RETURN
      END
