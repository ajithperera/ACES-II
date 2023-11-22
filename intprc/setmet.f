      SUBROUTINE SETMET
      IMPLICIT INTEGER(A-Z)
      CHARACTER*7 NAME
      CHARACTER*3 NAME0
      CHARACTER*11 NAME1
      CHARACTER*4 NAME2
      CHARACTER*5 NAME3
      CHARACTER*6 NAME4
      LOGICAL MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD
      COMMON/METH/MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD
      COMMON/ABCD/ABCDTYPE
      COMMON/FLAGS/IFLAGS(100)
C
      EQUIVALENCE(METHOD,IFLAGS(2))
C
2998  FORMAT('  ','The HF wave function is a ROHF function.')
2999  FORMAT('  ','Processing integrals for ',A3,' second derivatives.')
3000  FORMAT('  ','Processing integrals for ',A7,' calculation.')
3001  FORMAT('  ','Processing integrals for ',A11,' calculation.')
3002  FORMAT('  ','Processing integrals for ',A3,' calculation.')
3003  FORMAT('  ','Processing integrals for ',A4,' calculation.')
3004  FORMAT('  ','Processing integrals for ',A5,' calculation.')
C
      MBPT2=.FALSE.
      MBPT3=.FALSE.
      M4DQ=.FALSE.
      M4SDQ=.FALSE.
      M4SDTQ=.FALSE.
      CCD=.FALSE.
      QCISD=.FALSE.
      CCSD=.FALSE.
C
      IF(METHOD.EQ.0.AND.IFLAGS(3).EQ.2) THEN
       NAME0='SCF'
       write(6,2999) NAME0
       IF(IFLAGS(11).EQ.2) THEN 
        write(6,2998)
       ENDIF
      ELSE IF(METHOD.EQ.1) THEN
        MBPT2=.TRUE.
        NAME='MBPT(2)'
       write(6,3000)NAME
      ELSE IF(METHOD.EQ.2) THEN
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
       write(6,3000)NAME
      ELSE IF(METHOD.EQ.8) THEN
       CCD=.TRUE.
       NAME2='CCD'
       write(6,3002)NAME2
      ELSE IF(METHOD.EQ.10) THEN
       CCSD=.TRUE.
       NAME3='CCSD'
       write(6,3003)NAME3
      ELSE IF(METHOD.EQ.23) THEN
       QCISD=.TRUE.
       NAME4='QCISD'
       write(6,3004)NAME4
      ENDIF
C
      ABCDTYPE=IFLAGS(93)
C
      RETURN
      END
