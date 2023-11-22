      SUBROUTINE SETMET_QF(pCCD,pCCDS,pCCDTS)
      IMPLICIT INTEGER(A-Z)
      CHARACTER*9 NAME
      CHARACTER*11 NAME1
      CHARACTER*12 NAME11
      CHARACTER*4  NAME2
      CHARACTER*5  NAME3
      CHARACTER*6  NAME4
      CHARACTER*8  NAME5
      CHARACTER*8  NAME8
      CHARACTER*9  NAME9
      CHARACTER*10 NAME10
      CHARACTER*12 NAME12
      CHARACTER*13 NAME13
      CHARACTER*14 NAME14
      CHARACTER*15 NAME15
      LOGICAL MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC
      LOGICAL NONHF,LINCC,CICALC,TRIPNI,TRIPIT,TRIPNI1,T3STOR
      LOGICAL CIS,EOM,CC2,ADC2
      LOGICAL PCCD,PCCDS,PCCDTS
      COMMON/METH/MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC,
     &            CC2
      COMMON /FLAGS/   IFLAGS(100)
      COMMON /FLAGS2/  IFLAGS2(500)
      COMMON /MACHSP/  IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /LINEAR/  LINCC,CICALC
      COMMON /NHFREF/  NONHF
      COMMON /TRIPLES/ TRIPNI,TRIPNI1,TRIPIT,T3STOR
      COMMON /EXCITE/  CIS,EOM
C
      EQUIVALENCE(METHOD,IFLAGS(2))
C
 2098 FORMAT('   ',A4,' energy will be calculated.')
 2099 FORMAT('   ',A4,' energy will be calculated.')
 3000 FORMAT('   ',A,' energy will be calculated.')
 3001 FORMAT('   ',A11,' energy will be calculated.')
 3011 FORMAT('   ',A12,' energy will be calculated.')
 3002 FORMAT('   ',A3,' energy will be calculated.')
 3003 FORMAT('   ',A4,' energy will be calculated.')
 3004 FORMAT('   ',A5,' energy will be calculated.')
 3010 FORMAT('   ',A8,'  energy will be calculated.')
 3005 FORMAT('   The reference state is a ROHF wave function.')
 3006 FORMAT('   The reference state is a QRHF wave function.')
 3007 FORMAT('   The reference state is a general non-HF wave',
     &       ' function.')
 3008 FORMAT('   The reference state is a Brueckner wave function.')
 3009 FORMAT('   Sorry, not yet implemented !!!!.')
 2008 FORMAT('   ',A8 ,' energy will be calculated.')
 2009 FORMAT('   ',A9 ,' energy will be calculated.')
 2010 FORMAT('   ',A10,' energy will be calculated.')
 2011 FORMAT('   ',A11,' energy will be calculated.')
 2012 FORMAT('   ',A12,' energy will be calculated.')
 2013 FORMAT('   ',A13,' energy will be calculated.')
 2014 FORMAT('   ',A14,' energy will be calculated.')
 2015 FORMAT('   ',A15,' energy will be calculated.')
C
      MBPT2=.FALSE.
      MBPT3=.FALSE.
      M4DQ=.FALSE.
      M4SDQ=.FALSE.
      M4SDTQ=.FALSE.
      CCD=.FALSE.
      CC2=.FALSE.
      QCISD=.FALSE.
      CCSD=.FALSE.
      UCC=.FALSE.
      NONHF=.FALSE.
      TRIPIT=.FALSE.
      TRIPNI=.FALSE.
      TRIPNI1=.FALSE.
C
      IF(IFLAGS2(117).EQ. 2) METHOD = 1
C
      LINCC=.FALSE.
      CICALC=.FALSE.
C
      IF(METHOD.EQ.1) THEN
        IF (IFLAGS2(117).EQ.4) THEN
        NAME ='RPA'
        Write(6,2098) NAME 
        ELSE IF  (IFLAGS2(117).EQ.5) THEN
        NAME ='DRPA'
        Write(6,2099) NAME 
        ELSE
        MBPT2=.TRUE.
        NAME='MBPT(2)'
        write(6,3000)NAME
       ENDIF
       IF (IFLAGS2(117) .EQ. 2) CCSD = .TRUE.
      ELSE IF(METHOD.EQ.2) THEN
       MBPT3=.TRUE.
       NAME='MBPT(3)'
       write(6,3000)NAME
      ELSE IF(METHOD.EQ.3) THEN
       M4SDQ=.TRUE.
       NAME1='SDQ-MBPT(4)'
       write(6,3001) NAME1
      ELSE IF(METHOD.EQ.4) THEN
       M4SDQ=.TRUE.
       M4SDTQ=.TRUE.
       TRIPNI=.TRUE.
       NAME='MBPT(4)'
       write(6,3000)NAME
      ELSE IF(METHOD.EQ.5) THEN
       CCD=.TRUE.
       NAME3='LCCD'
       LINCC=.TRUE.
       write(6,3003)NAME3
      ELSE IF(METHOD.EQ.6) THEN
       CCSD=.TRUE.
       NAME5='LCCSD'
       LINCC=.TRUE.
       write(6,3004)NAME5
      ELSE IF(METHOD.EQ.7) THEN
       UCC=.TRUE.
       NAME='UCCSD(4)'
       write(6,3000)NAME
      ELSE IF(METHOD.EQ.8) THEN
       CCD=.TRUE.
       NAME2='CCD'
       write(6,3002)NAME2
      ELSE IF(METHOD.EQ.9) THEN
       UCC=.TRUE.
       TRIPIT=.TRUE.
       NAME='UCC(4)'
       write(6,3000)NAME
      ELSE IF(METHOD.EQ.10) THEN
        CCSD=.TRUE.
        NAME3='CCSD'
        write(6,3000)NAME3
      ELSE IF(METHOD.EQ.11) THEN
       CCSD=.TRUE.
       NAME11='CCSD+T(CCSD)'
       TRIPNI=.TRUE.
       write(6,3011) NAME11
      ELSE IF(METHOD.EQ.12) THEN
       CCSD = .TRUE.
       NAME14 ='CCSD+TQ*(CCSD)'
       TRIPNI  = .TRUE.
       TRIPNI1 = .TRUE.
       write(6,2014) NAME14
      ELSE IF(METHOD.EQ.13) THEN
       CCSD=.TRUE.
       NAME='CCSDT-1'
       TRIPIT=.TRUE.
       write(6,3000) NAME
      ELSE IF(METHOD.EQ.14) THEN
       CCSD=.TRUE.
       NAME='CCSDT-1b'
       TRIPIT=.TRUE.
       write(6,3000) NAME
      ELSE IF(METHOD.EQ.15) THEN
       CCSD=.TRUE.
       NAME='CCSDT-2'
       TRIPIT=.TRUE.
       write(6,3000) NAME
      ELSE IF(METHOD.EQ.16) THEN
       CCSD=.TRUE.
       NAME='CCSDT-3'
       TRIPIT=.TRUE.
       write(6,3000) NAME
      ELSE IF(METHOD.EQ.17) THEN
       CCSD=.TRUE.
       NAME='CCSDT-4'
       TRIPIT=.TRUE.
       write(6,3000) NAME
      ELSE IF(METHOD.EQ.18) THEN
      IF (PCCD)  THEN
          WRITE(6,3000) "pCCDQF"
          CCSD=.TRUE.
          TRIPIT=.TRUE.
       ELSE
       CCSD=.TRUE.
       NAME='CCSDT  '
       TRIPIT=.TRUE.
       write(6,3000) NAME
       ENDIF 
      ELSE IF(METHOD.EQ.21) THEN
       QCISD=.TRUE.
       NAME5='QCISD(T)'
       TRIPNI=.TRUE.
       TRIPNI1=.TRUE.
       write(6,3010) NAME5
      ELSE IF(METHOD.EQ.22) THEN
       CCSD=.TRUE.
       NAME='CCSD(T)'
       TRIPNI=.TRUE.
        TRIPNI1=.TRUE.
       write(6,3000) NAME
      ELSE IF(METHOD.EQ.23) THEN
       QCISD=.TRUE.
       NAME4='QCISD'
       write(6,3004)NAME4
      ELSEIF(METHOD.EQ.24)THEN
       CICALC=.TRUE.
       CCD=.TRUE.
       NAME2='CID'
       LINCC=.TRUE.
       write(6,3002)NAME2
       METHOD=5
      ELSEIF(METHOD.EQ.25)THEN
       CICALC=.TRUE.
       CCSD=.TRUE.
       NAME3='CISD'
       LINCC=.TRUE.
       write(6,3003)NAME3
       METHOD=6
      ELSEIF(METHOD.EQ.26)THEN
       QCISD=.TRUE.
       NAME9='QCISD(TQ)'
       TRIPNI =.TRUE.
       TRIPNI1=.TRUE.
       WRITE(6,2009) NAME9
      ELSEIF(METHOD.EQ.27)THEN
       CCSD=.TRUE.
       NAME8='CCSD(TQ)'
       TRIPNI =.TRUE.
       TRIPNI1=.TRUE.
       WRITE(6,2008) NAME8
      ELSEIF(METHOD.EQ.28)THEN
       CCSD=.TRUE.
       NAME13='CCSD+TQ(CCSD)'
       TRIPNI =.TRUE.
       TRIPNI1=.TRUE.
       WRITE(6,2013) NAME13
      ELSEIF(METHOD.EQ.29)THEN
       CCSD=.TRUE.
       NAME15='CCSDT+Q*(CCSDT)'
       TRIPIT=.TRUE.
       WRITE(6,2015) NAME15
      ELSEIF(METHOD.EQ.30)THEN
       CCSD=.TRUE.
       NAME14='CCSDT+Q(CCSDT)'
       TRIPIT=.TRUE.
       WRITE(6,2014) NAME14
      ELSEIF(METHOD.EQ.31)THEN
       CCSD=.TRUE.
       NAME8='CC5SD[T]'
       TRIPNI =.TRUE.
       TRIPNI1=.TRUE.
       WRITE(6,2008) NAME8
      ELSEIF(METHOD.EQ.32)THEN
       CCSD=.TRUE.
       NAME='CCSD-T '
       TRIPNI =.TRUE.
       TRIPNI1=.TRUE.
       WRITE(6,3000) NAME
      ELSEIF(METHOD.EQ.33)THEN
       CCSD=.TRUE.
       TRIPIT=.TRUE.
       NAME3='CC3 '
       WRITE(6,3003) NAME3
      ELSEIF(METHOD.EQ.34)THEN
       CCSD=.TRUE.
       TRIPIT=.TRUE.
       NAME10='CCSDT-T1T2'
       WRITE(6,2010) NAME10
      ELSEIF(METHOD.EQ.42)THEN
       CCSD=.TRUE.
       NAME8='aCCSD(T)'
       WRITE(6,2008) NAME8
      ELSEIF(METHOD.EQ.47) THEN
       CCSD=.TRUE.
       CC2 =.TRUE.
       NAME2="CC2 "
       write(6,3002)NAME2
      ELSEIF(METHOD.EQ.50)THEN
       CCD=.TRUE.
       NAME='CCD'
       TRIPNI=.TRUE.
       TRIPNI1=.TRUE.
       write(6,3000) NAME
      ELSEIF(METHOD.EQ.51)THEN
       CCSD=.TRUE.
       NAME="CCSD(T)"
       TRIPNI=.TRUE.
       TRIPNI1=.TRUE.
       write(6,3000) NAME
      ENDIF 
C
      IF(IFLAGS(11).EQ.2) THEN
       NONHF=.TRUE.
       write(6,3005)
      ENDIF
      IF(IFLAGS(22).NE.0) THEN
       NONHF=.TRUE.
       write(6,3008)
      ENDIF
      IF(IFLAGS(77).NE.0)THEN
       NONHF=.TRUE.
       write(6,3006)
      ENDIF
      IF(IFLAGS(38).NE.0) THEN
       NONHF=.TRUE.
       write(6,3007)
      ENDIF
CJDW/SRG 9/26/96.
C    Initialize /EXCITE/. The purpose of this common block is to pass
C    information to V1X2AA and V1X2AB. By doing this it is possible to
C    have only one version of V1X2AA and V1X2AB and to keep them in libr2.
C    This is useful for two reasons : (a) it reduces the number of
C    almost identical routines with the same name; (b) it removes a subtle
C    building difficulty on crunch (and perhaps machines we have not tried).
C    Note that xvcc only calls V1X2AA and V1X2AB in UCC calculations (EUCC
C    calls V1INX2 which calls V1X2AA, V1X2AB).
      CIS = .FALSE.
      EOM = .FALSE.

C This is warning message (just in case)
      ADC2 = .FALSE.
      ADC2 =  (IFLAGS(2) .EQ. 1 .AND. IFLAGS(87) .EQ. 3 .AND.
     &         IFLAGS2(117) .EQ. 10)
  
      IF (ADC2) THEN
          Write(*,*)
          Write(6,"(a,a)") "   The ADC(2) (EOM-MBPT(2)+Symm. SS block)",
     &                     " excitation energies are computed."
          Write(*,*)
      ENDIF
C
      RETURN
      END
