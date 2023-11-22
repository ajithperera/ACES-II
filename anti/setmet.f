      
      SUBROUTINE SETMET
      IMPLICIT INTEGER(A-Z)
      CHARACTER*7 NAME
      CHARACTER*11 NAME1
      CHARACTER*4 NAME2
      CHARACTER*5 NAME3
      CHARACTER*6 NAME4
      LOGICAL MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD
      LOGICAL DENS,GRAD,NONHF,ROHF,MOEQAO
      LOGICAL TDA,EOM,GABCD,DABCD
C	
      COMMON/EXCITE/TDA,EOM
      COMMON/METH/MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD
      COMMON/ABCD/GABCD,DABCD
      COMMON/IABCD/IOFFABCD
      COMMON/DERIV/DENS,GRAD,NONHF
      COMMON/REF/ROHF
      COMMON/FLAGS/IFLAGS(100)
      COMMON/FLAGS2/IFLAGS2(500)
      COMMON/SIZES/MOEQAO
      COMMON/STATSYM/IRREPX
C
      EQUIVALENCE(METHOD,IFLAGS(2))
C
      IONE=1
C
3000  FORMAT(' ',A7,' MO gammas will be sorted to Mulliken order.')
3001  FORMAT(' ',A11,' MO gammas will be sorted to Mulliken order.')
3002  FORMAT(' ',A3,' MO gammas will be sorted to Mulliken order.')
3003  FORMAT(' ',A4,' MO gammas will be sorted to Mulliken order.')
3004  FORMAT(' ',A5,' MO gammas will be sorted to Mulliken order.')
3005  FORMAT(' The reference state is a non HF wave functions.')
3020  FORMAT('  EOM-',A4,' density and intermediates are calculated.')
 3021 FORMAT('  ACES3 density and intermediates are calculated.')
C
      MBPT2=.FALSE.
      MBPT3=.FALSE.
      M4DQ=.FALSE.
      M4SDQ=.FALSE.
      M4SDTQ=.FALSE.
      CCD=.FALSE.
      QCISD=.FALSE.
      CCSD=.FALSE.
      NONHF=.FALSE.
      ROHF=.FALSE.
      TDA=.FALSE.
      EOM=.FALSE.
      GABCD=.FALSE.
      DABCD=.FALSE.
C
      IF (IFLAGS2(132).EQ.3 ) THEN
        CCSD=.TRUE.
        write(6,3021)
        EOM=.TRUE.
      ELSEIF(METHOD.EQ.0) THEN
       IF(IFLAGS(87).EQ.1)THEN
        TDA=.TRUE.
        NAME2='TDA'
        write(6,3002)NAME2
       ELSE
        CALL ERREX
       ENDIF
      ELSEIF (IFLAGS(87).EQ.3 .OR. IFLAGS(87).EQ.7) THEN
        CCSD=.TRUE.
        NAME3='CCSD'
        write(6,3020)NAME3
        EOM=.TRUE.
      ELSEIF(METHOD.EQ.1) THEN
        call aces_fin 
        STOP
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
      ELSE IF(METHOD.EQ.8.OR.METHOD.EQ.5) THEN
       CCD=.TRUE.
       NAME2='CCD'
       write(6,3002)NAME2
C If we don't include CCSD(T) the ABCD_TYPE=AOBASIS CCSD(T)
C gradients incorrect, 09/2019, Ajith Perera 
      ELSE IF(METHOD.EQ.10 .OR. METHOD .EQ.22) THEN
       CCSD=.TRUE.
       NAME3='CCSD'
       write(6,3003)NAME3
      ELSE IF(METHOD.EQ.23) THEN
       QCISD=.TRUE.
       NAME4='QCISD'
       write(6,3004)NAME4
      ENDIF
C
C SET FLAGS FOR DERIVATIVE CALCULATION
C
C  CALCULATE DENSITY MATRIX
C
      DENS=.TRUE.
C
C  CALCULATE GRADIENTS WITh RESPECT TO NUCLEAR COORDINATES
C  (THIS MEANS THE INTERMEDIATE I IS REQUIRED
C
      GRAD=.TRUE.
C
C CHECK IF WE ARE DEALING WITH A QRHF CASE
C
      IF(IFLAGS(77).NE.0)THEN
C  
       NONHF=.TRUE.
C
       WRITE(6,3005)
C
      ENDIF
C
C CHECK WHETHER ROHF REFERENCE FUNCTION
C
      IF(IFLAGS(11).EQ.2) THEN
C
       ROHF=.TRUE.
C
C ROHF MBPT3 REQUIRED FOURTH-ORDER LOGIC
C
       IF(MBPT3) THEN
        MBPT3=.FALSE.
        M4SDQ=.TRUE.
       ENDIF
C
      ENDIF
C
C CHECK TO SEE IF AO AND MO SIZES ARE THE SAME
C
      CALL GETREC(20,'JOBARC','NBASTOT ',IONE,NMO)
      CALL GETREC(20,'JOBARC','NAOBASFN',IONE,NAO)
      MOEQAO=NMO.EQ.NAO
      IF(.NOT.MOEQAO)THEN
       WRITE(6,4000)
4000   FORMAT(T3,'@SETMET-I, Unequal number of MO and AO orbitals.')
      ENDIF
C
C CHECK WHETHER G(Ab,CD) HAS BEEN CALCULATED OR IS RECALCULATED ON
C THE FLY HERE IN ANTI
C
C IFLAGS(100) corresponds to gamma_abcd with options {disk,direct}
C The default is set to disk and GABCD is true unless the user
C choose direct option. If the user choose direct option then
C G(AB,CD) is not on the disk and has to be made on the fly.
C This is controlled by DABCD (TRUE) and accomplished in sg1aaaa.f
C and sg1aaab.f. Ajith Perera,  05/2005.
C
      GABCD=IFLAGS(100).EQ.0
      IF(.NOT.GABCD)THEN
       DABCD=.TRUE.
       IF(MBPT3.OR.M4SDQ.OR.M4SDTQ) CALL ERREX
      ELSE
       IOFFABCD=0
      ENDIF 
C
      IF(EOM.AND.DABCD) THEN
C
C GET FOR THE DIRECT EVALUATION OF G(AB,CD) SYMMETRY OF EXCITED STATE
C
       IONE=1
       CALL GETREC(20,'JOBARC','STATESYM',IONE,IRREPX)
C
      
      ENDIF
C
      RETURN
      END
