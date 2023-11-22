      SUBROUTINE GETDENMETHDS(IUHF)
C
      IMPLICIT INTEGER(A-Z)
C
      CHARACTER*20 NAME
      CHARACTER*11 NAME1
      CHARACTER*3 NAME2
      CHARACTER*4 NAME3
      CHARACTER*5 NAME4
      CHARACTER*8 NAME5
      CHARACTER*12 NAME6
      CHARACTER*6 NAME7
      LOGICAL MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC
      LOGICAL DENS,GRAD,QRHF,NONHF,ROHF,SEMI,CANON,TRIP1,TRIP2
      LOGICAL CIS
C
      COMMON/EXCITE/CIS
      COMMON/METH/MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC
      COMMON/DERIV/DENS,GRAD,QRHF,NONHF,ROHF,SEMI,CANON,TRIP1,
     &             TRIP2
      COMMON/FLAGS/IFLAGS(100)
      COMMON/FLAGS2/IFLAGS2(500)
C
      EQUIVALENCE(METHOD,IFLAGS(2))
C
3000  FORMAT('  ',A7,' density and intermediates are calculated.')
3001  FORMAT('  ',A11,' density and intermediates are calculated.')
3002  FORMAT('  ',A,' density and intermediates are calculated.')
3003  FORMAT('  ',A4,' density and intermediates are calculated.')
3004  FORMAT('  ',A5,' density and intermediates are calculated.')
3010  FORMAT('  ',A8,' density and intermediates are calculated.')
3011  FORMAT('  ',A12,' density and intermediates are calculated.')
3012  FORMAT('  ',A6,' density and intermediates are calculated.')
3005  FORMAT('  The reference state is a QRHF wave functions.')
3006  FORMAT('  The reference state is a non HF wave functions.')
3007  FORMAT('  The reference state is a ROHF wave function.')
3008  FORMAT('  Semi-canonical orbitals are used.')
3009  FORMAT('  The perturbed orbitals are chosen canonical.')
3013  FORMAT(' ',A4,' unrelaxed density is evaluated.') 
C     
      MBPT2  = .FALSE.
      MBPT3  = .FALSE.
      M4DQ   = .FALSE.
      M4SDQ  = .FALSE.
      M4SDTQ = .FALSE.
      CCD    = .FALSE.
      QCISD  = .FALSE.
      CCSD   = .FALSE.
      UCC    = .FALSE.
      NONHF  = .FALSE.
      ROHF   = .FALSE.
      SEMI   = .FALSE.
      CANON  = .FALSE.
      TRIP1  = .FALSE.
      TRIP2  = .FALSE.
C     
      IF (METHOD .EQ. 0) THEN
         IF(IFLAGS(87) .EQ. 1)THEN
            CIS   = .TRUE.
            NAME2 = 'TDA'
            write (6,3002) NAME2
         ENDIF
C
      ELSEIF (METHOD .EQ. 1) THEN
         MBPT2 = .TRUE.
         NAME  = 'MBPT(2)'
         write(6,3000) NAME
C
      ELSE IF(METHOD .EQ. 2) THEN
         MBPT3 = .TRUE.
         NAME  = 'MBPT(3)'
         write(6,3000)NAME
C
      ELSE IF(METHOD .EQ. 3) THEN
         M4SDQ = .TRUE.
         NAME1 = 'SDQ-MBPT(4)'
         write(6,3001) NAME1   
C
      ELSE IF(METHOD .EQ. 4) THEN
         M4SDTQ = .TRUE.
         TRIP1  = .TRUE.
         NAME = 'MBPT(4)'
         write(6,3000)NAME
C
      ELSE IF(METHOD .EQ. 7) THEN
         UCC = .TRUE.
         NAME5 = 'UCCSD(4)'
         write(6,3010) NAME5
C
      ELSE IF(METHOD .EQ. 8) THEN
         CCD   = .TRUE.
         NAME2 = 'CCD'
         write(6,3002)NAME2
C
      ELSE IF(METHOD .EQ. 9) THEN
         UCC   = .TRUE.
         TRIP1 = .TRUE.
         NAME7 = 'UCC(4)'
         write(6,3012)NAME7
C
      ELSE IF (METHOD .EQ. 10) THEN
         CCSD = .TRUE.
         NAME3 = 'CCSD'
         write(6,3013)NAME3
C
      ELSE IF (METHOD .EQ. 11) THEN
         CCSD  = .TRUE.
         NAME6 = 'CCSD+T(CCSDT)'
         TRIP1 = .TRUE.
         write(6,3011) NAME6
C
      ELSE IF(METHOD.EQ.21) THEN
         QCISD = .TRUE.
         NAME5 = 'QCISD(T)'
         TRIP1 = .TRUE.
         TRIP2 = .TRUE.
         write(6,3010) NAME5
C
      ELSE IF(METHOD.EQ.22) THEN
         CCSD  = .TRUE.
         NAME  = 'CCSD(T)'
         TRIP1 = .TRUE.
         TRIP2 = .TRUE.
         write(6,3000) NAME
C
      ELSE IF(METHOD.EQ.23) THEN
         QCISD = .TRUE.
         NAME4 = 'QCISD'
         write(6,3004)NAME4
      ENDIF
C     
C Sets flag for derivative calculations
C     
C     
      DENS = .TRUE.
C     
C Calcualte Gradients with respect to nuclear coordinates.
C I intermediates are required.
C     
      GRAD = .TRUE.
C
      IF (IFLAGS(77) .NE. 0) THEN
C  
         QRHF = .TRUE.
         CALL QRHFSET
C     
         WRITE(6,3005)
C
      ENDIF
C     
C Check we are delaing with non-HF cases.
C
      IF (IFLAGS(38) .NE.0. AND. IFLAGS(11) .NE. 2) THEN
C     
         NONHF = .TRUE.
C     
         WRITE(6,3006)
C     
      ENDIF
C
      IF (IFLAGS(11).EQ.2) THEN
C     
         ROHF = .TRUE.
         WRITE(6,3007)
         IUHF = 1
C
         IF (IFLAGS(39) .EQ. 1) THEN
            SEMI = .TRUE.
            WRITE(6,3008)
         ENDIF
C
      ENDIF
C     
C Check if perturbed canonical orbitals are used.
C     
      IF (IFLAGS(64) .EQ. 1) THEN
         CANON = .TRUE.
         WRITE(6,3009)
      ENDIF
C     
      RETURN
      END
