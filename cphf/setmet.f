      SUBROUTINE SETMET
C
C THIS ROUTINE FILLS THE COMMON BLOCK /METHOD/ WHICH 
C CONTAINS INFORMATION ABOUT THE METHODS FOR WHICH 
C GRADIENTS ARE CALCULATED
C
C  IUHF ..... RHF/UHF FLAG
C  SCF ...... SCF/CORRELATED METHOD FLAG
C  NONHF .... NONHF FLAG FOR CORRELATED METHODS
C
CEND
      IMPLICIT INTEGER(A-Z) 
      LOGICAL SCF,NONHF,QRHFP,ROHF,SEMI,FIELD,GEOM,THIRD,MAGN,SPIN,
     &        TRAINV,SWITCH,TDHF,JSO,JFC,JSD,MSZ,GEERTSEN
      COMMON/METHOD/IUHF,SCF,NONHF
      COMMON/QRHF/QRHFP
      COMMON/OPENSH/ROHF,SEMI
      COMMON/INVAR/TRAINV
      COMMON/PERT2/FIELD,GEOM,THIRD,MAGN,SPIN,TDHF
      COMMON/PERT3/MSZ
      COMMON/DIAMAG/GEERTSEN 
      COMMON/MODE/SWITCH
      COMMON/JNMR/JSO,JFC,JSD
      COMMON/FLAGS/IFLAGS(100)
C
      SWITCH=.FALSE.
      FIELD=.FALSE.
      GEOM=.FALSE.
      THIRD=.FALSE.
      MAGN=.FALSE.
      SPIN=.FALSE.
      ROHF=.FALSE.
      SEMI=.FALSE.
      TRAINV=.TRUE.
      TDHF=.FALSE.
      JSO=.FALSE.
      JFC=.FALSE.
      JSD=.FALSE.
      MSZ=.FALSE.
      GEERTSEN=.FALSE.
c
c When PROPS=EOM_NLO (IFLAGS(18).EQ.1) we should be able to do 
c HF-SCF polarizabilities. 03/2006. Ajith Perera. 
c
      IF(IFLAGS(18).EQ.2 .OR. IFLAGS(18).EQ.1 .OR. 
     &   IFLAGS(18).EQ.11) FIELD=.TRUE. 
      IF(IFLAGS(18).EQ.3) MAGN=.TRUE.
      IF(IFLAGS(18).EQ.4) MAGN=.TRUE.
      IF(IFLAGS(18).EQ.4) SWITCH=.TRUE.
      IF(IFLAGS(18).EQ.5) MAGN=.TRUE.
      IF(IFLAGS(18).EQ.5) MSZ=.TRUE.
      IF(IFLAGS(18).EQ.6) SPIN=.TRUE.
      IF(IFLAGS(18).EQ.7) TDHF=.TRUE.
      IF(IFLAGS(18).EQ.8) JSO=.TRUE.
      IF(IFLAGS(18).EQ.9) JFC=.TRUE.
      IF(IFLAGS(18).EQ.10) JSD=.TRUE.
      IF(IFLAGS(18).EQ.12) MAGN=.TRUE.
      IF(IFLAGS(18).EQ.12) MSZ=.TRUE.
      IF(IFLAGS(18).EQ.12) GEERTSEN=.TRUE.
      IF((IFLAGS(73).EQ.1).OR.MAGN.OR.JFC.OR.JSD.OR.JSO) TRAINV=.FALSE.
      GEOM=(.NOT.FIELD).AND.(.NOT.MAGN).AND.(.NOT.JFC)
     &     .AND.(.NOT.JSD).AND.(.NOT.JSO)
C
C  DERIVATIVE LEVEL
C
      DERIV=IFLAGS(3)
      write(6,100)
C
C SET QRHF FLAG
C
C      QRHFP=.FALSE.
C      IF(IFLAGS(77).NE.0) QRHFP=.TRUE.
C
C  SET IUHF FLAG
C
      IUHF=IFLAGS(11)
C
C  SET SCF FLAG
C
      SCF=IFLAGS(2).EQ.0
C
C  SET NONHF FLAG
C
      NONHF=(IFLAGS(77).NE.0)
C
      IF(SCF) THEN
       IF(IUHF.EQ.0) THEN
        IF(DERIV.EQ.1) call errex
        IF(DERIV.EQ.2.AND.(GEOM)) write(6,1010)
        IF(DERIV.EQ.2.AND.FIELD) write(6,1012)
c        IF(DERIV.EQ.2.AND.THIRD) write(6,1012)
        IF(DERIV.EQ.2.AND.MAGN) write(6,1013)
        IF(DERIV.EQ.2.AND.MAGN.AND.SWITCH) write(6,7000)
        IF(DERIV.EQ.2.AND.JFC) write(6,7001)
        IF(DERIV.EQ.2.AND.JSD) write(6,7002)
        IF(DERIV.EQ.2.AND.JSO) write(6,7003)
       ELSE IF(IUHF.EQ.1) THEN
        IF(DERIV.EQ.1) call errex
        IF(DERIV.EQ.2.AND.(GEOM)) write(6,2010)
        IF(DERIV.EQ.2.AND.FIELD) write(6,2012)
c        IF(DERIV.EQ.2.AND.THIRD) write(6,2012)
        IF(DERIV.EQ.2.AND.SPIN) write(6,2014)
        IF(DERIV.EQ.2.AND.JFC) write(6,7001)
        IF(DERIV.EQ.2.AND.JSD) write(6,7002)
        IF(DERIV.EQ.2.AND.JSO) write(6,7003)
       ELSE IF(IUHF.EQ.2) THEN
        ROHF=.TRUE.
        IF(DERIV.EQ.1) call errex
        IF(DERIV.EQ.2.AND.(GEOM)) write(6,2110)
        IF(DERIV.EQ.2.AND.FIELD) write(6,2111)
        IF(DERIV.EQ.2.AND.SPIN) write(6,2113)
        IF(IFLAGS(39).NE.0) THEN
         SEMI=.TRUE.
         write(6,2212)
         call errex
        ELSE
         write(6,2213)
        ENDIF
        IUHF=1
       ELSE
        CALL ERREX
       ENDIF 
      ELSE
       IF(NONHF) THEN
        write(6,3000)
       ELSE
        IF(IUHF.EQ.0) THEN
         IF(MAGN) THEN
          write(6,4001)
         ELSE IF(JSO) THEN
          write(6,4004)
         ELSE IF(JSD) THEN
          write(6,4003)
         ELSE IF(JFC) THEN
          write(6,4002)
         ELSE
          write(6,4000)
         ENDIF
        ELSE IF(IUHF.EQ.1) THEN
         IF(MAGN) THEN
          write(6,5001)
         ELSE IF(JSO) THEN
          write(6,5004)
         ELSE IF(JSD) THEN
          write(6,5003)
         ELSE IF(JFC) THEN
          write(6,5002)
         ELSE
          write(6,5000)
         ENDIF
        ELSE IF(IUHF.EQ.2) THEN
         write(6,6000)
         ROHF=.TRUE.
         IF(IFLAGS(39).NE.0) THEN
          SEMI=.TRUE.
          write(6,2212)
          call errex
         ELSE
          write(6,2213)
         ENDIF
         IUHF=1
        ENDIF
       ENDIF
      ENDIF
100   FORMAT(' Coupled-perturbed HF (CPHF) equations ')
1010  FORMAT(' are solved for RHF hessian and dipole derivatives.')
1011  FORMAT(' are solved for RHF polarizability.')
1012  FORMAT(' are solved for RHF polarizability and ',
     &       'hyperpolarizability.')
1013  FORMAT(' are solved for RHF magnetic susceptibilities and',
     &       ' NMR shifts.')
2010  FORMAT(' are solved for UHF hessian and dipole derivatives.')
2011  FORMAT(' are solved for UHF polarizability.')
2012  FORMAT(' are solved for UHF polarizability and ',
     &       'hyperpolarizability.')
2014  FORMAT(' are solved for UHF hessian, dipole derivatives, ',
     &       ' and spin density derivatives.')
2110  FORMAT(' are solved for ROHF hessian and dipole',
     &       ' derivatives.')
2111  FORMAT(' are solved for ROHF polarizability and ',
     &       'hyperpolarizability.')
2113  FORMAT(' are solved for ROHF spin density derivatives.')
2212  FORMAT(' Semicanonical orbitals are used.')
2213  FORMAT(' Standard orbitals are used.')
3000  FORMAT(' are solved for QRHF-CCSD hessian',
     &       ' and dipole derivatives.')
4000  FORMAT(' are solved for RHF-CC/MBPT hessian',
     &       ' and dipole derivatives.')
4001  FORMAT(' are solved for RHF-CC/MBPT magnetic',
     &       ' susceptibilities and NMR shifts.')
4002  FORMAT(' are solved for the RHF-CC/MBPT',
     &       ' FC contributions to J.')
4003  FORMAT(' are solved for the RHF-CC/MBPT',
     &       ' SD contributions to J.')
4004  FORMAT(' are solved for RHF-CC/MBPT',
     &       ' SO contributions to J.')
5000  FORMAT(' are solved for UHF-CC/MBPT hessian',
     &       ' and dipole derivatives.')
5001  FORMAT(' are solved for UHF-CC/MBPT magnetic',
     &       ' susceptibilities and NMR shifts.')
5002  FORMAT(' are solved for the UHF-CC/MBPT',
     &       ' FC contributions to J.')
5003  FORMAT(' are solved for the UHF-CC/MBPT',
     &       ' SD contributions to J.')
5004  FORMAT(' are solved for the UHF-CC/MBPT',
     &       ' SO contributions to J.')
6000  FORMAT(' are solved for ROHF-CC/MBPT hessian',
     &       ' and dipole derivatives.')
7000  FORMAT(' CPHF equations are solved for the nuclear magnetic',
     &       ' moments.')
7001  FORMAT(' are solved for the Fermi-contact hamiltonian')
7002  FORMAT(' are solved for the spin-dipolar hamiltonian')
7003  FORMAT(' are solved for the spin-orbit hamiltonian')
      RETURN
      END
