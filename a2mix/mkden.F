      SUBROUTINE MKDEN(COEF, DENS, SCFDEN, RELDEN, FOCKMAT, SPH2CART, 
     &                 TMP, NBAS, NBASP, NOCC, IFLG, ONEP_ENERGY)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION COEF(NBASP*NBASP), DENS(NBASP*NBASP),SCFDEN(NBAS*NBAS), 
     &          NOCC(16), RELDEN(NBAS*NBAS), FOCKMAT(NBASP, NBASP),
     &          SPH2CART(NBAS*NBASP), TMP(NBAS*NBAS)
C
      LOGICAL SPHERICAL
      CHARACTER*8 LABELPSCF, LABELDSCF, LABELDENS, LABELSDEN
      CHARACTER*4 LABELF
      CHARACTER*1 SPIN(2)
C
      COMMON /FLAGS/ IFLAGS(100)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /IPAR/ LUOUT
      COMMON /ISYMINF/ NIRREP,NSOIRP(8)
C
      DATA SPIN /'A', 'B'/
C
      SPHERICAL = (IFLAGS(62) .EQ. 1)
      NBAS2=NBASP*NBASP
C
      LABELF   ='FOCK'
      LABELPSCF='HFDENSTY'
      LABELDSCF='HDDENSTY'
      LABELDENS='TDENSITY'
      LABELSDEN='DDENSITY'
C
      CALL ZERO(DENS,   NBASP*NBASP)
      CALL ZERO(SCFDEN, NBAS*NBAS)
      CALL ZERO(RELDEN, NBAS*NBAS)
C     
      IF(IFLG.EQ.1) THEN
C
C     Get the matrix of alpha MO's
      CALL GETREC(20,'JOBARC','SCFEVCA0',IINTFP*NBAS2,COEF)
      ISOFF=0
      ONEP_ENERGY = 0.0D0
C
      ELSE
C     Get the matrix of alpha MO's in terms of AO's.
      CALL GETREC(20,'JOBARC','SCFEVCB0',IINTFP*NBAS2,COEF)
      ISOFF=8
C
      ENDIF
C
      IOFF=0
      IOFF2=0
      DO 10 ISYM=1,NIRREP
         NOCCT=NOCC(ISOFF+ISYM)
         DO 20 IOCC=1,NOCCT
            DO 30 I=1,NSOIRP(ISYM)
               DO 40 J=1,NSOIRP(ISYM)
C
C  Total density
C
      DENS((IOFF+I-1)*NBASP+IOFF+J)=DENS((IOFF+I-1)*NBASP+IOFF+J)+
     &                              COEF((IOFF+IOCC-1)*NBASP+IOFF+I)*
     &                              COEF((IOFF+IOCC-1)*NBASP+IOFF+J)
C
   40          CONTINUE
 30         CONTINUE
   20    CONTINUE
         IOFF=IOFF+NSOIRP(ISYM)
         IOFF2=IOFF2+NOCCT*NBAS2
   10 CONTINUE
#ifdef _DEBUG_LVLM
      Write(6,"(a)") "The HF density in MO basis"
      Write(6,"(5(F12.6))")(Dens(i), i=1,NBAS2)
#endif 
C
C Get the SCF and correlated (relaxed) density matirices in AO basis. For 
C open shell calculations, we also read the density diferences. Note that
C the "relaxed" density read here only contain the correlated contributions.
C
      CALL GETREC(-20, 'JOBARC', LABELPSCF, NBAS2*IINTFP, SCFDEN)
#ifdef _DEBUG_LVLM
      Write(6,"(a)") "The HF density"
      CALL OUTPUT(SCFDEN, 1, NBAS, 1, NBASP, NBAS, NBASP, 1)
#endif 

      CALL GETREC(0, 'JOBARC', LABELDENS, LENGTH, RELDEN)
      IF (LENGTH .GT. 0) THEN
          CALL GETREC(-20, 'JOBARC', LABELDENS, NBAS2*IINTFP, RELDEN)
      ENDIF 

#ifdef _DEBUG_LVLM
      if (lenght .gt.0) then
      Write(6,"(a)") "The total (SCF+Corr.) density"
      CALL OUTPUT(RELDEN, 1, NBAS, 1, NBASP, NBAS, NBASP, 1)
      endif 
#endif
      CALL GETREC(-20, 'JOBARC', LABELF//SPIN(IFLG)//'   ', 
     &            NBAS2*IINTFP, FOCKMAT)
C      
      IF (IFLG .EQ. 2) THEN
         CALL GETREC(-20, 'JOBARC', LABELDSCF, NBAS2*IINTFP, SCFDEN)
         CALL GETREC(0, 'JOBARC', LABELSDEN, LENGTH, RELDEN)
         IF (LENGTH .GT. 0) THEN
            CALL GETREC(-20, 'JOBARC', LABELSDEN, NBAS2*IINTFP, 
     &                   RELDEN)
         ENDIF 
         CALL GETREC(-20, 'JOBARC', LABELF//SPIN(IFLG)//'   ', 
     &               NBAS2*IINTFP, FOCKMAT)
      ENDIF
C      
      CALL SAXPY(NBAS2, -1.0D0, SCFDEN, 1, RELDEN, 1)
      ONEP_ENERGY= ONEP_ENERGY + SDOT(NBAS2, RELDEN, 1, FOCKMAT, 1)
C
C$$$      CALL GETREC(20, 'JOBARC', 'ONEHAO  ', NBAS2*IINTFP, FOCKMAT, 1)
C$$$      Write(6,*) "The One Hamiltonian matrix in AO basis"
C$$$      CALL OUTPUT(FOCKMAT, 1, NBAS, 1, NBASP, NBAS, NBASP, 1)

#ifdef _DEBUG_LVLM
      Write(6,*) "The Density matrix in AO basis"
      CALL OUTPUT(RELDEN, 1, NBAS, 1, NBASP, NBAS, NBASP, 1)
#endif 

C$$$
C$$$      ONEHAMIL =  SDOT(NBAS2, RELDEN, 1, FOCKMAT, 1)
C$$$      Write(6,*) 'one hamil energy =', onehamil
C
C Get back the total density matrix into RELDEN array
C      
      IF (LENGTH .GT. 0) THEN
         CALL GETREC(-20, 'JOBARC', LABELDENS, NBAS2*IINTFP, RELDEN)
      ENDIF 
C
C At this point we need to transform the density matrices to Cartesian
C basis in order to directly contract with the Cartesian AO integrals.
C 
C$$$      IF (SPHERICAL) THEN
C$$$         CALL GETREC(20, 'JOBARC','CMP2CART', NBAS*NBASP*IINTFP, 
C$$$     &               SPH2CART)

C$$$      Write(6,*) "The spherical to Cartesian transformation matrix"
C$$$      CALL OUTPUT(SPH2CART, 1, NBAS, 1, NBASP, NBAS, NBASP, 1)
C
C$$$         CALL TRANS_SPH2CART(SCFDEN, SPH2CART, TMP, NBAS, NBASP)
C$$$         CALL TRANS_SPH2CART(RELDEN, SPH2CART, TMP, NBAS, NBASP)
C$$$      ENDIF
C
C
      RETURN
      END


