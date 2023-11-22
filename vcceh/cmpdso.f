C
C This subroutine drives the calculation of diamagnetic spin-orbit
C contribution to the NMR coupling tensor based on HF and correlated
C response density matrix. The diamagnetic contribution is second-order
C in the Hamiltonian and only requires zeroth order density. 
C
      SUBROUTINE CMPDSO(ICORE, MAXCOR, IUHF)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      LOGICAL ISOTOPES_PRESENT
      CHARACTER*80 FNAME
C
      DIMENSION LENVV(2), LENOO(2), LENVO(2)
      DIMENSION DMOM(3), DSO(3), ICORE(MAXCOR)
C
      COMMON /MACHSP/ IINTLN, IFLTLN, IINTFP, IALONE, IBITWD
      COMMON/FLAGS/IFLAGS(100)
      COMMON /INFO /  NOCCO(2),NVRTO(2)
      COMMON /FILES/ LUOUT, MOINTS
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
C
      IONE = 1
      ONE  = 1.00D+00
      ONEM = -1.00D+00

      CALL GFNAME("iSOTOPES",FNAME,LENTH)
      INQUIRE(FILE=FNAME(1:LENTH),EXIST=iSOTOPES_PRESENT)
      IUNITIS=1
      IF (iSOTOPES_PRESENT) THEN
         OPEN(UNIT=IUNITIS,FILE=FNAME(1:LENTH),FORM="FORMATTED")
      ENDIF
C
      CALL IZERO (LENOO, 2)
      CALL IZERO (LENVV, 2)
      CALL IZERO (LENVO, 2)
      NMO = NOCCO(1) + NVRTO(1)
      FACT = 2.00D+00 - DFLOAT(IUHF)
C
      CALL GETREC (20, 'JOBARC', 'NBASTOT', IONE, NAO)
      CALL GETREC (20, 'JOBARC', 'NATOMS  ', IONE, NATOM)
      CALL GETREC (20, 'JOBARC', 'NREALATM', IONE, NCENTR)
C
      DO 10 ISPIN = 1, (IUHF + 1)
         LENOO(ISPIN) = IRPDPD(1, 20+ISPIN)
         LENVV(ISPIN) = IRPDPD(1, 18+ISPIN)
         LENVO(ISPIN) = IRPDPD(1,  8+ISPIN)
 10   CONTINUE
C
C Allocate memory for ALPHA density
C
      I000 = IONE
      I010 = I000 + LENOO(1)*IINTFP
      I020 = I010 + LENVV(1)*IINTFP
      I030 = I020 + LENVO(1)*IINTFP
      I040 = I030 + LENVO(1)*IINTFP
      I050 = I040 + NAO*NAO*IINTFP
      I060 = I050 + NAO*NAO*IINTFP
      I070 = I060 + NAO*NAO*IINTFP
      I080 = I070 + NAO*NAO*IINTFP
      I090 = I080 + NAO*NAO*IINTFP
      I100 = I090 + NAO*NAO*IINTFP
      I110 = I100 + 3*NATOM*IINTFP
      I120 = I110 + 9*NCENTR*NCENTR*IINTFP
      IF (I120 .GT. MAXCOR) CALL INSMEM("cmpdso", I200, MAXCOR)
C
C Allocate memory for BETA density (UHF only)
C
      IF (IUHF .NE. 0) THEN 
         I130 = I120 + LENOO(2)*IINTFP
         I140 = I130 + LENVV(2)*IINTFP
         I150 = I140 + LENVO(2)*IINTFP
         I160 = I150 + LENVO(2)*IINTFP
         I170 = I160 + NAO*NAO*IINTFP
         I180 = I170 + NAO*NAO*IINTFP
         I190 = I180 + NAO*NAO*IINTFP
         I200 = I190 + NAO*NAO*IINTFP
         IF (I200 .GT. MAXCOR) CALL INSMEM("cmpdso", I200, MAXCOR)
      ENDIF
C    
C Read the OO, VV, OV blocks of the non-realx density
C form GAMLAM file
C
      CALL READDG (ICORE(I000), ICORE(I010), ICORE(I020),
     &             ICORE(I030), LENOO(1), LENVV(1), LENVO(1), 1)
      CALL EXPDEN (ICORE(I000), ICORE(I010), ICORE(I020),
     &             ICORE(I030), ICORE(I040), NMO, 1, 1, .FALSE.)
      CALL SYMMET2 (ICORE(I040), NMO)
      CALL MO2AO3 (ICORE(I040), ICORE(I050), ICORE(I060), ICORE(I070),
     &             NAO, NMO, 1)
      CALL GETREC (20, 'JOBARC', 'SCFDENSA', NAO*NAO*IINTFP, 
     &             ICORE(I060))
      CALL SAXPY (NAO*NAO, ONE/FACT, ICORE(I060), 1, ICORE(I050), 1)
C
      CALL IZERO (ICORE(I070), IINTFP*NAO*NAO)
      CALL SAXPY (NAO*NAO, ONE/FACT, ICORE(I060), 1, ICORE(I070), 1)
c YAU : old
c     CALL ICOPY(IINTFP*NAO*NAO,ICORE(I050),1,ICORE(I040),1)
c     CALL ICOPY(IINTFP*NAO*NAO,ICORE(I070),1,ICORE(I060),1)
c YAU : new
      CALL DCOPY(NAO*NAO,ICORE(I050),1,ICORE(I040),1)
      CALL DCOPY(NAO*NAO,ICORE(I070),1,ICORE(I060),1)
c YAU : end
C
C Put ALPHA + BETA density into ICORE(I050) and  ALPHA - BETA
C to ICORE(I040) for spin-density calculations
C
      If (IUHF .NE. 0) THEN
         CALL READDG (ICORE(I120), ICORE(I130), ICORE(I140),
     &                ICORE(I150), LENOO(2), LENVV(2), LENVO(2), 2)
         CALL EXPDEN (ICORE(I120), ICORE(I130), ICORE(I140),
     &                ICORE(I150), ICORE(I160), NMO, 1, 2,
     &               .FALSE.)
         CALL SYMMET2 (ICORE(I160), NMO)
         CALL MO2AO3 (ICORE(I160), ICORE(I170), ICORE(I180), 
     &                ICORE(I190), NAO, NMO, 2)
         CALL GETREC (20, 'JOBARC', 'SCFDENSB', NAO*NAO*IINTFP, 
     &                ICORE(I180))
         CALL SAXPY (NAO*NAO, ONE/FACT, ICORE(I180), 1, ICORE(I170), 1)
C
C Scf density
C
         CALL SAXPY (NAO*NAO, ONE/FACT, ICORE(I180), 1, ICORE(I070), 1)
         CALL SAXPY (NAO*NAO, ONEM/FACT, ICORE(I180), 1,ICORE(I060), 1)
C 
C Correlated density
C
         CALL SAXPY (NAO*NAO, ONE, ICORE(I170), 1, ICORE(I050), 1)
         CALL SAXPY (NAO*NAO, ONEM, ICORE(I170), 1, ICORE(I040), 1)
      ENDIF
C
C Now carefully open the VPOUT file
C
      OPEN (UNIT=30, FILE='VPOUT',FORM='UNFORMATTED',STATUS='OLD')
C
      CALL HEADER('PROPERTIES FROM HARTREE-FOCK DENSITY MATRIX',
     &            -1, 6)
C
      CALL IZERO (ICORE(I110), IINTFP*NCENTR*NCENTR)
C
C Diamagnetic spin orbit contribution to the J coupling
C
      IF (IFLAGS(18) .EQ. 8) THEN
         CALL DRVDSO (ICORE(I070), ICORE(I090), ICORE(I110), NMO, NAO,
     &                NCENTR, IUHF)
         WRITE(LUOUT, *)
         CALL FACTOR (ICORE(I110), 3*NCENTR, 0, 0, 0, 1, 0, -1)
      ENDIF
C
      CALL HEADER('PROPERTIES FROM UNRELAXED CORRELATED DENSITY MATRIX',
     &            -1, 6)
C
      IF (IFLAGS(18) .EQ. 8) THEN
         CALL IZERO (ICORE(I110), IINTFP*NCENTR*NCENTR)
         CALL DRVDSO (ICORE(I050), ICORE(I090), ICORE(I110),  NMO, NAO, 
     &                NCENTR, IUHF)
         WRITE(LUOUT, *)
         IF (ISOTOPES_PRESENT) THEN
            CALL FACTOR_IS(ICORE(I110), 3*NCENTR, 0, 0, 0, 1, 0, 1,
     &                     IUNITIS)
         ELSE
            CALL FACTOR (ICORE(I110), 3*NCENTR, 0, 0, 0, 1, 0, 1)
         ENDIF 
      ENDIF
C
      CLOSE(30)
      CLOSE(IUNITIS)
      WRITE(6,*)
C
      RETURN
      END
