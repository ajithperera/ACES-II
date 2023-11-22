      SUBROUTINE RFIELD1
C
C This routine calculates parameters for SELF CONSISTENT 
C REACTION FIELD method.
c 
c Piotr Rozyczko, QTP Spring 1995
c
      IMPLICIT NONE
      CHARACTER*8 LABEL
      CHARACTER*32 JUNK
      CHARACTER*80 FNAME
      LOGICAL REXN_FLD_IN_EXIST
      INTEGER IONE,ITWO,IDIELC,IMCHRG,IZCHRG,IFLAGS,IFLAGS2,ILENGTH
      DOUBLE PRECISION RADIUS,VOLUME,GFACT,FFACT,UNIT,
     $   EBORN,CAVIT,REPULS,RADIUSA,E1,E2,WEIGHT,DENSITY
      PARAMETER (UNIT=0.52917724924D+00)
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /FLAGS2/ IFLAGS2(500)
      COMMON /SCRF / FFACT,GFACT,EBORN,IZCHRG,E1,E2,REPULS
      IONE=1
      ITWO=2
c
      WRITE(6,*)' *** Self Consistent Reaction Field Contribution'
     $,' *** '
C
c  if CAVITY radius supplied; read it from file 'radius'. Otherwise
c  read it from JOBARC.
c 'radius' must contain cavity radius in ANGSTROMS!
c
      INQUIRE(FILE='REXN_FLD_IN',EXIST=REXN_FLD_IN_EXIST)

      IF (REXN_FLD_IN_EXIST) THEN
         OPEN(23,FILE='REXN_FLD_IN',FORM='FORMATTED', STATUS="OLD")
         REWIND(23)
         READ(23,*) RADIUSA, DENSITY, WEIGHT, IDIELC
         IF (RADIUSA .NE. -1.0D0) THEN
             RADIUS=RADIUSA/UNIT
             WRITE(6,'(a,a)')' @RFIELD-I: CAVITY AND DIELECTRIC  ',
     $            'DATA TAKEN FROM FILE REXN_FLD_IN'
         ELSE IF (WEIGHT .NE. 0.0D0 .AND. DENSITY .NE. 0.0D0) THEN
             RADIUS=(((21.0D0*WEIGHT)/(88.0D0*DENSITY))*10.0D0)
             RADIUS= (RADIUS)**(1.0D0/3.0D0)
             RADIUS=RADIUS + 0.5D0
             WRITE(6,'(a,a,F10.5)')' @RFIELD-I: Radius+0.5, computed',
     &                        ' using molar mass and the density = ', 
     &                         RADIUS
         ENDIF
      ELSE
         WRITE(6, '(a)') "REXN_FLD_IN input file is needed!"
         WRITE(6,*)
         CALL ERREX
      ENDIF
C
C      IF (RADIUS.EQ.0) THEN 
C        CALL GETREC(20,'JOBARC','CAVITY',IONE,RADIUSA)
C        WRITE(6,*)'  @RFIELD-I: CAVITY DATA CALCULATED INTERNALLY'
C        radius=radiusa
C      ENDIF 
C
      CAVIT=RADIUS*RADIUS*RADIUS
c
c calculate g factor
c
      GFACT=2.0D0*(IDIELC-1)/((2*IDIELC+1)*CAVIT)
c
c  get overall charge and the total number of protons.
c
      CALL GETREC(20,'JOBARC','NMPROTON',1,IZCHRG)
      IMCHRG=IFLAGS(28)
C
      
      IF (IFLAGS(1) .GT. 10) WRITE(6,102) IMCHRG,IZCHRG
 102  FORMAT(/,T3,'IMCHRG = ',I3,'  IZCHRG = ',I3,/)
c
c calculate the f factor for ions
c
      FFACT=(IZCHRG-IMCHRG+1.0D0)/(IZCHRG+1.0D0)
      EBORN=(IONE-IDIELC)*IMCHRG*IMCHRG/(2.0D0*IDIELC*RADIUS)
      EBORN=EBORN
C
C      IF (IFLAGS(1) .GT. 10) THEN
          WRITE(6,100)IDIELC
          WRITE(6,101)GFACT,FFACT,RADIUS,CAVIT
          WRITE(6,103)EBORN
C      ENDIF
C
 100  FORMAT(T2,'@RFIELD-I: Adding reaction field perturbation to',
     $          'the  one-electron hamiltonian.',//,T14,
     $          'Dielectric constant = ',I3)
 101  FORMAT(T14,'G FACTOR = ',F12.6,/,T14,'F FACTOR = ',
     $          F12.6,/,T14,'CAVITY RADIUS = ',F12.6,' A.',
     $          /,T14,'CAVITY VOLUME = ',F12.6,' a.u.^3')
 103  FORMAT(T14,'BORN energy = ',F12.6,' a.u.')
 190  RETURN
      END
