      SUBROUTINE FINFLD(HCORE,BUF,IBUF,IMAP,REPULS,
     &                  NBAS,NIRREP,NBFIRR,ILNBUF)
C
C THIS ROUTINE ADDS A FINITE-FIELD CONTRIBUTION TO THE ONE-ELECTRON
C HAMILTONIAN.
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*8 LABEL
      CHARACTER*32 JUNK
      CHARACTER*80 FNAME
      DIMENSION HCORE(1),BUF(ILNBUF),IBUF(ILNBUF),IMAP(1)
      DIMENSION NBFIRR(NIRREP)
      COMMON /FLAGS/ IFLAGS(100)
C
      CALL GFNAME('VPOUT   ',FNAME,ILENGTH)
      OPEN(UNIT=30,FILE=FNAME(1:ILENGTH),FORM='UNFORMATTED',
     &     STATUS='OLD')
      REWIND(30)
      CALL MAKMAP(NBAS,NIRREP,NBFIRR,IMAP)
C  
3     I23=IFLAGS(23)
      I24=IFLAGS(24)
      I25=IFLAGS(25)
      IF(IFLAGS(23).NE.0)THEN
       LABEL='     X  '
       IND=IFLAGS(23)
       IFLAGS(23)=0
      ELSEIF(IFLAGS(24).NE.0)THEN
       LABEL='     Y  '
       IND=IFLAGS(24)
       IFLAGS(24)=0
      ELSEIF(IFLAGS(25).NE.0)THEN
       LABEL='     Z  '
       IND=IFLAGS(25)
       IFLAGS(25)=0
      ELSE
       GOTO 20
      ENDIF
      FIELD=DFLOAT(IND)*1.E-6
C   
C ADD THE NUCLEAR DIPOLE CORRECTION TO THE ENERGY TO THE "REPULSION"
C  CONTRIBUTION.  NOT CORRECT, BUT THE REPULSION ENERGY IS NOT WRITTEN
C  OUT AFTER THIS POINT, SO WHO CARES?
C
      CALL SEEKLB(LABEL,IERR,0)
      BACKSPACE(30)
      READ(30)JUNK,DIPNUC
      REPULS=REPULS+FIELD*DIPNUC
C
30    WRITE(6,100)LABEL(6:6),FIELD
100   FORMAT(T3,'@FINFLD-I, Adding electric field perturbation to ',
     &          'one-electron Hamiltonian.',/,T10,'Field direction',
     &          T30,': ',A1,/,T10,'Field strength',T30,': ',F10.7,
     &          ' a.u.')
C
10    READ(30)BUF,IBUF,NUT
      IF(NUT.EQ.-1)GOTO 3
      DO 15 INDEX=1,NUT
       IOFF=IMAP(IBUF(INDEX))
       IF(IOFF.NE.0)HCORE(IOFF)=HCORE(IOFF)+BUF(INDEX)*FIELD   
15    CONTINUE
      GOTO 10
C
20    CLOSE(30,STATUS='KEEP')
      IFLAGS(23)=I23
      IFLAGS(24)=I24
      IFLAGS(25)=I25
C
      RETURN
      END
