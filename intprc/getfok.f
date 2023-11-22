      
      SUBROUTINE GETFOK(FOCKA,FOCKB,NSIZA,NSIZB,IUHF)
C
C RETRIEVES OFF-DIAGONAL PART OF FOCK MATRIX FROM HF2 AND PUTS IT OUT
C  ON JOBARC.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*80 FNAME
      CHARACTER*8 HF2NAME
      DIMENSION FOCKA(NSIZA),FOCKB(NSIZB)
      IF(IUHF.EQ.0) THEN
       HF2NAME='HF2     '
      ELSE
       HF2NAME='HF2AB   '
      ENDIF
      CALL GFNAME(HF2NAME,FNAME,ILENGTH)
      OPEN(UNIT=20,FILE=FNAME(1:ILENGTH),FORM='UNFORMATTED',
     &     STATUS='OLD')
      REWIND(20)
      READ(20)JUNK
      READ(20)FOCKA
      IF(IUHF.EQ.1)READ(20)FOCKB
      CLOSE(UNIT=20,STATUS='KEEP')
      RETURN
      END 