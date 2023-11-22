
C THIS FUNCTION RETURNS THE RELATIVE POSITION IN A CIRCULAR BUFFER
C  OF LENGTH IORDER FOR ABSOLUTE INDEX INDABS, GIVEN THAT 
C  RELATIVE INDEX = 1  OCCUPIES THE POSITION GIVEN BY ISTART.
C
C FOR EXAMPLE, THE POSITION OF ABSOLUTE INDEX = 3 IN THE CIRCULAR
C  BUFFER (7 8 9 1 2 3 4 5 6) IS 6.

      INTEGER FUNCTION ICRCLC(INDABS,ISTART,IORDER)
      ICRCLC=MOD(INDABS+ISTART-1,IORDER)
      IF (ICRCLC.EQ.0) ICRCLC=IORDER
      RETURN
      END

