      

      SUBROUTINE SEEKLB(TARGET,IERR,IMODE)
C
C THIS ROUTINE FINDS A RECORD LABEL IN A VPROPS OUTPUT FILE.
C
      CHARACTER*8 DUM1,DUM2,DUM3,TITLE,TARGET
      IF(IMODE.EQ.0)THEN
       REWIND(30)
      ENDIF
      IERR=0
1     READ(30,END=300,ERR=301)DUM1,DUM2,DUM3,TITLE
      IF(TITLE.EQ.TARGET)THEN
       RETURN
      ELSE
       GOTO 1
      ENDIF
300   IERR=1
      RETURN
301   WRITE(6,2000)
2000  FORMAT(T3,'@SEEKLB-F, I/O error on file VPOUT.')
      CALL ERREX
      END
