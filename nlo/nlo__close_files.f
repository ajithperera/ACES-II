         SUBROUTINE  NLO__CLOSE_FILES
     +
     +                    ( UNITID,UNITDB )
     +
C------------------------------------------------------------------------
C  OPERATION   : NLO__OPEN_FILES
C  MODULE      : Natural Localized Orbitals
C  MODULE-ID   : NLO
C  DESCRIPTION : This operation closes the two sequential files:
C
C                         UNITID : printout file
C                         UNITDB : diagnostic/debug file
C
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C
C             ...include files and declare variables.
C
C
         LOGICAL    OPSTAT

         INTEGER    UNITID,UNITDB
C
C
C------------------------------------------------------------------------
C
C
C             ...inquire the open status of the printout file and
C                close it if open.
C
C
         INQUIRE  ( UNIT   = UNITID ,
     +              OPENED = OPSTAT )

         IF  ( OPSTAT )  THEN
               CLOSE ( UNIT = UNITID )
         ENDIF
C
C
C             ...inquire the open status of the diagnostic/debug file
C                and close it if open.
C
C
         INQUIRE  ( UNIT   = UNITDB ,
     +              OPENED = OPSTAT )

         IF  ( OPSTAT )  THEN
               CLOSE ( UNIT = UNITDB )
         ENDIF
C
C
C             ...ready!
C
C
         RETURN
         END
