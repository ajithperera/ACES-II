         SUBROUTINE  NLO__OPEN_FILES
     +
     +                    ( NAME1,NAME2,
     +                      FILEST,
     +
     +                              UNITID,
     +                              UNITDB )
     +
C------------------------------------------------------------------------
C  OPERATION   : NLO__OPEN_FILES
C  MODULE      : Natural Localized Orbitals
C  MODULE-ID   : NLO
C  DESCRIPTION : This operation opens the two sequential files:
C
C                         NAME1 : printout file
C                         NAME2 : diagnostic/debug file
C
C                and returns their unit numbers.
C
C                Processing is aborted, if one of the following cases
C                apply:
C
C                     1) Files are already open
C                     2) Files should have status old, but are
C                        not existing
C                     3) Files should have status new, but are
C                        existing
C                     4) Files should be scratch files, but a
C                        name was given for the file (not strikt
C                        Fortran Standard!)
C
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C
C             ...include files and declare variables.
C
C
         CHARACTER*(*)      FILEST
         CHARACTER*(*)      NAME1,NAME2

         LOGICAL            EXISTS
         LOGICAL            OPSTAT

         INTEGER            UNITID,UNITDB
C
C
C------------------------------------------------------------------------
C
C
C             ...inquire the open and existence status of the
C                printout file.
C
C
         INQUIRE  ( FILE   = NAME1,
     +              OPENED = OPSTAT,
     +              EXIST  = EXISTS )
C
C
C             ...printout file already open.
C
C
         IF  ( OPSTAT )  THEN
               WRITE (1,*) ' The following file is already open: '
               WRITE (1,*) NAME1
               WRITE (*,*) ' The following file is already open: '
               WRITE (*,*) NAME1
               STOP
         END IF
C
C
C             ...check printout file existence if status = old or new.
C
C
         IF  ( .NOT.EXISTS  .AND.  FILEST.EQ.'OLD' )  THEN
               WRITE (1,*) ' The following file does not exist: '
               WRITE (1,*) NAME1
               WRITE (*,*) ' The following file does not exist: '
               WRITE (*,*) NAME1
               STOP
         END IF

         IF  ( EXISTS  .AND.  FILEST.EQ.'NEW' )  THEN
               WRITE (1,*) ' The following file already exists: '
               WRITE (1,*) NAME1
               WRITE (*,*) ' The following file already exists: '
               WRITE (*,*) NAME1
               STOP
         END IF
C
C
C             ...inquire the open and existence status of the
C                diagnostic/debug file.
C
C
         INQUIRE  ( FILE   = NAME2,
     +              OPENED = OPSTAT,
     +              EXIST  = EXISTS )
C
C
C             ...diagnostic/debug file already open.
C
C
         IF  ( OPSTAT )  THEN
               WRITE (1,*) ' The following file is already open: '
               WRITE (1,*) NAME2
               WRITE (*,*) ' The following file is already open: '
               WRITE (*,*) NAME2
               STOP
         END IF
C
C
C             ...check diagnostic/debug file existence if status = old
C                or new.
C
C
         IF  ( .NOT.EXISTS  .AND.  FILEST.EQ.'OLD' )  THEN
               WRITE (1,*) ' The following file does not exist: '
               WRITE (1,*) NAME2
               WRITE (*,*) ' The following file does not exist: '
               WRITE (*,*) NAME2
               STOP
         END IF

         IF  ( EXISTS  .AND.  FILEST.EQ.'NEW' )  THEN
               WRITE (1,*) ' The following file already exists: '
               WRITE (1,*) NAME2
               WRITE (*,*) ' The following file already exists: '
               WRITE (*,*) NAME2
               STOP
         END IF
C
C
C             ...open the two files.
C
C
         UNITID = 2
         UNITDB = 1

         OPEN  ( UNIT   = UNITID,
     +           ACCESS = 'SEQUENTIAL',
     +           FILE   = NAME1,
     +           FORM   = 'FORMATTED',
     +           STATUS = FILEST )

         OPEN  ( UNIT   = UNITDB,
     +           ACCESS = 'SEQUENTIAL',
     +           FILE   = NAME2,
     +           FORM   = 'FORMATTED',
     +           STATUS = FILEST )
C
C
C             ...ready!
C
C
         RETURN
         END
