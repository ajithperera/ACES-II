      SUBROUTINE SYMMETW3(IRREP,LISTG,DISSYG,SCR,MAXSIZE)
C
C   THIS ROUTINE SYMMETRIZES A MATRIX A(N,N) :
C
C    B(I,J) = ( A(I,J) + A(J,I))
C
C   THIS VERSION HANDLES THE CASES WHERE THE MATRIX A DOES NOT FIT INTO
C   CORE
C
CEND
C 
C  CODED JG/SEPT/90
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER DISSYG,DISLEFT,DISREAD
      DIMENSION SCR(MAXSIZE)
C
C  CHECK IF A FITS INTO THE AVAILABLE CORE MEMORY
C
      IF(MAXSIZE.GE.DISSYG*DISSYG) THEN
C
C   DO IT IN CORE BY USING SYMMET2 
C
       CALL GETLST(SCR,1,DISSYG,1,IRREP,LISTG)
       CALL SYMMETW2(SCR,DISSYG) 
       CALL PUTLST(SCR,1,DISSYG,1,IRREP,LISTG)
C
      ELSE
C
C  WE HAVE TO DO IT OUT OF CORE
C
C  DETERMINE MAXIMUM NUMBER OF DISTRIBUTIONS WHICH CAN BE HELD IN
C  CORE
C
       MAXDIS=MAXSIZE/DISSYG
C
C  HOWEVER, ONE DISTRIBUTION IS REQUIRED AS BUFFER, SO SUBTRACT ONE
C
       MAXDIS=MAXDIS-1
       IF(MAXDIS.LE.0) STOP 'SYMMET3'
C
C  ALLOCATE MEMORY
C
       IBUF=1+MAXDIS*DISSYG
C
C  SET THE NUMBER OF DISTRIBUTIONS WHICH MUST BE TREATED
C
       DISREAD=DISSYG
C
C  SET THE NUMBER oF DISTRIBUTIONS LEFT ON DISK
C
       DISLEFT=DISSYG
C
C  SET OFFSET FOR LIST LISTG 
C 
       IOFFSET=1
C
10     CONTINUE
C
C  DETERMINE THE NUMBER OF DISTRIBUTIONS PROCESSED DURING THIS
C  PATH
C
        DISREAD=MIN(DISLEFT,MAXDIS)
        DISLEFT=DISLEFT-DISREAD
C 
C READ IN THE DISTRIBUTIONS WHICH ARE UPDATED
C
        CALL GETLST(SCR,IOFFSET,DISREAD,1,IRREP,LISTG)
C
C SET IFIRST (FIRST DISTRIBUTION IN CORE)
C
        IFIRST=IOFFSET
C
C UPDATE IOFFSET
C
        IOFFSET=IOFFSET+DISREAD
C
C SET ILAST (LAST DISTRIBUTION IN CORE)
C
        ILAST=IOFFSET-1
C
C  WE HAVE NOW THREE DIFFERENT CATEGORIES OF DISTRIBUTIONS
C
C   1)  THE DISTRIBUTIONS IN MEMORY
C   2)  THE DISTRIBUTIONS ON DISK WHICH ARE ALREADY SYMMETRIZED
C   3)  THE DISTRBIBUTIONS ON DISK WHICH ARE STILL NOT SYMMETRIZED
C
C WE HAVE ONLY TO LOOP OVER THE THIRS SET !
C
C   USE FIRST THE DISTRIBUTIONS IN CORE FOR SYMMETRIZATION
C
C   THIS IS DONE IN SYMMETW4
C
        CALL SYMMETW4(SCR(IFIRST),DISSYG,DISREAD)
C
C   NOW PROCESS THE THIRD SET OF DISTRIBUTIONS
C
C   LOOP OVER ALL REMAINING AND READ THEM
C
        DO 100 IDIS=IOFFSET,DISSYG
C
C  GET THE DISTRIBUTION IN THE BUFFER
C
         CALL GETLST(SCR(IBUF),IDIS,1,1,IRREP,LISTG)
C
C  NOW SYMMETRIZE THE PART WHICH CORRESPONDS TO THAT IN CORE
C
C  APPLY IVDEP DIRECTIVE IN ORDER TO VECTORIZE
C
*VOCL LOOP,NOVREC
CDIR$ IVDEP
         DO 150 IPOS=IFIRST,ILAST
          SCR(IBUF+IPOS-1)=(SCR(IBUF+IPOS-1)+ 
     &                      SCR(IDIS+(IPOS-IFIRST)*DISSYG))
150      CONTINUE
C
C  COPY THE SYMMETRIZED ELEMENTS TO THE UODATING AREA
C
C  APPLY IVDEP DIRECTIVE TO VECTORIZE
C
*VOCL LOOP,NOVREC
CDIR$ IVDEP
         DO 160 IPOS=IFIRST,ILAST
          SCR(IDIS+(IPOS-IFIRST)*DISSYG)=SCR(IBUF+IPOS-1)
160      CONTINUE
C
C  WRITE THE DISTRIBUTION IN THE BUFFER BACK TO DISK
C
         CALL PUTLST(SCR(IBUF),IDIS,1,1,IRREP,LISTG)
C
100     CONTINUE
C
C SAVE THE UPDATES DISTRIBUTIONS ON DISK
C
        CALL PUTLST(SCR,IFIRST,DISREAD,1,IRREP,LISTG)
C
C  IF NOT ALL DISTRIBUTIONS HAVE BEEN PROCESSED, GO BACK TO 10 
C
       IF(DISLEFT.NE.0) GO TO 10
C
C  ALL DONE, RETURN
C
      ENDIF
      RETURN
      END