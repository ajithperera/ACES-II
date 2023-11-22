      SUBROUTINE PUTADD(Z,DISSYZ,NUMSYZ,SCR,MAXSIZE,LISTZ)
C
C  THIS ROUTINE ADDS THE ARRAY Z TO AN ARRAY Z' WHICH RESIDES
C  ON DISK. IT IS ASSUMED THAT ONLY Z FITS INTO CORE AND THAT
C  ONLY PARTS OF Z' FIT INTO THE SCR ARRAY WITh SIZE MAXSIZE
C
C   ARGUMENTS :
C
C      Z ....... ARRAY Z WHICH SHOULD BE ADDED To LISTZ
C      DISSYZ .. DISTRIBUTION SIZE OF Z
C      NUMSYZ .. NUMBER OF DISTRIBUTION OF Z
C      SCR ..... SCRATCH ARRAY
C      MAXSIZE . SIZE OF THE SCRATCH ARRAY
C      LISTZ ... LIST NUMBER
C
C   NOTE SINCE, THIS IS A VERY FANCY ROUTINE PAY SOME ATTENTION
C   TO THE ALLOCATION. SINCE THE SCR ARRAY IS INCREASED IN EVERY PATH
C   BY THE AMOUNT OF MEMORY FREED BY ADDING UP PARTS OF Z. SCR MUST HAVE
C   A LOWER ADDRESS THAN Z
C   I.E. :
C             ISCR=1
C             IZ=1+MAXSIZE
C
C   NOTE THAT THIS WILL DECREASE THE NUMBER OF CALLS TO GETLST AND PUTLST !
C
C   FANCY STUFF !!!
C 
CEND
C
C CODED SEPTEMBER/90  JG
C   
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER  DISSYZ,DISLEFT,DISREAD,DISMAX
      DIMENSION Z(DISSYZ,NUMSYZ),SCR(MAXSIZE)
C
      DATA ONE /1.0D0/
C
C  DETERMINE MAXIMUM NUMBER OF DISTRIBUTIONS WHICH CAN BE KEPT IN SCR
C
      DISMAX=MAXSIZE/DISSYZ
C
C  SET NUMBER OF DISTRIBUTIONS TO READ
C
      DISLEFT=NUMSYZ
C
C  SET OFFSET FOR Z
C
      IOFFSET=1
C
10    CONTINUE 
C
C  DETERMINE NUMBER OF DISTRIBUTIONS READ DURING THIS PASS
C
       DISREAD=MIN(DISLEFT,DISMAX)
       DISLEFT=DISLEFT-DISREAD
C
C  GET THESE DISTRIBUTIONS FROM DISK
C
       CALL GETLST(SCR,IOFFSET,DISREAD,2,LISTZ)
C
C  ADD THE CORRESPONDING CONTRIBUTION TO SCR
C
       CALL SAXPY(DISREAD*DISSYZ,ONE,Z(DISSYZ,IOFFSET),1,SCR,1)
C
C  PUT THESE DISTRIBUTIONS BACK TO DISK
C
       CALL PUTLST(SCR,IOFFSET,DISREAD,2,LISTZ)
C
C UPDATE OFFSET 
C
       IOFFSET=IOFFSET+DISREAD
C
C INCREASE DISMAX
C
       DISMAX=DISMAX+DISREAD
C
C  IF ANY OF THE DISTRIBUTIONS HAVEN'T PROCESSED GO BACK TO 10
C
      IF(DISLEFT.NE.0) GO TO 10
C
C  ALL DONE
C
      RETURN
      END