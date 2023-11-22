      SUBROUTINE READ_OCCNUM(DOCC,MAXCOR,ACT_MIN_A,ACT_MAX_A,ACT_MIN_B,
     +                       ACT_MAX_B,NBASIS,IUHF)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      
      DIMENSION DOCC(MAXCOR), NSUM(2), NOC_ORB(2)
      CHARACTER*80 FNAME, Blank
      INTEGER ACT_MIN_A,ACT_MIN_B,ACT_MAX_A,ACT_MAX_B
#include "machsp.com"
C
      Write(6,"(a,a)") " The occupation numbers are read",
     &                " from OCCNUMS file"
      Write(6,*)

      IUNIT = 5
      OPEN(UNIT=IUNIT,FILE="OCCNUMS",FORM="FORMATTED")

      READ(IUNIT, "(80a)") Blank
      READ(IUNIT, "(80a)") Blank
C
      READ(IUNIT,10,END=19) NOC_ORB(1)
      READ(IUNIT,10,END=19) NOC_ORB(2)

      READ(IUNIT, "(80a)") Blank

      READ(IUNIT,10,END=19) NSUM(1)
      READ(IUNIT,10,END=19) NSUM(2)
C
#ifdef _DEBUG_LVLM
         Write(6,"(a)") "The Alpha and Beta NSUM vector"
         Write(6,"(4(1x,I4))") NSUM(1), NSUM(2), NOC_ORB(1),
     +                         NOC_ORB(2)
#endif
      READ(IUNIT, "(80a)") Blank

      DO ISPIN =1, (IUHF+1)
         ISTART = 0
         DO IBF = 1, NOC_ORB(ISPIN)
            ISTART = ISTART + 1
            READ(IUNIT,*,END=19) DOCC_NUM
            DOCC(ISTART) = DOCC_NUM 
         ENDDO

         DO IBF = 1, (NSUM(ISPIN)-NOC_ORB(ISPIN))
            ISTART = ISTART + 1
            DOCC(ISTART) = 0.0D0
         ENDDO
         READ(IUNIT, "(80a)") Blank
         IF (ISPIN .EQ. 1) THEN
            CALL PUTREC(20,"JOABRC","ORB_OCCA",NBASIS,DOCC)
#ifdef _DEBUG_LVL0
      Write(6,"(a)")" The Alpha and Beta occupation numbers"
      write(6,*)
      Write(6,"(6F10.5)")(DOCC(I), I=1, NBASIS)
#endif
         ELSE 
            CALL PUTREC(20,"JOABRC","ORB_OCCB",NBASIS,DOCC)
#ifdef _DEBUG_LVL0
      write(6,*)
      Write(6,"(6F10.5)")(DOCC(I), I=1, NBASIS)
      write(6,*)
#endif
         ENDIF
      ENDDO 

C Read the active info blocks.

      READ(IUNIT,*,END=19) ACT_MIN_A
      READ(IUNIT,*,END=19) ACT_MAX_A
      READ(IUNIT, "(80a)") Blank
      READ(IUNIT,*,END=19) ACT_MIN_B
      READ(IUNIT,*,END=19) ACT_MAX_B

 10   FORMAT(16I5)
 19   CLOSE(IUNIT)
C
C Copy the occupation in the case of RHF
C
      IF (IUHF .EQ. 0) 
     &    CALL PUTREC(20,"JOABRC","ORB_OCCB",NBASIS,DOCC)

      RETURN
      END
                
