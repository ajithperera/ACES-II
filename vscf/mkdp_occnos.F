      SUBROUTINE MKDP_OCCNOS(DOCC, NOCC, NSUM, NBF, IUHF, PRINT)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      
      DIMENSION NOCC(8, 2), DOCC(NBF*2), NSUM(8,2), NOC_ORB(8,2)
C
      CHARACTER*80 FNAME, Blank
      LOGICAL OCCNUMS_FILE, PRINT
C
#include "machsp.com"
#include "symm2.com"
C
      IUNIT = 5
      CALL DZERO(DOCC, NBF*2)
      CALL ICOPY(16, NOCC, 1, NSUM, 1)
#ifdef _DEBUG_LVLM
         Write(6,*)
         Write(6,"(a)") " The Alpha and Beta occupation vector"
         Write(6,*) (NOCC(i,1), i=1, NIRREP)
         Write(6,*)
         Write(6,*) (NOCC(i,2), i=1, NIRREP)
#endif
      DO ISPIN = 1, (IUHF+1)
         IROFF = (ISPIN - 1)*NBF
         DO IRREP = 1, NIRREP
            ISTART = IROFF 
            DO IBF = 1, NOCC(IRREP, ISPIN)

               ISTART = ISTART + 1
               DOCC(ISTART) = 1.0D0
C
            ENDDO
            IROFF = IROFF + NBFIRR(IRREP)
         ENDDO
      ENDDO

C
C If the occ. nos given in an external file, then override 
C what is read from the occupation number array.
C
      CALL GFNAME('OCCNUMS',FNAME,ILENGTH)
      INQUIRE(FILE=FNAME(1:7), EXIST=OCCNUMS_FILE)

      IF (OCCNUMS_FILE) THEN 
         CALL IZERO(NOCC, 16)
         IF (PRINT) THEN
          write(6,*)
          Write(6,"(a,a)") " The occupation numbers are read",
     &                     " from OCCNUMS file"
          write(6,*)
         ENDIF 
         OPEN(UNIT=IUNIT, FILE="OCCNUMS", FORM="FORMATTED")

         READ(IUNIT, "(80a)") Blank
         READ(IUNIT, "(80a)") Blank
C
         READ(IUNIT,10,END=19) (NOC_ORB(i, 1), i=1, NIRREP)
         READ(IUNIT,10,END=19) (NOC_ORB(i, 2), i=1, NIRREP)

#ifdef _DEBUG_LVLM
         Write(6,"(a)") "The Alpha and Beta no. of occupied orbitals"
         Write(6,*) (NOC_ORB(i,1), i=1, NIRREP)
         Write(6,*)
         Write(6,*) (NOC_ORB(i,2), i=1, NIRREP)
#endif
         READ(IUNIT, "(80a)") Blank
         READ(IUNIT, "(80a)") Blank

C These additional blank lines have a specific role in the fractional 
C occupation work (now abandoned). 07/2021,  Ajith Perera

         READ(IUNIT, "(80a)") Blank
         READ(IUNIT, "(80a)") Blank
C
         READ(IUNIT,10,END=19) (NSUM(i, 1), i=1, NIRREP)
         READ(IUNIT,10,END=19) (NSUM(i, 2), i=1, NIRREP)
C
#ifdef _DEBUG_LVLM
         Write(6,"(a)") "The Alpha and Beta NSUM vector"
         Write(6,*) (NSUM(i,1), i=1, NIRREP)
         Write(6,*)
         Write(6,*) (NSUM(i,2), i=1, NIRREP)
#endif
         READ(IUNIT, "(80a)") Blank

         DO ISPIN =1, (IUHF+1)
            IROFF = (ISPIN - 1)*NBF 
            DO IRREP = 1, NIRREP
               ISTART = IROFF
               DO IBF = 1, NOC_ORB(IRREP, ISPIN)
                  ISTART = ISTART + 1
CSSS                  READ(IUNIT, "(F4.2)",END=19) DOCC_NUM
                  READ(IUNIT, * ,END=19) DOCC_NUM
                  DOCC(ISTART) = DOCC_NUM 
                  IF (DOCC_NUM .NE. 0) NOCC(IRREP,ISPIN) = 
     &                                 NOCC(IRREP,ISPIN) + 1
               ENDDO
               DO IBF = 1, (NSUM(IRREP,ISPIN)-NOC_ORB(IRREP,ISPIN))
                  ISTART = ISTART + 1
                  DOCC(ISTART) = 0.0D0
               ENDDO
               IROFF = IROFF + NBFIRR(IRREP)
             ENDDO
         READ(IUNIT, "(80a)") Blank
         ENDDO 
      ENDIF
 10   FORMAT(16I5)
 19   CLOSE(IUNIT)
C
C Copy the occupation in the case of RHF
C
      IF (IUHF .EQ. 0) CALL DCOPY(NBF, DOCC, 1, DOCC(NBF+1), 1)
      IF (IUHF .EQ. 0) CALL ICOPY(8, NOCC(1,1), 1, NOCC(1,2),1)

      CALL PUTREC(20,"JOABRC","ORB_OCCA",NBF,DOCC)
      IF (IUHF .NE.0) CALL PUTREC(20,"JOABRC","ORB_OCCB",NBF,
     &                            DOCC(NBF+1))

      IF (PRINT) THEN
      Write(6,"(a)")" The Alpha and Beta double pre. occupation numbers"
      write(6,*)
      Write(6,"(6F10.5)")(DOCC(I), I=1, nbf)
      Write(6,*)
      Write(6,"(6F10.5)")(DOCC(nbf+I), I=1, nbf)
      write(6,*) 
      Write(6,"(2a)") "  The No. of Alpha and Beta occupied orbitals",
     &                " (can be greater than" 
      Write(6,"(2a)") "  no. electons if fractional occupations are",
     &                " used)"
      write(6,*)
      Write(6,"(5x,a,8i5)") "No. of alpha occupied orbitals:",
     &                      (NOCC(i,1), i=1, NIRREP)
      Write(6,"(5x,a,8i5)") "No. of beta occupied orbitasl :",
     &                      (NOCC(i,2), i=1, NIRREP)
      write(6,*)
      ENDIF 

      RETURN
      END
                
