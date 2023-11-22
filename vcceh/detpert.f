C
      SUBROUTINE DETPERT (MAXCENT, NUCIND, NSYMOP, IPTCNT, CSTRA)
C
C This subroutine assign a perturbation to each  basis function      
C     
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C    
      DIMENSION IPTCNT(6*MAXCENT, 0:7), CSTRA(6*MAXCENT, 6*MAXCENT)
C
      LOGICAL JFC, JPSO, JSD, NUCLEI
C
      COMMON /NMR/JFC, JPSO, JSD, NUCLEI
      COMMON /MACHSP/ IINTLN, IFLTLN, IINTFP, IALONE, IBITWD
      COMMON /FILES/ LUOUT, MOINTS
      COMMON /FLAGS/IFLAGS(100)
      COMMON /PERT/NTPERT, NPERT(8), IPERT(8)
C
      IBTAND(I, J) = AND(I, J)
      IBTOR(I, J)  = OR(I, J)
      IBTXOR(I, J) = XOR(I, J)
C
      MAXLOP = 2**NSYMOP - 1 
      NIRREP = MAXLOP + 1 
      NTPERT = 0
      CALL IZERO(NPERT, 8)
      CALL IZERO(IPERT, 8)
C
      IF (JFC) THEN
         DO 10 IRREP = 1, MAXLOP + 1
            DO 20 ICENT = 1, NUCIND
               IF (IPTCNT(ICENT, IRREP - 1 ) .NE. 0) THEN
                  NPERT(IRREP) = NPERT(IRREP) + 1
               ENDIF
 20         CONTINUE
            NTPERT = NTPERT + NPERT(IRREP)
 10      CONTINUE
         IPERT(1) = 0
         DO 30 IRREP = 1, NIRREP - 1
             IPERT(IRREP + 1) = IPERT(IRREP) + NPERT(IRREP)
 30      CONTINUE
C
      ELSE IF (JSD) THEN
         DO 40 IRREP = 1, MAXLOP + 1
            DO 50 ICENT = 1, 6*NUCIND
               IF (IPTCNT(ICENT, IRREP - 1) .NE. 0) THEN 
                  NPERT(IRREP) = NPERT(IRREP) + 1
               ENDIF
 50         CONTINUE
            NTPERT = NTPERT + NPERT(IRREP)
 40      CONTINUE
C
         IPERT(1)=0
         DO 60 IRREP = 1,NIRREP - 1 
            IPERT(IRREP + 1) = IPERT(IRREP) + NPERT(IRREP)
 60      CONTINUE   
C
      ELSE IF (JPSO) THEN            
         DO 70 IRREP = 1, MAXLOP + 1
            DO 80 ICENT = 1, 3*NUCIND
               IF (IPTCNT(ICENT, IRREP - 1) .NE. 0) THEN 
                  NPERT(IRREP) = NPERT(IRREP) + 1
               ENDIF
 80         CONTINUE
            NTPERT = NTPERT + NPERT(IRREP)
 70      CONTINUE
C
         IPERT(1)=0
         DO 90 IRREP = 1, NIRREP - 1 
            IPERT(IRREP + 1) = IPERT(IRREP) + NPERT(IRREP)
 90      CONTINUE   
C
      ELSE
         CALL ERREX 
         RETURN
      ENDIF
C
      RETURN
      END
