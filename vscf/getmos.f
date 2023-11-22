      SUBROUTINE GETMOS(A,N,ISPIN,IRREP,NIRREP,NBFIRR,LUGSS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(N,N),NBFIRR(8)
      DIMENSION NSKIP(8,2),NTIMES(8),NLAST(8)
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /POPUL/  NOCC(8,2)
C
C     Note that N should be equal to NBFIRR(IRREP)
C
C     Subroutine for reading a symmetry block of alpha or beta MOs from
C     a formatted file.
C     The file contains the MOs in the SO basis, with all alpha
C     symmetry blocks appearing before all the beta ones. Each
C     symmetry block is written as NTIMES(IRREP) blocks of NCOLS columns
C     of length NBFIRR(IRREP) and one block containing NLAST(IRREP)
C     columns of length NBFIRR(I). NSKIP(IRREP,ISPIN) is how many records
C     to be skipped before reading a symmetry block.
C
      IF(IFLAGS(1).GE.10)THEN
       WRITE(6,1000)
 1000  FORMAT(' @GETMOS-I, Getting the MOs ! ')
      ENDIF
C
      NCOLS = 4
C
      DO   10 I=1,NIRREP
      NTIMES(I) = NBFIRR(I) / NCOLS
      NLAST(I)  = NBFIRR(I) - NTIMES(I) * NCOLS
      IF(IFLAGS(1).GE.10)THEN
       WRITE(6,*) NTIMES(I),NLAST(I)
      ENDIF
   10 CONTINUE
C
      DO   20 I=1,NIRREP
      IF(I.EQ.1)THEN
      NSKIP(I,1) = 0
      ELSE
      NSKIP(I,1) = NSKIP(I-1,1) + 
     1             (NTIMES(I-1) + MIN(NLAST(I-1),1)) * NBFIRR(I-1)
      ENDIF
      IF(IFLAGS(1).GE.10)THEN
       WRITE(6,*) NSKIP(I,1)
      ENDIF
   20 CONTINUE
C
      DO   30 I=1,NIRREP
      IF(I.EQ.1)THEN
      NSKIP(I,2) = NSKIP(NIRREP,1) + 
     1             (NTIMES(NIRREP) + MIN(NLAST(NIRREP),1)) * 
     1             NBFIRR(NIRREP)
      ELSE
      NSKIP(I,2) = NSKIP(I-1,2) + 
     1             (NTIMES(I-1) + MIN(NLAST(I-1),1)) * NBFIRR(I-1)
      ENDIF
      IF(IFLAGS(1).GE.10)THEN
       WRITE(6,*) NSKIP(I,2)
      ENDIF
   30 CONTINUE
C
C Read the occupation vector
C  
      IF (IRREP .EQ. 1) THEN
         READ(71,"(8(1x,I5))") (NOCC(i,1), i=1, NIRREP)
         READ(71,"(8(1x,I5))") (NOCC(i,2), i=1, NIRREP)
      ENDIF 

      IF (IFLAGS(1).GE.10) THEN
         Write(6,"(8(1x,I5))") (NOCC(i,1), i=1, NIRREP)
         Write(6,"(8(1x,I5))") (NOCC(i,1), i=1, NIRREP)
      ENDIF 
C
C     Skip records if necessary.
C
      IF(NSKIP(IRREP,ISPIN).NE.0)THEN
      DO   40 ISKIP=1,NSKIP(IRREP,ISPIN)
c     READ(71,*) CRAP
      READ(LUGSS,*) CRAP
   40 CONTINUE
      ENDIF
C
C This is to skip the two new lines that contains the occupation numbers.
C
      IF (IRREP .NE. 1) THEN
        READ(LUGSS,*) CRAP
        READ(LUGSS,*) CRAP
      ENDIF 
C
C     Read all the groups of NCOLS vectors.
C
      IF(NTIMES(IRREP).NE.0)THEN
      DO   60 ITIMES=1,NTIMES(IRREP)
      DO   50 I     =1,N
c     READ(71,*) (A(I,(ITIMES-1)*NCOLS + J),J=1,NCOLS)
      READ(LUGSS,*) (A(I,(ITIMES-1)*NCOLS + J),J=1,NCOLS)
   50 CONTINUE
   60 CONTINUE
      ENDIF
C
C     Read the remainder.
C
      IF(NLAST(IRREP).NE.0)THEN
      DO   70 I     =1,N
c     READ(71,*) (A(I,NTIMES(IRREP)*NCOLS + J),J=1,NLAST(IRREP))
      READ(LUGSS,*) (A(I,NTIMES(IRREP)*NCOLS + J),J=1,NLAST(IRREP))
   70 CONTINUE
      ENDIF
C
      RETURN
      END
