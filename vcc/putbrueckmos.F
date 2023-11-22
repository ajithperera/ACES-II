      SUBROUTINE PUTBRUECKMOS(A,LDIM,IUHF,NIRREP,NBFIRR,NOCC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL YESNO
      DIMENSION A(LDIM * (IUHF+1)),NBFIRR(8)
      DIMENSION NSKIP(8,2),NTIMES(8),NLAST(8)
      DIMENSION NOCC(8,2)
      COMMON /FLAGS/ IFLAGS(100)
C
C     Subroutine for writing the natural orbitals to a formatted file. The
C     format is supposed to be such that these MOs can be retrieved in a
C     subsequent calculation by GETMOS. Note that
C     while GETMOS can read a symmetry block of alpha or beta orbitals,
C     PUTMOS writes all blocks of alpha and beta symmetry at once. As a
C     result, PUTMOS does not need the skipping parameters. The MOs, in
C     the matrix A, are organized in IUHF+1 sets of LDIM elements, all
C     alpha coefficients coming before all beta ones. The coefficients
C     are stored in consecutive symmetry blocks for a given spin.
C
C     Note that N should be equal to NBFIRR(IRREP)
C
C     The file contains the MOs in the SO basis, with all alpha
C     symmetry blocks appearing before all the beta ones. Each
C     symmetry block is written as NTIMES(IRREP) blocks of NCOLS columns
C     of length NBFIRR(IRREP) and one block containing NLAST(IRREP)
C     columns of length NBFIRR(I). NSKIP(IRREP,ISPIN) is how many records
C     to be skipped before reading a symmetry block.
C
      IF(IFLAGS(1).GE.10)THEN
       WRITE(6,1000)
 1000  FORMAT(' @PUTBRUECKMOS-I, Writing natural orbitals to NATORBS.')
      ENDIF
C
      NCOLS = 4
C
      DO   10 I=1,NIRREP
      NTIMES(I) = NBFIRR(I) / NCOLS
      NLAST(I)  = NBFIRR(I) - NTIMES(I) * NCOLS
      IF(IFLAGS(1).GE.10)THEN
       WRITE(6,1010) I,NTIMES(I),NLAST(I)
 1010  FORMAT(' @PUTMOS-I, Symmetry ',I3,' Full ',I3, ' Partial ',I3)
      ENDIF
   10 CONTINUE
C
C If BRUECKMOS already exists, delete it and open a new one.
C
      INQUIRE(FILE='BRUECKMOS',EXIST=YESNO)
      IF(YESNO)THEN
       OPEN(71,FILE='BRUECKMOS',STATUS='OLD',ACCESS='SEQUENTIAL',
     1      FORM='FORMATTED')
       CLOSE(UNIT=71,STATUS='DELETE')
       IF(IFLAGS(1).GE.10)THEN
        WRITE(6,1030)
       ENDIF
 1030  FORMAT("@-PUTBRUECKMOS BRUECKMOS already exists and will be ",
     &         "deleted.")
      ENDIF
      OPEN(71,FILE='BRUECKMOS',STATUS='NEW',ACCESS='SEQUENTIAL',
     1     FORM='FORMATTED')
C
      WRITE(71,"(8(1x,I5))") (NOCC(i,1), i=1, NIRREP)
      WRITE(71,"(8(1x,I5))") (NOCC(i,2), i=1, NIRREP)

      IOFF = 0
      DO  100 ISPIN=1,IUHF+1
C
      DO   90 IRREP=1,NIRREP
C
C     Write all the groups of NCOLS vectors.
C
      IF(NTIMES(IRREP).NE.0)THEN
      DO   60 ITIMES=1,NTIMES(IRREP)
      DO   50 I     =1,NBFIRR(IRREP)
      WRITE(71,1020) (A(IOFF + (ITIMES-1)*NCOLS*NBFIRR(IRREP) + 
     1                      (J-1)*NBFIRR(IRREP) + I),J=1,NCOLS)
 1020 FORMAT(4F20.10)
   50 CONTINUE
   60 CONTINUE
      ENDIF
C
C     Write the remainder.
C
      IF(NLAST(IRREP).NE.0)THEN
      DO   70 I     =1,NBFIRR(IRREP)
      WRITE(71,1020) (A(IOFF + NTIMES(IRREP)*NCOLS*NBFIRR(IRREP) + 
     1                      (J-1)*NBFIRR(IRREP) + I),J=1,NLAST(IRREP))
   70 CONTINUE
      ENDIF
C
      IOFF = IOFF + NBFIRR(IRREP) * NBFIRR(IRREP)
   90 CONTINUE
  100 CONTINUE
C
      CLOSE(71,STATUS='KEEP')
      RETURN
      END
