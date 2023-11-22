      SUBROUTINE GETLEN2(LEN, DISSYZ, NUM1, NUM2)
C
C GETS THE TOTAL LENGTH OF AN ARRAY WITH ROWS OF SIZE DISSYZ AND COLUMNS
C OF SIZE NUM1 * NUM2 FOR EACH IRREP.
C
      IMPLICIT INTEGER(A-Z)
      DIMENSION DISSYZ(8), NUM1(8), NUM2(8)
C
      COMMON/SYMINF/NSTART,NIRREP,IRREP0(255,2),DIRPRD(8,8)
C
      LEN = 0
      DO 10 IRREP = 1, NIRREP
         NUMDSS = 0
         DO 20 JIRREP = 1, NIRREP
            KIRREP = DIRPRD(IRREP,JIRREP)
            NUMDSS = NUMDSS + NUM1(JIRREP)*NUM2(KIRREP)
 20      CONTINUE
         LEN = LEN + DISSYZ(IRREP) * NUMDSS
 10   CONTINUE
C
      RETURN
      END