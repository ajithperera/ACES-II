C
      SUBROUTINE DRVSPNDIP (SPNDINT, MAXCOR, MAXCENT, NATOM, NAO, NMAX,
     &                      NUCDEG, SYMTRN, IRREPERT, IUHF)
C
C Symmetry adapt spin dipole integrals
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      INTEGER DIRPRD
C
      DIMENSION SPNDINT(MAXCOR), SYMTRN(6*MAXCENT, 6*MAXCENT),
     &          NUCDEG(MAXCENT), IRREPERT(6*MAXCENT)
C
      CHARACTER*8 LABELS(6)
      CHARACTER*80 FNAME
      LOGICAL JFC, JPSO, JSD, YESNO, NUCLEI
C
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /MACHSP/ IINTLN, IFLTLN, IINTFP, IALONE, IBITWD
      COMMON /FILES/ LUOUT, MOINTS
      COMMON /NMR/JFC, JPSO, JSD, NUCLEI
      COMMON /PERT/NTPERT, NPERT(8), IPERT(8)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),DIRPRD(8,8)
C
      DATA LABELS/'  SDXX  ','  SDXY  ','  SDXZ  ','  SDYY  ',
     &            '  SDYZ  ','  SDZZ  '/
C
      IONE   = 1
      ISTART = 1
      ISIX   = 6 
      I000   = 1
      IRWND  = 0
      IZER   = 0
      IPRINT = IFLAGS(1)
      NSIZE  = NAO*NAO
      NLTRN  = NAO*(NAO + 1)/2
C        
c      CALL GFNAME('DERINT  ', FNAME, ILENGTH)
c      INQUIRE(FILE = FNAME(1:ILENGTH), EXIST=YESNO)
c      IF (YESNO) THEN
         IENTER = 0
         IOFF   = 0
c      ELSE
c         IENTER = 1
c         IOFF   = 0
c      ENDIF
C     
      DO 5 IRREP = 1, NIRREP
C     
C Create the MOIO pointers for output integral lists
C
         INPERT = NPERT(IRREP)
C
         IF (INPERT .NE. 0) THEN
            CALL UPDMOI(INPERT, NLTRN, IRREP, 398, IENTER, IOFF)
            IENTER = 0
         ENDIF
C     
 5    CONTINUE
C
      CALL ZERO(SPNDINT, MAXCOR)
C
      DO 10 IATOM = 1, NMAX
C
         IDEGN = NUCDEG(IATOM)
         NTSIZE =  ISIX*IDEGN*NSIZE
C
         I010 = I000 + ISIX*IDEGN*NSIZE
         I020 = I010 + ISIX*IDEGN*NSIZE
         I030 = I020 + NLTRN
         IF (I030 .GT. MAXCOR) CALL ERREX       
C  
         ISTEP = IZER
C
         DO 20 JATOM = 1, NUCDEG(IATOM)
            DO 30 IIII = 1, 6
C
               CALL SEEKLB (LABELS(IIII), IERR, IRWND, 30) 
               IF (IERR .NE. 0) CALL ERREX
C
               I000 =  I000 + ISTEP*NSIZE
C
               CALL LOADINT (SPNDINT(I000), NATOM, NSIZE, NAO, IUHF)
C     
               IF (IPRINT .GE. 40) THEN
                  CALL HEADER ('SD INTEGRALS', -1, 6)
                  CALL TAB (LUOUT, SPNDINT(I000), NAO, NAO, NAO, NAO)
               ENDIF
C     
               IRWND = IONE 
               ISTEP = IONE
C
 30         CONTINUE
 20      CONTINUE
C
         I000 = 1
C
         CALL SDADAPT (SPNDINT(I000), SPNDINT(I010), SPNDINT(I020), 
     &                 SYMTRN, NATOM, NAO, MAXCENT, IDEGN, ISTART,
     &                 NSIZE, NTSIZE, IRREPERT, NLTRN)
C
         ISTART = ISTART + ISIX*IDEGN
C     
 10   CONTINUE
C
      RETURN
      END
