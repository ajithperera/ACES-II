C
      SUBROUTINE FERMIADAPT (AOINT, SOINT, SINTLTR, SYMTRN, NATOM, NAO,
     &                       MAXCENT, IDEGN, ISTART, NSIZE, NTSIZE,
     &                       IRREPERT, NLTRN)
C
C This routine symmetry adapt the Fermi-contact integrals. 
C   
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      INTEGER DIRPRD
      DIMENSION AOINT(NTSIZE), SOINT(NTSIZE), SINTLTR(NLTRN),
     &          SYMTRN(6*MAXCENT, 6*MAXCENT), IRREPERT(6*MAXCENT),
     &          IOFFSET(8)
     
      LOGICAL JFC, JPSO, JSD, NUCLEI
C
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /MACHSP/ IINTLN, IFLTLN, IINTFP, IALONE, IBITWD
      COMMON /FILES/ LUOUT, MOINTS
      COMMON /NMR/JFC, JPSO, JSD, NUCLEI
      COMMON /PERT/NTPERT, NPERT(8), IPERT(8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),DIRPRD(8,8)
C
      SAVE ICOUNT, IOFFSET
C
      DATA ICOUNT /1/
      DATA IOFFSET /8*0/ 
C
      IONE   = 1
      IFIRST = IONE
      KORB   = IONE
      LORB   = IONE
      KKK    = IONE
      JJJ    = IONE
      LLL    = IONE
      IPRINT = IFLAGS(1)
C
      IEND = ISTART + IDEGN - 1
C
      CALL ZERO (SOINT, NTSIZE)
C      
      DO 10 INDEX = ISTART, IEND
C
         LIMIT = JJJ*NAO*NAO
C
         DO 20 IORB = IFIRST, LIMIT
            DO 30 JNDEX = ISTART, IEND
C
               KORB = KKK + (LLL - 1)*NAO*NAO
               SOINT(IORB) = SOINT(IORB) + AOINT(KORB)*
     &                       SYMTRN(INDEX, JNDEX)
               LLL = LLL + 1 
C
 30         CONTINUE
C
            KKK = KKK + 1 
            LLL = IONE
C
 20      CONTINUE
C
         IFIRST = LIMIT + 1
         JJJ = JJJ + 1
         KKK = IONE
C
 10   CONTINUE
C     
      IF (IPRINT .GE. 40) THEN
         ITRACK = IONE
         DO 40 II = 1, IDEGN
            CALL HEADER ('SYMMETRY AD. FERMI-CONTACT INTEGRALS', -1, 6)
            CALL TAB (LUOUT, SOINT(ITRACK), NAO, NAO, NAO, NAO)
            ITRACK = ITRACK + NAO*NAO
 40      CONTINUE
      ENDIF
C      
      INDEX  = IONE
C
      DO 50 II = 1, IDEGN
C
         CALL ZERO(SINTLTR, NLTRN)
         CALL SQUEZ2 (SOINT(INDEX), SINTLTR, NAO)
C
         IRREPX = IRREPERT(ICOUNT)
         INPERT = NPERT(IRREPX)
C
         IF (INPERT .NE. 0) THEN
            IOFFSET(IRREPX) = IOFFSET(IRREPX) + 1
         ENDIF
C
C Write the spin-adapted Fermi-contact integral to 'DERINT' file
C The list number is 398.
C
         CALL PUTLST(SINTLTR, IOFFSET(IRREPX), 1, 1, IRREPX, 398)
C
         INDEX = INDEX + NAO*NAO
         ICOUNT = ICOUNT + 1 

 50   CONTINUE
C
      RETURN
      END
