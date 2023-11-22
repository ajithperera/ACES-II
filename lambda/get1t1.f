      SUBROUTINE GET1T1(ICORE,MAXCOR,inext, ISPIN, IUHF,IOFFT1, IERR)
C
C  This routine loads the requested t1 vector at the bottom of icore,
C and returns a pointer array giving the offsets where each
c  irrep of the given spin case begins.  For RHF, calling for ISPIN=2
C  will return the same 
C  IRREP OF THE SPIN CASE BEGINS. This works for RHF or UHF cases
C
C  PARAMS:   USE                                                        CHANGED
C  -------   ---------------------------------------------------        -------
C    ICORE - THE CORE VECTOR (T1 RETURNED AT BOTTOM)                        YES
C   MAXCOR - THE TOTAL CORE SIZE  AVAILABLE for T1                           NO
C   INEXT  - The next position available in icore AFTER the t1 is loaded    YES
C    ISPIN - The spin case to use                                            NO
C     IUHF - THE UHF/RHF FLAG                                                NO
C   IOFFT1 - A ONE DIMENSIONAL ARRAY GIVING THE ADDRESS OF                  YES
C             THE BEGINNING OF EACH IRREP IN THE T1 VECTOR.
C             FOR EXAMPLE, IOFFT1(3) GIVES THE ADDRESS OF
C             THE FIRST ELEMENT OF THE THIRD IRREP FOR the SPIN
C             CASE 
C
CEND
      IMPLICIT INTEGER (A-Z)
      DIMENSION ICORE(MAXCOR),IOFFT1(8)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NF1AA,NF2AA,NF1BB,NF2BB
      MXCOR=MAXCOR
C
C COMPUTE OFFSETS FOR BEGINNING OF T1AA AND T1BB (1000 LOOP) AND OFFSETS FOR
C  BEGINNING OF IRREPS (2000 LOOP).
C
      IERR = 0
      TLIST=90
      TLIST2=ISPIN
      IF(IUHF.EQ.0)TLIST2=1
      IS = ISPIN
      IF(IUHF.EQ.0)THEN
         TLIST=90
         IS = 1
      ENDIF
      T1SIZ=NT(IS)
C
      IF ( T1SIZ .GT. MXCOR/IINTFP ) THEN
         IERR = 1
         RETURN
      ENDIF
      inext = 1 + t1siz*iintfp
C
      IOFF = 1
      DO 2000 IRREP=1,NIRREP
         IOFFT1(IRREP)=IOFF
         IOFF = IOFF + VRT(IRREP,IS)*POP(IRREP,IS)*IINTFP
 2000 CONTINUE
      CALL GETLST(ICORE,1,1,1,TLIST2,TLIST)
      RETURN
      END 
