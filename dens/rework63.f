      SUBROUTINE REWORK63(ICORE, MAXCOR, IUHF)
C
C The list 63 have been temporarily used to store 
C the intermediate in drvaovv and need to be restored
C to the original form.
C
      IMPLICIT INTEGER (A-Z)
C
      DIMENSION ICORE(MAXCOR)
C
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMINF/NSTART,NIRREP,IRREPY(255,2),DIRPRD(8,8)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
C
      DO ISPIN = 3-2*IUHF, 3
C
         LISTD  = 47 + ISPIN
         LISTDI = 63 + ISPIN
         ISIZE  = ISYMSZ(ISYTYP(1,LISTD),ISYTYP(2,LISTD))
C
         I000 = 1
         I010 = I000 + IINTFP*ISIZE
         CALL GETALL(ICORE(I000), ISIZE, 1,LISTD)
         CALL INVERS(ICORE(I000), ISIZE) 
         CALL PUTALL(ICORE(I000),ISIZE,1,LISTDI)
C
      ENDDO
C
      RETURN
      END 
