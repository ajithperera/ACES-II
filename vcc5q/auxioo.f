      SUBROUTINE AUXIOO
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DISTSZ
C      
      COMMON /AUXIO/ DISTSZ(8,100),NDISTS(8,100),INIWRD(8,100),LNPHYR,
     1               LUAUX
C
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
C
      CLOSE(UNIT=LUAUX,STATUS='DELETE')
C
      RETURN
      END