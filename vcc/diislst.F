      SUBROUTINE DIISLST(IRREPX,IUHF,SINGLE)
C
C THIS ROUTINE CREATES THE ERROR VECTOR LISTS FOR DIIS EXTRAPOLATION
C
CEND
      IMPLICIT INTEGER (A-Z)
      LOGICAL SINGLE
      COMMON/FLAGS/IFLAGS(100)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NF1(2),NF2(2)
cYAU - This is the same in DODIIS0 and should probably be an include file.
#include "bwcc.com"
      integer ld_diis
      parameter (ld_diis=25)
      double precision r(ld_diis,ld_diis,maxref)
      common/diisdat1/r,idimdiis(maxref)
      save  /diisdat1/
cjp
C
      MAXDIM=IFLAGS(12)
      IF(SINGLE)THEN
       idimdiis(iref)=IRPDPD(IRREPX,9)+IUHF*IRPDPD(IRREPX,10)
      ELSE
       idimdiis(iref)=0
      ENDIF
      DO 110 ISPIN=3,1,IUHF-2
       idimdiis(iref)=idimdiis(iref)+IDSYMSZ(IRREPX,
     &          ISYTYP(1,60+ISPIN),ISYTYP(2,60+ISPIN))
110   CONTINUE
      if (maxdim.gt.ld_diis) then
         print *, '@DIISLST: Assertion failed.'
         print *, '          Requested ',maxdim,' iterates but only ',
     &            ld_diis,' are allowed.'
         call errex
      end if
      call zero(r,ld_diis*ld_diis*maxref)
cjp
      CALL UPDMOI(MAXDIM,idimdiis(iref),1,70,0,0)
      CALL UPDMOI(MAXDIM,idimdiis(iref),2,70,0,0)
cjp
      RETURN
      END
