      SUBROUTINE DENOWT(TBUF,E,I,J,NSIZE,NDSSIZ,ISPIN)
C
C THIS ROUTINE DENOMINATOR WEIGHTS A DISTRIBUTION OF T AMPLITUDES.
C
cjp this applies only to rle method?!
cjp but in BW we must take care of denominator shift and of internal amplitudes
cjp and of orbital index mapping
cjp this routine was not checked, since we did not use rle
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LINCC,CICALC
      INTEGER A,B,atrue,btrue
      DIMENSION TBUF(NDSSIZ),E(NSIZE)
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      COMMON /LINEAR/ LINCC,CICALC
      COMMON /CORENG/ ELAST
cjp
#include "bwcc.com"
cjp
      NORB=NOCCO(1)+NVRTO(1)
      FACT=0.0
      IF(CICALC)FACT=ELAST
cjp   
      if(isbwcc) fact=ecorrbw
      IF(ISPIN.LE.2)THEN
       IX=0
       DO 10 A=2,NVRTO(ISPIN)
        DO 20 B=1,A-1
         IX=IX+1
cjp
         if(isbwcc) then
                itrue=invhnum(i,iref,ispin)
                jtrue=invhnum(j,iref,ispin)
                atrue=invpnum(a-nocco(1),iref,ispin)
                btrue=invpnum(b-nocco(2),iref,ispin)
                if(isactive(itrue,ispin).and.isactive(jtrue,ispin)
     &     .and.isactive(atrue,ispin).and.isactive(btrue,ispin)) then
                denom=denomblow
                else
                denom=e(atrue)+e(btrue)-e(itrue)-e(jtrue)
                endif
         else
         DENOM=E(NOCCO(ISPIN)+A)+E(NOCCO(ISPIN)+B)-E(I)-E(J)
         endif
         TBUF(IX)=TBUF(IX)*(DENOM+FACT)
20      CONTINUE
10     CONTINUE
      ELSEIF(ISPIN.EQ.3)THEN
       IX=0
       DO 30 B=1,NVRTO(2)
        DO 40 A=1,NVRTO(1)
         IX=IX+1
cjp
cjp: b,j are beta
         if(isbwcc) then
                itrue=invhnum(i,iref,1)
                jtrue=invhnum(j,iref,2)
                atrue=invpnum(a-nocco(1),iref,1)
                btrue=invpnum(b-nocco(2),iref,2)
        if(isactive(itrue,1).and.isactive(jtrue,2).and.
     &           isactive(atrue,1).and.isactive(btrue,2)) then
                denom=denomblow
                else
        denom=e(atrue)+e(btrue+norb)-e(itrue)-e(jtrue+norb)
        endif
         else
         DENOM=E(NOCCO(1)+A)+E(NOCCO(2)+B+NORB)
     &         -E(I)-E(J+NORB)
         endif
         TBUF(IX)=TBUF(IX)*(DENOM+FACT)
40      CONTINUE
30     CONTINUE
      ENDIF
      RETURN
      END
