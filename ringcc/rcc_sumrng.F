      SUBROUTINE RCC_SUMRNG(T2AIBJ,T2AJBI,ISCR,NSZ1,NSZ2,ISCRSZ)

      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION T2AIBJ(NSZ1),T2AJBI(NSZ2),ONE
      DIMENSION ISCR(ISCRSZ)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),D(18)
      DATA ONE /1.0/

      CALL GETALL(T2AJBI,NSZ1,1,42)
      CALL DCOPY(NSZ1,T2AJBI,1,T2AIBJ,1)

C Thiss multiplication by minus was in the standard CC code since
C we compute the negative of the P(ab)P(ij) T2(im,ae)W(mb,ej). 
C But this is not the case for drCCD or rCCD.
CSSS      CALL VMINUS(T2AIBJ,NSZ1)

      CALL PUTALL(T2AIBJ,NSZ1,1,42)
      RETURN
      END
