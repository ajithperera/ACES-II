      Subroutine Make_coulomb_abij(W,Maxcor,ListC,ISPIN,Type)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      DIMENSION W(Maxcor)
      CHARACTER*4 TYPE

#include "machsp.com"
#include "sympop.com"
#include "syminf.com"
#include "sym.com"
#include "info.com"
#include "parcc.com"

      LEN_C = IDSYMSZ(1,ISYTYP(1,LISTC),ISYTYP(2,LISTC))
C
C Form <AB|IJ> from <AI|BJ> or retrive <AB|IJ> from the disk.
C
      If (Type .EQ. "AIBJ") Then

         NSCRSZ=NVRTO(ISPIN)*NVRTO(ISPIN)+NOCCO(ISPIN)*NOCCO(ISPIN)+
     &          NVRTO(ISPIN)*NOCCO(ISPIN)

         I000 = 1
         I010 = I000 + Len_c
         I020 = I010 + Len_c
         I030 = I020 + Nscrsz

         CALL GETALL(W(I000), LEN_C, 1, LISTC) 

         CALL SSTGEN(W(I000),W(I010),LEN_C,VRT(1,ISPIN),
     &               POP(1,ISPIN),VRT(1,ISPIN),POP(1,ISPIN),
     &               W(I020),1,"1324") 
         CALL DCOPY(Len_c, W(I010),W(I000))
      Else

         I000 = 1
         I010 = I000 + Len_c
         CALL GETALL(W(I000), LEN_C, 1, LISTC) 

      Endif 

#ifdef _DEBUG_LVL0
      call checksum("W<AB|IJ>",W(I000),LEN_C)
#endif 

      I020 = I010 + Len_C

C Take the transpose; <AB|IJ>-><IJ|AB>; WIN(LEN1,LEN2)->WOUT(LEN2,LEN1)

      DO IRREPR = 1, NIRREP
         IRREPL = DIRPRD(IRREPR,1)
         LEN1   = IRPDPD(IRREPR,20+ISPIN)
         LEN2   = IRPDPD(IRREPL,18+ISPIN)
         CALL TRANSP(W(I000+Incrm1),W(I010+Incrm),LEN2,LEN1)
         INCRM = IINCRM + LEN1 * LEN2
      ENDDO 

#ifdef _DEBUG_LVL0
      call checksum("W<IJ|AB>",W(I020),LEN_C)
#endif 

      I030 = I020 + Len_C

C Change <IJ|AB> to <I<J|AB>; W(NDSSIZ,NUMDISC) -> W(NDSSIZW,NUMDISC) 

      IOFF2 = I020 
      IOFF1 = I010

      DO IRREPR = 1, NIRREP
         IRREPL = DIRPRD(IRREPR,1)
         NUMDISC=IRPDPD(IRREPR,18+ISPIN)

         NDSSIZC=IRPDPD(IRREPL,20+ISPIN)
         NDSSIZW=IRPDPD(IRREPL,2+ISPIN)

         CALL SQSYM(IRREPL,POP(1,ISPIN),NDSSIZW,NDSSIZC,NUMDISC,
     &              W(IOFF2),W(IOFF1))

         IOFF1 = IOFF1 + NDSSIZC * NUMDISC
         IOFF2 = IOFF2 + NDSSIZW * NUMDISC 
      ENDDO 

#ifdef _DEBUG_LVL0
      length = idsymsz(1,2+ISPIN,18+ISPIN)
      call checksum("W<I<J|AB>",W(I020),length)
#endif 

C Take the transpose; <I<J|AB> -> <AB|I<J>; WIN(LEN1,LEN2)->WOUT(LEN2,LEN1)

      DO IRREPR = 1, NIRREP
         IRREPL = DIRPRD(IRREPR,1)
         LEN1   = IRPDPD(IRREPR,12+ISPIN)
         LEN2   = IRPDPD(IRREPL,18+ISPIN)
         CALL TRANSP(W(I020+Incrm1),W(I010+Incrm),LEN2,LEN1)
         INCRM = IINCRM + LEN1 * LEN2
      ENDDO 

#ifdef _DEBUG_LVL0
      call checksum("W<AB|I<J>",W(I010),length)
#endif 

C Change <IJ|AB> to <I<J|AB>; W(NDSSIZ,NUMDISC) -> W(NDSSIZW,NUMDISC) 

      IOFF1 = I010 
      IOFF2 = I020

      DO IRREPR = 1, NIRREP
         IRREPL = DIRPRD(IRREPR,1)

         NUMDISC=IRPDPD(IRREPR,2+ISPIN)
         NDSSIZC=IRPDPD(IRREPL,18+ISPIN)
         NDSSIZW=IRPDPD(IRREPL,0+ISPIN)

         CALL SQSYM(IRREPL,VRT(1,ISPIN),NDSSIZW,NDSSIZC,NUMDISC,
     &             W(IOFF2),W(IOFF1))

         IOFF1 = IOFF1 + NDSSIZC * NUMDISC
         IOFF2 = IOFF2 + NDSSIZW * NUMDISC 

      ENDDO 

#ifdef _DEBUG_LVL0
      length = idsymsz(1,ISPIN,2+ISPIN)
      call checksum("W<A<B|I<J>",W(I020),length)
#endif 

      RETURN
      END
