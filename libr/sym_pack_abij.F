      Subroutine SYM_PACK_ABIJ(W,Maxcor,Ispin,Irrepx,Type)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      DIMENSION W(Maxcor)
      Character*4 Type

#include "machsp.com"
#include "sympop.com"
#include "syminf.com"
#include "sym.com"
#include "info.com"
#include "parcc.com"

      LEN_C = IDSYMSZ(18+ISPIN,20+ISPIN)
C
C Form <AB|IJ> from <AI|BJ> or retrive <AB|IJ> from the disk.
C
#ifdef _DEBUG_LVL0
       Write(6,"(a,1x,I4,1x,a,4a)") "Read the lists",Listc, "with type",
     &                            Type
#endif 
      I000 = 1
      If (Type .EQ. "AIBJ") Then
 
        LEN_C = IDSYMSZ(9+ISPIN,9+ISPIN)

         NSCRSZ=NVRTO(ISPIN)*NVRTO(ISPIN)+NOCCO(ISPIN)*NOCCO(ISPIN)+
     &          NVRTO(ISPIN)*NOCCO(ISPIN)

         I000 = 1
         I010 = I000 + Len_c
         I020 = I010 + Len_c
         I030 = I020 + Nscrsz

         CALL SSTGEN(W(I000),W(I010),LEN_C,VRT(1,ISPIN),
     &               POP(1,ISPIN),VRT(1,ISPIN),POP(1,ISPIN),
     &               W(I020),1,"1324") 
         CALL DCOPY(Len_c,W(I010),1,W(I000),1)
      Endif 

#ifdef _DEBUG_LVL0
      call checksum("W<AB|IJ>",W(I000),LEN_C)
#endif 

      I010 = I000 + Len_c
      I020 = I010 + Len_C

C Take the transpose; <AB|IJ>-><IJ|AB>; WIN(LEN2,LEN1)->WOUT(LEN1,LEN2)
     
      INCRM  = 0
      DO IRREPR = 1, NIRREP
         IRREPL = DIRPRD(IRREPR,IRREPX)
         LEN1   = IRPDPD(IRREPR,20+ISPIN)
         LEN2   = IRPDPD(IRREPL,18+ISPIN)
         CALL TRANSP(W(I000+INCRM),W(I010+INCRM),LEN1,LEN2)
         INCRM = INCRM + LEN1 * LEN2
      ENDDO 

#ifdef _DEBUG_LVL0
      call checksum("W<IJ|AB>",W(I010),LEN_C)
#endif 

      I030 = I020 + Len_C

C Change <IJ|AB> to <I<J|AB>; W(NDSSIZC,NUMDISC) -> W(NDSSIZW,NUMDISC) 

      IOFF2 = I020 
      IOFF1 = I010

      DO IRREPR = 1, NIRREP
         IRREPL = DIRPRD(IRREPR,IRREPX)
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

      INCRM = 0
      DO IRREPR = 1, NIRREP
         IRREPL = DIRPRD(IRREPR,1)
         LEN1   = IRPDPD(IRREPR,2+ISPIN)
         LEN2   = IRPDPD(IRREPL,18+ISPIN)
         CALL TRANSP(W(I020+INCRM),W(I010+INCRM),LEN2,LEN1)
         INCRM = INCRM + LEN1 * LEN2
      ENDDO 

#ifdef _DEBUG_LVL0
      length = idsymsz(1,2+ISPIN,18+ISPIN)
      call checksum("W<AB|I<J>",W(I010),length)
#endif 

C Change <AB|I<J> to <A<B|I<J>; W(NDSSIZC,NUMDISC) -> W(NDSSIZW,NUMDISC) 

      IOFF1 = I010 
      IOFF2 = I020

      DO IRREPR = 1, NIRREP
         IRREPL = DIRPRD(IRREPR,IRREPX)

         NUMDISC=IRPDPD(IRREPR,2+ISPIN)
         NDSSIZC=IRPDPD(IRREPL,18+ISPIN)
         NDSSIZW=IRPDPD(IRREPL,ISPIN)

         CALL SQSYM(IRREPL,VRT(1,ISPIN),NDSSIZW,NDSSIZC,NUMDISC,
     &             W(IOFF2),W(IOFF1))
         IOFF1 = IOFF1 + NDSSIZC * NUMDISC
         IOFF2 = IOFF2 + NDSSIZW * NUMDISC 

      ENDDO 

#ifdef _DEBUG_LVL0
      length = idsymsz(1,ISPIN,2+ISPIN)
      write(6,*) length
      call checksum("W<A<B|I<J>",W(I020),length)
      CALL GETALL(W(I000), length, 1, 14) 
      Write(6,"(5(1x,F10.7))") (W(I000+i),i=1,length)
      e=0.0d0
      ioff = 0
      do irrep = 1, nirrep
         ndissizw = irpdpd(irrep,ispin)
         numdissc = irpdpd(irrep,2+ispin)
         e = e + ddot(ndissizw*numdissc,w(i000+ioff),1,w(i020+iff),1)
         ioff = Ioff + ndissizw * numdissc
      enddo 
      Write(6,"(a,1x,F10.8)") "The Pseudo energy = ", e 
#endif 
      CALL DCOPY(LEN_C,W(I010), 1, W(I000), 1)
C
      RETURN
      END 
