      Subroutine Scrnc_make_lists(Iuhf,Irrepx,List_aibj,List_xibj,
     +                            List_xxbj,List_xxxx,List_ijxx, 
     +                            List_abxx,List_aixx,List_iaxx)

      Implicit Double Precision (A-H,O-Z)

      COMMON/AOSYM/IAOPOP(8),IOFFAO(8),ioffv(8,2),ioffo(8,2),
     &             IRPDPDAO(8),IRPDPDAOMO_OCCBK(8,2),
     &             IRPDPDAOMO_VRTBK(8,2),IRPDPDAOMO_OCCKB(8,2),
     &             IRPDPDAOMO_VRTKB(8,2),
     &             IRPDPDAOS(8),
     &             ISTART(8,8),ISTARTMO(8,3)

#include "machsp.com"
#include "syminf.com"
#include "sympop.com"
#include "flags.com"
#include "parcc.com"
#include "info.com"


C Actuall lists is List_xxx+Ispin

      List_aibj = 4

      List_iaxx = 203
      List_aixx = 205
      List_abxx = 207
      List_ijxx = 209

      List_xibj = 214
      List_xxbj = 216

      List_xxxx = 220 
   
      Imod = 0
      Do Ispin = 1, Iuhf+1
         Call inipck4(Irrepx,IRPDPDAOMO_OCCKB(1,Ispin),
     +                IRPDPD(1,8+Ispin),List_xibj+ispin,Imod,0,0)

         Call inipck4(Irrepx,IRPDPDAO(1),IRPDPD(1,8+Ispin),
     +                List_xxbj+Ispin,Imod,0,0)


         Call inipck4(Irrepx,IRPDPD(1,8+Ispin),IRPDPDAO(1),
     +                List_aixx+Ispin,Imod,0,0)

         Call inipck4(Irrepx,IRPDPD(1,8+Ispin),IRPDPDAO(1),
     +                 List_iaxx+Ispin,Imod,0,1)
         Call inipck4(Irrepx,IRPDPD(1,18+Ispin),IRPDPDAO(1),
     +                List_abxx+Ispin,Imod,0,1)
         Call inipck4(Irrepx,IRPDPD(1,20+Ispin),IRPDPDAO(1),
     +                List_ijxx+Ispin,Imod,0,1)

      Enddo 

      Call inipck4(Irrepx,IRPDPDAO(1),IRPDPDAO(1),List_xxxx,
     +             Imod,0,0)
      Return 
      End

 
