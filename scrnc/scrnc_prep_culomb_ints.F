      Subroutine Scrnc_prep_culomb_ints(Coulmb_oo,Coulmb_vv,Coulmb_vo,
     +                                  Iuhf,Len_oo,Len_vv,Len_vo,
     +                                  Ipq_pair_4irep,Irrepx,Irrepr)

      Implicit Double Precision (A-H,O-Z)
      Integer Ipq_pair_4irep

      Dimension Coulmb_oo(Len_oo),Coulmb_vv(Len_vv),Coulmb_vo(Len_vo)

#include "machsp.com"
#include "syminf.com"
#include "sympop.com"
#include "parcc.com"
#include "flags.com"

      List_aixx   = 205
      List_abxx   = 207
      List_ijxx   = 209

      Ioff_oo = 1
      Ioff_vv = 1
      Ioff_vo = 1
#ifdef _DEBUG_LVL0
      Write(6,"(a,a,i5,i)") " Integrals are obtained for pair",
     +                      " and irrep:",
     +                     Ipq_pair_4irep,Irrepr 
#endif 
     
      Do Ispin = 1, Iuhf+1
         Call Getlst(Coulmb_oo(Ioff_oo),Ipq_pair_4irep,1,1,Irrepr,
     +               List_ijxx+Ispin)
         Call Getlst(Coulmb_vv(Ioff_vv),Ipq_pair_4irep,1,1,Irrepr,
     +               List_abxx+Ispin)
         Call Getlst(Coulmb_vo(Ioff_vo),Ipq_pair_4irep,1,1,Irrepr,
     +               List_aixx+Ispin)

         Ioff_oo = Ioff_oo + Irpdpd(Irrepx,21)
         Ioff_vv = Ioff_vv + Irpdpd(Irrepx,19)
         Ioff_vo = Ioff_vo + Irpdpd(Irrepx,9)
      Enddo 

#ifdef _DEBUG_LVL0
      Write(6,*)
      write(6,"(a)") " <oo|J|XX>,<vv|J|XX> and <vo|J|XX> integrals"
      Ioff_oo = 1
      Ioff_vv = 1
      Ioff_vo = 1
      Lenoo  = Irpdpd(Irrepx,21)
      Lenvv  = Irpdpd(Irrepx,19)
      Lenvo  = Irpdpd(Irrepx,9)
      Write(6,"(a,3(1x,i3))")" Lenoo,Lenvv,Lenvo:", Lenoo,Lenvv,Lenvo
      call checksum("<oo|j|xx> Alpha: ",Coulmb_oo(ioff_oo),Lenoo)
      call checksum("<vv|j|xx> Alpha: ",Coulmb_vv(ioff_vv),Lenvv)
      call checksum("<vo|j|xx>:Alpha: ",Coulmb_vo(ioff_vo),Lenvo)
      Ioffoo = Ioff_oo + Irpdpd(Irrepx,21)
      Ioffvv = Ioff_vv + Irpdpd(Irrepx,19)
      Ioffvo = Ioff_vo + Irpdpd(Irrepx,9)
      call checksum("<oo|j|xx> Beta : ",Coulmb_oo(ioff_oo),Lenoo)
      call checksum("<vv|j|xx> Beta : ",Coulmb_vv(ioff_vv),Lenvv)
      call checksum("<vo|j|xx>:Beta : ",Coulmb_vo(ioff_vo),Lenvo)
#endif 
      Return
      End
