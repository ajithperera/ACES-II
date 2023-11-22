#ifndef _ACTIVE_SPACE_COM_
#define _ACTIVE_SPACE_COM_
c active_space.com : begin
      Parameter(Max_xp=100)
      integer Active_oo(8,2),Active_vv(8,2)
      integer Pactive_oo(100,8,2),Pactive_vv(100,8,2)
      integer Ioff_active_oo(8,2),Ioff_active_vv(8,2)
      integer Pioff_active_oo(100,8,2),Pioff_active_vv(100,8,2)
      Double Precision Oo_threshold,Vv_threshold
      common /actvsp_info/Active_oo,Pactive_oo,Active_vv,Pactive_vv,
     +                    Ioff_active_oo,Pioff_active_oo,
     +                    Ioff_active_vv,Pioff_active_vv,
     +                    Oo_threshold,Poo_threshold,
     +                    Vv_threshold,Pvv_threshold,Eta_val(Max_xp),
     +                    E_k(Max_xp),E_ks(Max_xp)
c active_space.com: end
#endif /* _SYMOFF_COM_ */
