      Subroutine Pccd_form_htau_1d_pqrs(Htau_pq,Htau_qp,Hpq,Dpq,Hoo_pq,
     +                                  Hoo_qp,Hvv_pq,Hvv_qp,Hvo,Hov,
     +                                  
     +                                  Htau_oo,Hoo,Doo,Htau_vv,Hvv,Dvv,
     +                                  Htau_ov,Htau_vo,H_ov,H_vo,D_ov,
     +                                  D_vo,Work,Nocc,Nvrt,Maxcor,
     +                                  Nbas)

      Implicit Double Precision(A-H,O-Z)

      Dimension Htau_pq(Nbas,Nbas)
      Dimension Htau_qp(Nbas,Nbas)
      Dimension Hpq(Nbas,Nbas)
      Dimension Dpq(Nbas,Nbas)
      Dimension Hoo_pq(Nocc,Nocc)
      Dimension Hoo_qp(Nocc,Nocc)
      Dimension Hvv_pq(Nvrt,Nvrt)
      Dimension Hvv_qp(Nvrt,Nvrt)
      Dimension Hvo(Nvrt*Nocc)
      Dimension Hov(Nocc*Nvrt)
      Dimension Work(Maxcor)

      Dimension Htau_oo(Nocc,Nocc)
      Dimension H_oo(Nocc,Nocc)
      Dimension D_oo(Nocc,Nocc)
      Dimension Htau_vv(Nvrt,Nvrt)
      Dimension H_vv(Nvrt,Nvrt)
      Dimension D_vv(Nvrt,Nvrt)
      Dimension Htau_ov(Nocc,Nvrt)
      Dimension Htau_vo(Nvrt,Nocc)
      Dimension H_ov(Nocc,Nvrt)
      Dimension H_vo(Nvrt,Nocc)
      Dimension D_ov(Nocc,Nvrt)
      Dimension D_vo(Nvrt,Nocc)
      Dimension Ioffo(8)
      Dimension Ioffv(8)

#include "sym.com"
#include "sympop.com"
#include "syminf.com"

      Data Ione /1/

C Symmetry pack the Htau_pq,Hpq and Dpq matrices. 

      Call Pccd_symm_pack(Htau_qp,Htau_oo,Htau_vv,Htau_ov,Htau_vo,
     +                  Nocc,Nvrt,Nbas)
      Call Pccd_symm_pack(Hpq,H_oo,H_vv,H_ov,H_vo,Nocc,Nvrt,Nbas)
      Call Pccd_symm_pack(Dpq,D_oo,D_vv,D_ov,D_vo,Nocc,Nvrt,Nbas)

#ifdef _DEBUG_LVLM
      call output(Htau_oo,1,Nocc,1,Nocc,Nocc,Nocc,1)
      call output(Htau_vv,1,Nvrt,1,Nvrt,Nvrt,Nvrt,1)
      call output(H_oo,1,Nocc,1,Nocc,Nocc,Nocc,1)
      call output(D_oo,1,Nocc,1,Nocc,Nocc,Nocc,1)
      call output(H_vv,1,Nvrt,1,Nvrt,Nvrt,Nvrt,1)
      call output(D_vv,1,Nvrt,1,Nvrt,Nvrt,Nvrt,1)
      call output(H_vo,1,Nvrt,1,Nocc,Nvrt,Nocc,1)
      call output(H_ov,1,Nocc,1,Nvrt,Nocc,Nvrt,1)
      call checksum("Htau_oo :",Htau_oo,Nocc*Nocc)
      call checksum("Hoo     :",H_oo,Nocc*Nocc)
      call checksum("Doo     :",D_oo,Nocc*Nocc)
      call checksum("Htau_vv :",Htau_vv,Nvrt*Nvrt)
      call checksum("Hvv     :",H_vv,Nvrt*Nvrt)
      call checksum("Dvv     :",D_vv,Nvrt*Nvrt)
      call checksum("Htau_vo :",Htau_vo,Nvrt*Nocc)
      call checksum("Htau_ov :",Htau_ov,Nvrt*Nocc)
      call checksum("H_vo    :",H_vo,Nvrt*Nocc)
      call checksum("H_ov    :",H_ov,Nvrt*Nocc)
      call checksum("D_vo    :",D_vo,Nvrt*Nocc)
      call checksum("D_ov    :",D_ov,Nvrt*Nocc)
      call checkvv(Htau_qp,Nbas,Nvrt)
#endif 

C Form Htau_pq contributions Htau_pqrs (This is four terns like 
C Htau_pq*Delta(rs) and all possible permutations. 

      I000 = Ione
      I010 = I000 + Nocc*Nocc
      Iend = I010 + Nvrt*Nvrt
      If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_1d_pqrs",
     +                                   Iend,Maxcor)

      Maxcor = Maxcor - Iend 
      Call Pccd_form_htau_1d_pqrs_1(Htau_qp,Htau_oo,Htau_vv,Htau_ov,
     +                              Htau_vo,Work(I000),Work(I010),
     +                              Work(Iend),Ioffo,Ioffv,
     +                              Maxcor,Nocc,Nvrt,Nbas)
      Call Pccd_form_htau_1d_pqrs_2(Hpq,Dpq,H_oo,H_vv,H_ov,H_vo,D_oo,
     +                              Dvv,D_ov,D_vo,Work(Iend),
     +                              Ioffo,Ioffv,Maxcor,Nocc,Nvrt,Nbas)

      Return
      End 
   
   
