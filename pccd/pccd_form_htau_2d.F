      Subroutine Pccd_form_htau_2d(Htau_pq,Htau_qp,Hoo_pq,Hoo_qp,Hvv_pq,
     +                             Hvv_qp,Hvo,Hov,Dpq,Dhf,Dcc,Work,Nocc,
     +                             Nvrt,Maxcor,Nbas,E)
   
      Implicit Double Precision(A-H,O-Z)
      Character*4 Spin
      Logical pCCD,CCD,LCCD

      Dimension Htau_pq(Nbas,Nbas)
      Dimension Htau_qp(Nbas,Nbas)
      Dimension Dpq(Nbas,Nbas)
      Dimension Dhf(Nbas,Nbas)
      Dimension Dcc(Nbas,Nbas)
      Dimension Hoo_pq(Nocc,Nocc)
      Dimension Hoo_qp(Nocc,Nocc)
      Dimension Hvv_pq(Nvrt,Nvrt)
      Dimension Hvv_qp(Nvrt,Nvrt)
      Dimension Hvo(Nvrt*Nocc)
      Dimension Hov(Nocc*Nvrt)
      Dimension Work(Maxcor)
    
      Common /Calc/pCCD,CCD,LCCD

      Data One,Onem,Dnull,Half,Quart,Two /1.0D0,-1.0D0,0.0D0,0.50D0,
     +                                0.25D0,2.0D0/
#ifndef _HF_ONLY
C The OO block

      List_v = 13
      List_g = 113
      Call Pccd_form_htau_2d_hhhh(Hoo_pq,Hoo_qp,Work,Maxcor,
     +                            Nocc,Nbas,List_v,List_g)
      List_v = 16
      List_g = 116
      Call Pccd_form_htau_2d_hhpp_hh(Hoo_pq,Hoo_qp,Work,Maxcor,
     +                               Nocc,Nbas,List_v,List_g)

      List_v = 23
      List_g = 123
      Call Pccd_form_htau_2d_phph_hh(Hoo_pq,Hoo_qp,Work,Maxcor,
     +                               Nocc,Nbas,List_v,List_g)

      List_v = 25
      List_g = 125
      Call Pccd_form_htau_2d_phph_hh(Hoo_pq,Hoo_qp,Work,Maxcor,
     +                               Nocc,Nbas,List_v,List_g)

      List_v = 18
      List_g = 118
      Call Pccd_form_htau_2d_phph_hh(Hoo_pq,Hoo_qp,Work,Maxcor,
     +                              Nocc,Nbas,List_v,List_g)

C The VV block

      List_v = 233
      List_g = 133
      Call Pccd_form_htau_2d_pppp(Hvv_pq,Hvv_qp,Work,Maxcor,
     +                            Nvrt,Nbas,List_v,List_g)

      List_v = 16
      List_g = 116
      Call Pccd_form_htau_2d_hhpp_pp(Hvv_pq,Hvv_qp,Work,Maxcor,
     +                               Nvrt,Nbas,List_v,List_g)

      List_v = 23
      List_g = 123
      Call Pccd_form_htau_2d_phph_pp(Hvv_pq,Hvv_qp,Work,Maxcor,
     +                               Nvrt,Nbas,List_v,List_g)

      List_v = 25
      List_g = 125
      Call Pccd_form_htau_2d_phph_pp(Hvv_pq,Hvv_qp,Work,Maxcor,
     +                               Nvrt,Nbas,List_v,List_g)

      List_v = 18
      List_g = 118
      Call Pccd_form_htau_2d_phph_pp(Hvv_pq,Hvv_qp,Work,Maxcor,
     +                               Nvrt,Nbas,List_v,List_g)

C The OV and VO blocks

#ifdef _NO_NRML_ORD
      Fact    = Half
      List_g1 = 125
      List_g2 = 113 

      Call Pccd_form_htau_2d_hf(Hvo,Hov,Dpq,Dhf,Dcc,Work,Maxcor,
     +                          Nocc,Nvrt,Nbas,List_g1,List_g2,
     +                          Fact)
#endif 

      List_v1 = 30
      List_v2 = 10
      List_g  = 116
      Call Pccd_form_htau_2d_hhpp_ph(Hov,Hvo,Work,Maxcor,Nocc,
     +                               Nvrt,Nbas,List_v1,List_v2,
     +                               List_g)

      List_v = 10
      List_g = 113

      Call Pccd_form_htau_2d_hhhh_ph(Hvo,Work,Maxcor,Nocc,
     +                               Nvrt,Nbas,List_v,List_g)

      List_v  = 30
      List_g  = 133

      Call Pccd_form_htau_2d_pppp_ph(Hov,Work,Maxcor,Nocc,
     +                               Nvrt,Nbas,List_v,List_g)

C Resorting List=118,Gamma(Ia,Bj) stored as BJ,ai-> Bi,aJ
C Resorting List=123,Gamma(IA,JB) stored as AI,BI-> AJ,BI
C Resorting List=125,Gamma(iA,jB) stored as Ai,Bj-> Aj,Bi
CSSS      Call Pccd_check_ijka(Work,Maxcor,7,10)

      Call Pccd_sortgam(Work,Maxcor,0)

      List_v1 = 30
      List_v2 = 7
      List_g  = 123
      Fact    = One

      Call Pccd_form_htau_2d_phph_ph_0(Hov,Hvo,Work,Maxcor,Nocc,
     +                                 Nvrt,Nbas,List_v1,List_v2,
     +                                 List_g,Fact)

      List_v1 = 30
      List_v2 = 10
      List_g  = 125
      Fact    = One

      Call Pccd_form_htau_2d_phph_ph_1(Hov,Hvo,Work,Maxcor,Nocc,
     +                                 Nvrt,Nbas,List_v1,List_v2,
     +                                 List_g,Fact)
C This contribution is zero for pCCD.

      If (.Not. Pccd) Then 
         List_v1 = 30 
         List_v2 = 10
         List_g  = 118
         Fact    = One

         Call Pccd_form_htau_2d_phph_ph_2(Hov,Hvo,Work,Maxcor,Nocc,
     +                                    Nvrt,Nbas,List_v1,List_v2,
     +                                    List_g,Fact)
      Endif 

      Call Dscal(Nvrt*Nocc,Onem,Hov,1)
#else
      Call Dzero(Hvv_pq,Nvrt*Nvrt)
      Call Dzero(Hvv_qp,Nvrt*Nvrt)
      Call Dzero(Hov,Nvrt*Nocc)
      Call Dzero(Hvo,Nvrt*Nocc)
#endif 

#ifdef _NO_NRML_ORD
      List_v = 13
      Call Pccd_form_htau_2d_hf_hh_1(Hoo_qp,Dpq,Dhf,Dcc,Work,Maxcor,
     +                               Nocc,Nvrt,Nbas,List_v,E)
#ifndef _HF_ONLY
      Fact   = Quart
      List_v = 25
      Spin   = "ABAB"
      Call Pccd_form_htau_2d_hf_hh_2(Hoo_qp,Dpq,Dhf,Dcc,Work,Maxcor,
     +                               Nocc,Nvrt,Nbas,List_v,Fact,Spin)

      Call Pccd_form_htau_2d_hf_pp_2(Hvv_qp,Dpq,Dhf,Dcc,Work,Maxcor,
     +                               Nocc,Nvrt,Nbas,List_v,Fact,Spin)

C This contribution is zero.
CSSS      List_v = 233
CSSS      Call Pccd_form_htau_2d_hf_pp_1(Hvv_qp,Dpq,Dhf,Dcc,Work,Maxcor,
CSSS     +                               Nocc,Nvrt,Nbas,List_v,Fact)

      List_v = 23 
      Spin   = "ABAB"
      Call Pccd_form_htau_2d_hf_hh_2(Hoo_qp,Dpq,Dhf,Dcc,Work,Maxcor,
     +                               Nocc,Nvrt,Nbas,List_v,Fact,Spin)
      Call Pccd_form_htau_2d_hf_pp_2(Hvv_qp,Dpq,Dhf,Dcc,Work,Maxcor,
     +                               Nocc,Nvrt,Nbas,List_v,Fact,Spin)
#endif 
#endif 

#ifdef _DEBUG_LVLM
      Call Pccd_zero_htau_2d(Hoo_qp,Hoo_pq,Hvv_qp,Hvv_pq,Hov,Hvo,
     +                       Nocc,Nvrt,Nbas)
#endif 
      Call Pccd_form_htau(Htau_qp,Htau_pq,Hoo_qp,Hoo_pq,Hvv_qp,Hvv_pq,
     +                    Hvo,Hov,Work,Maxcor,Nocc,Nvrt,Nbas)
#ifdef _DEBUG_LVLM
      Write(6,"(2a)") "The orbtial rotation gradient matrix:",
     +                " post G(pppp)"
      call output(Htau_pq,1,Nbas,1,Nbas,Nbas,Nbas,1)
#endif
      Return
      End 
