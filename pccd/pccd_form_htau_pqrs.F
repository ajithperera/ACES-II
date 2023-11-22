      Subroutine Pccd_form_htau_pqrs(Work,Maxcor,Iuhf)

      Implicit Double Precision(A-H,O-Z)

      Dimension Work(Maxcor)

C VVVV block

      List_v1 = 233
      List_v2 = 16
      List_g1 = 133
      List_g2 = 116
      List_h  = 234
      Call Pccd_form_htau_pqrs_pppp_0(Work,Maxcor,List_v1,List_v2,
     +                                List_g1,List_g2,List_h)

      List_v = 233
      List_g = 133
      List_h = 234
      Call Pccd_form_htau_pqrs_pppp_1(Work,Maxcor,List_v,List_g,List_h)

      List_v1 = 23
      List_v2 = 25
      List_v3 = 18
      List_g1 = 123
      List_g2 = 125 
      List_g3 = 118
      List_h  = 234
      Call Pccd_form_htau_pqrs_pppp_2(Work,Maxcor,List_v1,List_v2,
     +                                List_v3,List_g1,List_g2,
     +                                List_g3,List_h)
      
C OOOO block

      List_v1 = 13
      List_v2 = 16
      List_g1 = 113
      List_g2 = 116
      List_h  = 213
      Call Pccd_form_htau_pqrs_hhhh_0(Work,Maxcor,List_v1,List_v2,
     +                                List_g1,List_g2,List_h)

      List_v = 13
      List_g = 113
      List_h = 213
      Call Pccd_form_htau_pqrs_hhhh_1(Work,Maxcor,List_v,List_g,List_h)

      List_v1 = 23
      List_v2 = 25
      List_v3 = 18
      List_g1 = 123
      List_g2 = 125
      List_g3 = 118
      List_h  = 213
      Call Pccd_form_htau_pqrs_hhhh_2(Work,Maxcor,List_v1,List_v2,
     +                                List_v3,List_g1,List_g2,
     +                                List_g3,List_h)

C VVOO and OOVV blocks

      List_v1 = 16
      List_v2 = 233
      List_g1 = 113
      List_g2 = 116
      List_h  = 216
      Call Pccd_form_htau_pqrs_pphh_0(Work,Maxcor,List_v1,List_v2,
     +                                List_g1,List_g2,List_h)

      List_v1 = 13
      List_v2 = 16
      List_g1 = 116
      List_g2 = 133
      List_h  = 217
      Call Pccd_form_htau_pqrs_hhpp_0(Work,Maxcor,List_v1,List_v2,
     +                                List_g1,List_g2,List_h)

C H(ab,ij)^(t) + H(ij,ab) and H(ij,ab)(t)+H(Hab,ij)

      List_h1 = 216
      List_h2 = 217 
      Call Pccd_symmetrize_hhpp_pphh(Work,Maxcor,List_h1,List_h2)

      List_v1 = 18
      List_v2 = 25
      List_g1 = 116
      List_h  = 216
      Call Pccd_form_htau_pqrs_pphh_1(Work,Maxcor,List_v1,List_v2,
     +                                List_g1,List_h)

      List_v1 = 18
      List_v2 = 19
      List_v3 = 21
      List_g1 = 118
      List_g2 = 123
      List_g3 = 125
      List_h  = 216
      Call Pccd_form_htau_pqrs_pphh_2(Work,Maxcor,List_v1,List_v2,
     +                                List_v3,List_g1,List_g2,
     +                                List_g3,List_h)

      List_v1 = 18
      List_v2 = 25
      List_g1 = 116
      List_h  = 217
      Call Pccd_form_htau_pqrs_hhpp_1(Work,Maxcor,List_v1,List_v2,
     +                                List_g1,List_h)

      List_v1 = 18
      List_v2 = 19
      List_v3 = 21
      List_g1 = 118
      List_g2 = 123 
      List_g3 = 125
      List_h  = 217
      Call Pccd_form_htau_pqrs_hhpp_2(Work,Maxcor,List_v1,List_v2,
     +                                List_v3,List_g1,List_g2,
     +                                List_g3,List_h)

C VOVO block
      List_v1 = 23
      List_v2 = 25
      List_v3 = 18
      List_g1 = 123
      List_g2 = 125
      List_g3 = 118
      List_h  = 225

C H(ai,bj) + H(ai,bj)^t
      Call Pccd_form_htau_pqrs_phph_0(Work,Maxcor,List_v1,List_v2,
     +                                List_v3,List_g1,List_g2,
     +                                List_g3,List_h)
C H(aj,bi) + H(aj,bi)^t
      Call Pccd_form_htau_pqrs_phph_1(Work,Maxcor,List_v1,List_v2,
     +                                List_v3,List_g1,List_g2,
     +                                List_g3,List_h)
C H(aj,bi) + H(aj,bi)^t
      List_v1 = 18
      List_g1 = 116
      Call Pccd_form_htau_pqrs_phph_2(Work,Maxcor,List_v1,List_g1,
     +                                List_h)

C OVVO block

      List_v1 = 18
      List_g1 = 125
      List_v2 = 25 
      List_g2 = 118
      List_h  = 218

C H(ia,bj) + H(ia,bj)^t
      Call Pccd_form_htau_pqrs_hpph_0(Work,Maxcor,List_v1,List_v2,
     +                                List_g1,List_g2,List_h) 

C Form the VOVO antisymmetric list by from lists 225 (Coulomb) and
C 218 (exchange).

      List_h1 = 223
      List_h2 = 225
      List_h3 = 218
      Call Pccd_phph_antisym(Work,Maxcor,List_h1,List_h2,List_h3)

#ifdef _NOSKIP
C These contributions do not need to be explicitly constructed. They
C can be obtained by symmetrizing the contributions computed in _phph_1
C phph_2
      List_v1 = 23
      List_v2 = 25
      List_v3 = 18
      List_g1 = 123
      List_g2 = 125
      List_g3 = 118

      Call Pccd_form_htau_pqrs_phph_3(Work,Maxcor,List_v1,List_v2,
     +                                List_v3,List_g1,List_g2,
     +                                List_g3,List_h)
      List_v1 = 18
      List_g1 = 116
      Call Pccd_form_htau_pqrs_phph_4(Work,Maxcor,List_v1,List_g1,
     +                                List_h)
#endif 

#ifdef _DEBUG_LVL0
      Call Pccd_check_hess(Work,Maxcor)
CSSS      Call Pccd_zero_hess(Work,Maxcor)
#endif 
      Return
      End 
   
   
