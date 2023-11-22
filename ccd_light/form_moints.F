      Subroutine Form_moints(W0_aa,W0_bb,W0_ab,W1_aa,W1_bb,W1_ab,
     +                       W2_aa,W2_bb,W2_ab,W3_aaaa,W3_bbbb,
     +                       W3_abab,W3_baba,W3_abba,W3_baab,
     +                       W_aa,W_bb,W_ab,Cmo_a,Cmo_b,Nocc_a,
     +                       Nocc_b,Nvrt_a,Nvrt_b,Nbasis,Naobfns,
     +                       Temp,Maxcor,Frac_occ)

      Implicit Double Precision(A-H,O-Z)

      Dimension Buf(600),Ibuf(600)
      Dimension Temp(Maxcor)
      Logical Frac_occ

      Dimension W_pxxx(Nbasis,Naobfns,Naobfns,Naobfns)
      Dimension W_pqxx(Nbasis,Nbasis,Naobfns,Naobfns)
      Dimension W_pqrx(Nbasis,Nbasis,Nbasis,Naobfns)
      Dimension W_pqrs_aaaa(Nbasis,Nbasis,Nbasis,Nbasis)
      Dimension W_pqrs_bbbb(Nbasis,Nbasis,Nbasis,Nbasis)
      Dimension W_pqrs_aabb(Nbasis,Nbasis,Nbasis,Nbasis)
      Dimension W_pqrs_bbaa(Nbasis,Nbasis,Nbasis,Nbasis)

      Dimension W0_aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension W0_bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension W0_ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Dimension W1_aa(Nocc_a,Nocc_a,Nocc_a,Nocc_a)
      Dimension W1_bb(Nocc_b,Nocc_b,Nocc_b,Nocc_b)
      Dimension W1_ab(Nocc_a,Nocc_b,Nocc_a,Nocc_b)

      Dimension W2_aa(Nvrt_a,Nvrt_a,Nvrt_a,Nvrt_a)
      Dimension W2_bb(Nvrt_b,Nvrt_b,Nvrt_b,Nvrt_b)
      Dimension W2_ab(Nvrt_a,Nvrt_b,Nvrt_a,Nvrt_b)

      Dimension W3_aaaa(Nocc_a,Nvrt_a,Nvrt_a,Nocc_a)
      Dimension W3_bbbb(Nocc_b,Nvrt_b,Nvrt_b,Nocc_b)
      Dimension W3_abab(Nocc_a,Nvrt_b,Nvrt_a,Nocc_b)
      Dimension W3_baba(Nocc_b,Nvrt_a,Nvrt_b,Nocc_a)
      Dimension W3_abba(Nocc_a,Nvrt_b,Nocc_a,Nvrt_b)
      Dimension W3_baab(Nocc_b,Nvrt_a,Nocc_b,Nvrt_a)

      Dimension W_aa(Nocc_a,Nocc_a,Nvrt_a,Nvrt_a)
      Dimension W_bb(Nocc_b,Nocc_b,Nvrt_b,Nvrt_b)
      Dimension W_ab(Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Dimension Wmnij_1(Nocc_a,Nocc_a,Nocc_a,Nocc_a)
      Dimension Wmnij_2(Nocc_b,Nocc_b,Nocc_b,Nocc_b)
      Dimension Wmnij_3(Nocc_a,Nocc_b,Nocc_a,Nocc_b)

      Dimension Wabef_1(Nvrt_a,Nvrt_a,Nvrt_a,Nvrt_a)
      Dimension Wabef_2(Nvrt_b,Nvrt_b,Nvrt_b,Nvrt_b)
      Dimension Wabef_3(Nvrt_a,Nvrt_b,Nvrt_a,Nvrt_b)

      Dimension Cmo_a(Naobfns,Nbasis)
      Dimension Cmo_b(Naobfns,Nbasis)

      Maxdim   = Max(Naobfns,Nbasis)
      Maxdim2  = Maxdim*Maxdim
      Maxdim4  = Maxdim2*Maxdim2
      Naobfns4 = Naobfns*Naobfns*Naobfns*Naobfns
      Nbasis4  = Nbasis*Nbasis*Nbasis*Nbasis

#ifdef _DEBUG_LVL0
      Write(6,"(a,7(1x,I5))") "@-Entry to Form_moints:",
     +                        Maxdim,Naobfns,Nbasis,
     +                       Nocc_a,Nocc_b,Nvrt_a,Nvrt_b
      Write(6,"(a,ES8.2E2)") "Maxcor :", Dble(Maxcor)
      Write(6,"(a,ES8.2E2)") "Maxdim4:", Dble(Maxdim4) 
#endif
      Iw_xxxx  = 1
      Iw_pxxx  = Iw_xxxx + Maxdim4
      Iw_pqxx  = Iw_pxxx + Maxdim4
      Iw_pqrx  = Iw_pqxx + Maxdim4
      Iend     = Iw_pqrx + Maxdim4

      Ipqrs_aaaa  = Iend 
      Ipqrs_bbbb  = Ipqrs_aaaa + Maxdim4
      Ipqrs_aabb  = Ipqrs_bbbb + Maxdim4
      Ipqrs_bbaa  = Ipqrs_aabb + Maxdim4
      Iend        = Ipqrs_bbaa 

CSSS      Ipqrs_abab  = Ipqrs_bbaa + Maxdim4  
CSSS      Ipqrs_baba  = Ipqrs_abab + Maxdim4
CSSS      Iend        = Ipqrs_baba + Maxdim4

      If (Iend .Ge. Maxcor) Call Insmem("Form_moints",Iend,Maxcor)

C Get I(mu,nu,lambda,sigma)

#ifdef _DEBUG_LVL0
      Write(6,"(a)") "Enter get2ints,Ipqrs-aaaa,-aabb,-bbaa"
      write(6,"(4(1x,ES8.2E2))")Dble(Ipqrs_aaaa),Dble(Ipqrs_bbbb), 
     +                          Dble(Ipqrs_aabb),Dble(Ipqrs_bbaa)
      Write(6,*) 
#endif
      Call Dzero(Temp(Iw_xxxx),Naobfns4)
      Call Get2ints(Temp(Iw_xxxx),Buf,Ibuf,Naobfns)

C I(mu,nu,lambda,sigma) -> I(p,q,rs) all allowed spin combinations.
C Mulliken ordered.

#ifdef _DEBUG_LVL0
      call checksum("Ao-ints   :",Temp(Iw_xxxx),Naobfns4)
#endif
      Call form_wpqrs_aaaa(Temp(Iw_xxxx),Temp(Iw_pxxx),Temp(Iw_pqxx),
     +                     Temp(Iw_pqrx),Temp(Ipqrs_aaaa),Cmo_a,
     +                     Nbasis,Naobfns) 
#ifdef _DEBUG_LVL0
      Write(6,"(a)") " Enter Form_moints"
      call checksum("Ipqrs_aaaa:",Temp(Ipqrs_aaaa),Naobfns4)
#endif
      Call form_wpqrs_aaaa(Temp(Iw_xxxx),Temp(Iw_pxxx),Temp(Iw_pqxx),
     +                     Temp(Iw_pqrx),Temp(Ipqrs_bbbb),Cmo_b,
     +                     Nbasis,Naobfns)
#ifdef _DEBUG_LVL0
      call checksum("Ipqrs_bbbb:",Temp(Ipqrs_bbbb),Naobfns4)
#endif
      Call form_wpqrs_aabb(Temp(Iw_xxxx),Temp(Iw_pxxx),Temp(Iw_pqxx),
     +                     Temp(Iw_pqrx),Temp(Ipqrs_aabb),Cmo_a,Cmo_b,
     +                     Nbasis,Naobfns)
#ifdef _DEBUG_LVL0
      call checksum("Ipqrs_aabb:",Temp(Ipqrs_aabb),Naobfns4)
#endif
CSSS      Call form_wpqrs_aabb(Temp(Iw_xxxx),Temp(Iw_pxxx),Temp(Iw_pqxx),
CSSS     +                     Temp(Iw_pqrx),Temp(Ipqrs_bbaa),Cmo_b,Cmo_a,
CSSS     +                     Nbasis,Naobfns)
#ifdef _DEBUG_LVL0
CSSS      call checksum("Ipqrs_bbaa:",Temp(Ipqrs_bbaa),Naobfns4)
#endif
      Call form_wpqrs_abab(Temp(Iw_xxxx),Temp(Iw_pxxx),Temp(Iw_pqxx),
     +                     Temp(Iw_pqrx),Temp(Ipqrs_abab),Cmo_a,Cmo_b,
     +                     Nbasis,Naobfns)
#ifdef _DEBUG_LVL0
      call checksum("Ipqrs_abab:",Temp(Ipqrs_abab),Naobfns4)
#endif 
CSSS      Call form_wpqrs_abab(Temp(Iw_xxxx),Temp(Iw_pxxx),Temp(Iw_pqxx),
CSSS     +                     Temp(Iw_pqrx),Temp(Ipqrs_baba),Cmo_b,Cmo_a,
CSSS     +                     Nbasis,Naobfns)
#ifdef _DEBUG_LVL0
CSSS      call checksum("Ipqrs_baba:",Temp(Ipqrs_baba),Naobfns4)
#endif
C I(p,q,r,s) -> I(pr,qs) where p,q,r,s can be oooo,vvvv,ovov, oovv
C In Dirac ordering and antisymmetrized when allowed. 

      Call form_wmnij(Temp(Ipqrs_aaaa),Temp(Ipqrs_bbbb),
     +                Temp(Ipqrs_aabb),W1_aa,W1_bb,W1_ab,
     +                Nbasis,Nocc_a,Nocc_b)

      If (Frac_occ) Then

C Wmnij->Wabij 

         Call Dcopy(Nbasis4,W1_aa,1,W0_aa,1)
         Call Dcopy(Nbasis4,W1_bb,1,W0_bb,1)
         Call Dcopy(Nbasis4,W1_ab,1,W0_ab,1)

C Wmnij->Wabef 
         Call Dcopy(Nbasis4,W1_aa,1,W2_aa,1)
         Call Dcopy(Nbasis4,W1_bb,1,W2_bb,1)
         Call Dcopy(Nbasis4,W1_ab,1,W2_ab,1)

C Wmnij->Wmnef
         Call Dcopy(Nbasis4,W1_aa,1,W_aa,1)
         Call Dcopy(Nbasis4,W1_bb,1,W_bb,1)
         Call Dcopy(Nbasis4,W1_ab,1,W_ab,1)

         Call form_wmbej_fo(Temp(Ipqrs_aaaa),Temp(Ipqrs_bbbb),
     +                      Temp(Ipqrs_aabb),W3_aaaa,W3_bbbb,
     +                      W3_abab,W3_baba,W3_abba,W3_baab,
     +                      Nbasis,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)
         Return 
      Endif 

      Call form_wabij(Temp(Ipqrs_aaaa),Temp(Ipqrs_bbbb),
     +                Temp(Ipqrs_aabb),W0_aa,W0_bb,W0_ab,
     +                Nbasis,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Call form_wabef(Temp(Ipqrs_aaaa),Temp(Ipqrs_bbbb),
     +                Temp(Ipqrs_aabb),W2_aa,W2_bb,W2_ab,
     +                Nbasis,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Call form_wmnef(Temp(Ipqrs_aaaa),Temp(Ipqrs_bbbb),
     +                Temp(Ipqrs_aabb),W_aa,W_bb,W_ab,
     +                Nbasis,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Call form_wmbej(Temp(Ipqrs_aaaa),Temp(Ipqrs_bbbb),
     +                Temp(Ipqrs_aabb),W3_aaaa,W3_bbbb,
     +                W3_abab,W3_baba,W3_abba,W3_baab,
     +                Nbasis,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)
      Return
      End 


