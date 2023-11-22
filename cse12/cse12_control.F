      Subroutine Cse12_control(Work,Maxcor,Iuhf)

      Implicit Double precision (A-H,O-Z)

      Dimension Work(Maxcor)
      Logical Frac_occ,Open_shell

#include "sym.com"
#include "syminf.com"
#include "info.com"
#include "cse12_vars.com"

      Data Ione /1/
 
      If (Nirrep .Gt. 1) Then
      Write(6,"(a,a)") " Abelian point group symmetry is not allowed",
     +                 " use symmetry=off instead."
      Call Errex
      Endif 
      If (Iuhf .EQ.0 ) Then
     
      Write(6,"(a,a)") " REF=RHF is not supported (spin-orbital",
     +                 " only).Use REF=UHF instead."
      Call Errex
      Endif 

      Nocc_a = Pop(1,1)
      Nocc_b = Pop(1,2)
      Nvrt_a = Vrt(1,1)
      Nvrt_b = Vrt(1,2)

      Nbasis_a   = Nocco(1) + Nvrto(1)
      Nbasis_b   = Nocco(2) + Nvrto(2)
      Nbasis     = Max(Nbasis_a,Nbasis_b)
      Call Getrec(20,"JOBARC","NBASTOT",Ione,Naobfns)

      If (Nbasis .Ne. Naobfns) Then 
         Write(6,"(a,a)") "Caution: Spherical harmonics are being",
     +                    " used."
      Endif 

      Call Set_prelims(Work,Maxcor,Iuhf,Nbasis,Nocc_a,Nocc_b,
     +                 Nvrt_a,Nvrt_b,Frac_occ)

      Open_shell = (Nocc_a .Ne. Nocc_b)
      If (Frac_occ) Then
         Nocc_a = Nbasis 
         Nocc_b = Nbasis 
         Nvrt_a = Nbasis
         Nvrt_b = Nbasis
      Endif 

      If (Active_space) Then
         Call set_active_space(Work,Maxcor,Nocc_a,Nocc_b,Nvrt_a,
     +                          Nvrt_b,Length) 
      Endif 

      Ndim_oo_aa = Nocc_a * Nocc_a
      Ndim_oo_bb = Nocc_b * Nocc_b
      Ndim_oo_ab = Nocc_a * Nocc_b
      Ndim_vv_aa = Nvrt_a * Nvrt_a
      Ndim_vv_bb = Nvrt_b * Nvrt_b
      Ndim_vv_ab = Nvrt_a * Nvrt_b
      Ndim_ov_aa = Nocc_a * Nvrt_a
      Ndim_ov_bb = Nocc_b * Nvrt_b
      Ndim_ov_ab = Nocc_a * Nvrt_b
      Ndim_ov_ba = Nocc_b * Nvrt_a

C <ab|ij> MO ints.

      IW0_aaaa = Length   + Ione 
      IW0_bbbb = IW0_aaaa + Ndim_oo_aa*Ndim_vv_aa
      IW0_abab = IW0_bbbb + Ndim_oo_bb*Ndim_vv_bb
      IW1_aaaa = IW0_abab + Ndim_oo_ab*Ndim_vv_ab

C <mn|ij> MO ints.

      Iw1_bbbb = IW1_aaaa + Ndim_oo_aa*Ndim_oo_aa
      Iw1_abab = Iw1_bbbb + Ndim_oo_bb*Ndim_oo_bb
      Iw2_aaaa = Iw1_abab + Ndim_oo_ab*Ndim_oo_ab
  
C <ab|ef> MO ints.

      IW2_bbbb = Iw2_aaaa + Ndim_vv_aa*Ndim_vv_aa
      IW2_abab = Iw2_bbbb + Ndim_vv_bb*Ndim_vv_bb
      IW3_aaaa = Iw2_abab + Ndim_vv_ab*Ndim_vv_ab

C <mb|ej> MO ints.

      IW3_bbbb = IW3_aaaa + Ndim_ov_aa*Ndim_ov_aa
      IW3_abab = IW3_bbbb + Ndim_ov_bb*Ndim_ov_bb
      IW3_baba = IW3_abab + Ndim_ov_ab*Ndim_ov_ba
      IW3_abba = IW3_baba + Ndim_ov_ba*Ndim_ov_ab
      IW3_baab = IW3_abba + Ndim_oo_aa*Ndim_vv_bb
      IW__aaaa = IW3_baab + Ndim_oo_bb*Ndim_vv_aa

C <mn|ef> MO ints.

      IW__bbbb = IW__aaaa + Ndim_oo_aa*Ndim_vv_aa
      IW__abab = IW__bbbb + Ndim_oo_bb*Ndim_vv_bb
      Iend     = IW__abab + Ndim_oo_ab*Ndim_vv_ab

C <ab|ci> MO ints 

      Iw4_aaaa = Iend
      Iw4_bbbb = Iw4_aaaa + Ndim_vv_aa*Ndim_ov_aa
      Iw4_abab = Iw4_bbbb + Ndim_vv_bb*Ndim_ov_bb
      Iw4_abba = Iw4_abab + Ndim_vv_ab*Ndim_ov_ba
      Iw5_aaaa = Iw4_abba + Ndim_vv_ab*Ndim_ov_ab

C <ij|ka> Mo ints

      Iw5_bbbb = Iw5_aaaa + Ndim_oo_aa*Ndim_ov_aa
      Iw5_abab = Iw5_bbbb + Ndim_oo_bb*Ndim_ov_bb
      Iw5_abba = Iw5_abab + Ndim_oo_ab*Ndim_ov_ab
      Iend     = Iw5_abba + Ndim_oo_ab*Ndim_ov_ba

C Cmo_a(mu,p) and Cmo_b(nu,q) MO ints.

      Icmo_a   = Iend 
      Icmo_b   = Icmo_a + Nbasis*Naobfns 
      Iend     = Icmo_b + Nbasis*Naobfns 
 
      If (Iend .Ge. Maxcor) Call Insmem("ccd_light_control",Iend,
     +                                   Maxcor)
      Call Getrec(20,"JOBARC","SCFEVCA0",Nbasis*Naobfns,
     +            Work(Icmo_a))
      Call Getrec(20,"JOBARC","SCFEVCB0",Nbasis*Naobfns,
     +            Work(Icmo_b))

      Memleft = Maxcor - Iend
#ifdef _DEBUG_LVLM
      Write(6,"(a)") " Enter Form_moints"
      call output(Work(Icmo_a),1,Naobfns,1,Nbasis,Naobfns,Nbasis,1)
      call output(Work(Icmo_b),1,Naobfns,1,Nbasis,Naobfns,Nbasis,1)
#endif

      Call Form_moints(Work(IW0_aaaa),Work(IW0_bbbb),
     +     Work(IW0_abab),Work(IW1_aaaa),Work(IW1_bbbb),
     +     Work(IW1_abab),Work(IW2_aaaa),Work(IW2_bbbb),
     +     Work(IW2_abab),Work(IW3_aaaa),Work(IW3_bbbb),
     +     Work(IW3_abab),Work(IW3_baba),Work(IW3_abba),
     +     Work(IW3_baab),Work(IW__aaaa),Work(IW__bbbb),
     +     Work(IW__abab),Work(Iw4_aaaa),Work(Iw4_bbbb),
     +     Work(Iw4_abab),Work(Iw4_abba),Work(Iw5_aaaa),
     +     Work(Iw5_bbbb),Work(Iw5_abab),Work(Iw5_abba),
     +     Work(Icmo_a),Work(Icmo_b),
     +     Nocc_a,Nocc_b,Nvrt_a,Nvrt_b,
     +     Nbasis,Naobfns,Work(Iend),Memleft,Frac_occ,Repuls)

C Fock_a(p,q),Fock_b(p,q),Fockoo_a(i,j),Fockoo_b(i,j),Fockvv_a(a,b)
C Fockvv_b(a,b),Fockoo_od_a(i,j),Fockoo_od_b(i,j),Fockvv_od_a(a,b)
C Fockcc_od_b(a,b).

      Ifock_a    = Iend
      Ifock_b    = Ifock_a + Naobfns*Naobfns
      Ifockoo_a  = Ifock_b + Naobfns*Naobfns  
      Ifockoo_b  = Ifockoo_a + Nocc_a*Nocc_a 
      Ifockvv_a  = Ifockoo_b + Nocc_b*Nocc_b
      Ifockvv_b  = Ifockvv_a + Nvrt_a*Nvrt_a 
      Ifockov_a  = Ifockvv_b + Nvrt_b*Nvrt_b
      Ifockov_b  = Ifockov_a + Nvrt_a*Nocc_a 
      Iend       = Ifockov_b + Nvrt_b*Nocc_b
     
      Ifockoo_od_a = Iend
      Ifockoo_od_b = Ifockoo_od_a + Nocc_a*Nocc_a
      Ifockvv_od_a = Ifockoo_od_b + Nocc_b*Nocc_b
      Ifockvv_od_b = Ifockvv_od_a + Nvrt_a*Nvrt_a
      Iend         = Ifockvv_od_b + Nvrt_b*Nvrt_b

      Memleft = Maxcor - Iend

      If (Iend .Ge. Maxcor) Call Insmem("ccd_light_control",Iend,
     +                                   Maxcor)

      Call Prep_fock(Work(Ifock_a),Work(Ifock_b),Work(Ifockoo_a),
     +               Work(Ifockoo_b),Work(Ifockvv_a),Work(Ifockvv_b),
     +               Work(Ifockov_a),Work(Ifockov_b), 
     +               Work(Ifockoo_od_a),Work(Ifockoo_od_b),
     +               Work(Ifockvv_od_a),Work(Ifockvv_od_b),
     +               Work(IW0_aaaa),Work(IW0_bbbb),Work(Iw0_abab),
     +               E1_aa,E1_bb,Nbasis,Naobfns,Nocc_a,Nocc_b,
     +               Nvrt_a,Nvrt_b,Work(Iend),Memleft,Frac_occ,
     +               Open_shell)

      If (Frac_occ) Then
          Call Get_ref_energy(Work(IW0_aaaa),Work(IW0_bbbb),
     +                        Work(IW0_abab),Work(Ifock_a),
     +                        Work(Ifock_b),E1_aa,E1_bb,Nbasis,
     +                        Naobfns,Nocc_a,Nocc_b,Nvrt_a,
     +                        Nvrt_b,Work(Iend),Memleft,
     +                        Repuls,Ref_energy)
      Endif 

C W(mn,ij) intermediate 
      
      Imnij_aaaa = Iend 
      Imnij_bbbb = Imnij_aaaa + Ndim_oo_aa*Ndim_oo_aa
      Imnij_abab = Imnij_bbbb + Ndim_oo_bb*Ndim_oo_bb
      Iabef_aaaa = Imnij_abab + Ndim_oo_ab*Ndim_oo_ab

C W(ab,ef) intermediate 

      Iabef_bbbb = Iabef_aaaa + Ndim_vv_aa*Ndim_vv_aa
      Iabef_abab = Iabef_bbbb + Ndim_vv_bb*Ndim_vv_bb
      Imbej_aaaa = Iabef_abab + Ndim_vv_ab*Ndim_vv_ab

C W(mb,ej) intermediate 

      Imbej_bbbb = Imbej_aaaa + Ndim_ov_aa*Ndim_ov_aa
      Imbej_abab = Imbej_bbbb + Ndim_ov_bb*Ndim_ov_bb
      Imbej_baba = Imbej_abab + Ndim_ov_ab*Ndim_ov_ba
      Imbej_abba = Imbej_baba + Ndim_ov_ba*Ndim_ov_ab
      Imbej_baab = Imbej_abba + Ndim_oo_aa*Ndim_vv_bb
      Iend       = Imbej_baab + Ndim_oo_bb*Ndim_vv_aa

C T2(ab,ij) 

      It2_aaaa   = Iend 
      It2_bbbb   = It2_aaaa + Ndim_oo_aa*Ndim_vv_aa
      It2_abab   = It2_bbbb + Ndim_oo_bb*Ndim_vv_bb
      Iend       = It2_abab + Ndim_oo_ab*Ndim_vv_ab

C Tau(ab,ij)

      Itu_aaaa   = Iend
      Itu_bbbb   = Itu_aaaa + Ndim_oo_aa*Ndim_vv_aa
      Itu_abab   = Itu_bbbb + Ndim_oo_bb*Ndim_vv_bb
      Itp_aaaa   = Itu_abab + Ndim_oo_ab*Ndim_vv_ab

C Tap(ab,ij)

      Itp_bbbb   = Itp_aaaa + Ndim_oo_aa*Ndim_vv_aa
      Itp_abab   = Itp_bbbb + Ndim_oo_bb*Ndim_vv_bb
      It1_aa     = Itp_abab + Ndim_oo_ab*Ndim_vv_ab

C T1(a,i)

      It1_bb    = It1_aa + Nvrt_a*Nocc_a
      Iend      = It1_bb + Nvrt_b*Nocc_b

C F(m,i), F(a,e), F(m,e)  intermediate 

      Ifmi_a = Iend
      Ifmi_b = Ifmi_a + Ndim_oo_aa
      Ifae_a = Ifmi_b + Ndim_oo_bb
      Ifae_b = Ifae_a + Ndim_vv_aa
      Ifme_a = Ifae_b + Ndim_vv_bb
      Ifme_b = Ifme_a + Ndim_ov_aa
      Iend   = Ifme_b + Ndim_ov_bb

      If (Iend .Ge. Maxcor) Call Insmem("ccd_light_control",Iend,
     +                                   Maxcor)
      Memleft = Maxcor - Iend 

      Call Cse12_driver(Work(IW0_aaaa),Work(IW0_bbbb),
     +     Work(IW0_abab),Work(IW1_aaaa),Work(IW1_bbbb),
     +     Work(IW1_abab),Work(IW2_aaaa),Work(IW2_bbbb),
     +     Work(IW2_abab),Work(IW3_aaaa),Work(IW3_bbbb),
     +     Work(IW3_abab),Work(IW3_baba),Work(IW3_abba),Work(IW3_baab),
     +     Work(IW__aaaa),Work(IW__bbbb),Work(IW__abab),
     +     Work(Iw4_aaaa),Work(Iw4_bbbb),Work(IW4_abab),
     +     Work(IW4_abba),Work(Iw5_aaaa),Work(Iw5_bbbb),
     +     Work(IW5_abab),Work(IW5_abba),
     +     Work(Imnij_aaaa),Work(Imnij_bbbb),Work(Imnij_abab),
     +     Work(Iabef_aaaa),Work(Iabef_bbbb),Work(Iabef_abab),
     +     Work(Imbej_aaaa),Work(Imbej_bbbb),Work(Imbej_abab),
     +     Work(Imbej_baba),Work(Imbej_abba),Work(Imbej_baab),
     +     Work(It2_aaaa),Work(It2_bbbb),Work(It2_abab),
     +     Work(Itu_aaaa),Work(Itu_bbbb),Work(Itu_abab),
     +     Work(Itp_aaaa),Work(Itp_bbbb),Work(Itp_abab),
     +     Work(It1_aa),Work(It1_bb),
     +     Work(Ifockoo_a),Work(Ifockoo_b),Work(Ifockvv_a),
     +     Work(Ifockvv_b),Work(Ifockov_a),Work(Ifockov_b),
     +     Work(Ifockoo_od_a),Work(Ifockoo_od_b),
     +     Work(Ifockvv_od_a),Work(Ifockvv_od_b),
     +     Work(Ifmi_a),Work(Ifmi_b),Work(Ifae_a),
     +     Work(Ifae_b),Work(Ifme_a),Work(Ifme_b),
     +     Nocc_a,Nocc_b,Nvrt_a,Nvrt_b,Work(Ione),
     +     Work(Iend),Memleft,Nbasis,Length,Ref_energy,
     +     Frac_occ,Open_shell)

      Return
      End

