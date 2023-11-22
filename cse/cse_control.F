      Subroutine Cse_control(Work,Maxcor,Iuhf)

      Implicit Double precision (A-H,O-Z)

      Dimension Work(Maxcor)
      Logical Frac_occ
      Logical Oneh,Fock

#include "sym.com"
#include "syminf.com"
#include "info.com"
#include "cse_vars.com"

      Data Ione /1/
 
      If (Nirrep .Gt. 1) Then
      Write(6,"(a,a)") " Abelian point group symmetry is not allowed",
     +                    " use symmetry=off instead."
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

      If (Frac_occ) Then
         Nocc_a = Nbasis 
         Nocc_b = Nbasis 
         Nvrt_a = Nbasis
         Nvrt_b = Nbasis
      Endif 

      Ndim_oo_aa = Nocc_a * Nocc_a
      Ndim_oo_bb = Nocc_b * Nocc_b
      Ndim_oo_ab = Nocc_a * Nocc_b
      Ndim_vv_aa = Nvrt_a * Nvrt_a
      Ndim_vv_bb = Nvrt_b * Nvrt_b
      Ndim_vv_ab = Nvrt_a * Nvrt_b
      Ndim_ov_aa = Nocc_a * Nvrt_a
      Ndim_ov_bb = Nocc_b * Nvrt_b

C <ab|ij> MO ints.

      IW0_aaaa = Ione 
      IW0_bbbb = IW0_aaaa + Ndim_oo_aa*Ndim_vv_aa
      IW0_abab = IW0_bbbb + Ndim_oo_bb*Ndim_vv_bb
      IW1_aaaa = IW0_abab + Ndim_oo_ab*Ndim_vv_ab

C <mn|ij> MO ints.

      Iw1_bbbb = IW1_aaaa + Ndim_oo_aa*Ndim_oo_aa
      Iw1_abab = Iw1_bbbb + Ndim_oo_bb*Ndim_oo_bb
      Iend     = Iw1_abab + Ndim_oo_ab*Ndim_oo_ab

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
     +     Work(IW1_abab),Work(Icmo_a),Work(Icmo_b),
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

      If (Iend .Ge. Maxcor) Call Insmem("cse_control",Iend,
     +                                   Maxcor)

      Call Prep_fock(Work(Ifock_a),Work(Ifock_b),Work(Ifockoo_a),
     +               Work(Ifockoo_b),Work(Ifockvv_a),Work(Ifockvv_b),
     +               Work(Ifockov_a),Work(Ifockov_b), 
     +               Work(Ifockoo_od_a),Work(Ifockoo_od_b),
     +               Work(Ifockvv_od_a),Work(Ifockvv_od_b),
     +               Nbasis,Naobfns,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b,
     +               Work(Iend),Memleft,Frac_occ)

      Ionehoo_a = Iend 
      Ionehoo_b = Ionehoo_a + Nocc_a*Nocc_a
      Ionehvv_a = Ionehoo_b + Nocc_b*Nocc_b
      Ionehvv_b = Ionehvv_a + Nvrt_a*Nvrt_a
      Ionehov_a = Ionehvv_b + Nvrt_b*Nvrt_b
      Ionehov_b = Ionehov_a + Nvrt_a*Nocc_a
      Iend      = Ionehov_b + Nvrt_b*Nocc_b

      Ionehoo_od_a = Iend
      Ionehoo_od_b = Ionehoo_od_a + Nocc_a*Nocc_a
      Ionehvv_od_a = Ionehoo_od_b + Nocc_b*Nocc_b
      Ionehvv_od_b = Ionehvv_od_a + Nvrt_a*Nvrt_a
      Iend         = Ionehvv_od_b + Nvrt_b*Nvrt_b
      
      Idensp_aa  = Iend
      Idensp_bb  = Idensp_aa + Nbasis*Nbasis 
      Idensm_aa  = Idensp_bb + Nbasis*Nbasis 
      Idensm_bb  = Idensm_aa + Nbasis*Nbasis 
      Iend       = Idensm_bb + Nbasis*Nbasis 

      Ionehao1  = Iend
      Ionehao2  = Ionehao1 + Naobfns*(Naobfns+1)/2
      Ionehmoa  = Ionehao2 + Naobfns*Naobfns 
      Ionehmob  = Ionehmoa + Naobfns*Naobfns
      Iend      = Ionehmob + Naobfns*Naobfns

      Memleft = Maxcor - Iend

      If (Iend .Ge. Maxcor) Call Insmem("cse_control",Iend,
     +                                   Maxcor)
      Fock = .False.
      Oneh = .True. 

      Call Prep_oneh(Work(Ionehao1),Work(Ionehao2),Work(Ionehmoa),
     +               Work(Ionehmob),Work(Ionehoo_a),Work(Ionehoo_b),
     +               Work(Ionehvv_a),Work(Ionehvv_b),Work(Ionehov_a),
     +               Work(Ionehov_b),Work(Ionehoo_od_a),
     +               Work(Ionehoo_od_b),Work(Ionehvv_od_a),
     +               Work(Ionehvv_od_b),Work(IW0_aaaa),Work(IW0_bbbb),
     +               Work(Iw0_abab),Nbasis,Naobfns,Nocc_a,Nocc_b,
     +               Nvrt_a,Nvrt_b,Work(Iend),Memleft,Frac_occ) 

      Call Diag_oneh(Work(Ionehmoa),Work(Ionehmob),Work(Idensp_aa),
     +               Work(Idensp_bb),Repuls,Nbasis,Naobfns,Nocc_a,
     +               Nocc_b,Nvrt_a,Nvrt_b,Work(Iend),Memleft)

      If (Fock) Then 
          Call Get_ref_energy(Work(IW0_aaaa),Work(IW0_bbbb),
     +                        Work(IW0_abab),Work(Ifock_a),
     +                        Work(Ifock_b),Nbasis,Nocc_a,
     +                        Nocc_b,Nvrt_a,Nvrt_b,Work(Iend),
     +                        Memleft,Repuls,Ref_energy)
      Elseif (Oneh) Then
 
          Call Get_ref_energy(Work(IW0_aaaa),Work(IW0_bbbb),
     +                        Work(IW0_abab),Work(Ionehoo_a),
     +                        Work(Ionehoo_b),Nbasis,Nocc_a,
     +                        Nocc_b,Nvrt_a,Nvrt_b,Work(Iend),
     +                        Memleft,Repuls,Ref_energy)
      Endif 

C T1(a,i) and T1(b,j)

      Iend      = Ionehao1 
      It1_aa    = Iend 
      It1_bb    = It1_aa + Nvrt_a*Nocc_a
      Iend      = It1_bb + Nvrt_b*Nocc_b

C F(m,i), F(a,e), F(m,e)  intermediate 

      Ifmi_a = Iend
      ifmi_b = Ifmi_a + Ndim_oo_aa
      Ifae_a = Ifmi_b + Ndim_oo_bb
      Ifae_b = Ifae_a + Ndim_vv_aa
      Ifme_a = Ifae_b + Ndim_vv_bb
      Ifme_b = Ifme_a + Ndim_ov_aa
      Idlt_a = Ifme_b + Ndim_ov_bb
      Idlt_b = Idlt_a + Ndim_ov_aa
      Ih1b_a = Idlt_b + Ndim_ov_bb
      Ih1b_b = Ih1b_a + Ndim_ov_aa*Ndim_ov_aa
      Iend   = Ih1b_b + Ndim_ov_aa*Ndim_ov_aa

      If (Iend .Ge. Maxcor) Call Insmem("ccd_light_control",Iend,
     +                                   Maxcor)
      Memleft = Maxcor - Iend 

      If (Fock) Then 
          Call Cse_driver(Work(It1_aa),Work(It1_bb),
     +         Work(IW0_aaaa),Work(IW0_bbbb),Work(IW0_abab),
     +         Work(Ifockoo_a),Work(Ifockoo_b),Work(Ifockvv_a),
     +         Work(Ifockvv_b),Work(Ifockov_a),Work(Ifockov_b),
     +         Work(Ifockoo_od_a),Work(Ifockoo_od_b),
     +         Work(Ifockvv_od_a),Work(Ifockvv_od_b),
     +         Work(Ifmi_a),Work(Ifmi_b),Work(Ifae_a),
     +         Work(Ifae_b),Work(Ifme_a),Work(Ifme_b),
     +         Work(Idlt_a),Work(Idlt_b),Work(Ih1b_a),Work(Ih1b_b),
     +         Work(Idens_aa),Work(Idens_bb),
     +         Work(Idensm_aa),Work(Idensm_bb),
     +         Nocc_a,Nocc_b,Nvrt_a,Nvrt_b,Nbasis,
     +         Work(Iend),Memleft,Repuls,Ref_energy,Frac_occ)
      Elseif (Oneh) Then

         Call Cse_driver(Work(It1_aa),Work(It1_bb),
     +         Work(IW0_aaaa),Work(IW0_bbbb),Work(IW0_abab),
     +         Work(Ionehoo_a),Work(Ionehoo_b),Work(Ionehvv_a),
     +         Work(Ionehvv_b),Work(Ionehov_a),Work(Ionehov_b),
     +         Work(Ionehoo_od_a),Work(Ionehoo_od_b),
     +         Work(Ionehvv_od_a),Work(Ionehvv_od_b),
     +         Work(Ifmi_a),Work(Ifmi_b),Work(Ifae_a),
     +         Work(Ifae_b),Work(Ifme_a),Work(Ifme_b),
     +         Work(Idlt_a),Work(Idlt_b),Work(Ih1b_a),Work(Ih1b_b),
     +         Work(Idensp_aa),Work(Idensp_bb),
     +         Work(Idensm_aa),Work(Idensm_bb),
     +         Nocc_a,Nocc_b,Nvrt_a,Nvrt_b,Nbasis,
     +         Work(Iend),Memleft,Repuls,Ref_energy,Frac_occ)
      Endif 

      Return
      End

