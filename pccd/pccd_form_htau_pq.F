      Subroutine Pccd_form_htau_pq(Work,Maxcor,Iuhf,Igrad_calc,
     +                             IHess_calc,Scale)
       
      Implicit Double Precision(A-H,O-Z)

      Logical Nonhf_ref,Symmetry
      Dimension Work(Maxcor)
      Dimension Nbfirr(8)

#include "info.com"
#include "machsp.com"
#include "sym.com"
#include "syminf.com"
#include "sympop.com"
#include "flags.com"

      Common/Nhfref/Nonhf_ref
      Common/Symm/Symmetry 

      Data Ione,Itwo,Dnull,Two /1,2,0.0D0,2.0D0/

      Nbas = Nocco(1) + Nvrto(1)

      Call Getrec(20,"JOBARC",'NDROPGEO',1,Idrop)

      If (Idrop .Ne. 0) Then
         Call Getrec(20,"JOBARC","NUMDROPA",1,Ndrop)
         Nbas = (Nbas - Ndrop) 
      Endif 

      Itrlen = Nbas*(Nbas+1)/Itwo 
      Isqlen = Nbas*Nbas 
      Nocc   = Nocco(1)
      Nvrt   = Nvrto(1)

      Noccsqf = Nocc*Nocc
      Nvrtsqf = Nvrt*Nvrt
      Nvrtocf = Nocc*Nvrt

      If (Symmetry) Then
         Irrepx = Ione
         Noccsqf = Irpdpd(Irrepx,21)
         Nvrtsqf = Irpdpd(Irrepx,19)
         Nvrtocf = Irpdpd(Irrepx,9)
      Endif 
   
      Ndim_oo = Nocc*Nocc
      Ndim_vv = Nvrt*Nvrt
      Ndim_vo = Nvrt*Nocc

      I000 = Ione
      I010 = I000 + Isqlen
      I020 = I010 + Isqlen
      I030 = I020 + Isqlen
      I040 = I030 + Isqlen
      I050 = I040 + Isqlen

      I060 = I050 + Ndim_oo
      I070 = I060 + Ndim_vv
      I080 = I070 + Isqlen
      I090 = I080 + Isqlen 
      Iend = I090 + Ndim_vo
      Maxcor = Maxcor - Iend 
      If (Iend .Gt. Maxcor) Call Insmem("Pccd_htau_pq",Iend,Maxcor)

#ifdef _NO_NRML_ORD
C This block for non-normal ordered optimizations
C Read triangularly packe MO basis one electron Hamiltonian and 
C expand it. This is H(p,q) and p is ordered as occ-vrt. 
   
      Call Getrec(20,"JOBARC","ONEHMOA",Itrlen*IIntfp,Work(I030))

C The last argument indicates whether the data is symmetry packed or
C simply the lower triangle.

      Call Pccd_mkfull(Work(I020),Work(I030),Work(Iend),Maxcor,Itrlen,
     +                 Nbas,.False.)
#ifdef _DEBUG_LVL0
      Write(6,"(a)") "The MO basis one electron Hamiltonian"
      call output(Work(I020),1,nbas,1,nbas,nbas,nbas,1)
#endif
#else
      Call Getrec(20,"JOBARC",'FOCKA   ',Isqlen*IIntfp,Work(I020))
      Call Ao2mo2(Work(I020),Work(I020),Work(I070),Work(I080),Nbas,
     +            Nbas,1)

#ifdef _DEBUG_LVL0
      Write(6,"(a)") "The MO basis FOCK matrix"
      call output(Work(I020),1,nbas,1,nbas,nbas,nbas,1)
#endif 
#endif 
C 
      Call Getrec(20,"JOBARC","DENSOO  ",Noccsqf*IIntfp,Work(I050))
      Call Getrec(20,"JOBARC","DENSVV  ",Nvrtsqf*IIntfp,Work(I060))
#ifdef _DEBUG_LVL0
      Write(6,"(a,1x,l)") "Is nonhf set?",nonhf_ref
#endif
      If (Nonhf_ref) Call Getrec(20,"JOBARC","DENSVO  ",
     +                           Nvrtocf*IIntfp,Work(I090))

      Call Pccd_form_dpq(Work(I030),Work(I040),Work(I050),Work(I060),
     +                   Work(I090),Work(Iend),Maxcor,Nocc,Nvrt,Nbas,
     +                   Iuhf,Nonhf_ref)

C Built H_tau(pq) = h(r,p)*D(q,r). Note that D(o,v)=D(v,o) = 0 (i.e.occ-vrt or 
C vrt-occ block of H_tau is zero). 

      Call Pccd_form_htau_1D_1(Work(I000),Work(I010),Work(Iend),
     +                         Maxcor,Work(I020),Work(I030),
     +                         Work(I040),Nocc,Nvrt,Nbas,Iuhf)
      If (Symmetry) Then
         Call Pccd_symexp(Work(I050),Work(I060),Work(I090),
     +                    Work(Iend),Maxcor,Nocc,Nvrt,Nbas)
      Endif 

      Call Pccd_form_htau_1D_2(Work(I000),Work(I010),Work(I050),
     +                         Work(I060),Work(I090),Work(I040),
     +                         Work(Iend),Maxcor,Nocc,Nvrt,Nbas,
     +                         Iuhf)
#ifdef _NO_NRML_ORD
      Call Dcopy(Nbas*Nbas,Work(I010),1,Work(I000),1)
#endif 

#ifdef _NOSKIP
      Call Pccd_form_htau_1D(Work(I000),Work(I010),Work(Iend),Maxcor,
     +                       Work(I020),Work(I030),Work(I040),
     +                       Work(I050),Work(I060),Work(I090),Nocc,
     +                       Nvrt,Nbas,Iuhf,E)
#endif

      I050 = I040 + Max(Noccsqf,Nvrtocf)
      I060 = I050 + Max(Noccsqf,Nvrtocf)
      I070 = I060 + Nvrtsqf
      I080 = I070 + Nvrtsqf
      I090 = I080 + Nvrtocf
      I100 = I090 + Nvrtocf
      I110 = I100 + Isqlen
      Iend = I110 + Isqlen
      Maxcor = Maxcor - Iend 
      If (Iend .Gt. Maxcor) Call Insmem("Pccd_htau_pq",Iend,Maxcor)

      Call Pccd_form_htau_oo(Work(I050),Work(I060),Work(I090),
     +                       Work(I040),Work(I050),Work(Iend),
     +                       Maxcor,Nocc,Nvrt,Noccsqf,
     +                       Nvrtsqf,Nvrtocf,Iuhf,Nonhf_ref)

      Call Pccd_form_htau_ov(Work(I050),Work(I060),Work(I090),
     +                       Work(I040),Work(I050),Work(Iend),
     +                       Maxcor,Nocc,Nvrt,Noccsqf,
     +                       Nvrtsqf,Nvrtocf,Iuhf,Nonhf_ref)

      Call Pccd_form_htau_2D(Work(I000),Work(I010),Work(I040),
     +                       Work(I050),Work(I060),Work(I070),
     +                       Work(I080),Work(I090),Work(I030),
     +                       Work(I100),Work(I110),Work(Iend),
     +                       Nocco(1),Nvrto(1),Maxcor,Nbas,E)

      Call Putrec(20,"JOBARC","OBRTGRDA",Isqlen*IIntfp,Work(I000))

#ifdef _DEBUG_LVL0
      Write(6,*)
      Call Pccd_sortgam(Work(Iend),Maxcor,Iuhf)
      Call Pccd_gfock_debug(Work(Iend),Maxcor,Iuhf)
#endif
      If (IHess_calc .EQ. 0) Return

C This is part of the Hessian that is computed from H*{p+q}. This
C is done here since all the pieces that is needed for this (Htau_pq 
C D_pq, f_pq) is available. 

      I100 = Iend
      I110 = I100 + Noccsq
      I120 = I110 + Noccsq
      I130 = I120 + Noccsq
      I140 = I130 + Nvrtsq
      I150 = I140 + Nvrtsq
      I160 = I150 + Nvrtsq
      I170 = I160 + Nvrtoc
      I180 = I170 + Nvrtoc
      I190 = I180 + Nvrtoc
      I200 = I190 + Nvrtoc
      I210 = I200 + Nvrtoc
      Iend = I210 + Nvrtoc
      Maxcor = Maxcor - Iend 
      If (Iend .Gt. Maxcor) Call Insmem("Pccd_htau_pq",Iend,Maxcor)

      Call Pccd_form_htau_1D_pqrs(Work(I000),Work(I010),
     +                            Work(I020),Work(I030),
     +                            Work(I040),Work(I050),
     +                            Work(I060),Work(I070),
     +                            Work(I080),Work(I090),
     +
     +                            Work(I100),Work(I010),
     +                            Work(I120),Work(I130),
     +                            Work(I140),Work(I150),
     +                            Work(I160),Work(I170),
     +                            Work(I180),Work(I190),
     +                            Work(I200),Work(I210),
     +                            Work(Iend),Nocco(1),
     +                            Nvrto(1),Maxcor,Nbas)

      Return 
      End 
    
     
