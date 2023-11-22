      Subroutine Pccd_form_uhtau(Work,Maxcor,Iuhf)

      Implicit Double Precision(A-H,O-Z)
      Logical Sym_packed 
      Logical Symmetry
      Logical OO_constr,Opt_orbs 
      Character*5 Spin(2)

      Dimension Work(Maxcor)

#include "syminf.com"
#include "sym.com"
#include "machsp.com"
#include "sympop.com"
#include "info.com"
      Common /Symm/Symmetry 
      Common /OO_info/Opt_orbs,OO_constr

      Data Spin /"Alpha", "Beta "/
      Data Ione,Onem,One,Dnull,Half,Two/1,-1.0D0,1.0D0,0.0D0,0.50D0,
     +                                  2.0D0/

      Irrepx = Ione

      Call Getrec(20,"JOBARC",'NDROPGEO',1,Idrop)
      If (Idrop .Ne. 0) Then
         Call Getrec(20,"JOBARC","NUMDROPA",1,Ndrop)
      Endif

C Form the orbital rotation gradients for UHF reference. This uses the 
C ACES II route (RHF route can be used for debugging purposes)

      If (Idrop .Ne. 0) Then 
          Write(6,*) "Entering drop block"
          Call Pccd_form_dropmo_urotgrads(Work,Maxcor,Iuhf)
      Else
          Call Pccd_form_urotgrads(Work,Maxcor,Iuhf)
      Endif 

      Nbas = Nocco(1) + Nvrto(1)

      Nocca   = Nocco(1)
      Nvrta   = Nvrto(1)
      Noccb   = Nocco(2)
      Nvrtb   = Nvrto(2)

      If (Symmetry) Then
         Ndim_ooa = Irpdpd(Irrepx,21)
         Ndim_vva = Irpdpd(Irrepx,19)
         Ndim_voa = Irpdpd(Irrepx,9)
         Ndim_oob = Irpdpd(Irrepx,22)
         Ndim_vvb = Irpdpd(Irrepx,19)
         Ndim_vob = Irpdpd(Irrepx,10)
      Else
         Ndim_ooa   = Nocca*Nocca
         Ndim_vva   = Nvrta*Nvrta
         Ndim_voa   = Nocca*Nvrta
         Ndim_oob   = Noccb*Noccb
         Ndim_vvb   = Nvrtb*Nvrtb
         Ndim_vob   = Noccb*Nvrtb
      Endif

      Ndim_oo = Ndim_ooa+Iuhf*Ndim_oob
      Ndim_vv = Ndim_vva+Iuhf*Ndim_vvb
      Ndim_vo = Ndim_voa+Iuhf*Ndim_vob
      Length  = Nbas*Nbas 

      I000 = Ione
      I010 = I000 + Ndim_oo
      I020 = I010 + Ndim_vv
      I030 = I020 + Ndim_vo
      I040 = I030 + Ndim_vo
      I050 = I040 + Length+Iuhf*Length
      Iend = I050 
      Maxcor = Maxcor - Iend 
      If (Iend.Gt.Maxcor) Call Insmem("Pccd_form_uhtau",Iend,Maxcor)

      Call Getrec(20,"JOBARC","ROTGRDOO",Ndim_oo*Iintfp,Work(I000))
      Call Getrec(20,"JOBARC","ROTGRDVV",Ndim_vv*Iintfp,Work(I010))
      Call Getrec(20,"JOBARC","ROTGRDVO",Ndim_vo*Iintfp,Work(I020))
      Call Getrec(20,"JOBARC","ROTGRDOV",Ndim_vo*Iintfp,Work(I030))
      Call Dscal(Ndim_vo,Onem,Work(I030),1)

#ifdef _DEBUG_LVL0
      write(6,*)
      call checksum("RTGRDOOA:",Work(I000),Ndim_ooa)
      call checksum("RTGRDOOB:",Work(I000+Ndim_ooa),Ndim_oob)
      call checksum("RTGRDVVA:",Work(I010),Ndim_vva)
      call checksum("RTGRDVVB:",Work(I010+Ndim_vva),Ndim_vvb)
      call checksum("RTGRDVOA:",Work(I020),Ndim_voa)
      call checksum("RTGRDVOB:",Work(I020+Ndim_voa),Ndim_vob)
      call checksum("RTGRDOVA:",Work(I030),Ndim_voa)
      call checksum("RTGRDOVB:",Work(I030+Ndim_voa),Ndim_vob)
#endif 
C
      If (Symmetry) Sym_packed = .True.
      Ioff    = I040
      Ioff_oo = I000
      Ioff_vv = I010
      Ioff_vo = I020
      Ioff_ov = I030

      Do Ispin = 1, Iuhf + 1 
         If (Ispin .Eq. 1) Then
            Nocc = Nocca
            Nvrt = Nvrta
         Else if (Ispin .Eq. 2) Then
            Nocc = Noccb
            Nvrt = Nvrtb
         Endif 
         Call Putrec(20,"JOABRC","SPNINDEX",Ione,Ispin)

         Call Pccd_frmful(Work(Ioff),Work(Ioff_oo),Work(Ioff_vv),
     +                    Work(Ioff_vo),Work(Ioff_ov),Work(Iend),
     +                    Maxcor,Nocc,Nvrt,Nbas,"Special",Sym_packed)
#ifdef _DEBUG_LVL0
         Write(6,*)
         Write(6,"(a,5a)") "The orbital rotation gradients for spin: ",
     &                      Spin(Ispin)
         Call output(Work(Ioff),1,Nbas,1,Nbas,Nbas,Nbas,1)
#endif 
         Call DScal(Nbas,Dnull,Work(Ioff),Nbas+1)
         Call Pccd_antisymmetrize(Work(Ioff),Work(Iend),Maxcor,Nbas,
     &                           Nocc,Nvrt)
         Call Pccd_htau_trnspose(Work(Ioff),Work(Iend),Nbas,Nocc,Nvrt)
 
        If (OO_constr) Call pccd_insist_no_oovv(Work(Ioff),Nocc,Nvrt,
     &                                          Nbas)

         If (Idrop .Ne. 0) Then
C Remove the rows and columns corresponding to the frozen orbitals. 

            Call Pccd_dropmo_fix(Work(Ioff),Work(Iend),Nbas,Ndrop)

            Write(6,"(2a,5a)") " The anti-symm. orbital rotation",
     +                         " gradient for spin: ",Spin(Ispin)
            Call output(Work(Ioff),1,Nbas,1,Nbas,Nbas,Nbas,1)
         Else
#ifdef _DEBUG_LVL0
            Write(6,"(2a,5a)") " The anti-symm. orbital rotation",
     +                         " gradient for spin: ",Spin(Ispin)
            Call output(Work(Ioff),1,Nbas,1,Nbas,Nbas,Nbas,1)
#endif 
         Endif 

         If (Ispin .Eq. 1) Call Putrec(20,"JOBARC","OBRTGRDA",
     +                                 Length*Iintfp,Work(Ioff))
         If (Ispin .Eq. 2) Call Putrec(20,"JOBARC","OBRTGRDB",
     +                                 Length*Iintfp,Work(Ioff))

         Ioff    = Ioff    + Nbas*Nbas
         Ioff_oo = Ioff_oo + Ndim_ooa
         Ioff_vv = Ioff_vv + Ndim_vva
         Ioff_vo = Ioff_vo + Ndim_voa
         Ioff_ov = Ioff_ov + Ndim_voa
      Enddo 

      Return
      End


