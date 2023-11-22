#include "flags.h"
      Subroutine Job_control(Iuhf)  
      Implicit Integer(A-Z)
c
      Logical Ca,Fc,Bca,Bfc,Occnums_file,Move
      Character*80 Fname 
      Dimension Idrop(1000,2)
      Dimension Ndrop(2)
      Logical Hf,Ks,Bk,No,Lo,Fo,pCCD,pCCDS,pCCSTD
      Logical pCCDQF,CCDQF
      Logical Xpccd 
      LOgical Act_amps_exsit  
      Logical Geom_opt,Analytical_gradient,Numerical_gradient
      Logical Vib_specs 
      Logical Density_mats,props
      Logical Ee,Ip,Ea
      Logical Oo,Drop
      Logical Drop_mo
      
#include "flags.com"
#include "flags2.com"

      Data Izz,Ione,Izero /100,1,0/

      Pccd_refr = Iflags2(h_IFLAGS2_pccd_reference)
      Pccd_calc = Iflags(h_IFLAGS_calc)
      Xpccd   = .False.
      Drop_mo = .False.
      If (Iflags2(h_IFLAGS2_program) .eq.4) Xpccd = .True.
      Call Aces_ja_init
      Call Getrec(20,"JOBARC",'NDROPGEO',1,Drop)
      Call Aces_ja_fin
      If (Drop .Gt. 0) Drop_mo = .True. 
      geom_opt=(iflags2(h_IFLAGS2_geom_opt).ne.0)
      analytical_gradient=(iflags2(h_IFLAGS2_grad_calc).eq.1)
      numerical_gradient =(iflags2(h_IFLAGS2_grad_calc).eq.2)
      vib_specs=(iflags(h_IFLAGS_vibration).ne.0)
      Density_mats=(iflags(h_IFLAGS_density).eq.1)
      props=(iflags(h_IFLAGS_props).ne.1)
      Ee=(iflags(h_IFLAGS_excite).gt.0)
      Ip=(iflags2(h_IFLAGS2_ip_calc).gt.0)
      Ea=(iflags2(h_IFLAGS2_ea_calc).gt.0)
      Oo=(iflags2(h_IFLAGS2_oo_method).gt.0)

      Hf = .False.
      Ks = .False.
      Bk = .False.
      No = .False.
      Lo = .False.
      Fo = .False.

      If (Pccd_refr .Eq. 0) Hf = .True.
      If (Pccd_refr .Eq. 1) Ks = .True.
      If (Pccd_refr .Eq. 2) Bk = .True.
      If (Pccd_refr .Eq. 3) No = .True.
      If (Pccd_refr .Eq. 4) Lo = .True.
      If (Pccd_refr .Eq. 5 .and. oo) Then
         Bk = .True.
      Elseif (Pccd_refr .Eq. 5) Then
         Write(6,*)
         Write(6,"(3a)")  " The choice PCCD_REFRENCE=OPTIMIZE need",
     +                    " also to select OO_METHOD (orbital ",
     +                    " optimization"
         Write(6,"(a)")   " algorithm)"
         Write(6,*)
         Call Errex
      Endif 

      pCCD    = .False.
      pCCDS   = .False.
      pCCDTS  = .False.
      pCCDTSD = .False.
      pCCDQF  = .False.
      CCDQF   = .False.

      If (Pccd_calc .Eq. 52)  pCCD    = .True.
      If (Pccd_calc .Eq. 55)  pCCDS   = .True.
      If (Pccd_calc .Eq. 56)  pCCDTS  = .True.
      If (Pccd_calc .Eq. 58)  pCCDTSD = .True.
      If (Pccd_calc .Eq. 57)  pCCDQF  = .True.
      If (Pccd_calc .Eq. 57)  CCDQF   = .True.

      If (geom_opt .or. vib_specs .or. Density_mats) Go to 999
      If (Ee .or. Ea .or. Ip) Then
         If (.Not. (pCCDTS .or. pCCDTSD)) Go to 9999
      Else
         Go to 9
      Endif

    9 Continue

      If (Bk) Then
         If (Drop_mo) Then
            Call Do_brueckner_drpmo(Xpccd)
         Else 
            Call Do_brueckner(Xpccd)
         Endif 
      Else if (No) Then
         Call runit("xvmol")
         Call runit("xvmol2ja")
         Call runit("xvscf")
         Call Do_natural_orbs()
      Else if (Lo) Then
         Call Do_local_orbs()
      Else if (Fo) Then
         Call do_oopccd(Iuhf)
      Else if (Hf .or. Ks) Then
         Call Runit("xvmol")
         Call Runit("xvmol2ja")
         Call Runit("xvscf")
         Call Runit("xvtran")
         Call Runit("xintprc")
         If (pCCDTS .or. pCCDTSD) Then

C The pCCTS (Tailored Pccd involve two steps. The first is true PCCD calculation.
C It produces ACT_AMPS_T file and then run a pCCDS calculation in the presence of
C ACT_AMPS_T file
            Call Aces_ja_init
            Iflags(h_IFLAGS_calc)=52
            Call putrec(20,'JOBARC','IFLAGS  ',100,iflags)
            Call Aces_ja_fin
            Call Runit("xvcct")
            Inquire(file='ACT_AMPS_T',exist=Act_amps_exist)
            If (Act_amps_exist) Then
               If (Iuhf .Eq. 0) Then
                  Call Runit("modify_act_rhf")
               Else
                  Call Runit("modify_act_uhf")
               Endif 
               Call Aces_ja_init
               If (pCCDTS)  iflags(h_IFLAGS_calc)=10
               If (pCCDTSD) iflags(h_IFLAGS_calc)=10
               Call Getrec(20,'JOBARC',"NUMDROPA",Ione,Ndrop(1))
               If (Ndrop(1) .Ne. Izero) Then
                 Call Getrec(20,'JOBARC',"MODROPAB",Ndrop(1),Idrop)
                 Call Putrec(20,'JOBARC',"MODROPA ",Ndrop(1),Idrop)
               Endif 
               Call putrec(20,'JOBARC','IFLAGS  ',Izz,iflags)
               Call Aces_ja_fin

C In strict sense, I do not need to rerun, vmol-xintprc because the 
C entire file set with the correct data structure are on the disk except
C for DIIS lists which are preset to dimensions without singles. In order 
C to make them flexible to handle methods that can have or have not singles
C requires ugly work (calling inipck2,newtyp2...). Let's skip these and 
C regenerate the file set (may be costly but simple, less error prone). 

               Call Runit("rm JOBARC")
               Call Runit("rm JAINDX")
               Call Runit("xjoda")
               Call Runit("xvmol")
               Call Runit("xvmol2ja")
               If (Ee .or. Ea .or. Ip) Call Runit("xvprops")
               Call Runit("xvscf")
               Call Runit("xvtran")
               Call Runit("xintprc")
               Call Runit("xvcct")
               If (Ee .or. Ea .or. Ip) Call Runit("xlambda")
               If (Ee) Call Runit("xvee")
               If (Ea) Call Runit("xvea")
               If (Ip) Call Runit("xvip")
            Else
               Write(6,"(2a)") " The Tailored pCCD calculations requir",
     +                         "es ACT_AMPS file from pCCD calculation"
               Call Runit("xvcc")
            Endif 
         Else
            If (pCCDQF .or. CCDQF) Then
              Call Runit("xccsdtqf")
            Else
              Call Runit("xvcc")
            Endif 

         Endif 
      Endif 
      Return

 999  Continue 

      If (.not. (hf .or .ks)) Then
         If (geom_opt) Then
         Write(6,"(2a)")" A pCCD geometry optimization can not procced",
     +                  " as a one step calculation for this choice of"
         Write(6,"(2a)")" reference state. Generate rhe reference state"
     +                  ," first then use scf_maxcyc=o use the dummy",
     +                  " pccd_reference=hf to continue."
         Call Errex
     
         Elseif (vib_specs) Then
         Write(6,"(2a)")" A pCCD vib. frequence calcs. can not procced",
     +                  " as a one step calculation for this choice of"
         Write(6,"(2a)")" reference state. Generate rhe reference state"
     +                  ," first then use scf_maxcyc=o use the dummy",
     +                   " pccd_reference=hf to continue."
         Call Errex
         Elseif (density_mats) Then
         Write(6,"(2a)")" A pCCD density matrix calcs. can not procced",
     +                  " as a one step calculation for this choice of"
         Write(6,"(2a)")" reference state. Generate rhe reference state"
     +                  ," first then use scf_maxcyc=o use the dummy",
     +                   " pccd_reference=hf to continue."
         Call Errex
         Endif 
      Endif 
      If (pCCDTS .or. pCCDQF) Then 
         Write(6,"(2a)")" There are  no analytical gradient for pCCDTS",
     +                  " or pCCDQF."
         Call Errex
      Endif

      If (geom_opt)  Call Run_opt()
      If (vib_specs) Call Run_freq()
      If (density_mats .or. props) Call Run_dens()
      Return

 9999 Continue 
     
      If (.not. (hf .or .ks)) Then
         Write(6,"(2a)")" A pCCD EE/IP/EA calcs. can not procced",
     +                  " as a one step calculation for this choice of"
         Write(6,"(2a)")" reference state. Generate rhe reference state"
     +                  ," first then use scf_maxcyc=o use the dummy",
     +                   " pccd_reference=hf to continue."
         Call Errex
      Endif 
      If (pCCDTS .or. pCCDQF) Then 
         Write(6,"(2a)")" There are no well defined E/IP/EA method",
     +                  " for pCCDTS or pCCDQF."
         Call Errex
      Endif

      Call Runit("xvmol")
      Call Runit("xvmol2ja")
      Call Runit("xvprops")
      Call Runit("xvscf")
      Call Runit("xvtran")
      Call Runit("xintprc")
      Call Runit("xvcc")
      Call Runit("xlambda")
      If (Ee) Call Runit("xvee")
      If (Ea) Call Runit("xvea")
      If (Ip) Call Runit("xvip")

      return
      End
