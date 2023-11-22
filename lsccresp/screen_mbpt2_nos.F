      Subroutine Screen_mbpt2_nos(Doo,Dvv,Coo,Cvv,Iwork,Imaxcor,IUhf,
     &                            Tau_o,Tau_v)

      Implicit Double Precision(A-H,O-Z)
      Double Precision Mone 

#include "sym.com"
#include "syminf.com"
#include "active_space.com"

      Dimension Doo(Nfmi(1)+Iuhf*Nfmi(2))
      Dimension Dvv(Nfea(1)+Iuhf*Nfea(2))
      Dimension Coo(Nfmi(1)+Iuhf*Nfmi(2))
      Dimension Cvv(Nfea(1)+Iuhf*Nfea(2))
      Dimension Iwork(Imaxcor)
      Dimension Nbf_irrep(8)
      Dimension Eps_o(3),EPs_v(3)
      character*80 wrk
      Logical Print_def 
      Logical Name_list_input

      Data Ione,Done,Mone,Dnull,Inull /1,1.0D0,-1.0D0,0.0D0,0/

C The active space thresholds are read from the Active_space
C namelist.

      Name_list_input = .False.

      If (Name_list_input) Then
 
          Print_def = .True. 
          Call Nl_int("ACTIVE_SPACE", Ierr, Print_def)
          Call Nl_int("OO_threshold", Dnull, OO_threshold)
          Call Nl_int("VV_threshold", Dnull, VV_threshold)

          If (Ierr .eq. Inull) Then
             Write(6,*)
             Write(6,"(a,a)") " The ACTIVE_SPACE namelist was",
     +                    " not found. The default values for the " 
             Write(6,"(a,F3.1,a)") " OO_threshold and VV_threshold ",
     +                           Dnull," is used."
             Write(6,*)
          Else 
             Write(6,*)
             Write(6,"(2(a,F15.6))") " The OO_threshold", OO_threshold,
     +                               " and the VV_threshold", 
     +                                                VV_threshold
             Write(6,"(a)") " are read from the ACTIVE_SPACE namelist"
          Endif 

      Else 
          OO_threshold =  Tau_o
CSSS          VV_threshold = -OO_Threshold
          VV_threshold = -Tau_v 
      Endif 

      Write(6,"(a)") 
      Write(6,"(a,F15.6)") " The OO_threshold: ", OO_threshold
      Write(6,"(a,F15.6)") " The VV_threshold: ", VV_threshold
      Write(6,"(a)") 

      Eps_o(1) = Dnull
      Eps_o(2) = Dnull
      Eps_o(3) = Dnull
      Eps_v(1) = Dnull
      Eps_v(2) = Dnull
      Eps_v(3) = Dnull

      Do Ispin = 1, (Iuhf+1)
         Ioff = (Ispin-1)*Nfmi(1) + Ione
         Joff = (Ispin-1)*Nfea(1) + Ione

         Do Irrep = 1, Nirrep
            Nd1 = pop(Irrep,Ispin)
            Nd2 = vrt(Irrep,Ispin)
   
            Write(6,"(3a,1x,I2,a,1x,I2,a)") "The occupied and virtual",
     +                                      " occupation vectors for",
     +                                      " irrep",Irrep," spin", 
     +                                      ispin,"."

            Write(6,"(6(1x,F12.6))") (Doo(Ioff+K-1),K=1,Nd1*Nd1,Nd1+1)
            Write(6,"(6(1x,F12.6))") (Dvv(Joff+K-1),K=1,Nd2*Nd2,Nd2+1)

#ifdef _DEBUG_LVLM
            Write(6,"(a,I2)") " The irrep:", Irrep
            If (Irrep .eq.1 .and. Ispin .eq. 1) Then
               Write(6,"(a)") " The alpha Doo and Coo"
            Elseif(Irrep .eq. 1 .and. Ispin .eq. 2) Then
               Write(6,"(a)") " The beta Doo and Coo"
            Endif 

            Call Output(Doo(Ioff),1,Nd1,1,Nd1,Nd1,Nd1,1) 
            Call Output(Coo(Ioff),1,Nd1,1,Nd1,Nd1,Nd1,1) 

            If (Irrep .eq.1 .and. Ispin .eq. 1) Then
               Write(6,*)
               Write(6,"(a)") " The alpha Dvv and Cvv"
            Elseif(Irrep .eq. 1 .and. Ispin .eq. 2) Then
               Write(6,"(a)") " The beta Dvv and Cvv"
            Endif 

            Call Output(Dvv(Joff),1,Nd2,1,Nd2,Nd2,Nd2,1)
            Call Output(Cvv(Joff),1,Nd2,1,Nd2,Nd2,Nd2,1)
#endif  
            Active_oo(Irrep,Ispin) = Inull
            Active_vv(Irrep,Ispin) = Inull
            Koff = Inull 
            Do I = 1, Nd1
               If (Doo(Ioff+Koff+I-1) .Lt. OO_threshold) Then
                  Active_oo(Irrep,Ispin) = Active_oo(Irrep,Ispin)+
     +                                     Ione  
                  Eps_o(1) = Eps_o(1) - Doo(Ioff+Koff+I-1) 
               Endif 
               Eps_o(2) = Eps_o(2) - Doo(Ioff+Koff+I-1) 
               Koff = Koff + Nd1
            Enddo 
            Koff = Inull 
            Do I = 1, Nd2
               If (Dvv(Joff+Koff+I-1) .Gt. VV_threshold) Then
                  Active_vv(Irrep,Ispin) = Active_vv(Irrep,Ispin)+
     +                                     Ione
                  Eps_v(1) = Eps_v(1) - Dvv(Joff+Koff+I-1) 
               Endif 
               Eps_v(2) = Eps_v(2) - Dvv(Joff+Koff+I-1) 
               Koff = Koff + Nd2
            Enddo 
            
            Ioff = Ioff + Nd1*Nd1
            Joff = Joff + Nd2*Nd2
         Enddo 
      Enddo

      Eps_o(1) = Eps_o(1)/Eps_o(2)
      Eps_o(3) = Done - Eps_o(1)
      Eps_v(1) = Eps_v(1)/Eps_v(2)
      Eps_v(3) = Done - Eps_v(1)

#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(2a,3(1x,F12.6))")" The occupied orbital energy",
     +                           " interval: ", (Eps_o(i),i=1,3)
      Write(6,"(2a,3(1x,F12.6))")" The virtual orbital energy", 
     +                          " interval : ", (Eps_v(i),i=1,3)
#endif 

      Ioff_active_oo(1,1) = Ione
      Ioff_active_vv(1,1) = Ione
      Noo_active_aa       = Inull
      Nvv_active_aa       = Inull

      Do irrep = 1, Nirrep-1
         Ioff_active_oo(Irrep+1,1) = Ioff_active_oo(Irrep,1) +
     +                               Active_oo(Irrep,1)
         Ioff_active_vv(Irrep+1,1) = Ioff_active_vv(Irrep,1) +
     +                               Active_vv(Irrep,1) 
         Noo_active_aa = Noo_active_aa + Active_oo(Irrep,1)
         Nvv_active_aa = Nvv_active_aa + Active_vv(Irrep,1)
      Enddo

      Ioff_active_oo(1,2) = Noo_active_aa
      Ioff_active_vv(1,2) = Nvv_active_aa
      Do irrep = 1, Nirrep-1
         Ioff_active_oo(Irrep+1,2) = Ioff_active_oo(Irrep,2) +
     +                               Active_oo(Irrep,2)
         Ioff_active_vv(Irrep+1,2) = Ioff_active_vv(Irrep,2) +
     +                               Active_vv(Irrep,2) 
      Enddo

      Call Getrec(20,"JOBARC","NUMBASIR",Nirrep,Nbf_irrep)

      Icount = Ione
      Do Ispin = 1, (Iuhf+1)
         Do Irrep = 1, Nirrep
            Do I = 1, Pop(Irrep,Ispin)
               If (I .le. Active_oo(Irrep,Ispin)) Then
                  Iwork(Icount) = Ione
                  Icount = Icount + 1
               Else
                  Iwork(Icount) = Inull
                  Icount = Icount + 1
               Endif 
            Enddo 
         Enddo
      Enddo 

      Do Ispin = 1, (Iuhf+1)
         Do Irrep = 1, Nirrep
            Do I = 1, Vrt(Irrep,Ispin)
               If (I .le. Active_vv(Irrep,Ispin)) Then
                  Iwork(Icount) = Ione
                  Icount = Icount + 1
               Else
                  Iwork(Icount) = Inull
                  Icount = Icount + 1
               Endif 
            Enddo 
         Enddo
      Enddo 
      If (Icount .Gt. Imaxcor) Call Insmem("screen_mbpt2_nos",Icount,
     +                                      Imaxcor)

#ifdef _DEBUG
C set the all active or inactive limits by force to check the limits
      Do Ispin = 1, Iuhf+1
      Do Irrep = 1, Nirrep
         Nd1 = pop(Irrep,Ispin)
         Nd2 = vrt(Irrep,Ispin)
         Active_oo(Irrep,Ispin) = Inull
         Active_vv(Irrep,Ispin) = Inull
         Do I = 1, Nd1
            Active_oo(Irrep,Ispin) = Active_oo(Irrep,
     +                                         Ispin) + Ione
         Enddo
         Do I = 1, Nd2
            Active_vv(Irrep,Ispin) = Active_vv(Irrep,
     +                                         Ispin) + Ione
        Enddo
      Enddo
      Enddo
#endif

      If (Iuhf .Eq. 0) Then
         Call Icopy(Nirrep,Active_oo(1,1),1,Active_oo(1,2),1)
         Call Icopy(Nirrep,Active_vv(1,1),1,Active_vv(1,2),1)
      Endif 

#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a)") " No. Alpha OO and VV active:"
      Write(6,"(8(1x,i3))") (Active_oo(i,1),i=1,Nirrep)
      Write(6,"(8(1x,i3))") (Active_vv(i,1),i=1,Nirrep)
      Write(6,*)
      Write(6,"(a)") " No. Beta OO and VV active:"
      Write(6,"(8(1x,i3))") (Active_oo(i,2),i=1,Nirrep)
      Write(6,"(8(1x,i3))") (Active_vv(i,2),i=1,Nirrep)
#endif

      Return
      End


