#include "flags.h"
      Subroutine Lsccresp_driver(Work,IWork,Maxcor,Imaxcor,Iuhf)

      Implicit Double Precision (A-H,O-Z)

#include "machsp.com"
#include "syminf.com"
#include "sym.com"
#include "flags.com"
#include "flags2.com"
#include "sympop.com"
#include "info.com"
#include "symoff.com"
#include "active_space.com"

      Dimension Work(Maxcor)
      Dimension Iwork(Imaxcor)

      Integer T2ln_aa,T2ln_bb,T2ln_ab,T2ln
      Integer T1ln_aa,T1ln_bb
      Integer T1off,T1ln,Tln 
      Logical CCt1,Cct2,Singles,Tau,Nonhf
      Logical Lambda
      Logical T1_2ndord 
      Character*3 Type
      Dimension E(3),ES(3)

      Data Ione,Inull,I30,Dnull,Two,Twl /1,0,30,0.0D0,2.0D0,12.0D0/

      T2ln_aa  = 0
      T2ln_bb  = 0
      T2ln_ab  = 0
      T1ln_aa  = 0
      T1ln_bb  = 0 
      Lambda   = (Iflags(h_IFLAGS_props) .Ne. Inull)

C Precompute the symmetry packed lengths of <AB|IJ> like distributions

      Do Ispin = 1, (Iuhf+1)
         If (Ispin .Eq. 1) Then
            Do irrep = 1, Nirrep
               Ncol = Irpdpd(Irrep,Ispin)
               Nrow = Irpdpd(Irrep,2+Ispin)
               T2ln_aa  = T2ln_aa + Nrow*Ncol
               T1ln_aa  = T1ln_aa + Irpdpd(Irrep,9)
            Enddo
         Elseif (Ispin .Eq. 2) Then
            Do irrep = 1, Nirrep
               Ncol = Irpdpd(Irrep,Ispin)
               Nrow = Irpdpd(Irrep,2+Ispin)
               T2ln_bb  = T2ln_bb + Nrow*Ncol
               T1ln_bb  = T1ln_bb + Irpdpd(Irrep,10)
            Enddo
         Endif
      Enddo

      Do irrep = 1, Nirrep
         Ncol = Irpdpd(Irrep,14)
         Nrow = Irpdpd(Irrep,13)
         T2ln_ab = T2ln_ab + Nrow*Ncol
      Enddo

C Precompute the irrep offsets occ-occ and vrt-vrt blocks

      Ioff_oo(1,1) = Ione
      Ioff_oo(1,2) = Nfmi(1) + Ione 
      Ioff_vv(1,1) = Ione
      Ioff_vv(1,2) = Nfea(1) + Ione
      Do Ispin = 1, 2 
         Do irrep = 1, Nirrep-1
            Ioff_oo(Irrep+1,Ispin) = Ioff_oo(Irrep,Ispin) +
     +                               Pop(Irrep,Ispin)*Pop(Irrep,Ispin)
            Ioff_vv(Irrep+1,Ispin) = Ioff_vv(Irrep,Ispin) +
     +                               Vrt(Irrep,Ispin)*Vrt(Irrep,Ispin)
        Enddo
      Enddo

      T1ln = T1ln_aa + Iuhf*T1ln_bb
      T2ln = T2ln_aa + Iuhf*T2ln_bb + T2ln_ab
      Tln  = T2ln + T1ln

      I000 = 1
      Iend = I000 + T2ln + T1ln 
      if (Iend .Gt. Maxcor) Call Insmem("lsccresp_driver",Iend,
     +                                   Maxcor)

C Second order T2 amplittudes 

      Maxcor = Maxcor - Iend
      Call Get_mbpt2_t2s(Work(I000),Work(Iend),Maxcor,T2ln,T2ln_aa,
     +                    T2ln_bb,T2ln_ab,Iuhf)

C Create second-order T1 amplitudes.

      T1_2ndord = .True.
      T1off     = I000 + T2ln
      If (T1_2ndord) Then
          Call Get_2ndord_t1s(Work(I000),Work(Iend),Maxcor,Tln,
     +                        T2ln,T1ln,T2ln_aa,T2ln_bb,T2ln_ab,
     +                        T1ln_aa,T1ln_bb,Iuhf)
      Else
          Call Dzero(Work(T1off),T1ln)
      Endif 

      I010 = Iend 
      I020 = I010 + Nfmi(1)+Iuhf*Nfmi(2)
      I030 = I020 + Nfea(1)+Iuhf*Nfea(2)
      Iend = I030 + Nt(1)+Iuhf*Nt(2)
      if (Iend .Gt. Maxcor) Call Insmem("lsccresp_driver",Iend,
     +                                   Maxcor)
      Maxcor = Maxcor - Iend

C 2nd order density matrices. There is no ov block for 2nd order. 
C Lets leave the room for ov blocks if  we  choose to additional work.

      Call Get_mbpt2_dens(Work(I000),Work(I010),Work(I020),Work(I030),
     +                    Work(Iend),Maxcor,T2ln,T2ln_aa,T2ln_bb,
     +                    T2ln_ab,Iuhf)
#ifdef _DEBUG_LVL0
      Write(6,*)
      Singles = .False.
      Cct2    = .False.
      Nonhf   =  (Iflags(h_IFLAGS_nonhf) .Eq. 1)
      Call Get_energy(Work(I000),Work(Iend),Maxcor,90,43,0,T2ln,
     +                T2ln_aa,T2ln_bb,Singles,Cct2,Nonhf,Iuhf,E,ES)
      Write(6,"(a,a,ES17.9E2)") " The total MBPT(2) correlation",
     +                          " energy      :",E(1)+E(2)+E(3)
#endif 
      I040 = Iend
      I050 = I040 + Nfmi(1)+Iuhf*Nfmi(2)
      Iend = I050 + Nfea(1)+Iuhf*Nfea(2)

      if (Iend .Gt. Maxcor) Call Insmem("lsccresp_driver",Iend,
     +                                   Maxcor)
      Maxcor = Maxcor - Iend

C 2nd order natural orbitals

      Call Get_mbpt2_nos(Work(I010),Work(I020),Work(I040),Work(I050),
     +                   Work(Iend),Maxcor,Iuhf)

C Setup a loop over the number of extrapolation points. Currently the
C number of extrapolation point set to 30

      N_xp            = I30
      Eta_val(1)      = Dnull 
      Eta_val(2)      = Two
      Eta_val(n_xp+2) = Twl

      Do I = 3, N_xp+1
         Eta_val(I) = Eta_val(I-1) + 0.125D0
      Enddo 

      Call Getrec(0,"JOBARC","ETA_4EXP",Irecord,Ijunk)
      If (Irecord .Gt. 0) Then 
          Write(6,"(2a)") " The Eta values are read from the JOBARC",
     +                    " file."
          N_xp = 0
      Endif 
#ifdef _DEBUG
          N_xp = 1
#endif 
      Do I_xp = 1, n_xp + 2
         If (Irecord .Gt. 0) Then
            Call Getrec(20,"JOBARC","ETA_4EXP",Ione,Eta_o)
#ifdef _DEBUG
            Eta_o= -0.0000001
            Eta_v = -0.01
#endif 
         Else
            Eta_o = - 10.0d0**(-Eta_val(I_xp))
         Endif 
         Eta_v = Eta_o 

         Write(6,"(a,I3,a,2ES17.9E2)") " The iteration ",I_xp,
     +                                 " Eta_o and Eta_v ",Eta_o,
     +                                   Eta_v
#ifdef _DEBUG
         If (I_xp .eq. 2) Then
           Eta_o= -0.0000001
           Eta_v = -0.01
         Endif 
#endif 
         Call Screen_mbpt2_nos(Work(I010),Work(I020),Work(I040),
     +                         Work(I050),Iwork,Imaxcor,Iuhf,
     +                         Eta_o,Eta_v)

C Transform 2nd order T2, CCSD, T1, T2 and L1 and L2 to natural orbital
C basis. First T1 and T2. The type C2N convey that this is a canonical to
C NO transformation. The list2_h (host list) has converged T2 and after 
C the transformation it is written to the target list2_t. Then the target
C list get overwritten by the CC and MBPT(2) amplitude mixture. The back
C transformation also happens with the target list. We can use lists 64,
C 65 and 66 for this purpose,

         Cct2  = .True.
         Cct1  = .True.

         Type  = "C2N"
         List1 = 90
         Call Transform_2index(Work(I000),Work(I040),Work(I050),
     +                         Work(Iend),Maxcor,List1,Iuhf,T1ln,
     +                         T1ln_aa,Cct1,Type)

         List2_h = 43
         List2_t = 63
         Call Transform_4index_aa(Work(I000),Work(I040),Work(I050),
     +                            Work(Iend),Maxcor,List2_h,List2_t,
     +                            T2ln,T2ln_aa,Iuhf,Cct2,Type)

         List2_h = 46
         List2_t = 66
         Call Transform_4index_ab(Work(I000),Work(I040),Work(I050),
     +                            Work(Iend),Maxcor,List2_h,List2_t,
     +                            T2ln,T2ln_aa,T2ln_bb,Iuhf,Cct2,Type)

C#ifdef _DEBUG
         List2 = 13
         Call Transform_4index_aa(Work(I000),Work(I040),Work(I050),
     +                            Work(Iend),Maxcor,List2,List2,T2ln,
     +                            T2ln_aa,Iuhf,Cct2,Type)

         List2 = 16
         Call Transform_4index_ab(Work(I000),Work(I040),Work(I050),
     +                            Work(Iend),Maxcor,List2,List2,T2ln,
     +                            T2ln_aa,T2ln_bb,Iuhf,Cct2,Type)
C#endif 
C Now L1 and L2.
         If (Lambda) Then
             List1 = 190
             Call Transform_2index(Work(I000),Work(I040),Work(I050),
     +                             Work(Iend),Maxcor,List1,Iuhf,T1ln,
     +                             T1ln_aa,Cct1,Type)
             List2 = 143
             Call Transform_4index_aa(Work(I000),Work(I040),Work(I050),
     +                                Work(Iend),Maxcor,List2,List2,
     +                                T2ln,T2ln_aa,Iuhf,Cct2,Type)
             List2 = 146
             Call Transform_4index_ab(Work(I000),Work(I040),Work(I050),
     +                                Work(Iend),Maxcor,List2,List2,
     +                                T2ln,T2ln_aa,T2ln_bb,Iuhf,Cct2,
     +                                Type)

         Endif

C MBPT(2) T2.

         Cct1  = .False.
         Cct2  = .False.

         Call Transform_2index(Work(T1off),Work(I040),Work(I050),
     +                         Work(Iend),Maxcor,List1,Iuhf,T1ln,
     +                         T1ln_aa,Cct1,Type)

         Call Transform_4index_aa(Work(I000),Work(I040),Work(I050),
     +                            Work(Iend),Maxcor,List2,List2,T2ln,
     +                            T2ln_aa,Iuhf,Cct2,Type)

         Call Transform_4index_ab(Work(I000),Work(I040),Work(I050),
     +                            Work(Iend),Maxcor,List2,List2,T2ln,
     +                            T2ln_aa,T2ln_bb,Iuhf,Cct2,Type)

         Call form_new_amps(Work(I000),Work(Iend),Iwork,Maxcor,
     +                      Imaxcor,Tln,T1ln,T2ln,T2ln_aa,T2ln_bb,
     +                      T2ln_ab,T1ln_aa,T1ln_bb,Iuhf)

C Transform all CCSD T1, T2, L1 and L2 and MBPT(2) T back to Canonical basis
C First CCSD T1 and T2.

         Type  = "N2C"
         Cct1  = .True.
         Cct2  = .True.
         List1 = 90
         Call Transform_2index(Work(I000),Work(I040),Work(I050),
     +                         Work(Iend),Maxcor,List1,Iuhf,T1ln,
     +                         T1ln_aa,Cct1,Type)

         List2 = 63
         Call Transform_4index_aa(Work(I000),Work(I040),Work(I050),
     +                            Work(Iend),Maxcor,List2,List2,T2ln,
     +                            T2ln_aa,Iuhf,Cct2,Type)

         List2 = 66
         Call Transform_4index_ab(Work(I000),Work(I040),Work(I050),
     +                            Work(Iend),Maxcor,List2,List2,T2ln,
     +                            T2ln_aa,T2ln_bb,Iuhf,Cct2,Type)

C Now L1 and L2.

         If (Lambda) Then
             List1 = 190
             Call Transform_2index(Work(I000),Work(I040),Work(I050),
     +                             Work(Iend),Maxcor,List1,Iuhf,
     +                             T1ln,T1ln_aa,CCt1,Type)
             List2 = 143
             Call Transform_4index_aa(Work(I000),Work(I040),Work(I050),
     +                                Work(Iend),Maxcor,List2,List2,
     +                                T2ln,T2ln_aa,Iuhf,Cct2,Type)
             List2 = 146
             Call Transform_4index_ab(Work(I000),Work(I040),Work(I050),
     +                                Work(Iend),Maxcor,List2,List2,
     +                                T2ln,T2ln_aa,T2ln_bb,Iuhf,Cct2,
     +                                Type)
         Endif 

C MBPT(2) T2.

         Cct1  = .False.
         Cct2  = .False.
         Call Transform_2index(Work(T1off),Work(I040),Work(I050),
     +                         Work(Iend),Maxcor,List1,Iuhf,T1ln,
     +                         T1ln_aa,Cct1,Type)

         Call Transform_4index_aa(Work(I000),Work(I040),Work(I050),
     +                            Work(Iend),Maxcor,List2,List2,T2ln,
     +                            T2ln_aa,Iuhf,Cct2,Type)

         Call Transform_4index_ab(Work(I000),Work(I040),Work(I050),
     +                            Work(Iend),Maxcor,List2,List2,T2ln,
     +                            T2ln_aa,T2ln_bb,Iuhf,Cct2,Type)

C#ifdef _DEBUG
         CCt2  = .True.
         List2 = 13
         Call Transform_4index_aa(Work(I000),Work(I040),Work(I050),
     +                            Work(Iend),Maxcor,List2,List2,T2ln,
     +                            T2ln_aa,Iuhf,Cct2,Type)

         List2 = 16
         Call Transform_4index_ab(Work(I000),Work(I040),Work(I050),
     +                            Work(Iend),Maxcor,List2,List2,T2ln,
     +                            T2ln_aa,T2ln_bb,Iuhf,Cct2,Type)

C#endif 
         If (I_xp .Eq. 1) Then
            Singles = .False.
         Else 
            Singles = .True. 
         Endif 
         Nonhf   =  (Iflags(h_IFLAGS_nonhf) .Eq. 1)
         Call Get_energy(Work(I000),Work(Iend),Maxcor,90,63,2,T2ln,
     +                   T2ln_aa,T2ln_bb,Singles,Cct2,Nonhf,Iuhf,E,ES)
         Write(6,*)
         Write(6,"(a,a,ES17.9E2)")" The total CCSD correlation energy ",
     +                            "        :",E(1)+E(2)+E(3)
         E_k(I_xp)  = E(1) + E(2) + E(3)
         E_kS(I_xp) = ES(1) + ES(2) + ES(3) - E_k(I_xp)

         Call Rebuilt_intms(Work(Iend),Maxcor,Iuhf) 

      Enddo 

      Call Getrec(20,"JOBARC","SCFENEG ",Ione,E_scf)

#ifdef _DEBUG_LVLM
      Do I_xp = 1, N_xp+2
      Write(6,"(a,ES17.9E2)") " The singles contribution: ",E_kS(I_xp)
      Enddo 
#endif 

      Write(6,"(3x,2a)")   "----------------------------------------",
     +                 "----------------------------"
      Write(6,"(4x,a,5x,a,17x,a,2x,a)") " Eta_k ", "Tau", 
     +                                  "Correlation energy", 
     +                                  "Total Energy"
      Write(6,"(3x,2a)")   "----------------------------------------",
     +                 "----------------------------"
      Do I_xp = 1, N_xp+2
         Write(6,"(1x,F10.5,3x,ES17.9E2,3x,ES17.9E2,3x,ES17.9E2)") 
     +                                                 Eta_val(I_xp),
     +                                                 10.0d0**
     +                                                 (-Eta_val(I_xp)),
     +                                                 E_k(I_xp),
     +                                                 E_k(I_xp)+E_scf
      Enddo 
      Write(6,"(3x,2a)")   "----------------------------------------",
     +                 "----------------------------"

      Return 
      End
