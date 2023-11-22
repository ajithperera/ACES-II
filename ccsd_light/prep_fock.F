      Subroutine Prep_fock(Fock_a,Fock_b,Fockoo_a,Fockoo_b,
     +                     Fockvv_a,Fockvv_b,Fockov_a,Fockov_b,
     +                     Fockoo_od_a,Fockoo_od_b,
     +                     Fockvv_od_a,Fockvv_Od_b,
     +                     W0_aa,W0_bb,W0_ab,E1_aa,E1_bb,
     +                     Nbasis,Naobfns,Nocc_a,Nocc_b,
     +                     Nvrt_a,Nvrt_b,Work,Maxcor,Frac_occ,
     +                     Open_shell)

      Implicit Double Precision(A-H,O-Z)
      Dimension Work(Maxcor)
      Logical Frac_occ,Do_Ao_fock,Do_mo_fock
      Logical Open_shell
      Integer Dbuf 
      
      Dimension Fock_a(Nbasis,Nbasis),Fock_b(Nbasis,Nbasis)

      Dimension Fockoo_a(Nocc_a,Nocc_a),Fockoo_b(Nocc_b,Nocc_b)
      Dimension Fockvv_a(Nvrt_a,Nvrt_a),Fockvv_b(Nvrt_b,Nvrt_b)
      Dimension Fockoo_od_a(Nocc_a,Nocc_a),Fockoo_od_b(Nocc_b,Nocc_b)
      Dimension Fockvv_od_a(Nvrt_a,Nvrt_a),Fockvv_Od_b(Nvrt_b,Nvrt_b)
      Dimension Fockov_a(Nocc_a,Nvrt_a),Fockov_b(Nocc_b,Nvrt_b)

      Dimension W0_aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension w0_bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension W0_ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)
      Dimension Nbfns_irrep(8) 

#include "sym.com"
#include "syminf.com"
#include "info.com"
#include "ccsdlight_vars.com"

      Data Ione,Done /1,1.0D0/

      Ldim    = Nbasis*(Nbasis+1)/2
      Nbasis2 = Nbasis*Nbasis 
      Lbuf    = 600
      Iunit   = Ione

      Do_ao_fock = .True.
      Do_mo_fock = .False.

      If (Frac_occ) Then 
         
         Ioneh_aa = 1
         Ioneh_bb = Ioneh_aa + Nbasis2
         Ioneh1   = Ioneh_bb + Nbasis2
         Ioneh2   = Ioneh1   + Ldim
         Iscr1    = Ioneh2   + Nbasis2
         Iscr2    = Iscr1    + Nbasis2
         Iend     = Iscr2    + Nbasis2
         If (Iend .GT. Maxcor) Call Insmem("Get_ref_energy",Iend,Maxcor)

         Call Get_oneh(Work(Ioneh1),Work(Ioneh2),Work(Ioneh_aa),
     +                 Work(Ioneh_bb),Work(Iscr1),Work(Iscr2),
     +                 E1_aa,E1_bb,Nbasis)

         If (Do_ao_fock) Then 
C Built the AO basis fock matrix (from the AO integrals and the
C density matirices.

         Idens_aa = Ioneh1 
         Idens_bb = Idens_aa + Nbasis2
         Idens_tt = Idens_bb + Nbasis2
         Ievec_sa = Idens_tt + Nbasis2
         Iscr1    = Ievec_sa + Nbasis2
         Iscr2    = Iscr1    + Nbasis2
         Iend     = Iscr2    + Nbasis2
         If (Iend .GT. Maxcor) Call Insmem("Get_ref_energy",Iend,Maxcor)

         Call Built_frac_dens(Work(Idens_aa),Work(Idens_bb),
     +                        Work(Idens_tt),Work(Ievec_sa),
     +                        Work(Iscr1),Work(Iscr2),Nbasis)
         Ifock_a = Ievec_sa
         Ifock_b = Ifock_a + Ldim
         Dbuf    = Ifock_b + Ldim
         Ibuf    = Dbuf    + Lbuf 
         Imap    = Ibuf    + Lbuf
         Iend    = Imap    + Naobfns*Naobfns 
         If (Iend .GT. Maxcor) Call Insmem("Get_ref_energy",Iend,Maxcor)

         Call Squeez(Work(Idens_aa),Nbasis,0)
         Call Squeez(Work(Idens_bb),Nbasis,0)
         Call Squeez(Work(Idens_tt),Nbasis,0)

         Call Getrec("JOBARC",20,"NUMBASIR",Nirrep,Nbfns_Irrep) 
         Ldim_c = 0
         Do I = 1, Nirrep
            Ldim_c = Ldim_c + Nbfns_irrep(I)*(Nbfns_irrep(I)+1)/2
         Enddo 

         If (Ldim_c .Ne. Ldim) Then
            Write(6,"(3a,2(1x,I4))") "There is a inconsistency in",
     +                      " triangular dimensions. There are two", 
     +                      " values", Ldim,Ldim_c 
            Call Errex
         Endif 

         Call Built_frac_fock_ao(Work(Ifock_a),Work(Ifock_b),
     +                           Work(Idens_aa),Work(Idens_bb),
     +                           Work(Idens_tt),Work(Dbuf),
     +                           Work(Ibuf),Work(Imap),
     +                           Naobfns,Nbfns_Irrep,Ldim_c,
     +                           Lbuf,Iunit)

         Call Expnd2(Work(Ifock_a),Fock_a,Nbasis)
         Call Expnd2(Work(Ifock_b),Fock_b,Nbasis)

         Iscr1 = Dbuf
         Iscr2 = Iscr1 + Nbasis2 
         Iscr3 = Iscr2 + Nbasis2 
         Iscr4 = Iscr3 + Nbasis
         Iend  = Iscr4 + Nbasis
         If (Iend .GT. Maxcor) Call Insmem("prep_fock",Iend,Maxcor)

         Call Form_mo_Fock(Fock_a,Fock_a,Work(Iscr1),Work(Iscr2),
     +                     Nbasis,Naobfns,1) 
         Call Form_mo_Fock(Fock_b,Fock_b,Work(Iscr1),Work(Iscr2),
     +                     Nbasis,Naobfns,2)

         Call Daxpy(Nbasis2,Done,Work(Ioneh_aa),1,Fock_a,1)
         Call Daxpy(Nbasis2,Done,Work(Ioneh_bb),1,Fock_b,1)

#ifdef _DEBUG_LVL0
      Write(6,"(a)") "@-prep_fock: MO basis fock matrices" 
      call output(Fock_a,1,Nbasis,1,Nbasis,Nbasis,Nbasis,1)
      call output(Fock_b,1,Nbasis,1,Nbasis,Nbasis,Nbasis,1)
#endif 
         Elseif (Do_mo_fock) Then

C Built the MO basis fock matrix from the transformed integrals (MO)
C and MO density matrices. 

         Call Built_frac_fock_mo(W0_aa,w0_bb,W0_ab,Work(Ioneh_aa),
     +                           Work(Ioneh_bb),Nocc_a,Nocc_b,
     +                           Nvrt_a,Nvrt_b,Nbasis)

         Call Dcopy(Nbasis2,Work(Ioneh_aa),1,Fock_a,1)
         Call Dcopy(Nbasis2,Work(Ioneh_bb),1,Fock_b,1)

#ifdef _DEBUG_LVLM
      Write(6,"(a)") "@-prep_fock: MO basis fock matrices" 
      call output(Fock_a,1,Nbasis,1,Nbasis,Nbasis,Nbasis,1)
      call output(Fock_b,1,Nbasis,1,Nbasis,Nbasis,Nbasis,1)
#endif 
         Endif 
      Else 
C Read the fock matrices built during SCF (in AO basis)
         Call Getrec(20,"JOBARC","FOCKA",Nbasis2,Fock_a)
         Call Getrec(20,"JOBARC","FOCKB",Nbasis2,Fock_b)

#ifdef _DEBUG_LVLM
      Write(6,"(a)") "@-prep_fock: AO basis fock matrices" 
      call output(Fock_a,1,Nbasis,1,Nbasis,Nbasis,Nbasis,1)
      call output(Fock_b,1,Nbasis,1,Nbasis,Nbasis,Nbasis,1)
#endif 
         Iscr1 = 1
         Iscr2 = Iscr1 + Nbasis2 
         Iscr3 = Iscr2 + Nbasis2 
         Iscr4 = Iscr3 + Nbasis
         Iend  = Iscr4 + Nbasis
         If (Iend .GT. Maxcor) Call Insmem("prep_fock",Iend,Maxcor)

         Call Form_mo_Fock(Fock_a,Fock_a,Work(Iscr1),Work(Iscr2),
     +                     Nbasis,Naobfns,1) 
         Call Form_mo_Fock(Fock_b,Fock_b,Work(Iscr1),Work(Iscr2),
     +                     Nbasis,Naobfns,2)
#ifdef _DEBUG_LVLM
      Write(6,"(a)") "@-prep_fock: MO basis fock matrices" 
      call output(Fock_a,1,Nbasis,1,Nbasis,Nbasis,Nbasis,1)
      call output(Fock_b,1,Nbasis,1,Nbasis,Nbasis,Nbasis,1)
#endif 
      Endif 

      Do I = 1, Nocc_a
      Do J = 1, Nocc_a
         Fockoo_a(J,I) = Fock_a(J,I)
         Fockoo_od_a(J,I) = Fock_a(J,I)
      Enddo
      Enddo 

      Do I = 1, Nocc_b
      Do J = 1, Nocc_b
         Fockoo_b(J,I) = Fock_b(J,I)
         Fockoo_Od_b(J,I) = Fock_b(J,I)
      Enddo
      Enddo 

#ifdef _DEBUG_LVLM
      call output(Fockoo_a,1,Nocc_a,1,Nocc_a,Nocc_a,Nocc_a,1)
      call output(Fockoo_b,1,Nocc_b,1,Nocc_b,Nocc_b,Nocc_b,1)
#endif       
      If (Frac_occ) Then
         Do I = 1, Nvrt_a
         Do J = 1, Nvrt_a
            Fockvv_a(J,I) = Fock_a(J,I)
            Fockvv_od_a(J,I) = Fock_a(J,I)
         Enddo
         Enddo

         Do I = 1, Nvrt_b
         Do J = 1, Nvrt_b
            Fockvv_b(J,I) = Fock_b(J,I)
            Fockvv_od_b(J,I) = Fock_b(J,I)
         Enddo
         Enddo
      Else 
         Do I = 1, Nvrt_a
         Do J = 1, Nvrt_a
            Fockvv_a(J,I) = Fock_a(J+Nocc_a,I+Nocc_a)
            Fockvv_od_a(J,I) = Fock_a(J+Nocc_a,I+Nocc_a)
         Enddo
         Enddo 
         Do I = 1, Nvrt_b
         Do J = 1, Nvrt_b
            Fockvv_b(J,I) = Fock_b(J+Nocc_b,I+Nocc_b)
            Fockvv_od_b(J,I) = Fock_b(J+Nocc_b,I+Nocc_b)
         Enddo
         Enddo 
      Endif 

#ifdef _DEBUG_LVLM
      call output(Fockvv_a,1,Nvrt_a,1,Nvrt_a,Nvrt_a,Nvrt_a,1)
      call output(Fockvv_b,1,Nvrt_b,1,Nvrt_b,Nvrt_b,Nvrt_b,1)
#endif       

      If (Frac_occ) Then
         Do I = 1, Nvrt_a
         Do J = 1, Nocc_a
            Fockov_a(J,I) = Fock_a(J,I)
         Enddo
         Enddo
         Do I = 1, Nvrt_b
         Do J = 1, Nocc_b
            Fockov_b(J,I) = Fock_b(J,I)
         Enddo
         Enddo
      Else
         Do I = 1, Nvrt_a
         Do J = 1, Nocc_a
            Fockov_a(J,I) = Fock_a(J,I+Nocc_a)
         Enddo
         Enddo 
         Do I = 1, Nvrt_b
         Do J = 1, Nocc_b
            Fockov_b(J,I) = Fock_b(J,I+Nocc_b)
         Enddo
         Enddo 
      Endif 

      DO I = 1, Nocc_a
         Fockoo_od_a(I,I) = 0.0D0
      ENDDO  

      DO I = 1, Nocc_b
         Fockoo_od_b(I,I) = 0.0D0
      ENDDO  

      DO I = 1, Nvrt_a 
         Fockvv_od_a(I,I) = 0.0D0
      ENDDO  

      DO I = 1, Nvrt_b
         Fockvv_od_b(I,I) = 0.0D0
      ENDDO  

#ifdef _DEBUG_LVL0
      Write(6,*)
      call checksum("Fockoo_a  :",Fockoo_a,Nocc_a*Nocc_a)
      call checksum("Fockoo_b  :",Fockoo_b,Nocc_b*Nocc_b)
      call checksum("Fockvv_a  :",Fockvv_a,Nvrt_a*Nvrt_a)
      call checksum("Fockvv_b  :",Fockvv_b,Nvrt_b*Nvrt_b)
      call checksum("Fockov_a  :",Fockov_a,Nocc_a*Nvrt_a)
      call checksum("Fockov_b  :",Fockov_b,Nocc_b*Nvrt_b)
      call checksum("Focko_od_a:",Fockoo_od_a,Nocc_a*Nocc_a)
      call checksum("Focko_od_b:",Fockoo_od_b,Nocc_b*Nocc_b)
      call checksum("Fockv_od_a:",Fockvv_od_a,Nvrt_a*Nvrt_a)
      call checksum("Fockv_od_b:",Fockvv_od_b,Nvrt_b*Nvrt_b)
#endif 
      
      Return
      End 
      
