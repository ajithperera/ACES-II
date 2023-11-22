      Subroutine Prep_oneh(Oneh1,Oneh2,Oneh_aa,Oneh_bb,Onehoo_a,
     +                     Onehoo_b,Onehvv_a,Onehvv_b,Onehov_a,
     +                     Onehov_b,Onehoo_od_a,Onehoo_od_b,
     +                     Onehvv_od_a,Onehvv_od_b,W0_aa,W0_bb,
     +                     W0_ab,Nbasis,Naobfns,Nocc_a,Nocc_b,
     +                     Nvrt_a,Nvrt_b,Work,Maxcor,Frac_occ)

      Implicit Double Precision(A-H,O-Z)
      Dimension Work(Maxcor)
      Logical Frac_occ,Make_fock

      Dimension Oneh1(Naobfns*(Naobfns+1)/2) 
      Dimension Oneh2(Naobfns,Naobfns) 
      Dimension Oneh_aa(Naobfns,Naobfns)
      Dimension Oneh_bb(Naobfns,Naobfns)

      Dimension W0_aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension w0_bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension W0_ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Dimension Onehoo_a(Nocc_a,Nocc_a),Onehoo_b(Nocc_b,Nocc_b)
      Dimension Onehvv_a(Nvrt_a,Nvrt_a),Onehvv_b(Nvrt_b,Nvrt_b)
      Dimension Onehoo_od_a(Nocc_a,Nocc_a),Onehoo_od_b(Nocc_b,Nocc_b)
      Dimension Onehvv_od_a(Nvrt_a,Nvrt_a),Onehvv_Od_b(Nvrt_b,Nvrt_b)
      Dimension Onehov_a(Nocc_a,Nvrt_a),Onehov_b(Nocc_b,Nvrt_b)

      Ldim = Naobfns*(Naobfns+1)/2
      Call Getrec(20,"JOBARC","ONEHAO  ",Ldim,Oneh1)
      Call Expnd2(Oneh1,Oneh2,Naobfns)

#ifdef _DEBUG_LVL0
      Write(6,"(a)") "The AO basis one electron operator"
      call output(Oneh2,1,Naobfns,1,Naobfns,Naobfns,Naobfns,1)
#endif 

      Make_fock = .False.
      Nbasis2 = Naobfns*Naobfns
      Iscr1 = 1
      Iscr2 = Iscr1 + Nbasis2 
      Iscr3 = Iscr2 + Nbasis2 
      Iscr4 = Iscr3 + Nbasis
      Iend  = Iscr4 + Nbasis

      If (Iend .GT. Maxcor) Call Insmem("prep_fock",Iend,Maxcor)

      Call Form_mo_fock(Oneh2,Oneh_aa,Work(Iscr1),Work(Iscr2),Nbasis,
     +                  Naobfns,1) 
      Call Form_mo_fock(Oneh2,Oneh_bb,Work(Iscr1),Work(Iscr2),Nbasis,
     +                  Naobfns,2)
      If (Make_fock) Then
          Call Built_frac_fock_mo(W0_aa,w0_bb,W0_ab,Oneh_aa,Oneh_bb,
     +                           Nocc_a,Nocc_b,Nvrt_a,Nvrt_b,Nbasis)
      Endif 

#ifdef _DEBUG_LVL0
      Write(6,"(a)") "The MO basis alpha one electron operator"
      call output(Oneh_aa,1,Nbasis,1,Nbasis,Nbasis,Nbasis,1)
      Write(6,"(a)") "The MO basis beta one electron operator"
      call output(Oneh_bb,1,Nbasis,1,Nbasis,Nbasis,Nbasis,1)
#endif       

      Do I = 1, Nocc_a
      Do J = 1, Nocc_a
         Onehoo_a(J,I) = Oneh_aa(J,I)
         Onehoo_od_a(J,I) = Oneh_aa(J,I)
      Enddo
      Enddo 

      Do I = 1, Nocc_b
      Do J = 1, Nocc_b
         Onehoo_b(J,I) = Oneh_bb(J,I)
         Onehoo_Od_b(J,I) = Oneh_bb(J,I)
      Enddo
      Enddo 

#ifdef _DEBUG_LVLM
      call output(Onehoo_a,1,Nocc_a,1,Nocc_a,Nocc_a,Nocc_a,1)
      call output(Onehoo_b,1,Nocc_b,1,Nocc_b,Nocc_b,Nocc_b,1)
#endif       
      If (Frac_occ) Then
         Do I = 1, Nvrt_a
         Do J = 1, Nvrt_a
            Onehvv_a(J,I) = Oneh_aa(J,I)
            Onehvv_od_a(J,I) = Oneh_aa(J,I)
         Enddo
         Enddo

         Do I = 1, Nvrt_b
         Do J = 1, Nvrt_b
            Onehvv_b(J,I) = Oneh_bb(J,I)
            Onehvv_od_b(J,I) = Oneh_bb(J,I)
         Enddo
         Enddo
      Else 
         Do I = 1, Nvrt_a
         Do J = 1, Nvrt_a
            Onehvv_a(J,I) = Oneh_aa(J+Nocc_a,I+Nocc_a)
            Onehvv_od_a(J,I) = Oneh_aa(J+Nocc_a,I+Nocc_a)
         Enddo
         Enddo 
         Do I = 1, Nvrt_b
         Do J = 1, Nvrt_b
            Onehvv_b(J,I) = Oneh_bb(J+Nocc_b,I+Nocc_b)
            Onehvv_od_b(J,I) = Oneh_bb(J+Nocc_b,I+Nocc_b)
         Enddo
         Enddo 
      Endif 

#ifdef _DEBUG_LVLM
      call output(Onehvv_od_a,1,Nvrt_a,1,Nvrt_a,Nvrt_a,Nvrt_a,1)
      call output(Onehvv_od_b,1,Nvrt_b,1,Nvrt_b,Nvrt_b,Nvrt_b,1)
#endif       

      If (Frac_occ) Then
         Do I = 1, Nvrt_a
         Do J = 1, Nocc_a
            Onehov_a(J,I) = Oneh_aa(J,I)
         Enddo
         Enddo
         Do I = 1, Nvrt_b
         Do J = 1, Nocc_b
            Onehov_b(J,I) = Oneh_bb(J,I)
         Enddo
         Enddo
      Else
         Do I = 1, Nvrt_a
         Do J = 1, Nocc_a
            Onehov_a(J,I) = Oneh_aa(J,I+Nocc_a)
         Enddo
         Enddo 
         Do I = 1, Nvrt_b
         Do J = 1, Nocc_b
            Onehov_b(J,I) = Oneh_bb(J,I+Nocc_b)
         Enddo
         Enddo 
      Endif 

      DO I = 1, Nocc_a
         Onehoo_od_a(I,I) = 0.0D0
      ENDDO  

      DO I = 1, Nocc_b
         Onehoo_od_b(I,I) = 0.0D0
      ENDDO  

      DO I = 1, Nvrt_a 
         Onehvv_od_a(I,I) = 0.0D0
      ENDDO  

      DO I = 1, Nvrt_b
         Onehvv_od_b(I,I) = 0.0D0
      ENDDO  

#ifdef _DEBUG_LVL0
      Write(6,*)
      call checksum("Onehoo_a  :",Onehoo_a,Nocc_a*Nocc_a)
      call checksum("Onehoo_b  :",Onehoo_b,Nocc_b*Nocc_b)
      call checksum("Onehvv_a  :",Onehvv_a,Nvrt_a*Nvrt_a)
      call checksum("Onehvv_b  :",Onehvv_b,Nvrt_b*Nvrt_b)
      call checksum("Onehov_a  :",Onehov_a,Nocc_a*Nvrt_a)
      call checksum("Onehov_b  :",Onehov_b,Nocc_b*Nvrt_b)
      call checksum("Oneho_od_a:",Onehoo_od_a,Nocc_a*Nocc_a)
      call checksum("Oneho_od_b:",Onehoo_od_b,Nocc_b*Nocc_b)
      call checksum("Onehv_od_a:",Onehvv_od_a,Nvrt_a*Nvrt_a)
      call checksum("Onehv_od_b:",Onehvv_od_b,Nvrt_b*Nvrt_b)
#endif 
      
      Return
      End 
      
