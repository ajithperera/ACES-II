      Subroutine Prep_fock(Fock_a,Fock_b,Fockoo_a,Fockoo_b,
     +                     Fockvv_a,Fockvv_b,Fockov_a,Fockov_b,
     +                     Fockoo_od_a,Fockoo_od_b,
     +                     Fockvv_od_a,Fockvv_Od_b,
     +                     Nbasis,Naobfns,Nocc_a,Nocc_b,Nvrt_a,
     +                     Nvrt_b,Work,Maxcor,Frac_occ)

      Implicit Double Precision(A-H,O-Z)
      Dimension Work(Maxcor)
      Logical Frac_occ 

      Dimension Fock_a(Nbasis,Nbasis),Fock_b(Nbasis,Nbasis)

      Dimension Fockoo_a(Nocc_a,Nocc_a),Fockoo_b(Nocc_b,Nocc_b)
      Dimension Fockvv_a(Nvrt_a,Nvrt_a),Fockvv_b(Nvrt_b,Nvrt_b)
      Dimension Fockoo_od_a(Nocc_a,Nocc_a),Fockoo_od_b(Nocc_b,Nocc_b)
      Dimension Fockvv_od_a(Nvrt_a,Nvrt_a),Fockvv_Od_b(Nvrt_b,Nvrt_b)
      Dimension Fockov_a(Nocc_a,Nvrt_a),Fockov_b(Nocc_b,Nvrt_b)

      Nbasis2 = Naobfns*Naobfns
      Call Getrec(20,"JOBARC","FOCKA",Nbasis2,Fock_a)
      Call Getrec(20,"JOBARC","FOCKB",Nbasis2,Fock_b)

#ifdef _DEBUG_LVLM
      call output(Fock_a,1,Nbasis,1,Nbasis,Nbasis,Nbasis,1)
      call output(Fock_b,1,Nbasis,1,Nbasis,Nbasis,Nbasis,1)
#endif 

      Iscr1 = 1
      Iscr2 = Iscr1 + Nbasis2 
      Iscr3 = Iscr2 + Nbasis2 
      Iscr4 = Iscr3 + Nbasis
      Iend  = Iscr4 + Nbasis
      If (Iend .GT. Maxcor) Call Insmem("prep_fock",Iend,Maxcor)

      Call Form_mo_Fock(Fock_a,Fock_a,Work(Iscr1),Work(Iscr2),Nbasis,
     +                  Naobfns,1) 
      Call Form_mo_Fock(Fock_b,Fock_b,Work(Iscr1),Work(Iscr2),Nbasis,
     +                  Naobfns,2)
#ifdef _DEBUG_LVL0
CSS      call output(Fock_a,1,Nbasis,1,Nbasis,Nbasis,Nbasis,1)
CSS      call output(Fock_b,1,Nbasis,1,Nbasis,Nbasis,Nbasis,1)
CSS      Call Getrec(20,"JOBARC","SCFEVLA0",Nbasis,Work(Iscr3))
CSS      Call Getrec(20,"JOBARC","SCFEVLB0",Nbasis,Work(Iscr4))
CSS      Call Dzero(Fock_a,Nbasis2)
CSS      Call Dzero(Fock_b,Nbasis2)
CSS      Do I = 1, Nbasis
CSS         Fock_a(I,I) = Work(Iscr3-1+I)
CSS      Enddo
CSS      Do I = 1, Nbasis
CSS         Fock_b(I,I) = Work(Iscr4-1+I)
CSS      Enddo
#endif 

#ifdef _DEBUG_LVLM
      call output(Fock_a,1,Nbasis,1,Nbasis,Nbasis,Nbasis,1)
      call output(Fock_b,1,Nbasis,1,Nbasis,Nbasis,Nbasis,1)
#endif       
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
      
