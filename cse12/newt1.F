      Subroutine Newt1(T1aa_new,T1bb_new,T1resid_aa,T1resid_bb,
     +                 fockoo_a,fockoo_b,fockvv_a,fockvv_b,
     +                 Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Implicit Double Precision(A-H,O-Z)
      
      Dimension T1aa_new(Nvrt_a,Nocc_a)
      Dimension T1bb_new(Nvrt_b,Nocc_b)

      Dimension T1resid_aa(Nvrt_a,Nocc_a)
      Dimension T1resid_bb(Nvrt_b,Nocc_b)

      Dimension fockoo_a(Nocc_a,Nocc_a)
      Dimension fockoo_b(Nocc_b,Nocc_b)
      Dimension fockvv_a(Nvrt_a,Nvrt_a)
      Dimension fockvv_b(Nvrt_b,Nvrt_b)

      Integer I,J,A,B

#include "cse12_vars.com"
      
C AA block 

      DO I = 1, Nocc_a
      DO A = 1, Nvrt_a
         E_a = (1.0D0-Ocn_va(A))
         E_i = Ocn_oa(I)
         C   = E_a*E_i
         Daa = (Fockoo_a(I,I)*E_i-Fockvv_a(A,A)*E_a)
         If (Dabs(Daa) .Gt. Denom_tol) Then  
             T1aa_new(A,I) = T1resid_aa(A,I)/Daa
         ELse
             T1aa_new(A,I) = 0.0D0
         Endif 
#ifdef _DEBUG_LVLM
         Write(6,"(a,F15.10,1x,F15.10,2(1x,I2),1x,F15.10)") 
     +                                   "Daa,T1resid,I,A,T1aa_new:", 
     +                                    Daa,T1resid_aa(A,I),A,I,
     +                                    T1aa_new(a,i)
#endif 
      ENDDO
      ENDDO

C BB block 

      DO I = 1, Nocc_b
      DO A = 1, Nvrt_b
         E_a = (1.0D0-Ocn_vb(a))
         E_i = Ocn_ob(i)
         C   = E_a*E_i
         Dbb = (Fockoo_b(i,i)*E_i-Fockvv_b(a,a)*E_a)
         If (Dabs(Dbb) .Gt. Denom_tol) Then  
            T1bb_new(a,i) = T1resid_bb(a,i)/Dbb
         Else 
            T1bb_new(a,i) = 0.0D0
         Endif 
      ENDDO
      ENDDO

#ifdef _DEBUG_LVL0
      L_aa = Nocc_a*Nvrt_a
      L_bb = Nocc_b*Nvrt_b
      call checksum("T1aa      :",T1aa_new,L_aa)
      call checksum("T1bb      :",T1bb_new,L_bb)
      do i=1,Nocc_a
      do a=1,Nvrt_a
      Write(6,"(2I2,1x,F15.10)") a,i,T1aa_new(a,i)
     +                                             
      enddo
      enddo
#endif

      Return
      End
