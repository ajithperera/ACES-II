      Subroutine Form_jacobi(Hbar_mia,Hbar_mib,Hbar_aea,Hbar_aeb,
     +                       Hbar_mea,Hbar_meb,Hbar_ijkl1,
     +                       Hbar_ijkl2,Hbar_ijkl3,Hbar_abef1,
     +                       Hbar_abef2,Hbar_abef3,Hbar_mbej1, 
     +                       Hbar_mbej2,Hbar_mbej3,Hbar_mbej4,
     +                       Hbar_mbej5,Hbar_mbej6,Hbar_abci1,
     +                       Hbar_abci2,Hbar_abci3,Hbar_abci4,
     +                       Hbar_aibc1,Hbar_aibc2,Hbar_aibc3,
     +                       Hbar_aibc4,Hbar_ijka1,Hbar_ijka2,
     +                       Hbar_ijka3,Hbar_ijka4,Hbar_iajk1,
     +                       Hbar_iajk2,Hbar_iajk3,Hbar_iajk4,
     +                       Hbar_ijab1,Hbar_ijab2,Hbar_ijab3,
     +                       T1resid_aa,T1resid_bb,T2resid_aa,
     +                       T2resid_bb,T2resid_ab,T1aa,T1bb,
     +                       T2aa,T2bb,T2ab,CC_jacobi,
     +                       Work,Maxcor,Ileng_cc_jacobi,
     +                       Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Implicit Double Precision(A-H,O-Z)

#include "ccsdlight_vars.com"

      Dimension Work(Maxcor)
      Dimension CC_jacobi(Ileng_cc_jacobi,Ileng_cc_jacobi)

      Dimension Hbar_abef1(Nvrt_a,Nvrt_a,Nvrt_a,Nvrt_a)
      Dimension Hbar_abef2(Nvrt_b,Nvrt_b,Nvrt_b,Nvrt_b)
      Dimension Hbar_abef3(Nvrt_a,Nvrt_b,Nvrt_a,Nvrt_b)

      Dimension Hbar_ijkl1(Nocc_a,Nocc_a,Nocc_a,Nocc_a)
      Dimension Hbar_ijkl2(Nocc_b,Nocc_b,Nocc_b,Nocc_b)
      Dimension Hbar_ijkl3(Nocc_a,Nocc_b,Nocc_a,Nocc_b)

      Dimension Hbar_ijab1(Nocc_a,Nocc_a,Nvrt_a,Nvrt_a)
      Dimension Hbar_ijab2(Nocc_b,Nocc_b,Nvrt_b,Nvrt_b)
      Dimension Hbar_ijab3(Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Dimension Hbar_ijka1(Nocc_a,Nocc_a,Nocc_a,Nvrt_a)
      Dimension Hbar_ijka2(Nocc_b,Nocc_b,Nocc_b,Nvrt_b)
      Dimension Hbar_ijka3(Nocc_a,Nocc_b,Nocc_a,Nvrt_b)
      Dimension Hbar_ijka4(Nocc_b,Nocc_a,Nocc_b,Nvrt_a)

      Dimension Hbar_iajk1(Nocc_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension Hbar_iajk2(Nocc_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension Hbar_iajk3(Nocc_a,Nvrt_b,Nocc_a,Nocc_b)
      Dimension Hbar_iajk4(Nocc_b,Nvrt_a,Nocc_b,Nocc_a)

      Dimension Hbar_abci1(Nvrt_a,Nvrt_a,Nvrt_a,Nocc_a)
      Dimension Hbar_abci2(Nvrt_b,Nvrt_b,Nvrt_b,Nocc_b)
      Dimension Hbar_abci3(Nvrt_a,Nvrt_b,Nvrt_a,Nocc_b)
      Dimension Hbar_abci4(Nvrt_b,Nvrt_a,Nvrt_b,Nocc_a)

      Dimension Hbar_aibc1(Nvrt_a,Nocc_a,Nvrt_a,Nvrt_a)
      Dimension Hbar_aibc2(Nvrt_b,Nocc_b,Nvrt_a,Nvrt_b)
      Dimension Hbar_aibc3(Nvrt_a,Nocc_b,Nvrt_a,Nvrt_b)
      Dimension Hbar_aibc4(Nvrt_b,Nocc_a,Nvrt_b,Nvrt_a)

      Dimension Hbar_mbej1(Nocc_a,Nvrt_a,Nvrt_a,Nocc_a)
      Dimension Hbar_mbej2(Nocc_b,Nvrt_b,Nvrt_b,Nocc_b)
      Dimension Hbar_mbej3(Nocc_a,Nvrt_b,Nvrt_a,Nocc_b)
      Dimension Hbar_mbej4(Nocc_b,Nvrt_a,Nvrt_b,Nocc_a)
      Dimension Hbar_mbej5(Nocc_a,Nvrt_b,Nocc_a,Nvrt_b)
      Dimension Hbar_mbej6(Nocc_b,Nvrt_a,Nocc_b,Nvrt_a)

      Dimension Hbar_mia(Nocc_a,Nocc_a)
      Dimension Hbar_mib(Nocc_b,Nocc_b)
      Dimension Hbar_aea(Nvrt_a,Nvrt_a)
      Dimension Hbar_aeb(Nvrt_b,Nvrt_b)
      Dimension Hbar_mea(Nocc_a,Nvrt_a)
      Dimension Hbar_meb(Nocc_b,Nvrt_b)

      Dimension T1aa(Nvrt_a,Nocc_a)
      Dimension T1bb(Nvrt_b,Nocc_b)
      Dimension T2aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2ab(Nvrt_a,Nvrt_B,Nocc_a,Nocc_b)

      Dimension T1resid_aa(Nvrt_a,Nocc_a)
      Dimension T1resid_bb(Nvrt_b,Nocc_b)
      Dimension T2resid_aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2resid_bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2resid_ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Data Ione,Onem,DOne,Dnull,Thres/1,-1.0D0,1.0D0,0.0D0,
     +                                1.0D-08/

      Laa   = Nvrt_a*Nocc_a
      Lbb   = Nvrt_b*Nocc_b
      Laaaa = Nvrt_a*(Nvrt_a-1)*Nocc_a*(Nocc_a-1)/4
      Lbbbb = Nvrt_b*(Nvrt_b-1)*Nocc_b*(Nocc_b-1)/4
      Labab = Nvrt_a*Nvrt_b*Nocc_a*Nocc_b

      Call Dscal(Laa,Onem,T1resid_aa,1)
      Call Dscal(Lbb,Onem,T1resid_bb,1)
      Call Dscal(Laaaa,Onem,T2resid_aa,1)
      Call Dscal(Lbbbb,Onem,T2resid_bb,1)
      Call Dscal(Labab,Onem,T2resid_ab,1)

      Ivec1_in_aa = Ione
      Ivec1_in_bb = Ivec1_in_aa + Laa
      Ivec2_in_aa = Ivec1_in_bb + Lbb   
      Ivec2_in_bb = Ivec2_in_aa + Laaaa
      Ivec2_in_ab = Ivec2_in_bb + Lbbbb
      Ivec1_ou_aa = Ivec2_in_ab + Labab
      Ivec1_ou_bb = Ivec1_ou_aa + Laa
      Ivec2_ou_aa = Ivec1_ou_bb + Lbb   
      Ivec2_ou_bb = Ivec2_ou_aa + Laaaa
      Ivec2_ou_ab = Ivec2_ou_bb + Lbbbb
      Iend        = Ivec2_ou_ab + Labab
      Memleft = Maxcor - Iend 
      If (Iend .GT. Maxcor) Call Insmem("Form_jacobi",Iend,
     +                                   Maxcor)

#ifdef _DEBUG_LVLM
         Call Dcopy(Laa,T1aa,1,Work(Ivec1_in_aa),1)
         Call Dcopy(Lbb,T1bb,1,Work(Ivec1_in_bb),1)
         Call Dcopy(Laaaa,T2aa,1,Work(Ivec2_in_aa),1)
         Call Dcopy(Lbbbb,T2bb,1,Work(Ivec2_in_bb),1)
         Call Dcopy(Labab,T2ab,1,Work(Ivec2_in_ab),1)

         Call Hbar_t_mult(Hbar_mia,Hbar_mib,Hbar_aea,Hbar_aeb,
     +                    Hbar_mea,Hbar_meb,Hbar_ijkl1,
     +                    Hbar_ijkl2,Hbar_ijkl3,Hbar_abef1,
     +                    Hbar_abef2,Hbar_abef3,Hbar_mbej1,
     +                    Hbar_mbej2,Hbar_mbej3,Hbar_mbej4,
     +                    Hbar_mbej5,Hbar_mbej6,Hbar_abci1,
     +                    Hbar_abci2,Hbar_abci3,Hbar_abci4,
     +                    Hbar_aibc1,Hbar_aibc2,Hbar_aibc3,
     +                    Hbar_aibc4,Hbar_ijka1,Hbar_ijka2,
     +                    Hbar_ijka3,Hbar_ijka4,Hbar_iajk1,
     +                    Hbar_iajk2,Hbar_iajk3,Hbar_iajk4,
     +                    Hbar_ijab1,Hbar_ijab2,Hbar_ijab3,
     +                    Work(Ivec1_ou_aa),Work(Ivec1_ou_bb),
     +                    Work(Ivec2_ou_aa),Work(Ivec2_ou_bb),
     +                    Work(Ivec2_ou_ab),Work(Ivec1_in_aa),
     +                    Work(Ivec1_in_bb),Work(Ivec2_in_aa),
     +                    Work(Ivec2_in_bb),Work(Ivec2_in_ab),
     +                    T1aa,T1bb,T2aa,T2bb,T2ab,Work(Iend),
     +                    Memleft,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)
#endif 

      Len0 = Laa 
      Len1 = Len0 + Lbb
      Len2 = Len1 + Laaaa
      Len3 = Len2 + Lbbbb
      Len4 = Len3 + Labab

      Call Dzero(CC_jacobi,Ileng_cc_jacobi*Ileng_cc_jacobi)

      Do I = 1, Ileng_cc_jacobi

         Call Dzero(Work(Ivec1_in_aa),Laa)
         Call Dzero(Work(Ivec1_in_bb),Lbb)
         Call Dzero(Work(Ivec2_in_aa),Laaaa)
         Call Dzero(Work(Ivec2_in_bb),Lbbbb)
         Call Dzero(Work(Ivec2_in_ab),Labab)

         If (I .Le. Len0) Then
             Ioff = I-1 
             Work(Ivec1_in_aa+Ioff) = Done
         Elseif (I .Gt. Len0 .And. I .Le. Len1) Then
             Ioff = I-Len0-1
             Work(Ivec1_in_bb+Ioff) = Done
         Elseif (I .Gt. Len1 .And. I .Le. Len2) Then
             Ioff = I-Len1-1
             Work(Ivec2_in_aa+Ioff) = Done 
CSSS             Call Fill_vec(Work(Ivec2_in_aa),Nocc_a,Nvrt_a,I,Len1)
         Elseif (I .Gt. Len2 .And. I .Le. Len3) THen
             Ioff = I-Len2-1
             Work(Ivec2_in_bb+Ioff) = Done
CSSS             Call Fill_vec(Work(Ivec2_in_bb),Nocc_b,Nvrt_b,I,Len2)
         Elseif (I .Gt. Len3 .And. I .Le. Len4) THen
             Ioff = I-Len3-1
             Work(Ivec2_in_ab+Ioff) = Done
         Endif 

         Call Hbar_t_mult(Hbar_mia,Hbar_mib,Hbar_aea,Hbar_aeb,
     +                    Hbar_mea,Hbar_meb,Hbar_ijkl1,
     +                    Hbar_ijkl2,Hbar_ijkl3,Hbar_abef1,
     +                    Hbar_abef2,Hbar_abef3,Hbar_mbej1,
     +                    Hbar_mbej2,Hbar_mbej3,Hbar_mbej4,
     +                    Hbar_mbej5,Hbar_mbej6,Hbar_abci1,
     +                    Hbar_abci2,Hbar_abci3,Hbar_abci4,
     +                    Hbar_aibc1,Hbar_aibc2,Hbar_aibc3,
     +                    Hbar_aibc4,Hbar_ijka1,Hbar_ijka2,
     +                    Hbar_ijka3,Hbar_ijka4,Hbar_iajk1,
     +                    Hbar_iajk2,Hbar_iajk3,Hbar_iajk4,
     +                    Hbar_ijab1,Hbar_ijab2,Hbar_ijab3,
     +                    Work(Ivec1_ou_aa),Work(Ivec1_ou_bb),
     +                    Work(Ivec2_ou_aa),Work(Ivec2_ou_bb),
     +                    Work(Ivec2_ou_ab),Work(Ivec1_in_aa),
     +                    Work(Ivec1_in_bb),Work(Ivec2_in_aa),
     +                    Work(Ivec2_in_bb),Work(Ivec2_in_ab),
     +                    T1aa,T1bb,T2aa,T2bb,T2ab,
     +                    Work(Iend),Memleft,Nocc_a,Nocc_b,
     +                    Nvrt_a,Nvrt_b)

         Do J = 1 ,Ileng_cc_jacobi
            If (J .Le. Len0) Then
               Joff = J-1
               CC_jacobi(J,I)  = Work(Ivec1_ou_aa+Joff)
            Elseif (J .Gt. Len0 .And. J .Le. Len1) Then
                Joff = J-Len0-1
                CC_jacobi(J,I) = Work(Ivec1_ou_bb+Joff)
            Elseif (J .Gt. Len1 .And. J .Le. Len2) Then
                Joff = J-Len1-1
                CC_jacobi(J,I) = Work(Ivec2_ou_aa+Joff)
            Elseif (J .Gt. Len2 .And. J .Le. Len3) THen
                Joff = J-Len2-1
                CC_jacobi(J,I) = Work(Ivec2_ou_bb+Joff)
            Elseif (J .Gt. Len3 .And. J .Le. Len4) THen
                Joff = J-Len3-1
                CC_jacobi(J,I) = Work(Ivec2_ou_ab+Joff)
            Endif 
        Enddo

      Enddo 

      Write(6,*)
C Diagonalize the Jacobi matrix. The eigenvector corresponding to the
C lowest eigenvalue is increment

      Nrow    = Ileng_cc_jacobi
      Ncol    = Nrow
      Nrow2   = Nrow*Ncol
      Ieval_r = Ione
      Ieval_l = Ieval_r + Nrow
      Ievec_r = Ieval_l + Nrow
      Ievec_l = Ievec_r + Nrow2
      Iscr1   = Ievec_l + Nrow2
      Iend    = Iscr1   + Nrow2

#ifdef _DEBUG_LVL0
      Write(6,"(a)") " The Jacobi matrix"
      call output(CC_jacobi,1,Nrow,1,Ncol,Nrow,Ncol,1)
#endif

      Call Dgeev("V","V",Nrow,CC_jacobi,Nrow,Work(Ieval_r),
     +            Work(Ieval_l),Work(Ievec_r),Nrow,Work(Ievec_l),
     +            Nrow,Work(Iscr1),4*Nrow,Ierr)

#ifdef _DEBUG_LVL0
      Write(6,"(a)") "The real  eigenvalues if CC Jacobi Matrix"
      Write(6,"(4(1x,F15.10))") (Work(Ieval_r+I-1),i=1,Nrow)
      Write(6,"(a)") "The imag. eigenvalues if CC Jacobi Matrix"
      Write(6,"(4(1x,F15.10))") (Work(Ieval_l+I-1),i=1,Nrow)
      Stop
#endif

      Call Set_lr_norm(Work(Ievec_l),Work(Ievec_r),Work(Iscr1),
     +                 Nrow,Ncol)
     
      Call Dzero(CC_jacobi,Nrow2)
      Do I = 1, Nrow
         If (Abs(Work(Ieval_r+I-1)) .Gt. Thres) Then
             CC_jacobi(I,I) = 1.0d0/Work(Ieval_r+I-1)
         Else
             CC_jacobi(I,I) = 0.0D0
         Endif
      Enddo

      Call Dgemm("N","N",Nrow,Ncol,Ncol,Done,Work(Ievec_l),Nrow,
     +            CC_jacobi,Nrow,Dnull,Work(Iscr1),Nrow)

      Call Dgemm("N","T",Nrow,Ncol,Ncol,Done,Work(Iscr1),Nrow,
     +            Work(Ievec_r),Nrow,Dnull,CC_jacobi,Nrow)

      Len0 = Laa 
      Len1 = Len0 + Lbb
      Len2 = Len1 + Laaaa
      Len3 = Len2 + Lbbbb
      Len4 = Len3 + Labab
      
      Ivec_in = Ione 
      Ivec_dt = Ivec_in + Nrow
      Iend    = Ivec_dt + Nrow 
      If (Iend .GT. Maxcor) Call Insmem("Form_jacobi",Iend,
     +                                   Maxcor)

      Call Dcopy(Laa,T1resid_aa,1,Work(Ivec_in),1)
      Call Dcopy(Lbb,T1resid_bb,1,Work(Ivec_in+Len0),1)
      Call Dcopy(Laaaa,T2resid_aa,1,Work(Ivec_in+Len1),1)
      Call Dcopy(Lbbbb,T2resid_bb,1,Work(Ivec_in+Len2),1)
      Call Dcopy(Labab,T2resid_ab,1,Work(Ivec_in+Len3),1)
      
      Call Dgemm("N","N",Nrow,1,Nrow,Done,CC_Jacobi,Nrow,Work(Ivec_in),
     +            Nrow,Dnull,Work(Ivec_dt),Nrow)

C Construct T(new) = T(old) + dT

      Call Daxpy(Laa,Done,Work(Ivec_dt),1,T1aa,1)
      Call Daxpy(Lbb,Done,Work(Ivec_dt+Len0),1,T1bb,1)
      Call Daxpy(Laaaa,Done,Work(Ivec_dt+Len1),1,T2aa,1)
      Call Daxpy(Lbbbb,Done,Work(Ivec_dt+Len2),1,T2bb,1)
      Call Daxpy(Labab,Done,Work(Ivec_dt+Len3),1,T2ab,1)

#ifdef _DEBUG_LVL0
      write(6,"(a)") " The T(new) = T(old) + dT"
      call checksum("Newt1_aa  :",T1aa,Laa)
      call checksum("Newt1_bb  :",T1bb,Lbb)
      call checksum("Newt2_aa  :",T2aa,Laaaa)
      call checksum("Newt2_bb  :",T2bb,Lbbbb)
      call checksum("Newt2_ab  :",T2ab,Labab)
#endif  

      Return
      End 
