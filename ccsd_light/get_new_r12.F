      Subroutine Get_new_r12(Hbar_mia,Hbar_mib,Hbar_aea,Hbar_aeb,
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
     +                       T2resid_bb,T2resid_ab,D1hbar_aa,
     +                       D1hbar_bb,D2hbar_aa,D2hbar_bb,
     +                       D2hbar_ab,Asmall,T1aa,T1bb,T2aa,
     +                       T2bb,T2ab,Work,Maxcor,Nocc_a,Nocc_b,
     +                       Nvrt_a,Nvrt_b)

      Implicit Double Precision(A-H,O-Z)

#include "ccsdlight_vars.com"

      Dimension Work(Maxcor)
      Dimension Asmall(Lineq_mxcyc,Lineq_mxcyc)

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

      Dimension D1hbar_aa(Nvrt_a,Nocc_a)
      Dimension D1hbar_bb(Nvrt_b,Nocc_b)
      Dimension D2hbar_aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension D2hbar_bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension D2hbar_ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Data  Ione,Onem,One /1,-1.0D0,1.0D0/

C Lets try simple Jacobi iterations.  They are given by Hbar*T(k+1) =
C      -1         -1
C -Dhbar *R(t) - (Dhbar*Hbar*T(k) - T(k)) where R(t) is the current 
C residual vector (as an input), T(k) is the T(t) at kth iteration and
C Dhbar is the diagonal of Hbar (Hbar = Dhbar + Odhabr) (OD =off diagonals
C

C Form the Dhbar^(-1) *R(t) and keep it (this does not change from iteration
C to iteration (micro iterations).

      Laa   = Nvrt_a*Nocc_a
      Lbb   = Nvrt_b*Nocc_b
      Laaaa = Nvrt_a*Nvrt_a*Nocc_a*Nocc_a
      Lbbbb = Nvrt_b*Nvrt_b*Nocc_b*Nocc_b
      Labab = Nvrt_a*Nvrt_b*Nocc_a*Nocc_b

      Call Dscal(Laa,Onem,T1resid_aa,1)
      Call Dscal(Lbb,Onem,T1resid_bb,1)
      Call Dscal(Laaaa,Onem,T2resid_aa,1)
      Call Dscal(Lbbbb,Onem,T2resid_bb,1)
      Call Dscal(Labab,Onem,T2resid_ab,1)

      Id1invt_aa = Ione
      Id1invt_bb = Id1invt_aa + Laa
      Id2invt_aa = Id1invt_bb + Lbb
      Id2invt_bb = Id2invt_aa + Laaaa
      Id2invt_ab = Id2invt_bb + Lbbbb
      Iend       = Id2invt_ab + Labab

      If (Iend .GT. Maxcor) Call Insmem("get_new_r12",Iend,
     +                                   Maxcor)
      Memleft = Maxcor - Iend

      Call Form_dinvt(D1hbar_aa,D1hbar_bb,D2hbar_aa,D2hbar_bb,
     +                D2hbar_ab,Work(Id1invt_aa),Work(Id1invt_bb),
     +                Work(Id2invt_aa),Work(Id2invt_bb),
     +                Work(Id2invt_ab),T1resid_aa,T1resid_bb,
     +                T2resid_aa,T2resid_bb,T2resid_ab,
     +                Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

C Replace the locations of Hbar Diagonals with their inverse.

      Call Form_dinv(D1hbar_aa,D1hbar_bb,D2hbar_aa,D2hbar_bb,
     +               D2hbar_ab,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Inewt1_0_aa = Iend 
      Inewt1_0_bb = Inewt1_0_aa + Laa
      Inewt2_0_aa = Inewt1_0_bb + Lbb
      Inewt2_0_bb = Inewt2_0_aa + Laaaa
      Inewt2_0_ab = Inewt2_0_bb + Lbbbb
      Inewt1_1_aa = Inewt2_0_ab + Labab
      Inewt1_1_bb = Inewt1_1_aa + Laa
      Inewt2_1_aa = Inewt1_1_bb + Lbb
      Inewt2_1_bb = Inewt2_1_aa + Laaaa
      Inewt2_1_ab = Inewt2_1_bb + Lbbbb
      Inewt1_aa   = Inewt2_1_ab + Labab
      Inewt1_bb   = Inewt1_aa   + Laa 
      Inewt2_aa   = Inewt1_bb   + Lbb
      Inewt2_bb   = Inewt2_aa   + Laaaa
      Inewt2_ab   = Inewt2_bb   + Lbbbb
      Iend        = Inewt2_ab   + Labab

      If (Iend .GT. Maxcor) Call Insmem("get_new_r12",Iend,
     +                                   Maxcor)
      Memleft = Maxcor - Iend

      Call Form_newvec_0(Work(Id1invt_aa),Work(Id1invt_bb),
     +                   Work(Id2invt_aa),Work(Id2invt_bb),
     +                   Work(Id2invt_ab),Work(Inewt1_0_aa),
     +                   Work(Inewt1_0_bb),Work(Inewt2_0_aa),
     +                   Work(Inewt2_0_bb),Work(Inewt2_0_ab),
     +                   Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Icycle = 0
      Do K = 1, Lineq_mxcyc

         Icycle = Icycle + 1 
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
     +                    Work(Inewt1_1_aa),Work(Inewt1_1_bb),
     +                    Work(Inewt2_1_aa),Work(Inewt2_1_bb),
     +                    Work(Inewt2_1_ab),Work(Inewt1_0_aa),
     +                    Work(Inewt1_0_bb),Work(Inewt2_0_aa),
     +                    Work(Inewt2_0_bb),Work(Inewt2_0_ab),
     +                    T1aa,T1bb,T2aa,T2bb,T2ab,Work(Iend),
     +                    Memleft,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

         Call Form_newvec_1(D1hbar_aa,D1hbar_bb,D2hbar_aa,D2hbar_bb,
     +                      D2hbar_ab,T1resid_aa,T1resid_bb,
     +                      T2resid_aa,T2resid_bb,T2resid_ab,
     +                      Work(Inewt1_1_aa),
     +                      Work(Inewt1_1_bb),Work(Inewt2_1_aa),
     +                      Work(Inewt2_1_bb),Work(Inewt2_1_ab),
     +                      Work(Id1invt_aa),Work(Id1invt_bb),
     +                      Work(Id2invt_aa),Work(Id2invt_bb),
     +                      Work(Id2invt_ab),Work(Inewt1_0_aa),
     +                      Work(Inewt1_0_bb),Work(Inewt2_0_aa),
     +                      Work(Inewt2_0_bb),Work(Inewt2_0_ab),
     +                      Work(Inewt1_aa),Work(Inewt1_bb),
     +                      Work(Inewt2_aa),Work(Inewt2_bb),
     +                      Work(Inewt2_ab),Nocc_a,Nocc_b,
     +                      Nvrt_a,Nvrt_b)

        Call Daxpy(Laa,Onem,Work(Inewt1_aa),1,Work(Inewt1_0_aa),1)
        Call Daxpy(Lbb,Onem,Work(Inewt1_bb),1,Work(Inewt1_0_bb),1)
        Call Daxpy(Laaaa,Onem,Work(Inewt2_aa),1,Work(Inewt2_0_aa),1)
        Call Daxpy(Lbbbb,Onem,Work(Inewt2_bb),1,Work(Inewt2_0_bb),1)
        Call Daxpy(Labab,Onem,Work(Inewt2_ab),1,Work(Inewt2_0_ab),1)

        Dt1_aa = Ddot(Laa,Work(Inewt1_0_aa),1,Work(Inewt1_0_aa),1)
        Dt1_bb = Ddot(Lbb,Work(Inewt1_0_bb),1,Work(Inewt1_0_bb),1)
        Dt2_aa = Ddot(Laaaa,Work(Inewt2_0_aa),1,Work(Inewt2_0_aa),1)
        Dt2_bb = Ddot(Lbbbb,Work(Inewt2_0_bb),1,Work(Inewt2_0_bb),1)
        Dt2_ab = Ddot(Labab,Work(Inewt2_0_ab),1,Work(Inewt2_0_ab),1)
        
        DT_max = Max(Dt1_aa,Dt1_bb,Dt2_aa,Dt2_bb,Dt2_ab)

#ifdef _DEBUG_LVL0
        Write(6,"(5(1x,F15.9))") Dt1_aa,Dt1_bb,Dt2_aa,Dt2_bb,Dt2_ab
        Write(6,"(a,1x,F15.9))") "DT_max :", DT_max 
#endif 
        If (Dt_max .Lt. Lineq_tol) Go to 10

        Call Dcopy(Laa,Work(Inewt1_aa),1,Work(Inewt1_0_aa),1)
        Call Dcopy(Lbb,Work(Inewt1_bb),1,Work(Inewt1_0_bb),1)
        Call Dcopy(Laaaa,Work(Inewt2_aa),1,Work(Inewt2_0_aa),1)
        Call Dcopy(Lbbbb,Work(Inewt2_bb),1,Work(Inewt2_0_bb),1)
        Call Dcopy(Labab,Work(Inewt2_ab),1,Work(Inewt2_0_ab),1)

      Enddo 

 10   Continue
  
      Write(6,*)
      Write(6,"(2a,1x,I2,1x,a)")" Linear equation Hbar(Delta T) = R(t)",
     +                         " is converged in", Icycle, "iterations."

      Write(6,*)
C The solutions to the above linear equation gives the Delta(T). Add 
C this increment to the incomming T1 and T2 vectors and return

#ifdef _DEBUG_LVL0
      call checksum("Deltat1_aa:",Work(Inewt1_aa),Laa)
      call checksum("Deltat1_bb:",Work(Inewt1_bb),Lbb)
      call checksum("Deltat2_aa:",Work(Inewt2_aa),Laaaa)
      call checksum("Deltat2_bb:",Work(Inewt2_bb),Lbbbb)
      call checksum("Deltat2_ab:",Work(Inewt2_ab),Labab)
#endif
      
      Call Daxpy(Laa,One,Work(Inewt1_aa),1,T1aa,1)
      Call Daxpy(Lbb,One,Work(Inewt1_bb),1,T1bb,1)
      Call Daxpy(Laaaa,One,Work(Inewt2_aa),1,T2aa,1)
      Call Daxpy(Lbbbb,One,Work(Inewt2_bb),1,T2bb,1)
      Call Daxpy(Labab,One,Work(Inewt2_ab),1,T2ab,1)

#ifdef _DEBUG_LVL0
      call checksum("Newt1_aa  :",T1aa,Laa)
      call checksum("Newt1_bb  :",T1bb,Lbb)
      call checksum("Newt2_aa  :",T2aa,Laaaa)
      call checksum("Newt2_bb  :",T2bb,Lbbbb)
      call checksum("Newt2_ab  :",T2ab,Labab)
#endif
    

      Return
      End 
