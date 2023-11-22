      Subroutine Form_vvoo_aa(T2cc,T2mp,T2,Iwork,Imaxcor,Nij,Nab,
     +                        Irrepij,Irrepab,Ispin)

      Implicit Double Precision(A-H,O-Z)
      Integer A,B,Aend,Bend

#include "syminf.com"
#include "machsp.com"
#include "syminf.com"
#include "sym.com"
#include "symoff.com"
#include "sympop.com"
#include "active_space.com"

      Dimension T2cc(Nab,Nij),T2mp(Nab,Nij),T2(Nab,Nij)

      Data Ione,Inull /1,0/
      Data Done,Dnull /1.0D0,0.0D0/

C Note the both type of T2 amplitudes are stored as AB,IJ form 

      Ioff  = Ione 
      Nij_b = Nij

      Call Dcopy(Nab*Nij,T2mp,1,T2,1)

#ifdef _DEBUG_LVLM
      Write(6,"(a)") "T2-in "
      call output(T2,1,Nab,1,Nij,Nab,Nij,1)
#endif 

      Do Irrep_r = Nirrep, 1, -1

         Irrep_j = Irrep_r
         Irrep_i = Dirprd(Irrep_j,Irrepij)

         Ncol_j  = Pop(Irrep_j,Ispin)
         Nrow_i  = Pop(Irrep_i,Ispin)
        
         Iend    = Nrow_i - Active_oo(Irrep_i,Ispin) + Ione
         Jend    = Ncol_j - ActiVe_oo(Irrep_j,Ispin) + Ione

C Select active ij pairs

         Do J = Ncol_j, Jend, -1
            Do I = Nrow_i, Iend, -1 
          
               Nij_a =  (I-Nrow_i)+(J-Ncol_j)*Nrow_i
               Nij_c  = Nij_b + Nij_a

               Call Form_vv_aa(T2cc(Ioff,Nij_c),T2mp(Ioff,Nij_c),
     +                         T2(Ioff,Nij_c),Iwork,Imaxcor,Nab,
     +                         Irrepab,Ispin)
            Enddo
         Enddo
         Nij_b = Nij_b - Ncol_j*Nrow_i 
      Enddo

      Call Dcopy(Nab*Nij,T2,1,T2cc,1)

#ifdef _DEBUG_LVLM
      Write(6,"(a)") "T2-out"
      call output(T2cc,1,Nab,1,Nij,Nab,Nij,1)
#endif 
      Return
      End
