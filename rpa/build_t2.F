      Subroutine Build_t2(T2_in,T2_out,Nrow_s,Ncol_s,Nrow_t,Ncol_t,
     +                    Ioff_row,Ioff_col,Irrepx)

      Implicit Double Precision (A-H,O-Z)
      Dimension Ioff_row(8),Ioff_col(8)

      Dimension T2_in(Nrow_s*Ncol_s),T2_out(Nrow_t,Ncol_t)


      Icoff = Ioff_col(Irrepx)
      Iroff = Ioff_row(Irrepx)

      Do Icol = 1, Ncol_s
         Call Dcopy(Nrow_s,T2_in((Icol-1)*Nrow_s+1),1,
     +                     T2_out(Iroff,Icoff+Icol-1),1)
      Enddo
      
      Return 
      End 
