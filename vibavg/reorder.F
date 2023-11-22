      Subroutine Reorder(W_in,Tmp,T,Ndims)

      Implicit Double Precision(A-H,O-Z)

      Dimension W_in(Ndims,Ndims)
      Dimension Tmp(Ndims,Ndims)
      Dimension T(Ndims,Ndims)

      Data Done,Dnull /1.0D0,0.0D0/

      Call Dgemm("N","T",Ndims,Ndims,Ndims,Done,W_in,Ndims,T,Ndims,
     +            Dnull,Tmp,Ndims)

      Call Dgemm("N","N",Ndims,Ndims,Ndims,Done,T,Ndims,Tmp,Ndims,
     +            Dnull,W_in,Ndims)

      Return
      End 
