      Subroutine check_evals(E,Nbas)

      Implicit double precision(a-h,o-z)

      Dimension E(Nbas)


      Write(6,"(6(1x,f15.7))") (E(i),i=1,nbas)

      Return
      End 
