      Subroutine  Check_xocc(Xocc,Nbasp)

      Implicit Double Precision(A-H,O-Z)

      Dimension Xocc(Nbasp,2)

      Write(6,"(6(1x,F12.5))") (Xocc(i,1),i=1,Nbasp)
      Write(6,*)
      Write(6,"(6(1x,F12.5))") (Xocc(i,2),i=1,Nbasp)
     
      Return
      End 
