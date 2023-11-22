      Subroutine Reord_2zmat_is(Center_in,Imap,Label_in,Iatmchrg_in,
     +                          Atmmass_in,Label_out,Center_out,
     +                          Iatmchrg_out,Atmmass_out,Natoms,
     +                          Maxatms)

      Implicit Double Precision(A-H,O-Z)
     
      Dimension Imap(Natoms)
      Dimension Center_in(Maxatms,3)
      Dimension Center_out(Maxatms,3)
      Dimension Iatmchrg_in(Maxatms)
      Dimension Iatmchrg_out(Maxatms)
      Dimension Atmmass_in(Maxatms)
      Dimension Atmmass_out(Maxatms)

      Character*6 Label_in(Maxatms)
      Character*6 Label_Out(Maxatms)

      Data Dnull,Izero /0.0D0,0/

      Do Jatmvml = 1, Natoms
         Jatmzmt = Imap(Jatmvml)
         If (Jatmzmt .Ne. Izero) Then 
            Do Ixyz = 1, 3
               Center_out(Jatmzmt,Ixyz) = Center_in(Jatmvml,Ixyz)
            Enddo 
         Endif 
      Enddo 

      Do Iatmvml = 1, Natoms 
         Iatmzmt = Imap(Iatmvml)
         If (Iatmzmt .Ne. Izero) Then 
             Label_Out(Iatmzmt)    = Label_in(Iatmvml)
             Iatmchrg_out(Iatmzmt) = Iatmchrg_in(Iatmvml)
             Atmmass_out(Iatmzmt)  = Atmmass_in(Iatmvml)
         Endif 
      Enddo

      Return
      End 


