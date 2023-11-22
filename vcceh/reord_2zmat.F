      Subroutine Reord_2zmat(JCC_in,Center_in,Imap,Label_in,JCC_out,
     +                       Label_out,Center_out,NATOMS,MXATOMS)

      Implicit Double Precision(A-H,O-Z)
      Double Precision JCC_in,Jcc_out
     
      Dimension JCC_in(MXATOMS,MXATOMS)
      Dimension JCc_out(MXATOMS,MXATOMS)
      Dimension Center_in(MXATOMS,3)
      Dimension Center_out(MXATOMS,3)
      Dimension Imap(NATOMS)
      Character*6 Label_in(MXATOMS)
      Character*6 Label_Out(MXATOMS)
      Data Dnull,Izero /0.0D0,0/

      Call Zero(Jcc_out,Natoms*Natoms)

      Do Jatms =1, Natoms
         Do Iatms = 1, Natoms 
            Jcc_in(Jatms,Iatms) = Jcc_in(Iatms,Jatms)
         Enddo
      Enddo 
      
      Do Jatmvml = 1, Natoms
         Jatmzmt = Imap(Jatmvml)
         Do Iatmvml = 1, Natoms
            Iatmzmt = Imap(Iatmvml)
            If (Jatmzmt .Ne. Izero .And. Iatmzmt .Ne. Izero) Then
               Jcc_out(Iatmzmt,Jatmzmt)  = Jcc_in(Iatmvml,Jatmvml)
            Endif 
         Enddo
      Enddo
 
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
            Label_Out(Iatmzmt) = Label_in(Iatmvml)
         Endif 
      Enddo

      Return
      End 


