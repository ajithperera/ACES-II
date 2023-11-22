










      Subroutine Orbdens(Evecs,Dens,Nocc,Nbfns,Ispin)

      Implicit Double Precision(A-H,O-Z)

      Dimension Evecs(Nbfns,Nbfns), Dens(Nbfns*Nbfns,Nocc)

      Data Done,Dnull /1.0D0,0.0D0/

      Call Dzero(Dens,Nbfns*Nbfns*Nocc)

      Do Iocc = 1, Nocc

         Call Mkden(Evecs(1,Iocc),Dens(1,Iocc),Nbfns)
     
        Write(6,*) 
        If (Ispin .EQ. 1) Then
        write(6,"(a)") "The alpha orbital density matrices"
        call output(Dens(1,Iocc),1,Nbfns,1,Nbfns,Nbfns,Nbfns,1)
        else 
        write(6,"(a)") "The beta orbital density matrices"
        call output(Dens(1,Iocc),1,Nbfns,1,Nbfns,Nbfns,Nbfns,1)
        endif 

      Enddo 

      Return
      End 


 
