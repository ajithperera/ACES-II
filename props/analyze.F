      Subroutine Analyze(Dens_a,Dens_b,Work,Maxcor,Nbfns,Norbs)

      Implicit Double Precision(A-H,O-Z)

      Dimension Dens_a(Nbfns*Nbfns,Norbs)
      Dimension Dens_b(Nbfns*Nbfns,Norbs)
      Dimension Work(Maxcor)

      Data Done /1.0D0/
      
      Nbfns2 = Nbfns*Nbfns 

      Idens_a = I000
      Idens_b = Idens_a + Nbfns2 
      Iend    = Idens_b + Nbfns2 

      Call Dzero(Work(Idens_a),Nbfns2)
      Call Dzero(Work(Idens_b),Nbfns2)

      Do Iorb = 1, Norbs
         Call Daxpy(Nbfns2,Done,Dens_a(1,Iorb),1,Work(Idens_a),1)
      Enddo 
      Do Iorb = 1, Nocc_b
         Call Daxpy(Nbfns2,Done,Dens_b(1,Iorb),1,Work(Idens_b),1)
      Enddo 

#ifdef _DEBUG_LVL0
        write(6,"(a)") "The alpha density matrix"
        call output(Work(Idens_a),1,Nbfns,1,Nbfns,Nbfns,Nbfns,1)
        write(6,"(a)") "The beta density matrix"
        call output(Work(Idens_b),1,Nbfns,1,Nbfns,Nbfns,Nbfns,1)
#endif

      Return
      End 


 
