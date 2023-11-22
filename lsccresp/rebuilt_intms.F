      Subroutine Rebuilt_intms(Work,Maxcor,Iuhf)

      Implicit Double Precision(A-H,O-Z)

      Dimension Work(Maxcor)
      Logical Hbar_4lccsd

C We need to rebuilt all the intermediate with new Ts so that
C new Hbar can be formed.

C Both of the following variable says that we are not doing linear
C CCSD.

      Intype      = 2 
      Hbar_4lccsd = .False.

      Call Genint(Work,Maxcor,Iuhf,Intype,Hbar_4lccsd)
      Call Feacont(Work,Maxcor,Iuhf)
      Call Fmicont(Work,Maxcor,Iuhf)
      Call Fmecont(Work,Maxcor,Iuhf,1)
      If (iuhf .Ne. 0) Call Fmecont(Work,Maxcor,Iuhf,2)

      Return
      End 
  

    
