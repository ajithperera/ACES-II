      SUbroutine Tdee_return_mubar_dot(Work,Memleft,Iuhf,Irrepx,Iside)

      Implicit Double Precision(A-H,O-Z)
      Logical Source, Target

      Dimension Work(Memleft)

#ifdef _DEBUG_LVLM
      Write(6,"(a)") "Entered - tdee_return_mubar_dot"
      Write(6,"(a,i10)") "The memory left: ",Memleft
      Write(6,"(a)") "At entry"
      Call check_source_lists(Work,Memleft,Irrepx,Iuhf,Iside)
#endif 
      Source = .False.
      Target = .True. 
      Call Tdee_init_lists(Work,Memleft,Iuhf,Irrepx,Source,
     +                     Target,Iside)

      Call Tdee_hbar_mult(Work,Memleft,Iuhf,Irrepx,Iside)

#ifdef _DEBUG_LVLM
      Write(6,*)
      Write(6,"(a)") "At Exit"
      Call check_target_lists(Work,Memleft,Irrepx,Iuhf,Iside)
#endif 

      Return
      End
