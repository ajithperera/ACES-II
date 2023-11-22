      Subroutine Scrnc_rcc_hbar_inv(work,Memleft,Iuhf,Imult,Iloc)
  
      Implicit Double Precision (A-H,O-Z)

      Dimension Work(Memleft) 
      Dimension Iloc(8)

#ifdef _DEBUG_LVL0
      Write(6,"(1x,a)") "-----Entering scnc_rcc_hbar_inv------"
#endif
      Call Scrnc_form_rcc_hbar(Work,Memleft,Iuhf,Imult)
      Call Scrnc_form_rcc_hbar_inv(Work,Memleft,Iuhf,Iloc) 
      
      Return
      End

