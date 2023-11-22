      Subroutine Get_flags_c4(FlgACESElec, FlgACESGrad, FlgACESHess,
     &                        FlgACESForc)

      Implicit Double Precision(A-H, O-Z)
C  
      Logical FlgACESElec, FlgACESGrad, FlgACESHess, FlgACESForc
C
      Call Getrec_c4(0,'JOBARC', 'TOTENERG', JLength, A2Ener)
      If (JLength .GT. 0) FlgACESElec = .TRUE.
      Call Getrec_c4(0,'JOBARC', 'GRADIENT', JLength, VGrad)
      If (JLength .GT. 0) FlgACESGrad = .TRUE.
      Call Getrec_c4(0,'JOBARC', 'CART_HES', JLength, VHess)
      If (JLength .GT. 0) FlgACESHess = .TRUE.
      Call Getrec_c4(0, "JOBARC", 'FORCECON', JLength, Omega)
      If (JLength .GT. 0) FlgACESForc = .TRUE.

      Return
      End

