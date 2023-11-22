      Subroutine Assign_spin_state(Coefs, Nsize, Spin_state,
     &                             Nroot,Irrep) 

      Implicit Double Precision (A-H, O-Z)
      
      Dimension Coefs(2*Nsize,2*Nsize),Nroot(8)
      Character*7 Spin_state(40)

      Thres1  = 1.0D-05
      Thres2  = 0.05D0

#ifdef _DEBUG_LVLM
      call output(coefs, 1, 2*nsize, 1, 2*nsize, 2*nsize,2*nsize,1)
#endif 

      If (Irrep .EQ. 1) Then
         Do Iroot = 1, Nroot(Irrep)
            Do I = 1, Nsize

               If (Dabs(Coefs(I,Iroot)) .GT. Thres2  .AND. 
     &             Dabs(Coefs(I+Nsize,Iroot)) .GT. Thres2) Then

                    Diff = Coefs(I,Iroot)/Coefs(I+Nsize,Iroot)
                  
                    If (Dabs(1.0D0+Diff) .lt. Thres1) Then
                        Spin_state(Iroot) = "Triplet"
                    Else if (Dabs(1.0D0-Diff) .lt. Thres1) then
                        Spin_state(Iroot) = "Singlet"
                    Else 
                        Spin_state(Iroot) = "??????"
                    Endif
               Endif

               If (Spin_state(Iroot) .EQ. "Singlet" .OR.
     &             Spin_state(Iroot) .EQ. "Triplet")
     &             Go to 10

               Spin_state(Iroot) = "??????"
            Enddo
 10      Continue
         Enddo
      Else
         Do Iroot = 1, Nroot(Irrep)
            Do I = 1, Nsize

              If (Dabs(Coefs(I,Iroot)) .NE. Thres2 .AND.
     &            Dabs(Coefs(I+Nsize,Iroot)) .NE. Thres2) Then

                    Diff = Coefs(I,Iroot)/Coefs(I+Nsize,Iroot)

                    If (Dabs(1.0D0+Diff) .lt. Thres1) Then
                        Spin_state(Iroot) = "Triplet"
                    Else if (Dabs(1.0D0-Diff) .lt. Thres1) then
                        Spin_state(Iroot) = "Singlet"
                    Else 
                        Spin_state(Iroot) = "??????"
                    Endif 

                    If (Spin_state(Iroot) .EQ. "Singlet" .OR.
     &                  Spin_state(Iroot) .EQ. "Triplet")
     &                  Go to 20

               Endif
               Spin_state(Iroot) = "??????"
            Enddo
 20      Continue
         Enddo
      Endif 

      Return
      End 
