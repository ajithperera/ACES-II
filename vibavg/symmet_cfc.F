      Subroutine Symmet_cfc(Cfc,Nmodes)

      Implicit Double Precision(A-H,O-Z)
      Dimension Cfc(Nmodes,Nmodes,Nmodes)

      Data One,Three,Diffmax_t/1.0D0,3.0D0,1.0D-04/

      Do I = 1, Nmodes
         Do J = 1, Nmodes
            Do K = 1, Nmodes

               Diff1   = Cfc(I,J,K) - Cfc(J,I,K)
               Diff2   = Cfc(I,J,K) - Cfc(J,K,I)
               Diff3   = Cfc(I,J,K) - Cfc(K,J,I)
               Diff4   = Cfc(I,J,K) - Cfc(K,I,J)
               Diff5   = Cfc(I,J,K) - Cfc(I,K,J)
               DiffMax = Max(Dabs(Diff1),Dabs(Diff2),Dabs(Diff3),
     +                       Dabs(Diff4),Dabs(Diff5))

               If (DiffMax .Gt. Diffmax_t) Then
                  Write(6,*) 
                  Write(6,10) "Warning: The maximum deviation",
     +                         DiffMax, " is greater than ",
     +                         Diffmax_t
                  Write(6,"(2a,3(1x,I2),2a,ES12.5E2,a))") 
     +                                       " for the vibrational",
     +                                       " modes", I,J,K,": ",
     +                                       "cubic force constant = ",
     +                                        Cfc(I,J,K)," a.u."
               Endif 

            Enddo
         Enddo
      Enddo

  10  Format(1x,a,ES12.5E2,a,ES8.2E2)

#ifdef _SYMMETRIZE
      Cfc(I,J,K) = Cfc(I,J,K) + Cfc(I,K,J) + Cfc(K,J,I)
      Nmodes3  = Nmodes**3 
      Onethird = One/Three 
      Call Dscal(Nmodes3,Onethird,Cfc,1)
#endif

      Return
      ENd
