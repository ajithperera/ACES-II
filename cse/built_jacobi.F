      Subroutine Built_Jacobi(H1bar,Foo,Fvv,Nocc,Nvrt)

      Implicit Double Precision(A-H,O-Z)

      Dimension H1bar(Nvrt,Nocc,Nvrt,Nocc)
      Dimension Foo(Nocc,Nocc),Fvv(Nvrt,Nvrt)

      Integer I,J,A,B
 
      Do J = 1, Nocc
      Do B = 1, Nvrt 
         Do I = 1, Nocc
         Do A = 1, Nvrt 
            If (I .Eq. J) Then
               H1bar(A,I,B,J) = H1bar(A,I,B,J) + Fvv(A,B) 
            Endif 
         Enddo
         Enddo
      Enddo
      Enddo 

      Do J = 1, Nocc
      Do B = 1, Nvrt
         Do I = 1, Nocc
         Do A = 1, Nvrt
            If (A .Eq. B) Then
               H1bar(A,I,B,J) = H1bar(A,I,B,J) - Foo(I,J)
            Endif
         Enddo
         Enddo
      Enddo
      Enddo

      Return
      End 
   
