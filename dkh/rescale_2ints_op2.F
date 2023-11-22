      Subroutine rescale_2ints_op2(SOI, ACES_SOI, Nbasis, ERD_index,
     &                         ERD_scale,Symmetry)

      Implicit double precision (a-h,o-z)

      Dimension SOI(Nbasis,Nbasis), ERD_scale(Nbasis)
      Dimension ACES_SOI(Nbasis,Nbasis)
      Integer ERD_index(Nbasis)
      Logical Spherical,Symmetry 

#ifdef  _DEBUG_LVL0
      Write(*,"(a)") "Reorder Index array"
      Write(*,"(10I5)") (ERD_Index(i), i=1, Nbasis)
      Write(*,"(a)") "The scale values"
      Write(*,"(1x,10F5.2)") (ERD_scale(i), i=1, Nbasis)
#endif
      If (Symmetry) Then
         Do J = 1, Nbasis
            DO I = 1, Nbasis

               II = ERD_index(I)
               JJ = ERD_index(J)
                    
               C1 = ERD_scale(I)
               C2 = ERD_scale(J)
                  
               SCALE = C1*C2
               ACES_SOI(I,J) = SOI(I,J)*SCALE
            ENDDO
         ENDDo

      Else

         Do J = 1, Nbasis
            DO I = 1, Nbasis

               II = ERD_index(I)
               JJ = ERD_index(J)

               C1 = ERD_scale(II)
               C2 = ERD_scale(JJ)

               SCALE = C1*C2
              ACES_SOI(II,JJ) = SOI(I,J)*SCALE
            ENDDO
         ENDDO

      Endif 

      RETURN
      END
