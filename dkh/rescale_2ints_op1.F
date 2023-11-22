      Subroutine rescale_2ints_op1(SOI, ACES_SOI, Nbasis, ERD_index,
     &                             ERD_scale,Symmetry)

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

C If symmetry is used, we do not reorder here, but just do the scalling, 
C becuse crazy symmetry adaptation requires that the integrals are in the
C order the OED compute (i.e ZMAT centered s,s,px,py,pz,px,py,pz...) insted 
C of ZMAT centered s,s,px,px,py,py,pz,pz...order the SCF operate. The symmetry
C adapatation does tne reordering too.

      If (symmetry) Then
         Do J = 1, Nbasis
            DO I = 1, Nbasis

               II = ERD_index(I)
               JJ = ERD_index(J)

               C1 = ERD_scale(I)
               C2 = ERD_scale(J)

               SCALE = C1*C2
               ACES_SOI(I,J) = SOI(I,J)*SCALE
            ENDDO
         ENDDO
      Else 

C Since there is no symmetry adaptation we need to do both rescaling 
C reodering (see above). 

         DO J = 1, Nbasis
            DO I = 1, Nbasis

               II = ERD_index(I)
               JJ = ERD_index(J)

               C1 = ERD_scale(I)
               C2 = ERD_scale(J)

               SCALE = C1*C2
               ACES_SOI(II,JJ) = SOI(I,J)*SCALE
            ENDDO
         ENDDO

      Endif 

      RETURN
      END
