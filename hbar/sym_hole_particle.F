      Subroutine Sym_hole_particle(Win_L, Win_R, Wout_L, Wout_R, Nrows,
     &                             Ncolumns)

      Implicit Double Precision (A-H, O-Z)

      Dimension Win_L(Nrows,Ncolumns), Win_R(Ncolumns,Nrows)
      Dimension Wout_L(Nrows,Ncolumns), Wout_R(Ncolumns,Nrows)

      Do I = 1, Ncolumns
         Do J = 1, Nrows

            Wout_L(J,I) = 0.5D0 * (Win_L(J, I) + Win_R(I, J))

         Enddo
      Enddo  

      Do I = 1, Ncolumns
         Do J = 1, Nrows

            Wout_R(I,J) = 0.5D0 * (Win_R(I, J) + Win_L(J, I))

         Enddo
      Enddo 
         
      Return
      End
