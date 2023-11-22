      Subroutine Get_massw_coords(Coords,AtmMass,Nreals)

      Implicit Double Precision(A-H,O-Z)

      Double Precision Coords(3,Nreals)
      Dimension AtmMass(Nreals)

      Data Done,Thres /1.0D0,1.0D-09/

      Do J = 1, Nreals
         SqrtMass    = Dsqrt(AtmMass(J))
         DO I = 1, 3
            Coords(I,J) = Coords(I,J)*SqrtMass
         Enddo
      Enddo 

#ifdef _DEBUG_LVL0
      Write(6,*) "Mass-weighted reference geometry"
      do i=1,Nreals
      Write(6, "(3F10.5)") (Coords(J,I),j=1,3)
      enddo 
      Write(6,*)
#endif 
      Return
      End 


