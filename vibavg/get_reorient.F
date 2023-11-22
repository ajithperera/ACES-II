      Subroutine Get_reorient(R,Iunitor)

      Implicit Double Precision(A-H,O-Z)

      Dimension R(3,3), D(3)
      Character*80 Wrk

      Read(Iunitor,*) Wrk
      Print*, Wrk

      Do I = 1, 3
         Read(Iunitor,*) D
         Do J =1, 3
            R(I,J) = D(J)
         Enddo
      Enddo

#ifdef _DEBUG_LVL0
      write(6,"(a)") "The reorientation matrix"
      call output(R,1,3,1,3,3,1)
#endif

      Return
      End 


