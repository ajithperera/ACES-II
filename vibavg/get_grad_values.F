      Subroutine Get_grad_values(Data,Ndims,Iunit)

      Implicit Double Precision(A-H,O-Z)
     
      Dimension Data(Ndims,3)
      Logical NMR_spncc,NMR_shifts
      Character*80 Wrk

      Read(Iunit,*) Wrk
      Do I = 1, Ndims
         Read(Iunit,*) (Data(I,J),J=1,3)
      Enddo  

#ifdef _DEBUG_LVL0
         Write(6,*)
         write(6,"(a,a)") "Gradient vectors from: ",Wrk
         Do I = 1, Ndims
         Write(6,"(5(1x,F15.10))") (Data(I,J),J=1,3)
         Enddo 
#endif 

      Return
      End


