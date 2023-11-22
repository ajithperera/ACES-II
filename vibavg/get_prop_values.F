      Subroutine Get_prop_values(Data,NMR_spncc,NMR_shifts,Ndims,
     +                           Iunit)

      Implicit Double Precision(A-H,O-Z)
     
      Dimension Data(Ndims,Ndims)
      Logical NMR_spncc,NMR_shifts
      Character*80 Wrk

      Read(Iunit,*) Wrk
      Print*, Wrk

      If (NMR_spncc) Then
          M  = Mod(Ndims,5)
          N  = Ndims/5
          Jc = 1
          Do K = 1, N
             Ib = (K-1)*5 + 1
             Do I = Ib, Ndims
                Read(Iunit,"(5F16.8)")(Data(I,J),J=Jc,Min(I-1,Jc+4))
             Enddo   
             Jc = Jc + 5
          Enddo 
          If (M .Ne. 0) Then
             Ib = N*5 + 1
             Jc = Ib
             Do I = Ib, Ndims 
                Read(Iunit,"(5F16.8)")(Data(I,J),J=Jc,Min(I-1,Jc+4))
             Enddo 
          Endif 
          
      Elseif (NMR_shifts) Then
          Read(Iunit,"(a)") Wrk
          Do I = 1, Ndims
             Read(Iunit,"(F15.8)") Data(I,1)
          Enddo  
      Endif 

#ifdef _DEBUG_LVL0
      If (NMR_spncc) Then
         write(6,"(a)") " Coupling constant matrices"
         call output(Data,1,Ndims,1,Ndims,Ndims,Ndims,1)
      Elseif (NMR_shifts) Then
         write(6,"(a)") " NMR shift vectors "
         Write(6,"(6(1x,F15.8))") (Data(I,1),I=1,Ndims)
      Endif 
#endif 

      Return
      End


