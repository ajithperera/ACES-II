        subroutine Get_RMS_MAX(inputVar, inputVarDim, RMS, MAX_elem)
        integer, intent(in)::inputVarDim
        double precision, intent(in)::inputVar(inputVarDim)
        double precision, intent(inout)::RMS, MAX_elem

        integer::i
        RMS=0.0d0
        MAX_elem=0.0d0
        do i=1,inputVarDim
          if (abs(inputVar(i)).gt.MAX_elem) then
                MAX_elem=abs(inputVar(i))
          endif
          RMS=RMS+inputVar(i)**2
        enddo 
        RMS=sqrt(RMS/inputVarDim)        

        end subroutine
