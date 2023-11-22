      Subroutine Reorient_fc(Fc,R,Work,Maxcor,Ndims,Nreals)

      Implicit Double Precision(A-H,O-Z)

      Dimension Work(Maxcor)
      Dimension Fc(Ndims,Ndims),R(3,3)

#include "machsp.com"

      Data Ione /1/

      I000 = Ione 
      I010 = I000 + Ndims*Ndims 
      Iend = I010 + Ndims*Ndims 

#ifdef _DEBUG_LVLM
      Write(6,"(a)") "The incoming force constant matrix"
      call output(Fc,1,Ndims,1,Ndims,Ndims,Ndims,1)
#endif

      Call Big_reori_mat(Work(I000),R,Nreals)

#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a)") "The reorientation matrix"
      call output(R,1,3,1,3,3,3,1)
      Write(6,"(a)") "The big reorientation matrix"
      call output(Work(I000),1,Ndims,1,Ndims,Ndims,Ndims,1)
#endif

      Call Reorder(Fc,Work(I010),Work(I000),Ndims)

#ifdef _DEBUG_LVL0
      Write(6,"(a)") "The reoriented force constant matrix"
      call output(Fc,1,Ndims,1,Ndims,Ndims,Ndims,1)
#endif 

      Return
      End
