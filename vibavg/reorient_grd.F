      Subroutine Reorient_grd(Grd,R,Work,Maxcor,Nreals)

      Implicit Double Precision(A-H,O-Z)

      Dimension Work(Maxcor)
      Dimension Grd(Nreals,3),R(3,3)

#include "machsp.com"

      Data Ione,Done,Dnull /1,1.0D0,0.0D0/

      I000 = Ione 
      Iend = I000 + Ndims*3

#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a)") "The reorientation matrix"
      call output(R,1,3,1,3,3,3,1)
#endif 

      Call Dgemm("N","N",Nreals,3,3,Done,Grd,Nreals,R,3,Dnull,
     +           Work(I000),Nreals)
      Call Dcopy(Nreals*3,Work(I000),1,Grd,1)

#ifdef _DEBUG_LVL0
       Write(6,*)
       Write(6,"(a)") "The reoriented gradients"
       Do I = 1, Nreals 
       Write(6,"(6(1x,F15.6))") (Grd(I,J),J=1,3)
       Enddo
#endif

      Return
      End
