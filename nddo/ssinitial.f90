subroutine ssinitial(coef,expa,xyz,ssout,n1,n2)
use constants
implicit double precision (a-h,o-z)
double precision,dimension(:,:),intent(in)::coef,expa,xyz
integer,intent(in)::n1,n2
double precision,intent(inout)::ssout
dist=(xyz(1,1)-xyz(1,2))**2+(xyz(2,1)-xyz(2,2))**2+(xyz(3,1)-xyz(3,2))**2
         w=0.0
         do 51 it=1,n1
            do 52 jt=1,n2
               r=0.0D0
               t=0.0D0
               r=expa(it,1)
               t=expa(jt,2)
      w=w+((pi/(r+t))**(1.5D0))*dexp((-r*t/(r+t))*dist) &
     *coef(it,1)*coef(jt,2)
 52   continue
 51   continue
      ssout=w
end subroutine ssinitial

