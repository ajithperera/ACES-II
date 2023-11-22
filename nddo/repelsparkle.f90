subroutine repelsparkle(nuclear,r,iatom,jatom,ssss)
use constants
use tables
use control
implicit double precision (a-h,o-z)
!**************
interface
subroutine repam1(iatom,jatom,term,r)
double precision,intent(inout)::term
double precision,intent(in)::r
integer,intent(in)::iatom,jatom
end subroutine repam1
end interface
!**************


double precision,intent(inout)::nuclear,r
double precision,intent(in)::ssss
integer,intent(in)::iatom,jatom
integer,dimension(9)::isite


!nuclear=eff_core(zeff(iatom))*eff_core(zeff(jatom))*27.21d0/(r/.529167d0)
!return




nuclear=zero
nuclear=one+dexp(-alpha(species(iatom))*r)+dexp(-alphasp(zeffsp(jatom))*r)
nuclear=eff_core(zeff(iatom))*spcharge(zeffsp(jatom))*ssss*nuclear
!return
! if(isite(iatom).eq.isite(jatom))then

   if(method=='AM1'  .or.  core=='GAUSSIAN'  .or. method=='PM3')then
      term=zero
do k=1,4
   exp1=g2(k,species(iatom))*(r-g3(k,species(iatom)))**2!(rij-g3(k,species(iatom)))
   term=term+g1(k,species(iatom))*exp(-exp1)
end do
                   
term=term*eff_core(zeff(iatom))*spcharge(zeffsp(jatom))/r
      nuclear=nuclear+term
   end if






end subroutine repelsparkle

