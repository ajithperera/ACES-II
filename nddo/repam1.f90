subroutine repam1(iatom,jatom,term,rij)
use tables
use constants

implicit double precision (a-h,o-z)
double precision,intent(inout)::term
double precision,intent(in)::rij
integer,intent(in)::iatom,jatom

term=zero
do k=1,4
! take out redundant computation    
   exp1=g2(k,species(iatom))*(rij-g3(k,species(iatom)))**2!(rij-g3(k,species(iatom)))
   exp2=g2(k,species(jatom))*(rij-g3(k,species(jatom)))**2!(rij-g3(k,species(jatom)))
   term=term+g1(k,species(iatom))*exp(-exp1)+g1(k,species(jatom))*exp(-exp2)
end do
                   
term=term*eff_core(zeff(iatom))*eff_core(zeff(jatom))/rij

end subroutine repam1
