       subroutine ppp(nbas,fock)
        implicit none
        integer j, i,nbas
         integer fock(nbas)

CSSS           do j=1,nbas
CSSS         write(*,*)j,    fock(j) 
CSSS           end do
           return
          end  
