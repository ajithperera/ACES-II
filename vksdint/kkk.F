       subroutine kkk(nbas,fock)
        implicit none
        integer j, i,nbas
         double precision fock(nbas)

           do j=1,nbas
CSSS            write(*,*) j, Fock(j)
           end do
           return
          end  
