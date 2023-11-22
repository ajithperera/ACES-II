       subroutine gradsymgen(iorbit,ref,gradE,cord,pop,ncount,equ,
     &                       pmapSA)

       implicit none
#include <aces.h>
#include <mol.com>

        M_REAL
     &  CORD(3,natoms),zero,gradE(natoms,3,natoms),xtemp

         integer xyz,ncount,iatom,equ(natoms,natoms)
         integer r,ri,i,ref,pop(ncount),iorbit,jatom

         integer pmapSA(natoms,natoms,natoms)

          zero=0.0d0

c         do i=1,natoms
c          write(*,100) cord(1,ref),cord(2,ref),cord(3,ref) 
c           write(*,100) gradE(i,1,ref),gradE(i,2,ref),gradE(i,3,ref)
c         end do
  
c        write(*,*) 'k=',equ(iorbit,1),equ(iorbit,pop(iorbit)) 
c        write(*,*) 'p=',pop(iorbit) 
         do i=1,pop(iorbit)-1
c          write(*,*) 'k=',equ(ref,ref),equ(ref,i+ref)

            r=equ(ref,ref)
            ri=equ(ref,ref+i)

            do xyz=1,3
        if(cord(xyz,ri) .eq. zero .and. cord(xyz,r) .eq. zero) then
                   xtemp=1.d0  
            else if( cord(xyz,ri) .eq. zero) then
                    xtemp=zero
               else
                    xtemp=cord(xyz,r)/cord(xyz,ri)
               end if
               do iatom=1,natoms     
                  if( iatom .ne. r   .and.  iatom  .ne. ri) then
                       if(pmapSA(iatom,r,ri) .ne. 0) then
                          jatom=pmapSA(iatom,r,ri)
                       else
                          jatom=iatom
                       end if
c                      write(*,*) ' jatom=',jatom,iatom,r,ri
                     gradE(Jatom,xyz,ri)=gradE(iatom,xyz,r)*
     &                                      xtemp 
                  else if ( iatom .eq. r) then
                     gradE(iatom,xyz,ri)=gradE(ri,xyz,r)*
     &                                      xtemp 
                  else if( iatom .eq. ri) then
                    gradE(iatom,xyz,ri)=gradE(r,xyz,r)*
     &                                     xtemp 
                 end if
            end do
           end do
        end do
c       write(*,*) 'ref ref+2 =',ref,ref+2
c          do i=1,natoms
c       write(*,100) cord(1,ref+1),cord(2,ref+1),cord(3,ref+1)
c       write(*,100) gradE(i,1,ref+2),gradE(i,2,ref+2),gradE(i,3,ref+2)
c         end do
 100  format (3f16.10)

c         write(*,*) pmapSA(5,2,4)


        return
        end 






