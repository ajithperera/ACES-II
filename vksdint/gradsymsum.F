       subroutine gradsymsum(ncount,gradE,gradF,cord,pop,equ,
     &                       pmapatom,gradS)

       implicit none
#include <aces.h>
#include <mol.com>
#include <grid.com>
#include <kscore.com>

        M_REAL
     &  CORD(3,natoms),zero,gradE(natoms,3,natoms),xtemp,
     &  gradF(natoms,3,natoms),gradS(natoms,3)


         integer igrid,xyz,ncount,iatom,equ(natoms,natoms)
         integer r,ri,i,ref,pop(ncount),iorbit,pmapatom(natoms)

          zero=0.0d0
       if ( natoms .ne. ncount) then
        do iorbit=1,ncount
         do i=0,pop(iorbit)-1
            ref=kscore(polist+iorbit-1)
            r=equ(ref,ref)
            ri=equ(ref,ref+i)
             do xyz=1,3
               gradF(r,xyz,r)=gradF(r,xyz,r)+gradE(r,xyz,ri)
             end do
         end do
        end do

c          write(*,*) ' grad='
c             do i=1,natoms
c           write(*,100) gradF(i,1,i),gradF(i,2,i),gradF(i,3,i)
c         end do
 100    format (3f16.10)

       do iorbit=1,ncount
         do i=1,pop(iorbit)-1
            ref=kscore(polist+iorbit-1)
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
            gradF(ri,xyz,ri)=gradF(r,xyz,r)*xtemp            
           end do
           end do
        end do

c       write(*,*) ' full grad'
c       do i=1,natoms
c           write(*,100) gradF(i,1,i),gradF(i,2,i),gradF(i,3,i)
c         end do

         do iatom=1,natoms
           do iorbit=1,ncount
              do i=1,pop(iorbit)-1
                ref=kscore(polist+iorbit-1)
                r=equ(ref,ref)
                ri=equ(ref,ref+i)
                if(iatom .eq. r  .OR.  iatom .eq. ri) then
                   pmapatom(iatom)=iorbit 
                end if
              end do
           end do
         end do

        do iatom=1,natoms
           do igrid=1,natoms
c              write(*,*) pmapatom(igrid),igrid,pmapatom(iatom),iatom
              if(pmapatom(igrid) .ne. pmapatom(iatom)) then
c              write(*,*) 'ig=',igrid,'ia=',iatom
                do xyz=1,3
                gradS(iatom,xyz)=gradS(iatom,xyz)+gradE(iatom,xyz,igrid) 
                end do
              end if
           end do
         end do

c          write(*,*)'chek'
c           do i=1,natoms
c           write(*,100) gradS(i,1),grads(i,2),gradS(i,3)
c         end do

        do xyz=1,3
        do iatom=1,natoms
           gradS(iatom,xyz)=gradS(iatom,xyz)+gradF(iatom,xyz,iatom)
        end do
        end do

        else
           do iatom =1 ,natoms
             do igrid=1,natoms
               do xyz=1,3
                gradS(iatom,xyz)=gradS(iatom,xyz)+gradE(iatom,xyz,igrid)
               end do
             end do
           end do
        end if 
 
        return
        end 






