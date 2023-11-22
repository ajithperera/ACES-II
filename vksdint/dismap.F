       subroutine dismap(pop,equ,ncount,dmap,cord)

       implicit none
#include <aces.h>
#include <iuhf.com>
#include <mol.com>
#include <grid.com>
#include <kscore.com>

        M_REAL
     &  CORD(3,natoms),dmap(natoms,natoms)

       integer equ(natoms,natoms),pop(natoms)

       integer ncount
       integer io,i,ref,mem1,jo,j,reff,mem2
       do mem1=1,natoms
       do mem2=1,natoms  
c      do io=1,ncount
c         do i=0,pop(io)-1
c           ref=kscore(polist+io-1)
c           mem1=equ(ref,ref+i)
c                do jo=1,ncount
c                   if( io .ne. jo) then
c                     do j=0,pop(jo)-1
c                      reff=kscore(polist+jo-1)
c                      mem2=equ(reff,reff+j)
      dmap(mem1,mem2)=dsqrt(
     &                (cord(1,mem1)-cord(1,mem2))**2.d0+
     &                (cord(2,mem1)-cord(2,mem2))**2.d0+
     &                (cord(3,mem1)-cord(3,mem2))**2.d0
     &                      )
c                      end do
c                    end if
c                end do     
c         end do
c      end do     

        end do
        end do
      return
      end
 
