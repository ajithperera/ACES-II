       subroutine mkdenpra(CA,CC,O,D,CB,order)

       implicit none
#include <aces.h>
#include <iuhf.com>
#include <mol.com>

        M_REAL
     &  CA(naobasfn,nbastot),cc(naobasfn,naobasfn,iuhf+1),
     &   O(nbastot,2),D(naobasfn,nbastot,iuhf+1),
     &  CB(naobasfn,nbastot)

        integer is,ispin,i,j
        integer ji, ii,jj

        integer order(naobasfn)
        
        
         
         do J=1,nbastot
            do i=1,naobasfn
                D(i,j,1)=CA(i,j)*O(j,1)
            end do
         end do
    
        if(iuhf .eq.1) then
            do J=1,nbastot
            do i=1,naobasfn
                D(i,j,2)=CB(i,j)*O(j,2)
            end do
         end do
         end if 
c
       if ( iuhf .eq. 0) then
       call xgemm('N','T',naobasfn,naobasfn,nbastot,
     &             2.0d0,D(1,1,1),naobasfn,
     &                  D(1,1,1),naobasfn,
     &             0.0d0,CC(1,1,1),naobasfn)
      else
       do is=1,iuhf+1 
       call xgemm('N','T',naobasfn,naobasfn,nbastot,
     &             1.0d0,D(1,1,is),naobasfn,
     &                  D(1,1,is),naobasfn,
     &             0.0d0,CC(1,1,is),naobasfn)
       end do  
      end if

      return
      end
 
