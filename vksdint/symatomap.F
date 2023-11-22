       subroutine symatomap(dmap,pop,equ,ncount,pmapSA,ppjj,
     &                     cord,cosym,pk,xtmp)

       implicit none
#include <aces.h>
#include <iuhf.com>
#include <mol.com>
#include <grid.com>
#include <kscore.com>

        M_REAL
     &  dmap(natoms,natoms),distmp,xtmp(3),cord(3,natoms)

c       logical T(natoms)
       integer equ(natoms,natoms),pop(natoms)
       integer ix, pk(natoms),cosym(3),zz
       integer ncount,pmapSA(natoms,natoms,natoms)

       integer ii,io,i,ref,mem1,jo,j,jj,reff,k,pp,mem2,memA,memB
       integer ikk,ij,memC
       integer kk,ppjj(natoms,natoms)
       pp=0 
c      write(*,*) ' naount=',ncount
c       do i=1,natoms
c       write(*,*)i,cord(1,i),cord(2,i),cord(3,i)
c       end do

      do io=1,ncount
         do i=1,pop(io)-1
c          write(*,*) 'p=',pop(io)-1
           ref=kscore(polist+io-1)
           mem1=equ(ref,ref)
           mem2=equ(ref,ref+i)
c           write(*,*) 'mem1=',mem1,mem2
           do ix=1,3
         if(cord(ix,mem1) .eq. 0.d0 .and. cord(ix,mem2) .eq. 0.d0) then
                   xtmp(ix)=1.d0
            else if( cord(ix,mem2) .eq. 0.d0) then
                    xtmp(ix)=0.d0
               else
                    xtmp(ix)=cord(ix,mem1)/cord(ix,mem2)
               end if
          end do
               do jo=1,ncount
                  if( io .ne. jo) then
                   do 50 j=0,pop(jo)-1
                     reff=kscore(polist+jo-1)
                     memA=equ(reff,reff+j)
c                  write(*,*) 'memA=',memA
                        kk=0
                     do  k=0,pop(jo)-1 
                         memB=equ(reff,reff+k)
c                  write(*,*) ' memB=',memB
          distmp=dmap(memA,mem1)-dmap(memB,mem2)
c      write(*,*) 'di=',distmp,dmap(memA,mem1),dmap(memB,mem2)

       if (abs(distmp).le. 1.0D-10 .and. 
     &     memA .eq. memB         )then
           pmapSA(memA,mem1,mem2)=memB
c       write(*,*)'pma1=',pmapSA(memA,mem1,mem2),memA,mem1,mem2
                              pp=memb
                              goto 50
       else if( pp .ne. memB .AND. 
     &   abs(distmp) .le. 1.0d-10 )then
           ppjj(memA,kk)=memB
                   Kk=kk+1
                   pk(memA)=kk
c                   write(*,*) 'kk=',kk,pk(memA)
         end if
                  end do 
         cosym(1)=xtmp(1)*SIGN(1,int(cord(1,memA)))
         cosym(2)=xtmp(2)*SIGN(1,int(cord(2,memA)))
         cosym(3)=xtmp(3)*SIGN(1,int(cord(3,memA)))
         if ( pk(memA) .eq. 1) then
            pmapSA(memA,mem1,mem2)=ppjj(memA,0)
c         write(*,*)'pma2=',pmapSA(memA,mem1,mem2),memA,mem1,mem2
               goto 50
         else
        do zz=0,pk(memA)-1
          if( cosym(1) .eq. SIGN(1,int(cord(1,ppjj(memA,zz))))
     &    .AND. cosym(2) .eq. SIGN(1,int(cord(2,ppjj(memA,zz))))
     &    .AND. cosym(3) .eq. SIGN(1,int(cord(3,ppjj(memA,zz))))
     &     )  then
           pmapSA(memA,mem1,mem2)=ppjj(memA,zz)
c        write(*,*)'pma3=',pmapSA(memA,mem1,mem2),memA,mem1,mem2
        goto 50
           end if
        end do
       end if
  50                  continue
          

                   end if
                end do     
         end do
      end do     
      return
      end
 
