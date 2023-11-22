         subroutine findmo(                             
     &   pcoeffA,pcoeffB,valmo,gradmo,xnat,
     &   xocc,
     &   ncount,iradpt,
     &   groupoint,dens,graddens,valao,
     &   valgradao,intnumradpts,max_angpts,ncnt)


      implicit none
#include <aces.h>
#include <maxangshell.par>
#include <iuhf.com>
#include <intgrtflags.com>
#include <mol.com>
#include <grid.com>
#include <int.com>
#include <dftfunc.com>


      integer groupoint,ncount,iradpt,ncnt

      integer intnumradpts,max_angpts
      M_REAL
     &    pcoeffA(naobasfn,nbastot),
     &    pcoeffB(naobasfn,nbastot),
     &    valmo(nbastot,groupoint,2),
     & gradmo(nbastot,groupoint,3,2),xnat(nbastot,groupoint),
     &    xocc(nbastot,2)

      M_REAL
     &   valao(naobasfn,max_angpts,intnumradpts,ncount),
     &   valgradao(naobasfn,max_angpts,intnumradpts,ncount,3)


      M_REAL
     &   dens(1,groupoint,2),graddens(1,groupoint,3,2)
      integer
     &    ispin,xyz,ip,ik,il


      call callstack_push('FINDMO')
      
        call xgemm('t','n',nbastot,groupoint,naobasfn,
     &              1.d0, pcoeffA,naobasfn,
     &       valao(1,1,iradpt,ncnt),naobasfn,
     &              0.d0,valmo(1,1,1), nbastot)

      if(iuhf .eq. 1) then
       call xgemm('t','n',nbastot,groupoint,naobasfn,
     &              1.d0, pcoeffB,naobasfn,
     &       valao(1,1,iradpt,ncnt),naobasfn,
     &              0.d0,valmo(1,1,2), nbastot)
      end if 
       do ispin=1,iuhf+1   
        call vecprd_pra(valmo(1,1,ispin),valmo(1,1,ispin),
     &          xnat,nbastot,groupoint)
 
       call xgemm('t','n',1,groupoint,nbastot,
     &            1.d0,xocc(1,ispin),nbastot,
     &            xnat,nbastot,
     &            0.d0,dens(1,1,ispin),1)
       end do
          
            do xyz=1,3
        call xgemm('t','n',nbastot,groupoint,naobasfn,
     &                    1.d0, pcoeffA,  naobasfn,
     &        valgradao(1,1,iradpt,ncnt,xyz), naobasfn,
     &                    0.d0,gradmo(1,1,xyz,1),nbastot)

        call vecprd_pra(valmo(1,1,1),gradmo(1,1,xyz,1),
     &                     xnat,nbastot,groupoint)

         call xgemm('t','n',1,groupoint,nbastot,
     &            2.d0,xocc(1,1),nbastot,
     &            xnat,nbastot,
     &           0.d0,graddens(1,1,xyz,1),1)


       if(iuhf .eq.1 ) then
       call xgemm('t','n',nbastot,groupoint,naobasfn,
     &                    1.d0, pcoeffB,  naobasfn,
     &        valgradao(1,1,iradpt,ncnt,xyz), naobasfn,
     &                    0.d0,gradmo(1,1,xyz,2),nbastot)
       
       call vecprd_pra(valmo(1,1,2),gradmo(1,1,xyz,2),
     &                     xnat,nbastot,groupoint)
          
       call xgemm('t','n',1,groupoint,nbastot,
     &            2.d0,xocc(1,2),nbastot,
     &            xnat,nbastot,
     &           0.d0,graddens(1,1,xyz,2),1)

      end if
            end do

       if(iuhf .eq.0)then
         do ik=1,groupoint
           dens(1,ik,2)=dens(1,ik,1)
          end do
            do il=1,3
            do ip=1,groupoint
              graddens(1,ip,il,2)=graddens(1,ip,il,1)
            end do
         end do
       end if  
   
      call callstack_pop
      return
      end
