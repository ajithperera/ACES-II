          subroutine  setnumint(valao,valgradao,itotwt,ncount)

          implicit none

#include "aces.h"
#include "iuhf.com"
#include "machsp.com"
#include "kscore.com"
#include "mol.com"
#include "grid.com"
#include "int.com"
#include "intgrtflags.com"

        INTEGER  znull,pnull,setptr,zcdnt,zvalprim,zgradprim,grid,
     & zradpt,zrwt,zwtintr,zsize,zrelsiz,zrtmp,zrint,
     & ncount,ncnt,nnn,ifctr,iradpt,iangpt,zrsqrd,
     & zirrtmp

      M_REAL
     &  valao(naobasfn,maxangpts,int_numradpts,ncount),
     &  valgradao(naobasfn,maxangpts,int_numradpts,ncount,3),
     &  itotwt(ncount,int_numradpts,maxangpts)
   
       integer ivalprim,igradprim

       integer pprimao,pprimaol

        M_REAL
     &  factor,times(20),totwt

       LOGICAL evalpt
c zcdnt : the coordinate of the current point with respect to each atom

        call callstack_push('SETNUMINT')


        znull= setptr(1,F_REAL,1)
        pnull= setptr(1,F_INTEGER,1)

        zcdnt     = setptr(1,F_REAL, natoms*3)
        zrsqrd    = setptr(1,F_REAL, natoms)

       ivalprim=setptr(1,F_REAL,totprim)
       igradprim=setptr(1,F_REAL,totprim*3)

        zirrtmp = setptr(1,F_REAL, natoms)
      zradpt  = setptr(1,F_REAL, int_numradpts)
      zrwt    = setptr(1,F_REAL, int_numradpts*maxangpts)
       zwtintr = setptr(1,F_REAL, natoms)
      zsize   = setptr(1,F_REAL, natoms*natoms) 
      zrelsiz = setptr(1,F_REAL, natoms*natoms)
      zrtmp   = setptr(1,F_REAL, natoms)
      zrint   = setptr(1,F_REAL, natoms*maxangpts) 

      pprimao=setptr(1,F_INTEGER,naobasfn)
      pprimaol=setptr(1,F_INTEGER,naobasfn)

       call getrec(1,'JOBARC','NPRIMST',naobasfn,
     &              kscore(pprimao))

       call getrec(1,'JOBARC','NPRIMLT',naobasfn,
     &              kscore(pprimaol))

       call setintpath(dcore(zrelsiz),dcore(zbslrd))

c Set up ordered list of unique atoms

       do ncnt=1,ncount
         nnn=ncnt

c prakash start
c not going for unique atoms

         icntr=kscore(polist+ncnt-1)
cx         write(*,*) 'icntr=',icntr
         ifctr=kscore(pcomppopv+ncnt-1)
         factor=dble(ifctr)
         if (ncount.eq.natoms) then
            icntr=ncnt
            ifctr=1
            factor=1.d0
         end if

c prakash end

       call oct(kscore(patomchrg),dcore(zatmvc),
     &             dcore(zrij), kscore(pradgrid),
     &            dcore(zrelsiz),  dcore(zrtmp),
     &            dcore(zrint),    dcore(zradpt),
     &            dcore(zrwt),     dcore(zbslrd),
     &            dcore(zgridxyz), kscore(pgrdangpts))

       do iradpt=1,int_numradpts
            grid=kscore(pradgrid+iradpt-1)
c            write(*,*) 'grid=',kscore(pradgrid+iradpt-1)
c            write(*,*) 'iangt=',kscore(pgrdangpts+grid-1)

            do iangpt=1,kscore(pgrdangpts+grid-1)
               evalpt=.true.


          
        call symoct(kscore(patomchrg), dcore(zatmvc),
     &                     dcore(zrij),      dcore(zaij),
     &                     dcore(zcdnt),     dcore(zrsqrd),
     &                     dcore(zirrtmp),   dcore(zwtintr),
     &                     totwt,            kscore(pradgrid),
     &                     dcore(zradpt),    dcore(zrwt),
     &                     kscore(pintegaxis),dcore(zgridxyz),
     &                     dcore(zgridwt),   kscore(pgrdangpts),
     &                     iradpt,           iangpt,
     &                     grid,             evalpt)
      call abcd(itotwt,iradpt,iangpt,int_numradpts,
     &          maxangpts,totwt,ncount,ncnt)
            if (evalpt) then
      call setinteg(kscore(pnshellatom),kscore(pnshellprim),
     &                       dcore(zcdnt),      dcore(zalpha),
     &                       dcore(zrsqrd),     dcore(ivalprim),
     &                        dcore(igradprim),
     &                       ifctr,
     &                       dcore(zcoord),
     &                       kscore(pcompmemb),  kscore(pangfct),
     &                       iradpt,            iangpt,
     &                       times,ncount,
     &                       valao,valgradao,
     &                       dcore(zprimcoef),kscore(pprimao),
     &                       kscore(pprimaol),ncnt)

               end if
               end do
               end do
             end do



c               call relptr(1,F_REAL,znull)
c               call relptr(1,F_INTEGER,pnull)
   
             call callstack_pop
             return
             end 














