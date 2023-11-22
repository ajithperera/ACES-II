      subroutine ordlst(olist,compmemb,comppopv,coord,integaxis)

      implicit none
#include <aces.h>
#include <mol.com>
      integer
     &    olist(compnorb),compmemb(natoms),comppopv(compnorb),
     &    integaxis(3,natoms)
      M_REAL
     &    coord(3,natoms)

      integer iatom,ixyz,icnt,i,j,ioct,ioct1

      call callstack_push('ORDLST')

      do iatom=1,natoms
        do ixyz=1,3
          if (abs(coord(ixyz,iatom)).lt.1.d-12)
     &        coord(ixyz,iatom)=0.d0
        end do
      end do

      icnt=1
      do 10 i=1,compnorb
        ioct1=9
        do 20 j=1,comppopv(i)

          call symcntr(compmemb(icnt),coord,integaxis,ioct)

          if(ioct.lt.ioct1) then
            ioct1=ioct
            olist(i)=compmemb(icnt)
          endif

          icnt=icnt+1
   20   continue
   10 continue

      call callstack_pop
      return
      end
