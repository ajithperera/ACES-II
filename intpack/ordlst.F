      subroutine ordlst(olist,fullmemb,fullpopv,coord,integaxis)

#include <aces.h>
#include <const.par>
#include <mol.com>
      integer
     &    olist(fullnorb),fullmemb(natoms),fullpopv(fullnorb),
     &    integaxis(3,natoms)
      M_REAL
     &    coord(3,natoms)

      integer iatom,ixyz,icnt,i,j,ioct,ioct1

      call callstack_push('ORDLST')

      do iatom=1,natoms
        do ixyz=1,3
          if (abs(coord(ixyz,iatom)).lt.ten**(-12))
     &        coord(ixyz,iatom)=zero
        end do
      end do

      icnt=1
      do 10 i=1,fullnorb
        ioct1=9
        do 20 j=1,fullpopv(i)

          call symcntr(fullmemb(icnt),coord,integaxis,ioct)

          if(ioct.lt.ioct1) then
            ioct1=ioct
            olist(i)=fullmemb(icnt)
          endif

          icnt=icnt+1
   20   continue
   10 continue

      call callstack_pop
      return
      end
