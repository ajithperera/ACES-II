
      subroutine chklst(scr,maxcor,irrep,list,ispin)
      implicit none

c ARGUMENT LIST
      integer maxcor, irrep, list, ispin
      double precision scr(maxcor)

c INTERNAL VARIABLES
      integer nsize
      character*7 spcase(3)

c EXTERNAL FUNCTIONS
      integer idsymsz
      double precision ddot

c COMMON BLOCKS
#include "machsp.com"
#include "info.com"
#include "syminf.com"
#include "sym.com"
#include "sympop.com"
#include "flags.com"
#include "flags2.com"
#include "istart.com"
#include "files.com"
#include "iopos.com"
#include "lists.com"

c DATA STATEMENTS
      data spcase /'AAAA = ','BBBB = ','ABAB = '/

c --------------------------------------------------------------

C   o get the size of the list
      nsize = idsymsz(irrep,isytyp(1,list),isytyp(2,list))

c   o get the amplitudes
      call getall(scr,nsize,irrep,list)

C   o print the checksum/dot product
      write(*,*) spcase(ispin),ddot(nsize,scr,1,scr,1)   

      return
      end
