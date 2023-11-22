
      program ja
      implicit none

      integer icore, iuhf, i,j
      character czInput
      logical bDone

#include "syminf.com"
#include "sympop.com"
#include "disttype.h"

      call aces_init(icore,1,0,iuhf,.false.)
      print *, 'dirprd(1,*) = ',(dirprd(1,i),i=1,8)
      print *, 'dirprd(2,*) = ',(dirprd(2,i),i=1,8)
      print *, 'dirprd(3,*) = ',(dirprd(3,i),i=1,8)
      print *, 'dirprd(4,*) = ',(dirprd(4,i),i=1,8)
      print *, 'dirprd(5,*) = ',(dirprd(5,i),i=1,8)
      print *, 'dirprd(6,*) = ',(dirprd(6,i),i=1,8)
      print *, 'dirprd(7,*) = ',(dirprd(7,i),i=1,8)
      print *, 'dirprd(8,*) = ',(dirprd(8,i),i=1,8)
      print *, 'irpdpd(*,AltB) = ',(irpdpd(i,h_DistType_AltB),i=1,8)
      print *, 'irpdpd(*,altb) = ',(irpdpd(i,h_DistType_altb),i=1,8)
      print *, 'irpdpd(*,Ab)   = ',(irpdpd(i,h_DistType_Ab  ),i=1,8)
      call aces_fin

      end

