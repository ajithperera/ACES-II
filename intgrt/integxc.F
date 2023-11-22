      subroutine integxc(v,z1,screx,scr,scfkslastiter,valao,
     &     valgradao,totwt,intnumradpts,max_angpts,ncount,kshf)

#include <aces.h>
      implicit none
#include <iuhf.com>
#include <mol.com>
c#include <scfflags.com>
#include "int.com"
#include "intgrtflags.com"
#include "jodaflags.com"
#include <kscore.com>

C prakash ebd

       integer ncount,max_angpts,intnumradpts
       logical kshf

      M_REAL
     &     valao(naobasfn,max_angpts,intnumradpts,ncount),
     &   valgradao(naobasfn,max_angpts,intnumradpts,ncount,3),
     &    totwt(ncount,intnumradpts,max_angpts)

      integer 
     &    z1(naobasfn,2)
      M_REAL
     &    v(naobasfn,naobasfn,iuhf+1),
     &    screx(naobasfn,naobasfn),
     &    scr(naobasfn,naobasfn,iuhf+1)

       logical scfkslastiter
      call callstack_push('INTEGXC')

        int_ks=.true.

          int_ks_finaliter=scfkslastiter

      memknown=-1

      call anlytc
         
      idns=0
      call primvec
      call numinteff(valao,valgradao,totwt,.false.,
     &   intnumradpts,max_angpts,ncount,scr,kshf)
      call fixv(v,z1,screx,scr)
      call callstack_pop
      return
      end
