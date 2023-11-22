       subroutine invert_vector(length,orig,power,res_mat,nsing)

#include <aces.h>
#include <const.par>
#include <machsp.com> 
#include <kscore.com>
#include <mol.com>

      integer length,i,j,idamax,nsing

      M_REAL
     &    orig(length),power,res_mat(length,length),
     &    scalefac,
     &    rcond
      parameter(rcond=1.d-8)
     
c This program returns the inverse of orig and 
c raises it to 1/power.
c It is assumed that the sign of the elements and power are
c considered before calling the subroutine.

      call callstack_push('INVERT_VECTOR')

      call dzero(res_mat,length*length)

      i=idamax(length,orig,1)

c Quick return
      if ((dabs(orig(i))).gt.(1.d-14)) then
        scalefac=one/(orig(i))
        nsing=0
        do 10 j=1,length
           if ((dabs(scalefac*orig(j))).lt.rcond) then
              res_mat(j,j)=zero
              nsing=nsing+1
           else
              res_mat(j,j)=one/((orig(j))**(one/power))
           end if
  10    continue
      end if
      call callstack_pop 
      return
      end
