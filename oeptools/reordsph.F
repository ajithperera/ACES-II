      SUBROUTINE REORDSPH(ncen,nbas,angmom,swap,nbasis,spherical)
      IMPLICIT NONE
#include "maxangshell.par"      
C Reorders spherical gaussians from xxx,yyy,zzz to xyz,xyz,xyz
C Input variables
      integer nbasis,ncen
      logical spherical
C Output variables
      double precision swap(nbasis,nbasis)
C Pre-allocated Local variables
      integer nbas(ncen),angmom(nbasis)
C Local variables
      integer cen,pp,off,ncont,funcs,ang(maxangshell),cont,ll
      double precision one
      data one /1.0d0/
C - - - - - - - - - - - -- - -- - - - - - - - - - - -- - -  - -- - - - 
c      call getrec(0,'JOBARC','NBASATOM',ncen,1)
      call getrec(20,'JOBARC','NBASATOM',ncen,nbas)
      call getrec(20,'JOBARC','ANGMOSPH',nbasis,angmom)
      call dzero(swap,nbasis*nbasis)
      off=0
      do cen=1,ncen
        do ll=1,6
          ang(ll)=0
        end do
        do pp=off+1,off+nbas(cen)
          ang(angmom(pp)+1)=ang(angmom(pp)+1)+1
        end do
        do ll=1,6
C Determine number of copies of each angular momentum
          if (ang(ll).ne.0) then
            if (spherical) then
              funcs=2*ll-1
            else
              funcs=ll*(ll+1)/2
            endif
            ncont=ang(ll)/funcs
            write(6,*) 'f=',funcs,' n=',ncont,' a=',ang(ll),' l=',ll
            if ((ncont*funcs).ne.ang(ll)) call errex
            do cont=1,ncont
              do pp=1,funcs              
                swap(off+(cont-1)*funcs+pp,off+cont+(pp-1)*ncont)=one
              end do
            end do
            off=off+ang(ll)
          endif
        end do
      end do
      return
      end
      