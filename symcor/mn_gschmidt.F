      subroutine mn_gschmidt(vec,vorth,nsize,ndim,tmp)
c
c this projects out all parts of an input vector (vec)
c which lie in the space spanned by the orthogonal basis
c vorth.
c
c   |v'> = |v> - sum <i|v> |i>
c                 i 
c
c where the |i> are normalized basis vectors for the space vorth
c
c input:
c       vec : the vector which is to be orthogonalized to
c             the existing basis.
c     vorth : the basis vectors for the existing orthogonal
c             basis
c     nsize : the length of the basis vectors
c      ndim : the dimension of the orthogonal space
c       tmp : a scratch vector of length nsize
c     resid : the norm of vec, after orthogonalization and
c             before normalization
c
      implicit none
c
      integer nsize, ndim
      double precision
     &   vec(nsize),vorth(nsize,ndim),tmp(nsize)
c
c local variables
c
      double precision tol, x, fact, resid, DDOT
      integer i
      data tol /1.d-14/
c
c normalize input vector
c
c
c make sure vector is not identically zero
c
c      vec(1) = vec(1) + 1.0d-12
       resid=sqrt(DDOT(nsize,vec,1,vec,1))
       if(resid.gt.tol)then
         x = 1.0d0 / resid
         call DSCAL(nsize,x,vec,1)
       else
         write(6,*) ' norm input vector in gschmidt close to zero'
         call errex
       endif
       if (ndim .eq. 0) goto 999
c
c orthogonalize on vorth
c
      call DCOPY(nsize,vec,1,tmp,1)
      do 10 i = ndim ,1 ,-1
       fact=  DDOT(nsize,vorth(1,i),1,vec,1)
       call DAXPY(nsize,- fact,vorth(1,i),1,tmp,1)
10    continue
      call DCOPY(nsize,tmp,1,vec,1)
c
c renormalize the residual
c
       resid=sqrt(DDOT(nsize,vec,1,vec,1))
      if(resid.gt.tol)then
       x=1.0d0/resid
       call DSCAL(nsize,x,vec,1)
      else
        write(6,*) ' residual very small in mn_gschmidt', resid
        call errex
      endif 
c
c orthogonalize on vorth now in reverse order
c
      call DCOPY(nsize,vec,1,tmp,1)
      do i = 1, ndim 
       fact=  DDOT(nsize,vorth(1,i),1,vec,1)
       call DAXPY(nsize,- fact,vorth(1,i),1,tmp,1)
      enddo
      call DCOPY(nsize,tmp,1,vec,1)
c
c renormalize the residual
c
       resid=sqrt(DDOT(nsize,vec,1,vec,1))
      if(resid.gt.tol)then
       x=1.0d0/resid
       call DSCAL(nsize,x,vec,1)
      else
        write(6,*) ' residual very small in gschmidt', resid
        call errex
      endif 
c
 999  return
      end

c Local Variables: c
c mode: fortran c
c fortran-do-indent: 2 c
c fortran-if-indent: 2 c
c fortran-continuation-indent: 4 c
c fortran-comment-indent-style: nil c
c End: c
