
#ifndef _DIM1
#   define _DIM1 10
#   define _DIM2 100
#endif

      program main
      implicit none

      integer dim1, dim2, workdim
      parameter (dim1=_DIM1)
      parameter (dim2=_DIM2)
      parameter (workdim=dim1*dim2)

      double precision dmat(dim1*dim2), dmat2(dim2*dim1), dtmp
      integer iwork(workdim)
      logical lwork(workdim)
      integer m, n

      dtmp = 1.
      do n = 0, dim2-1
         do m = 0, dim1-1
            dmat(1+n*dim1+m) = dtmp
            dtmp = dtmp + 1.
         end do
      end do
      dtmp = 1.
      do n = 0, dim1-1
         do m = 0, dim2-1
            dmat2(1+n*dim2+m) = dtmp + ( dim1 * m )
         end do
         dtmp = dtmp + 1.
      end do
      m = 0
c      call dmat_prt(dmat ,dim1,dim2,dim1)
c      call dmat_prt(dmat2,dim2,dim1,dim2)
c      call dmat_xpose(dmat,dim1,dim2)
c      call dmat_xpose(dmat,dim1,dim2,dim1*dim2,iwork,workdim,m)
      call trans380(dmat,dim1,dim2,dim1*dim2,iwork,workdim,m)
c      call dxpose467(dmat,dim1,dim2,dim1*dim2,lwork,workdim)
c      call dxpose467(dmat,dim1,dim2,dim1*dim2,lwork,1)
c      call dxpose_467(dmat,dim1,dim2)
      if (m.ne.0) stop 1
c      print *, 'Checking ',dim1,' by ',dim2,':'
      n = 0
      do m = 1, (dim1*dim2)
         if (dmat(m).ne.dmat2(m)) then
            write(*,*) 'TRANSPOSE ERROR: ',dmat(m),' != ',dmat2(m)
            n = 1
         end if
      end do
      if (n.ne.0) then
         write(*,*) 'Failed test case: ',dim1,' x ',dim2
         call dmat_prt(dmat ,dim2,dim1,dim2)
         call dmat_prt(dmat2,dim2,dim1,dim2)
         write(*,*)
      end if
c      call dmat_prt(dmat ,dim2,dim1,dim2)
c      call dmat_prt(dmat2,dim2,dim1,dim2)

      end

