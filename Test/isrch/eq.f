
      program main
      implicit none
      integer ia(10), i, target, isrcheq
      parameter (target=0)
      do i = 1, 10
         ia(i) = i
      end do
      i = isrcheq(10,ia(1),-1,target)
      print *, target,' is at ',i,' in ',ia
      end

