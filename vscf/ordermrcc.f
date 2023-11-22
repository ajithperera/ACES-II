      subroutine ordermrcc(p1, p2, n, nbas, eigvec, scr, ip1, ip2,
     $    reverse, ierror, np10, np20, eigval, scr2)
c     
c     p1 contains the active character of type 1 for all orbitals in a
c     certain symmetry block. Likewise for p2. 
c     There are various possibilities
c     
c     A. reverse i false (for virtual orbitals)
c     1.    set p1 in set p2, p2 is larger  -> order orbitals as p1, p2, rest
c     2.    set p2 in set p1, p1 is larger  -> order orbitals as p2, p1, rest
c     
c     B. Reverse is true (for occupied orbitals)
c     1.    p1 in p2   -> rest, p2, p1
c     2.    p2 in p1   -> rest, p1, p2
c     
      implicit none
c     
      integer n, ierror, nbas
      double precision p1(n), p2(n), eigvec(nbas,n), scr(nbas,n),
     $    eigval(n), scr2(n)
      integer ip1(n), ip2(n)
      logical reverse, print
c     
      double precision threshl, threshs
      integer np1, np2, i, j, np10, np20
c     
      ierror = 0
c     
c     first determine which orbitals active, from projected norms in p1, p2
c     
      threshl = 0.8d0
      threshs = 0.2d0
      np1 = 0
      np2 = 0
c     
      do i = 1, n
        if (p1(i) .gt. threshl) then
          ip1(i) = 1
          np1 = np1 + 1
        elseif (abs(p1(i)) .lt. threshs) then
          ip1(i) = 0
        else
          write(6,*) ' character of orbitals unclear'
          ierror = 1
        endif
      enddo
c     
      do i = 1, n
        if (p2(i) .gt. threshl) then
          ip2(i) = 1
          np2 = np2 + 1
        elseif (abs(p2(i)) .lt. threshs) then
          ip2(i) = 0
        else
          write(6,*) ' character of orbitals unclear'
          ierror = 1
        endif
      enddo
c
      np10 = np1
      np20 = np2
c
      print = .false.
      if (print) then
        write(6,*) '@ordermrcc, np1, np2 ', np1, np2
cmn        call output(p1, 1, 1, 1, n, 1, n, 1)
        write(6,999) (ip1(i), i=1,n)
cmn        call output(p2, 1, n, 1, 1, n, 1, 1)
        write(6,999) (ip2(i), i=1,n)
      endif
 999  format(10i4)
c     
c     Let us define that p1 is always a subset of p2. 
c     
      if (np2 .lt. np1) then
c     
c     interchange meaning of p1 and p2
c     
        call icopy(n, ip1, 1, scr, 1)
        call icopy(n, ip2, 1, ip1, 1)
        call icopy(n, scr, 1, ip2, 1)
        i = np1
        np1=np2
        np2=i
      endif
c     
c     check that p1 is indeed a subset
c     
      do i = 1, n
        if (ip1(i) .gt. ip2(i)) then
          write(6,*) ' something wrong in ordermrcc'
          write(6,*) ' active orbitals do not form subsets '
          ierror = 1
        endif
      enddo
c     
c     ready to reorder
c     
      call scopy(nbas*n, eigvec, 1, scr, 1)
      call scopy(n, eigval, 1, scr2, 1)
c     
c     j is the counter for eigvec
c     
      if (reverse) then
        j = n+1
        do i = n, 1, -1
          if (ip1(i) .eq. 1) then
            ip2(i) = -1         ! meaning: don't treat this orbital anymore
            j = j - 1
            call scopy(nbas, scr(1,i), 1, eigvec(1,j), 1)
            eigval(j) = scr2(i)
          endif
        enddo
        do i = n, 1, -1
          if (ip2(i) .eq. 1) then
            ip2(i) = -1
            j = j - 1
            call scopy(nbas, scr(1,i), 1, eigvec(1,j), 1)            
            eigval(j) = scr2(i)
          endif
        enddo
        do i = n, 1, -1
          if (ip2(i) .eq. 0) then
            j = j - 1
            call scopy(nbas, scr(1,i), 1, eigvec(1,j), 1)            
            eigval(j) = scr2(i)
          endif
        enddo
        if (j .ne. 1) then
          write(6,*) ' something wrong in ordermrcc'
          call output(p1, 1, n, 1, n, n, n, 1)
          call output(p2, 1, n, 1, n, n, n, 1)
          ierror = 1
        endif
c     
      else
c     
c     not reverse
c     
        j = 0
        do i = 1, n
          if (ip1(i) .eq. 1) then
            ip2(i) = -1
c     
c     meaning: don't treat this orbital anymore
c     
            j = j + 1
            call scopy(nbas, scr(1,i), 1, eigvec(1,j), 1)
            eigval(j) = scr2(i)
          endif
        enddo
        do i = 1, n
          if (ip2(i) .eq. 1) then
            ip2(i) = -1
            j = j + 1
            call scopy(nbas, scr(1,i), 1, eigvec(1,j), 1)            
            eigval(j) = scr2(i)
          endif
        enddo
        do i = 1, n
          if (ip2(i) .eq. 0) then
            j = j + 1
            call scopy(nbas, scr(1,i), 1, eigvec(1,j), 1)            
            eigval(j) = scr2(i)
          endif
        enddo
        if (j .ne. n) then
          write(6,*) ' something wrong in ordermrcc'
          call output(p1, 1, n, 1, n, n, n, 1)
          call output(p2, 1, n, 1, n, n, n, 1)
          ierror = 1
        endif
c     
      endif
c     
c     This symmetry block has been reordered.
c     
      if (print) then
        write(6,*) 'original eigenvalues'
        call output(scr2, 1, 1, 1, n, 1, n, 1)
        write(6,*) 'resorted eigenvalues'
        call output(eigval, 1, 1, 1, n, 1, n, 1)
      endif
c
      return
      end
c
c Local Variables: c
c mode: fortran c
c fortran-do-indent: 2 c
c fortran-if-indent: 2 c
c fortran-continuation-indent: 4 c
c fortran-comment-indent-style: nil c
c End: c
