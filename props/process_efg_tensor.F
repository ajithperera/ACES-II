      subroutine process_efg_tensor(efg_tensor, evec, eval,
     $     efg_max, efg_asym)
c
      double precision efg_tensor(3,3), evec(3,3), eval(3),
     $     efg_max, efg_asym, aamax
      integer iimax, ii, i
      logical print
      
c
      print = .false.
      if (print) then
         write(6,*) ' input efg tensor '
         call output(efg_tensor, 1, 3, 1, 3, 3, 3, 1)
      endif
      call eig(efg_tensor, evec, 1, 3, 1)
      if (print) then
         write(6,*) ' diagonal efg tensor '
         call output(efg_tensor, 1, 3, 1, 3, 3, 3, 1)
      endif
c
      do i = 1, 3
         eval(i) = efg_tensor(i,i)
      enddo
      if (print) then
         write(6,*) ' Eigenvalues in process_efg_tensor '
         call output(eval, 1, 1, 1, 3, 1, 3, 1)
      endif
c
c organize eigenvalues: 3: abs max value, 1, abs min value, 2 intermediate
c
      aamax = abs(efg_tensor(1,1))
      iimax = 1
      do i = 2, 3
         if (abs(efg_tensor(i,i)) .gt. aamax) then
            aamax = abs(efg_tensor(i,i))
            iimax = i
         endif
      enddo
c
      eval(3) = efg_tensor(iimax, iimax)
      ii = 0
      do i = 1, 3
         if (i .ne. iimax) then
            ii = ii + 1
            eval(ii) = efg_tensor(i,i)
         endif
      enddo
c
      if (print) then
         write(6,*) ' Sorted eigenvalues in process_efg_tensor '
         call output(eval, 1, 1, 1, 3, 1, 3, 1)
      endif
c
      efg_max = eval(3)
      if (abs(efg_max) .lt. 1.0d-7) then
         efg_asym = 0.0d0
      else
         efg_asym = abs((eval(1) - eval(2)) / efg_max)
      endif
c
      return
      end
