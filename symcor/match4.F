      function match4(string1, string2, n0)
c     
c     Apparently the strings on input can be offset by a space !? 
c     This is a stupid fix to compare the two. At least n characters are the same.
c     If there is no starting space all 4 characters should match
c
      implicit none

      character*4 string1, string2
      integer i, n, icount, jcount, n0
      logical match, match4, found, print
c
      print = .false.
      if (n0 .le. 0 .or. n0 .gt. 4) then
         write(6,*) ' Wrong arguments in match4 ', n0
         call aces_exit(1)
      endif
c     
      if (print) then
         write(6,*) ' @match4 string info '
         write(6,*) ' string1 ', string1
         write(6,*) ' string2 ', string2
         do i = 1, 4
            write(6,*) i, string1(i:i), string2(i:i),
     $           string1(i:i) .eq. string2(i:i)
         enddo
      endif
c     
      icount = 0
      found = .false.
      do i = 1, 4
         if (string1(i:i) .eq. ' ' .and. .not. found) then
            icount = icount + 1
         else
            found = .true.
         endif
      enddo
      if (print)
     $     write(6,*) ' number of leading blanks string1 ', icount
      icount = icount + 1

      jcount = 0
      found = .false.
      do i = 1, 4
         if (string2(i:i) .eq. ' ' .and. .not. found) then
            jcount = jcount + 1
         else
            found = .true.
         endif
      enddo
      if (print)
     $     write(6,*) ' number of leading blanks string2 ', jcount
      jcount = jcount + 1
c
      n= min(4-icount+1, 4-jcount+1)
      n = min(n, n0)
      if (print) write(6,*) ' matching length ', n
      match = .true.
      if (icount .eq. 0 .and. jcount .eq. 0) then
         do i = 1, 4
            match = match .and. string1(i:i) .eq. string2(i:i)
         enddo
      else
         do i = 1, n
            match = match .and. string1(icount:icount)
     $           .eq. string2(jcount:jcount)
            icount = icount + 1
            jcount = jcount + 1
         enddo
      endif
c     
      match4 = match
      if (match4 .and. print) write(6,*) ' Strings MATCH '
      return
      end
