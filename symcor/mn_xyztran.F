      subroutine mn_xyztran(natom, coord, refcoord,
     $     dummy_coord, masses, 
     $     xyz_tran, scr, mxcor)
c     
c     In this subroutine we determine a transformation that
c     amounts to a permutation/rotation of the axis such that 
c     the current geometry (in COORD) lines up (best) with the coordinates
c     in the reference geometry
c     
      implicit none

      integer natom, nsize, pass1, i, icase, i0, 
     $     imap(3,3), k, ipiv(3), mxcor, j, icount
      double precision coord(3, natom), refcoord(3, natom),
     $     dummy_coord(3, natom), masses(natom),
     $     xyz_tran(3,3), overlap(3),
     $     DDOT, diff, x, y, Amat(3,3), Amatinv(3,3),
     $     wmat(3,3), wtran(3,3), scr(mxcor), diff2, xx, det3,
     $     mn_det3x3
      character*3 xyzcase(3)
      logical ok, print, print0, least_squares

c     
      nsize = 3 * natom
      print = .false.
      print0 = .false.
c     
      xyzcase(1) = 'XYZ'
      xyzcase(2) = 'ZXY'
      xyzcase(3) = 'YZX'
c     
c     When comparing with REFCOORD, we rotate XYZ to be in best agreement with 
c     the reference geometry.
c     Alternatively, when doing a scan, we might compare with nextgeom
c     This does not work at all when normal_curvy is switched on, because poor cartesians.
c     Let us try REFCOORD therefore.
c     
c     
      if (print) then
         write(6,*) ' Reference coordinates '
         call output(refcoord, 1, 3, 1, natom, 3, natom, 3)
         write(6,*) ' Current Coordinates '
         call output(coord, 1, 3, 1, natom, 3, natom, 3)
      endif
c     
c     calculate norm of coordinate and refcoordinate for good measure
c     
      x = SQRT(DDOT(3*natom, coord, 1, coord, 1))
      y = SQRT(DDOT(3*natom, refcoord, 1, refcoord, 1))
c
      if (.false.) then
c
      diff = DDOT(3*natom, coord, 1, refcoord, 1) / (x*y)
c     
c     perform better check
c     
      diff2 = -1.0d0
      do j = 1, natom
         do i = 1, 3
            xx = coord(i,j) - refcoord(i,j)
            diff2 = max(diff2, abs(xx))
         enddo
      enddo
c
      write(6,886) 1.0d0-diff, diff2
 886  format(' @xyz_tran initial: overlap deviation:', E12.4,
     $     ' Max deviation: ', E12.4)
      endif
c     
c     
c     try to use permutation algorithm first
c     
c     consider 'XYZ', 'ZXY' and 'YZX' definition of axis
c
      if (.false.) then
c
c no longer useful / needed.
c
      do icase = 1, 3
         if (icase .eq. 1) then
            imap(1, icase) = 1
            imap(2, icase) = 2
            imap(3, icase) = 3
         elseif (icase.eq. 2) then
            imap(1, icase) = 2
            imap(2, icase) = 3
            imap(3, icase) = 1
         elseif (icase.eq. 3) then
            imap(1, icase) = 3
            imap(2, icase) = 1
            imap(3, icase) = 2
         endif
      enddo
c     
      do icase = 1, 3
c     
c     permute coordinates and evaluate overlap with ref coordinates '
c     
         do i = 1, natom
            do k = 1, 3
               dummy_coord(imap(k, icase), i) = coord(k,i)
            enddo
         enddo
         overlap(icase) = DDOT(3*natom, dummy_coord,
     $        1, refcoord, 1)  / (x*y)
      enddo
c     
c     find the best match
c     
      call mn_fndmaxd(3, overlap, diff, icase, 1.0d-8)
c     
      do i = 1, natom
         do k = 1, 3
            dummy_coord(imap(k, icase), i) = coord(k,i)
         enddo
      enddo
c     
c     next determine transformation matrix and perform check
c     
      call zero(xyz_tran, 9)
      do i = 1, 3
         xyz_tran(i, imap(i, icase)) = 1.0d0
      enddo
c     
      call DGEMM('T', 'N', 3, natom, 3, 1.0d0, xyz_tran, 3,
     $     coord, 3, 0.0d0, dummy_coord, 3)     
c     
      diff = DDOT(3*natom, dummy_coord, 1, refcoord, 1)
     $     / (x*y)
      if (diff .lt. -0.9d0) then
c     
c     create a sign change of the overall transformation
c     
         write(6,*) ' @xyztran unexpected overlap REF / COORD',
     $        diff
c     
c     the suggested multiplication changes handedness. Not a good idea ...
c     
c     mn            call B_SCAL(9, -1.0d0, xyz_tran, 1)
c     
      endif
c     
      call DGEMM('T', 'N', 3, natom, 3, 1.0d0, xyz_tran, 3,
     $     coord, 3, 0.0d0, dummy_coord, 3)     
c     
      if (print) then
         write(6,*) ' Results from Permutation procedure '
         write(6,*) ' Transformation matrix '
         call output(xyz_tran, 1, 3, 1, 3, 3, 3, 1)
         write(6,*) ' Transformed coordinates '
         call output(dummy_coord, 1, 3, 1, natom, 3, natom, 1)
      endif
c     
      diff = DDOT(3*natom, dummy_coord, 1, refcoord, 1) / (x*y)
c     
c     perform better check
c     
      diff2 = -1.0d0
      do j = 1, natom
         do i = 1, 3
            xx = dummy_coord(i,j) - refcoord(i,j)
            diff2 = max(diff2, abs(xx))
         enddo
      enddo
      least_squares = .false.
c     
      if (abs(diff-1.0d0) .gt. 1.0d-14 .or. diff2 .gt. 1.0d-12) then
         least_squares = .true.
      endif
c
      else
         least_squares = .true.
      endif

      if (least_squares) then
         if (print0) write(6,*) ' @xyztran: use least squares '
c     
c     use a least square procedure to align coordinates
c     
         call DGEMM('N', 'T', 3, 3, natom, 1.0d0, refcoord, 3,
     $        refcoord, 3, 0.0d0, amat, 3)
c     
         call DGEMM('N', 'T', 3, 3, natom, 1.0d0, coord, 3,
     $        refcoord, 3, 0.0d0, wmat, 3)
c     
         if (print) then
            write(6,*) ' Original (Reference) Metric matrix '
            call output(amat, 1, 3, 1, 3, 3, 3, 1)
            write(6,*) ' Original Rotation matrix '
            call output(wmat, 1, 3, 1, 3, 3, 3, 1)
         endif
c     
c     account for zero coordinates 
c     
         do i = 1, 3
            if (abs(amat(i,i)) .lt. 1.0d-12) then
               amat(i,i) = 1.0d0
            endif
         enddo
c     
         if (print) then
            write(6,*) ' Adjusted Metric matrix '
            call output(amat, 1, 3, 1, 3, 3, 3, 1)
         endif
c     
         call mn_invert(amat, amatinv,3, ipiv)
c     
         if (print) then
            write(6,*) ' Inverse Metric matrix '
            call output(amatinv, 1, 3, 1, 3, 3, 3, 1)
         endif
c     
c     Wtran =  wmat * amatinv
c     
         call DGEMM('N', 'N', 3, 3, 3, 1.0d0,Wmat, 3,
     $        amatinv, 3, 0.0d0, xyz_tran, 3)
c     
         if (print) then
            write(6,*) 'xyz_tran before orthonormalization '
            call output(xyz_tran, 1, 3, 1, 3, 3, 3, 1)
         endif
c     
c     account for zero coordinates ... once again
c     
         icount = 0
         do i = 1, 3
            diff = DDOT(3, xyz_tran(1,i), 1, xyz_tran(1,i), 1)
            if (diff .gt. 1.0d-8) then
               icount = icount + 1
               imap(icount,1) = i
               call DCOPY(3, xyz_tran(1, i), 1, wmat(1,icount), 1)
            endif
         enddo
c     
         if (print) write(6,*) ' Number of independent vectors ',
     $        icount

         call mn_lowdin(wmat, amat, 3, icount, 0, scr, mxcor)
c     
         if (print) then
            write(6,*) ' Independent orthogonal vectors '
            call output(wmat, 1, 3, 1, icount, 3, icount, 1)
         endif
c     
         if (icount .eq. 3) then
            do j = 1, icount
               i = imap(j, 1)
               call DCOPY(3, wmat(1,j), 1, xyz_tran(1,i), 1)
            enddo
         elseif (icount .eq. 2) then
c     
c     molecule spans a plane. obtain third vector as orthonormal one
c     Take 1, 1, 1 direction as arbitrary direction to orthogonalize to
c     
            do j = 1, 3
               scr(j) = 1.0d0 / sqrt(3.0d0)
            enddo
            call mn_gschmidt(scr, wmat, 3, icount, amatinv)
            if (print) then
               write(6,*) ' Third vector '
               call output(scr, 1, 3, 1, 1, 3, 1, 1)
            endif
            do j = 1, icount
               i = imap(j, 1)
               
               call DCOPY(3, wmat(1,j), 1, xyz_tran(1,i), 1)
            enddo
            do i = 1, 3
               diff = DDOT(3, xyz_tran(1,i), 1, xyz_tran(1,i), 1)
               if (diff .lt. 1.0d-8) then
                  call DCOPY(3, scr, 1, xyz_tran(1,i), 1)
               endif
            enddo
         elseif (icount .eq. 1) then
c     
c     use original permutation algorithm, which will work for this case.
c     
c     consider 'XYZ', 'ZXY' and 'YZX' definition of axis
c     
            do icase = 1, 3
               if (icase .eq. 1) then
                  imap(1, icase) = 1
                  imap(2, icase) = 2
                  imap(3, icase) = 3
               elseif (icase.eq. 2) then
                  imap(1, icase) = 2
                  imap(2, icase) = 3
                  imap(3, icase) = 1
               elseif (icase.eq. 3) then
                  imap(1, icase) = 3
                  imap(2, icase) = 1
                  imap(3, icase) = 2
               endif
            enddo
c     
            do icase = 1, 3
c     
c     permute coordinates and evaluate overlap with ref coordinates '
c     
               do i = 1, natom
                  do k = 1, 3
                     dummy_coord(imap(k, icase), i) = coord(k,i)
                  enddo
               enddo
               overlap(icase) = DDOT(3*natom, dummy_coord,
     $              1, refcoord, 1)  / (x*y)
            enddo
c     
c     find the best match
c     
            call mn_fndmaxd(3, overlap, diff, icase, 1.0d-8)
c     
            do i = 1, natom
               do k = 1, 3
                  dummy_coord(imap(k, icase), i) = coord(k,i)
               enddo
            enddo
c     
c     next determine transformation matrix and perform check
c     
            call zero(xyz_tran, 9)
            do i = 1, 3
               xyz_tran(i, imap(i, icase)) = 1.0d0
            enddo
c     
            call DGEMM('T', 'N', 3, natom, 3, 1.0d0, xyz_tran, 3,
     $           coord, 3, 0.0d0, dummy_coord, 3)     
c     
            diff = DDOT(3*natom, dummy_coord, 1, refcoord, 1)
     $           / (x*y)
            if (diff .lt. -0.9d0) then
c     
c     create a sign change of the overall transformation
c     
               write(6,*) ' @xyztran unexpected overlap REF / COORD',
     $              diff   
c     mn                  call B_SCAL(9, -1.0d0, xyz_tran, 1)
c     
            endif
c     
         endif                  ! permutation algorithm
c     
c     check that W W^T = unit matrix
c     
         call DGEMM('N', 'T', 3, 3, 3, 1.0d0, xyz_tran, 3,
     $        xyz_tran, 3, 0.0d0, wmat, 3)
c     
         if (print) then
            write(6,*) ' @xyz_tran: check unitarity W W^T'
            call output(wmat, 1, 3, 1, 3, 3, 3, 1)
         endif
c     
         do i = 1, 3
            wmat(i,i) = wmat(i,i) - 1.0d0
         enddo
         call mn_fndmaxd(9, wtran, diff, icase, 1.0d-8)
         if (diff .gt. 1.0d-10) then
            write(6,*) ' @xyztran: Something wrong unitarity '
            write(6,*) ' Current transformation matrix '
            call output(xyz_tran, 1, 3, 1, 3, 3, 3, 1)
            write(6,*) ' Current W W^T - 1 '
            call output(wmat, 1, 3, 1, 3, 3, 3, 1)
            call errex
         endif
c     
         call DGEMM('T', 'N', 3, natom, 3, 1.0d0, xyz_tran, 3,
     $        coord, 3, 0.0d0, dummy_coord, 3)     
c     
         det3 = mn_det3x3(xyz_tran)
         if (abs(det3 -1.0d0) .gt. 0.1d0) then
            write(6,*) ' @xyztran: Determinant deviates from unity ',
     $           det3
            write(6,*) ' try to fix ...'
c     
c     Flip the sign of the x, y or z axis. Use the axis with poorest overlap ...
c     
            call DGEMM('N', 'T', 3, 3, natom, 1.0d0, dummy_coord, 3,
     $           dummy_coord, 3, 0.0d0, amat, 3)
            do i = 1, 3
               if (abs(amat(i,i)) .lt. 1.0d-12) then
                  amat(i,i) = 1.0d0
               endif
            enddo
            call mn_invert(amat, amatinv,3, ipiv)


            call DGEMM('N', 'T', 3, 3, natom, 1.0d0, dummy_coord, 3,
     $           refcoord, 3, 0.0d0, wmat, 3)
c     
            call DGEMM('N', 'N', 3, 3, 3, 1.0d0, amatinv, 3,
     $           wmat, 3, 0.0d0, amat, 3)
c     
            if (print) then
               write(6,*) ' Matrix used for overlap '
               call output(amat, 1, 3, 1, 3, 3, 3, 1)
            endif

            call DCOPY(3, amat, 4, overlap, 1)
            diff2 = overlap(1)
            i0 = 1
            do i = 2, 3
               if (overlap(i) .lt. diff2) then
                  diff2 = overlap(i)
                  i0 = i
               endif
            enddo
c     
            write(6,*) ' Switch sign of poorest overlap ', i0
            do i = 1, 3
               xyz_tran(i,i0) = - xyz_tran(i,i0)
            enddo
            call DGEMM('T', 'N', 3, natom, 3, 1.0d0, xyz_tran, 3,
     $           coord, 3, 0.0d0, dummy_coord, 3)
c     
            det3 = mn_det3x3(xyz_tran)
            if (abs(det3 - 1.0d0) .gt. 0.1d0) then
               write(6,*) ' determinant still wrong ', det3
               write(6,*) ' Fix did not work or some error '
               write(6,*) ' Once more back to xyztran ...'
               call errex
            endif
c     
         endif

         if (print) then
            write(6,*) ' Results from least squares procedure '
            write(6,*) ' Transformation matrix '
            call output(xyz_tran, 1, 3, 1, 3, 3, 3, 1)
            write(6,*) ' Transformed coordinates '
            call output(dummy_coord, 1, 3, 1, natom, 3, natom, 3)
         endif
c     
      endif                     ! permutation vs. least squares '
c     
      det3 = mn_det3x3(xyz_tran)
      if (print0) then
         if (least_squares) then
            write(6,*) ' @xyztran: Results from least squares '
         else
            write(6,*) ' @xyztran: Results from permutation  '
         endif
         write(6,*) ' Transformation matrix '
         call output(xyz_tran, 1, 3, 1, 3, 3, 3, 1)
         write(6,*) ' Determinant transformation ', det3
         write(6,*) ' Reference coordinates '
         call output(refcoord, 1, 3, 1, natom, 3, natom, 3)
         write(6,*) ' Transformed coordinates '
         call output(dummy_coord, 1, 3, 1, natom, 3, natom, 3)
      endif
      diff = DDOT(3*natom, dummy_coord, 1, refcoord, 1) / (x*y)
c     
c     perform another / better check
c     
      diff2 = -1.0d0
      do j = 1, natom
         do i = 1, 3
            xx = dummy_coord(i,j) - refcoord(i,j)
            diff2 = max(diff2, abs(xx))
         enddo
      enddo
c     
      write(6,881) 1.0d0-diff, diff2
 881  format(' @xyz_tran: overlap deviation:', E12.4,
     $     ' Max deviation: ', E12.4)
c     
      if (diff .lt. 0.95 .or. diff2 .gt. 0.2d0 .and. .false.) then
         write(6,*) ' Something suspicious in xyztran '
         if (.not. print0) then
            write(6,*) ' Reference coordinates '
            call output(refcoord, 1, 3, 1, natom, 3, natom, 3)
            write(6,*) ' Results from least squares procedure '
            write(6,*) ' Transformation matrix '
            call output(xyz_tran, 1, 3, 1, 3, 3, 3, 1)
            write(6,*) ' Transformed coordinates '
            call output(dummy_coord, 1, 3, 1, natom, 3, natom, 3)
         endif
         if (print) then
            write(6,*) ' Current input Coordinates '
            call output(coord, 1, 3, 1, natom, 3, natom, 3)
         endif
c     
c     if we are not doing a scan, we do not expect this much difference. Bomb out and analyse!
c     
c     bomb out unless scan_normal or geo_sample calculation
c     
c     if (iflags2(158) .eq. 0 .and. iflags2(160) .eq. 0)
c     $           call mn_errex
      endif
c     
      return
      end
         
         
         
