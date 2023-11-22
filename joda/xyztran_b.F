      subroutine xyztran_b(natom, coord, refcoord,
     $     dummy_coord, masses, refmasses,
     $     xyz_tran, scr, mxcor)
c     
c     In this subroutine we determine a transformation that
c     amounts to a permutation of the axis such that 
c     the current geometry (in COORD) lines up (best) with the coordinates
c     in the reference geometry.
c     This subroutine is hooked up from a more general subroutine in MRCC
c     It is used when NOREORI=ON and we wish to keep the axes same as some reference geometry.
c     Perhaps this should be used only when doing a 'findif' or prep_vib calculation.
c     
      implicit none

#include <machsp.com>
#include <flags.com>      
#include <flags2.com>      

      integer natom, nsize, pass1, i, icase, i0, 
     $     imap(3,3), k, ipiv(3), mxcor, j, icount
      double precision
     $     coord(3, natom), refcoord(3, natom),
     $     dummy_coord(3, natom), masses(natom),
     $     refmasses(natom), xyz_tran(3,3), overlap(3),
     $     XDOT, diff, x, y, Amat(3,3), Amatinv(3,3),
     $     wmat(3,3), wtran(3,3), scr(mxcor), diff2, xx, det3,
     $     det3x3
      character*3 xyzcase(3)
      logical ok, print, print0, least_squares

c     
      CALL GETREC(-1,'JOBARC','PASS1',1,PASS1)
      CALL GETREC(1,'JOBARC', 'NATOMS',1,NATOM)
      nsize = 3 * natom
      print = .true.
      print0 = .true.
c     
      xyzcase(1) = 'XYZ'
      xyzcase(2) = 'ZXY'
      xyzcase(3) = 'YZX'
c     
      if (pass1 .eq. -1) then
c     
c     the initialization of the reference geometry
c     
         write(6,*) ' Retrieve initial geometry from Jobarc '
         CALL GETREC(1,'JOBARC', 'COORD',NSIZE*IINTFP,coord)
         write(6,*) ' Initial coordinates '
         call output(coord, 1, 3, 1, natom, 3, natom, 1)
         CALL PUTREC(1,'JOBARC', 'REFCOORD',NSIZE*IINTFP,coord)
         CALL GETREC(1,'JOBARC', 'ATOMMASS',NATOM*IINTFP,masses)
         CALL PUTREC(1,'JOBARC', 'REFMASS',NATOM*IINTFP,masses)
         call dzero(xyz_tran, 9)
         do i = 1, 3
            xyz_tran(i, i) = 1.0d0
         enddo
         call PUTREC(1,'JOBARC', 'XYZ_TRAN', 9*IINTFP, xyz_tran)
      else
         CALL GETREC(1,'JOBARC', 'COORD',NSIZE*IINTFP,coord)
c     
c     When comparing with REFCOORD, we rotate XYZ to be in best agreement with 
c     the reference geometry.
c
         CALL GETREC(1,'JOBARC', 'REFCOORD',NSIZE*IINTFP,refcoord)
         CALL GETREC(1,'JOBARC', 'ATOMMASS',NATOM*IINTFP,masses)
         CALL GETREC(1,'JOBARC','REFMASS',NATOM*IINTFP,refmasses)
c     
c     first check that masses and refmasses are equal
c     
         ok = .true.
         do i = 1, natom
            ok = ok .and.
     $           (abs(masses(i) - refmasses(i)) .lt. 1.0d-6)
         enddo
         if (.not. ok) then
            write(6,*) ' @xyztran: unexpected reshuffle of atoms '
            write(6,*) ' Reference masses '
            call output(refmasses, 1, 1, 1, natom, 1, natom, 1)
            write(6,*) ' Current masses '
            call output(masses, 1, 1, 1, natom, 1, natom, 1)
            call errex
         endif
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
         x = SQRT(XDOT(3*natom, coord, 1, coord, 1))
         y = SQRT(XDOT(3*natom, refcoord, 1, refcoord, 1))
c     
c     try to use permutation algorithm first
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
            overlap(icase) = XDOT(3*natom, dummy_coord,
     $           1, refcoord, 1)  / (x*y)
         enddo
c     
c     find the best match
c     
         call fndmaxd(3, overlap, diff, icase, 1.0d-8)
c     
         do i = 1, natom
            do k = 1, 3
               dummy_coord(imap(k, icase), i) = coord(k,i)
            enddo
         enddo
c     
c     next determine transformation matrix and perform check
c     
         call dzero(xyz_tran, 9)
         do i = 1, 3
            xyz_tran(i, imap(i, icase)) = 1.0d0
         enddo
c     
         call XGEMM('T', 'N', 3, natom, 3, 1.0d0, xyz_tran, 3,
     $        coord, 3, 0.0d0, dummy_coord, 3)     
c     
         diff = XDOT(3*natom, dummy_coord, 1, refcoord, 1)
     $        / (x*y)
         if (diff .lt. -0.9d0) then
c     
c     create a sign change of the overall transformation
c     
            write(6,*) ' @xyztran unexpected overlap REF / COORD',
     $           diff
c     
         endif
c     
         call XGEMM('T', 'N', 3, natom, 3, 1.0d0, xyz_tran, 3,
     $        coord, 3, 0.0d0, dummy_coord, 3)     
c     
         if (print) then
            write(6,*) ' Results from Permutation procedure '
            write(6,*) ' Transformation matrix '
            call output(xyz_tran, 1, 3, 1, 3, 3, 3, 1)
            write(6,*) ' Transformed coordinates '
            call output(dummy_coord, 1, 3, 1, natom, 3, natom, 3)
         endif
c     
         diff = XDOT(3*natom, dummy_coord, 1, refcoord, 1) / (x*y)
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

         det3 = det3x3(xyz_tran)
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
         diff = XDOT(3*natom, dummy_coord, 1, refcoord, 1) / (x*y)
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
 881     format(' @xyz_tran: overlap deviation:', E12.4,
     $        ' Max deviation: ', E12.4)
c     
         if (diff .lt. 0.95 .or. diff2 .gt. 0.2d0) then
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
            write(6,*) ' Current input Coordinates '
            call output(coord, 1, 3, 1, natom, 3, natom, 3)
c     
         endif
c     

         call PUTREC(1,'JOBARC','XYZ_TRAN', 9*IINTFP, xyz_tran)
c     
      endif                     ! pass1
c     
      return
      end
         
         
         
