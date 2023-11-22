      subroutine finite_field()
c     
c     Calculate numerical derivatives of electric moments.
c     This allows numerical evaluation of static multipole polarizabilities
c     dipole, quadrupole, octupole
c     
      implicit none
c     
#include "aces2.h"
#include "flags.h"
#include "machsp.com"
c     
      integer maxprop
      parameter (maxprop = 19)  ! 19 = 3 + 6 + 10
      integer maxscr, ifield
      parameter (maxscr = 10 000 000)
      double precision scr(maxscr), ddot
c     
      character*8 PERTSTR,STRING
      integer nopert, ialpha_order, ipert, jpert, nao2, icase,
     $     iunit, nao
      double precision ffield, xprop(maxprop,maxprop), 
     $     xplus(maxprop), xminus(maxprop), x,
     $     xdummy(maxprop), xprop0(maxprop), xnorm
c     
      dimension pertstr(maxprop)
      character*3 Label(19)
      data label /'  X', '  Y', '  Z',
     $     ' XX', ' YY', ' ZZ', ' XY', ' XZ', ' YZ',
     $     'XXX', 'YYY', 'ZZZ',
     $     'XXY', 'XXZ', 'XYY', 'YYZ', 'XZZ', 'YZZ',
     $     'XYZ' /
      integer i, j, ilow(3), ihigh(3), ib, jb, i1, i2
c     
      if (.not.direct) call runit(integral_package)
      call runit('xvmol2ja')
      call runit('xvprops')
c     
      ialpha_order = iflags2(h_IFLAGS2_order_alpha)
      ffield = iflags2(h_IFLAGS2_ffield) * 10.0d-6
      ifield = iflags2(h_IFLAGS2_ffield)
      call a2getrec(20, 'JOBARC', 'NBASTOT ', 1, nao)
      nao2 = (nao*(nao+1)) / 2
c     
      i1 = 1
      i2 = 1 + nao2
c     call check_multipole(nao, nao2, scr(i1), scr(i2))
c     
      write(6,*) ' Number of property integrals ', nao2
c     
      if (ialpha_order .eq. 0) then
         NOPERT=3
         PERTSTR(1)='DIPOLE_X'
         PERTSTR(2)='DIPOLE_Y'
         PERTSTR(3)='DIPOLE_Z'
      elseif (ialpha_order .eq. 1) then
         NOPERT=9
         PERTSTR(1)='DIPOLE_X'
         PERTSTR(2)='DIPOLE_Y'
         PERTSTR(3)='DIPOLE_Z'
c     
         PERTSTR(4)='QUAD_XX '
         PERTSTR(5)='QUAD_YY'
         PERTSTR(6)='QUAD_ZZ '
         PERTSTR(7)='QUAD_XY '
         PERTSTR(8)='QUAD_XZ '
         PERTSTR(9)='QUAD_YZ '
c     PERTSTR(9)='2NDMO_YZ'
      elseif (ialpha_order .eq. 2) then
         NOPERT=19
         PERTSTR(1)='DIPOLE_X'
         PERTSTR(2)='DIPOLE_Y'
         PERTSTR(3)='DIPOLE_Z'
c     
         PERTSTR(4)='QUAD_XX '
         PERTSTR(5)='QUAD_YY'
         PERTSTR(6)='QUAD_ZZ '
         PERTSTR(7)='QUAD_XY '
         PERTSTR(8)='QUAD_XZ '
         PERTSTR(9)='QUAD_YZ '
c     
         PERTSTR(10)='OCTUPXXX'
         PERTSTR(11)='OCTUPYYY'
         PERTSTR(12)='OCTUPZZZ'
         PERTSTR(13)='OCTUPXXY'
         PERTSTR(14)='OCTUPXXZ'
         PERTSTR(15)='OCTUPXYY'
         PERTSTR(16)='OCTUPYYZ'
         PERTSTR(17)='OCTUPXZZ'
         PERTSTR(18)='OCTUPYZZ'
         PERTSTR(19)='OCTUPXYZ'
      endif
c     
      do ipert = 1, nopert
         do jpert = 1, nopert
            xprop(ipert, jpert) = 0.0d0
         enddo
      enddo
c     
      do i = 1, maxprop
         xminus(i) = 0.0d0
      enddo
      call a2putrec(20,'JOBARC','MULTMOM ',
     $     maxprop*iintfp, xminus)
      call a2putrec(20,'JOBARC','MULTMOM2',
     $     maxprop*iintfp, xminus)
c     
c     calculate zerofield static properties first
c     
      call DSCAL(nao2, 0.0d0, scr, 1)
      call a2putrec(20,'JOBARC','FFIELD  ',
     $     nao2*iintfp, scr)
c     
      write(6,*) ' @finite_field: run zero field reference '
c     
c     run single point calculation including evaluation of correlated properties 
c     
      if (hf_scf) then
         if (nddo_guess) call runit('xnddo')
         If (dkh .and. contract) Call Runit("xdkh")
         call runit('xvscf')
      else
         If (dkh .and. contract) Call Runit("xdkh")
         call runit('xvscf')
         call runit('xintgrt')
      endif
c     
      if (fno) then
         call runfno
      end if
      Call Runit('xvtran')
      call runit('xintprc')
      Call Runit('xvcc')
      Call Post_scf_anlytcl_grad(.false.)
      Call Runit("xprops")
c     
c     Static electric multipoles are calculated
c     
      call a2getrec(20,'JOBARC','MULTMOM ',
     $     maxprop*iintfp, xprop0)
c
c ifield = -1: only calculate static multipole moments 
c
      if (ifield .eq. -1) go to 900
c     
c     now calculate static polarizabilities
c     
      do ipert = 1, nopert
         do icase = -1, 1, 2
            string = pertstr(ipert)
            call a2getrec(20,'JOBARC',string, nao2*iintfp, scr)
c     
c     adjust ffield based on norm perturbation
c     
            xnorm = DDOT(nao2, scr, 1, scr, 1)
            if (xnorm .gt. 1.0d0) then
               xnorm = 1.0d0 / sqrt(xnorm)
            else
               xnorm = 1.0d0
            endif
            call DSCAL(nao2, icase*ffield*xnorm, scr, 1)
            call a2putrec(20,'JOBARC','FFIELD  ',
     $           nao2*iintfp, scr)
c     
            write(6,700) ipert, icase*ffield*xnorm,
     $           label(ipert), nao2
 700        format(' Include perturbation ', i4, E14.4, 3x, A3, i8)
c     
c     run single point calculation including evaluation of correlated properties 
c     
            if (hf_scf) then
               if (nddo_guess) call runit('xnddo')  
               If (dkh .and. contract) Call Runit("xdkh")
               call runit('xvscf')
            else 
               If (dkh .and. contract) Call Runit("xdkh")
               call runit('xvscf')
               call runit('xintgrt')
            endif
c     
            if (fno) then
               call runfno
            end if
            Call Runit('xvtran')
            call runit('xintprc')
            Call Runit('xvcc')
c     Call Runit("xlambda")     
            Call Post_scf_anlytcl_grad(.false.)
            Call Runit("xprops")
c     
c     properties are calculated
c     
            if (icase .eq. -1) then
               call a2getrec(20,'JOBARC','MULTMOM2',
     $              maxprop*iintfp, xminus)
               if (.false.) then
c     
c     the difference between multmom and multmom2 is nuclear contribution, which cancels in finite difference.
c     
                  call a2getrec(20,'JOBARC','MULTMOM2',
     $                 maxprop*iintfp, xdummy)
                  call DAXPY(maxprop, -1.0d0, xminus, 1, xdummy, 1)
                  x = ddot(maxprop, xdummy, 1, xdummy,1)
                  if (x .gt. 1.0d-12) then
                     write(6,*) ' Problem case ', ipert, icase, x
                     write(6,*) ' MULTMOM '
                     call output(xminus, 1, 1, 1, maxprop,
     $                    1, maxprop, 1)
                     call a2getrec(20,'JOBARC','MULTMOM2',
     $                    maxprop*iintfp, xdummy)
                     write(6,*) ' MULTMOM2'
                     call output(xdummy, 1, 1, 1, maxprop,
     $                    1, maxprop, 1)
                  endif
               endif
            elseif (icase .eq. 1) then
               call a2getrec(20,'JOBARC','MULTMOM2',
     $              maxprop*iintfp, xplus)
               if (.false.) then
                  call a2getrec(20,'JOBARC','MULTMOM2',
     $                 maxprop*iintfp, xdummy)
                  call DAXPY(maxprop, -1.0d0, xplus, 1, xdummy, 1)
                  x = ddot(maxprop, xdummy, 1, xdummy,1)
                  if (x .gt. 1.0d-12) then
                     write(6,*) ' Problem case ', ipert, icase, x
                  endif
               endif
            endif

         enddo                  ! icase = +/-
c     
         do jpert = 1, nopert
            xprop(ipert, jpert) = (xplus(jpert) - xminus(jpert) )
     $           / (2.0d0 * ffield * xnorm)
         enddo
c     
      enddo                     ! ipert
c     
      if (.false.) then
         do ipert = 1, nopert
            do jpert = 1, ipert-1
               x = xprop(ipert, jpert) + xprop(jpert, ipert)
               x = x /2.0d0
               xprop(ipert,jpert) = x
               xprop(jpert,ipert) = x
            enddo
         enddo
      endif
c
 900  continue
c     
c     print a nice output file following xvcceh  write_alpha.f
c     
      iunit = 80
      open(unit=iunit, file='alpha_static', status = 'UNKNOWN')
      rewind(iunit)
c     
      write(iunit,*) ' Static Electric multipoles '
      write(iunit,*) ' n_moments ', nopert
      do i = 1, nopert
         write(iunit, 705) i, label(i), xprop0(i)
      enddo
 705  format(3x,I4, 3x, '<<', A4, ' >>', 5x, E16.8)
      write(iunit,*)
c
      if (ifield .eq. -1) then
         write(iunit,701) ifield
         go to 9010
      endif
c
      write(iunit,701)  ialpha_order + 1
 701  format(' Static multipole polarizabilies up to order: ', i4)
      write(iunit, *)
c     
c     now print coefficients in tabular format
c     
      ilow(1) = 1
      ihigh(1) = 3
      ilow(2) = 4
      ihigh(2) = 9
      ilow(3) = 10
      ihigh(3) = 19
c     
c     write diagonal blocks, (ialpha+1) indicates highest multipole considered.
c     
      do ib = 1, ialpha_order+1
         jb = ib
         do i = ilow(ib), ihigh(ib)
            do j = ilow(ib), i
               write(iunit,800) i, j, label(i), label(j), xprop(i,j),
     $              xprop(j,i), abs(xprop(i,j)-xprop(j,i))
            enddo
         enddo
         write(iunit,*)
      enddo
c     
c     
c     write off-diagonal blocks
c     
      do ib = 1, ialpha_order + 1
         do jb = 1, ib-1
            do i = ilow(ib), ihigh(ib)
               do j = ilow(jb), ihigh(jb)
                  write(iunit,800) i, j, label(i), label(j), xprop(i,j),
     $                 xprop(j,i), abs(xprop(i,j)-xprop(j,i))
               enddo
            enddo
            write(iunit,*)
         enddo
      enddo
c     
 800  format(3x,2I4, 3x, '<<', A3, ';', A3, '>>', 5x, E16.8, 8x,
     $     E16.8, 4x, E12.4)
c
 9010 continue
      close(iunit,status='KEEP')
c
      return
      end
