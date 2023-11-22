      subroutine write_alpha(iunit0, ialpha, xfreq, imag_freq, nfreq,
     $     nopert, xprop)
c     
c     this subroutine writes multipole polarizabilities in convenient table.
c     
      implicit none
c     
      integer iunit, ialpha, nfreq, nopert, iunit0
      double precision xfreq(nfreq), xprop(nopert, nopert, nfreq)
      integer imag_freq(nfreq)
      character*3 Label(19)
      data label /'  X', '  Y', '  Z',
     $     ' XX', ' YY', ' ZZ', ' XY', ' XZ', ' YZ',
     $     'XXX', 'YYY', 'ZZZ',
     $     'XXY', 'XXZ', 'XYY', 'YYZ', 'XZZ', 'YZZ',
     $     'XYZ' /
      integer i, j, ifreq, icount, ilow(3), ihigh(3), ib, jb, n
      double precision freq, xvec(100), x
C
      iunit = iunit0
      rewind(iunit)
      n = (nopert * (nopert+1) ) / 2
c     
c     read in data from original file
c
      read(iunit,*)
      do ifreq = 1, nfreq
         read(iunit,*) freq
         read(iunit,*) (xvec(j), j=1,n)
c     
         icount = 0
         do i = 1, nopert
            do j = 1, i
               icount = icount + 1
               x = xvec(icount)
               xprop(i,j,ifreq) = x
               xprop(j,i,ifreq) = x
            enddo
         enddo
c     
      enddo
c     
      rewind(iunit)
c
c      iunit=6
c
      write(iunit,702) nfreq
 702  format(' Numerical Frequencies: ', i4)
      do i = 1, nfreq
         write(iunit,700) xfreq(i), imag_freq(i)
      enddo
 700  format(F18.10, I4)
c     
      write(iunit,*)
      write(iunit,701)  ialpha + 1
 701  format(' Multipole polarizabilies up to order: ', i4)
      write(iunit, *)
c     
c     now print coefficients in nice tabular format
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
      do ib = 1, ialpha+1
         jb = ib
         do i = ilow(ib), ihigh(ib)
            do j = ilow(ib), i
               write(iunit,800) i, j, label(i), label(j)
               do ifreq = 1, nfreq
                  xvec(ifreq) = xprop(i,j,ifreq)
               enddo
               write(iunit,801) (xvec(ifreq), ifreq=1,nfreq)
            enddo
         enddo
         write(iunit,*)
      enddo
c     
c     write off-diagonal blocks
c     
      do ib = 1, ialpha+1
         do jb = 1, ib-1
            do i = ilow(ib), ihigh(ib)
               do j = ilow(jb), ihigh(jb)
                  write(iunit,800) i, j, label(i), label(j)
                  do ifreq = 1, nfreq
                     xvec(ifreq) = xprop(i,j,ifreq)
                  enddo
                  write(iunit,801) (xvec(ifreq), ifreq=1,nfreq)
               enddo
            enddo
            write(iunit,*)
         enddo
      enddo
c     
 800  format(3x,2I4, 3x, '<<', A3, ';', A3, '>>')
 801  format(100E16.8)
c     
c      close(iunit, status='KEEP')
c     
      return
      end
