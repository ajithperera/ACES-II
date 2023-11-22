      subroutine droplist(nocc,nvir,nbastot,trans,occthresh,virthresh,
     & vec,droparray)
      implicit none

      integer ppp, qqq, qqqmax
 
      integer nocc, nvir, nbastot, ndropocc, ndropvir

      integer droparray(nbastot)
 
      double precision occthresh, virthresh

      double precision trans(nbastot,nbastot), vec(nbastot)

      double precision transmax, transtemp

      if (.false.) then
      open(900,file='dropvec')
      do qqq = 1, nbastot
         transmax = 0.0D0
         do ppp = 1, nbastot
            if (abs(trans(ppp,qqq)).ge.transmax) then
               transmax = abs(trans(ppp,qqq))
            end if
         end do
         if ((qqq.le.nocc).and.(transmax.ge.occthresh)) then
            write(*,10) qqq
            write(900,20) qqq
         end if
         if ((qqq.gt.nocc).and.(transmax.ge.virthresh)) then
            write(*,10) qqq
            write(900,20) qqq
         end if
      end do
      close(900)
      end if

      do qqq = 1, nbastot
         transmax = 0.0D0
         do ppp = 1, nbastot
            if (abs(trans(ppp,qqq)).ge.transmax) then
               transmax = abs(trans(ppp,qqq))
            end if
         end do
         vec(qqq) = transmax
      end do

      ndropocc = 0
      ndropvir = 0

      do ppp = 1, nbastot

      transmax = 0.0D0

      do qqq = 1, nbastot
         transtemp = vec(qqq)
         if (transtemp.ge.transmax) then
            transmax = transtemp
            qqqmax = qqq
         end if
      end do

      vec(qqqmax) = -1.0D0

      if (qqqmax.le.nocc) then
         if (dint(((dble(ndropocc)/dble(nocc))*100.0D0)).le.
     & occthresh) then
            ndropocc = ndropocc + 1
            droparray(qqqmax) = qqqmax
         end if
      end if

      if (qqqmax.gt.nocc) then
         if (dint(((dble(ndropvir)/dble(nvir))*100.0D0)).le.
     & virthresh) then
            ndropvir = ndropvir + 1
            droparray(qqqmax) = qqqmax
         end if
      end if

      end do

      open(900,file='dropvec')
      do ppp = 1, nbastot
         if (droparray(ppp).gt.0) then
            write(*,10) droparray(ppp)
            write(900,20) droparray(ppp)
         end if
      end do
      close(900)

 10   format('Drop orbital ',I6,' for NLS-EOM')
 20   format(I6)

      return
      end
 
