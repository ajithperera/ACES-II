
      subroutine conv_set(icore,maxcor,uhf,ncycle)
      implicit none
C Common blocks
      integer iflags(100)
      common/flags/iflags
      integer iintln,ifltln,iintfp
      common/machsp/iintln,ifltln,iintfp
      integer nstart,nirrep
      common/syminf/nstart,nirrep
C Input variables
      integer maxcor,uhf
C Input/output variable
      integer ncycle
C Pre-allocated Local variables
      integer icore(maxcor)
C Local variables
      integer maxiter,ii,i2,i3,nmo,totlen,mxirr2,ndrop


      call getrec(20,'JOBARC','OOMBPTCY',1,maxiter)
      write(6,*) '@OOMBPT2: The OO-MBPT(2) equations converged'
      call putrec(20,'JOBARC','OOMBPTST',1,1)

      call putrec(20,'JOBARC','OOMBPTIT',1,ncycle+1)

      call getrec(20,'JOBARC','IFLAGS',100,iflags)
C Reset parameters for running first natural determinant
      iflags(16)=0
      if (iflags(11).gt.0) iflags(11)=1
      iflags(38)=1
      iflags(39)=1
      iflags(34)=0
      iflags(77)=0
      call putrec(20,'JOBARC','IFLAGS',100,iflags)
      call form_dropmo_length(nmo,totlen,mxirr2)
      call getrec(20,'JOBARC','NUMDROPA',1,ndrop)
C Reset some vector lengths
      if (ndrop .gt. 0) then
        call getrec(20,'JOBARC','MODROPA0',ndrop,icore)
        call putrec(20,'JOBARC','MODROPA ',ndrop,icore)
        call getrec(20,'JOBARC','REORDEA0',nmo,icore)
        call putrec(20,'JOBARC','REORDERA',nmo,icore)
        if (uhf .gt. 0) then
          call getrec(20,'JOBARC','MODROPB0',ndrop,icore)
          call putrec(20,'JOBARC','MODROPB ',ndrop,icore)
          call getrec(20,'JOBARC','REORDEB0',nmo,icore)
          call putrec(20,'JOBARC','REORDERB',nmo,icore)
        endif
      else
         call getrec(20,'JOBARC','OCCUPYA ',nirrep,icore)
         call putrec(20,'JOBARC','OCCUPYA0',nirrep,icore)
         if (uhf .gt. 0) then
           call getrec(20,'JOBARC','OCCUPYB ',nirrep,icore)
           call putrec(20,'JOBARC','OCCUPYB0',nirrep,icore)
         endif
      endif

      call getrec(20,'JOBARC','OOMBITER',(maxiter+1)*iintfp,icore)
      i2 = 1 + (maxiter+1)*iintfp
      call getrec(20,'JOBARC','OOMBITR2',(maxiter+1)*iintfp,icore(i2))
      if (ncycle .eq. 0) ncycle = ncycle+1
      call print_cycle(icore,icore(i2),ncycle)
      return
      end

      subroutine print_cycle(iterations,it2,cycles)
      implicit none
      integer cycles,unit
      double precision iterations(cycles),it2(cycles)
      integer ii
      write(6,*) '  OO-MBPT(2) Iterations  '
      write(6,*) ' Cycle            norm(Via)           norm(Fia)     '
      write(6,*) '--------      ------------------  ------------------'
      do ii = 1,cycles+1
        write(6,1000) ii-1,iterations(ii),it2(ii)
      end do

      unit=10
      open(file='FIAOUT',unit=10,status='REPLACE')
      do ii = 1,cycles+1
        write(10,1001) ii,iterations(ii),it2(ii)
      end do
      close(10)

      return
 1000 format(T3,I6,8X,ES16.10,4X,ES16.10)
 1001 format(I6,',',ES16.10,',',ES16.10)
      end
