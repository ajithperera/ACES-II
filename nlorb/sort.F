      subroutine resort(evals,esize,iintfp,beta)
      integer esize,iintfp,beta
      double precision evals(esize,esize)
      integer place,i,j
      logical, dimension(esize):: mk
      integer resortedArr(esize)
      double precision HFevec(esize,esize),tmp(esize),val
      double precision resortedEvecs(esize,esize)
      double precision resortedEvals(esize)
      print*,'INSIDE RESORT.F'
      tmp=0.0d0
      do i=1,esize
      tmp(i)=evals(i,i)
      enddo

      mk=.true.
      resortedArr=0
      do i=1,esize
      val=minval(tmp,mk)
      resortedArr(i)=minloc(tmp,1,mk)!minval(tmp,mk)!minloc(tmp,val)

      print*, 'MO',i,'becomes NLMO',val
      print*, 'MO',i,'resorts to',resortedArr(i)
      print*
      mk(minloc(tmp,mk))=.False.
      enddo


      HFevec=0.0d0
      resortedEvecs=0.0d0
      if (beta.eq.1) then
        call getrec(20,'JOBARC','SCFEVCA0',esize*esize*iintfp,
     & HFevec) 
      else
         call getrec(20,'JOBARC','SCFEVCB0',esize*esize*iintfp,
     & HFevec)
      endif

      resortedEvals=0.0d0
      do i=1,esize
        place=resortedArr(i)
        do j=1,esize
         resortedEvecs(j,i)=HFevec(j,place)
        enddo
      enddo

      call prnmat(" ZWW NLMO Evecs b4 ", HFevec, esize,esize,.true.)
      call prnmat(" ZWW NLMO Evecs after ", resortedEvecs, esize,
     &                esize,.true.)

      if (beta.eq.1) then
        call putrec(20,'JOBARC','SCFEVCA0',esize*esize*iintfp,
     & resortedEvecs)
      else
        call putrec(20,'JOBARC','SCFEVCB0',esize*esize*iintfp,
     & resortedEvecs)
      endif

      endsubroutine 
