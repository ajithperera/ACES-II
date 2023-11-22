      subroutine sortEvals(evals,esize,iintfp,beta)
      integer esize,iintfp,beta
      double precision evals(esize)
      integer place,i
      logical, dimension(esize):: mk
      integer resortedArr(esize)
      double precision HFevec(esize,esize),tmp(esize),val
      double precision resortedEvecs(esize,esize)
      double precision resortedEvals(esize)
      print*,'INSIDE RESORT.F -- resorting eVals'
      tmp=evals
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

      resortedEvals=0.0d0
      do i=1,esize
        place=resortedArr(i)
        resortedEvals(i)=evals(place)
        print*,'MO',i,'now has energy',evals(place)
      enddo

      if (beta.eq.1) then
        call putrec(20,'JOBARC','SCFEVLA0',esize*iintfp,
     & resortedEvals)
      else
        call putrec(20,'JOBARC','SCFEVLB0',esize*iintfp,
     & resortedEvals)
      endif

      endsubroutine 
