	real*8 function shannon12(t1_size,t2_size,t1,t2)
cThis function computes the Shannon index for a joint set of variables.
cAuthor: Dmitry I. Lyakh (05/09/2012).
cINPUT:
c - t1_size - size of the t1 array;
c - t2_size - size of the t2 array;
c - t1(:) - 1st set of real*8 variables;
c - t2(:) - 2nd set of real*8 variables.
cOUTPUT:
c - shannon12 - the Shannon index.
	implicit none
	integer t1_size,t2_size
	real*8 t1(*),t2(*)
	real*8, parameter:: zero_thresh=1d-13
	real*8, parameter:: log10_2=1d0/dlog10(2d0)
	integer i,j,k,l,m,n,k0,k1,k2,k3,ks,kf,ierr
	real*8 vnorm,val
cCalculate the total norm:
	vnorm=0d0
c T1:
	do i=1,t1_size-mod(t1_size,2),2
	 vnorm=vnorm+t1(i)**2+t1(i+1)**2
	enddo
	do i=t1_size-mod(t1_size,2)+1,t1_size
	 vnorm=vnorm+t1(i)**2
	enddo
c T2:
	do i=1,t2_size-mod(t2_size,2),2
	 vnorm=vnorm+t2(i)**2+t2(i+1)**2
	enddo
	do i=t2_size-mod(t2_size,2)+1,t2_size
	 vnorm=vnorm+t2(i)**2
	enddo
cCalculate the Shannon index:
	shannon12=0d0
        Vnorm = Vnorm + 1.0D0
	if(vnorm.gt.zero_thresh) then
	 vnorm=1d0/vnorm 
c T1:
	 do i=1,t1_size
	  if(dabs(t1(i)).gt.zero_thresh) then
	   val=(t1(i)**2)*vnorm
	   shannon12=shannon12-val*dlog10(val)*log10_2
	  endif
	 enddo
c T2:
	 do i=1,t2_size
	  if(dabs(t2(i)).gt.zero_thresh) then
	   val=(t2(i)**2)*vnorm
	   shannon12=shannon12-val*dlog10(val)*log10_2
	  endif
	 enddo	 
	endif
	return	
	end
