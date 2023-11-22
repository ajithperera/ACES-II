subroutine definePQ(CISmat,NSize,Nblock,P,Q,projP,projQ,sVecAlgo)
  integer,intent(in)::NSize,Nblock
  logical,intent(in)::sVecAlgo
  double precision,intent(in)::CISmat(NSize,NSize)
  double precision,intent(inout)::P(NSize,Nblock),Q(Nsize,NSize-Nblock)
  double precision, intent(inout)::projP(NSize,NSize),projQ(NSize,NSize)
  double precision::diagCIS(NSize),tmp,cmp
  integer::i,indx,j
  
  integer::MATDIM,resortedArr(NSize)
  logical,dimension(NSize)::mk






  P=0.0d0
  Q=0.0d0
  projP=0.0d0
  projQ=0.0d0
  resortedArr=0
  mk=.true.
  do i=1,NSize
!   if (abs(CISmat(i,i)).lt.0.1d-4) then
!    print*,'modifying CISmat(i,i)'
!    diagCIS(i)=100
!   else
    diagCIS(i)=CISmat(i,i)
!   endif
  enddo

  print*,'DiagCIS:',diagCIS
    
  do i=1,NSize
    resortedArr(i)=minloc(diagCIS,1,mk)
    mk(minloc(diagCIS,mk))=.false.
  enddo
  if (sVecAlgo) then
    P(resortedArr(1),1)=1.0d0
    do i=2,NSize
      Q(resortedArr(i),i-1)=1.0d0
    enddo
!    do i=1,NSize
!      if (i.ne.resortedArr(1)) then
!        Q(i,1)=1.0d0
!      endif
!    enddo
!    tmp=NSize**(1.0d0/2.0d0)
!    tmp=1.0d0/tmp
!    call dscal(NSize,tmp,Q,1)
!    Q=(1.0d0/dsqrt(NSize))*Q
    projP=matmul(P,transpose(P))
    projQ=matmul(Q,transpose(Q))

  else
    do i=1,Nblock
      P(resortedArr(i),i)=1.0d0
    enddo

    do i=NBlock+1,Nsize
      Q(resortedArr(i),i-NBlock)=1.0d0
    enddo

    print*,'lowest element of diagonal:',resortedArr(1)
    print*,'the rest of sorted arr',resortedArr(2:)

!    do i=2,NSize
!      Q(resortedArr(i),i-1)=1.0d0
!    enddo

    MATDIM=NSize
    print*,'p'

  print*,p
  print*
  print*,'q'
  print*,q
  print*
!  call output(P,1,MATDIM,1,MATDIM,MATDIM,MATDIM,1)
!  call output(Q,1,MATDIM,1,MATDIM,MATDIM,MATDIM,1)
    projP=matmul(P,transpose(P))
    projQ=matmul(Q,transpose(Q))
 

    print*,'projP'

  print*
    print*,'projQ' 
!  call output(projP,1,MATDIM,1,MATDIM,MATDIM,MATDIM,1)
!  call output(projQ,1,MATDIM,1,MATDIM,MATDIM,MATDIM,1)

   endif


end subroutine
