subroutine extrapBlockSD(revR,oldR,Tvecs,thetaOld, &
                     Nblocks,NSize,itr,vecNum)
  integer, intent(in)::Nblocks,NSize,itr,vecNum
  double precision,intent(inout)::revR(NSize,itr*Nblocks)
  double precision, intent(inout)::oldR(NSize,itr*Nblocks)
  double precision,intent(inout)::Tvecs(Nblocks*itr,Nblocks*itr)
  double precision,intent(in)::thetaOld(Nblocks*itr)

  integer::counter
  double precision::diff(Nblocks)
  double precision::HbarDiag(Nsize),r1(Nsize)
  double precision::SCR(Nsize,itr*Nblocks)
  double precision::factor
  double precision::revRSCR(Nsize,Nblocks)
! calculate (Hbar*C - eig*C)*Tvec
  print*,'Inside ExtrapBlock'
  counter=Nsize*(Nblocks*itr)
  HbarDiag=0.0d0
  call getlst(HbarDiag,1,1,1,1,472)
call checksum("HbarDiag:",HbarDiag,NSize,s)
  revRSCR=0.0d0
  r1=0.0d0
  do j=1,Nblocks
    print*

    SCR=revR
    call daxpy(counter,-thetaOld(j),oldR,1,SCR,1)
    call checksum("Tvecs:",Tvecs(:,j),itr*Nblocks,s)
    r1=matmul(SCR,Tvecs(:,j))! r1
    call checksum("r1:",r1,NSize,s)
    print*,thetaOld(j)-HbarDiag(vecNum)
    factor=1.0d0/(thetaOld(j)-HbarDiag(vecNum))
    !print*,"HbarDiag:", HbarDiag(j)
    print*,'Factor:',thetaOld(j)-HbarDiag(vecNum)

    call dscal(Nsize,factor,r1,1)
    call checksum("q1:",r1,NSize,s)
    print*,'Deposing on rec',Nblocks*itr+j
    call putlst(r1,Nblocks*itr+j,1,1,1,497)

    print*
  enddo


end subroutine
