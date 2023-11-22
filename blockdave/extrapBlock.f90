subroutine extrapBlock(revR,oldR,Tvecs,thetaOld,TOL, &
       Nblocks,NSize,itr,totVecs,nroots,rootIndx,unconvergedRoots)
  integer, intent(in)::Nblocks,NSize,itr
  integer,intent(inout)::totVecs
  integer,intent(inout)::nroots
  double precision, intent(in)::TOL
  double precision,intent(inout)::revR(NSize,totVecs)
  double precision, intent(inout)::oldR(NSize,totVecs)
  double precision,intent(inout)::Tvecs(totVecs,totVecs)
  double precision,intent(in)::thetaOld(Nblocks)
  integer,intent(inout)::rootIndx(NBlocks)
  logical,intent(inout)::unconvergedRoots(Nblocks)


  integer:: zz,bkuprootIndx(Nblocks)
  integer::vecCount,counter,tmproot,j
  double precision::diff(Nblocks)
  double precision::HbarDiag(Nsize),r1(Nsize)
  double precision::SCR(Nsize,totVecs)
  double precision::ratio,factor
  double precision::r1norm,deltanorm,deltapnorm
! revR - Hbar*C
! oldR - C
! r1   - residual for that iteration 
! r1norm - convergence criteria for a particular vector
! ***PURPOSE: Iterate over all remaining unconverged roots in the block
!             to determine the next iteration's guess vector.
! ***GOAL:  calculates next iterations guess vector r1=(Hbar*C - eig*C)*Tvec
!           ,adds it to the subspace, & orthogonalizes it iff ||r1||_2 < TOL
!
! ***CAVEAT: Once a vector r1 is converged, we stop expanding the subspace using
!            that r1.
  print*,'Inside ExtrapBlock',nroots
  bkuprootIndx=rootIndx
  counter=Nsize*(totVecs)
  HbarDiag=0.0d0
  call getlst(HbarDiag,1,1,1,1,472)
!call checksum("HbarDiag:",HbarDiag,NSize,s)
  r1=0.0d0
  tmproot=nroots
  vecCount=0
! source:
! https://riptutorial.com/fortran/example/8789/using-pack-to-select-elements-meeting-a-condition
  rootIndx=PACK([(zz,zz=1,Nblocks)],unconvergedRoots)
  print*,'** Unconverged roots prior to iteration:'
  print*,rootIndx(1:nroots)
  do j=1,nroots
    print*

    SCR=revR
    call daxpy(counter,-thetaOld(rootIndx(j)),oldR,1,SCR,1)!Hbar*C - eig*C
    !call checksum("Tvecs:",Tvecs(:,j),totVecs,s)

    r1=matmul(SCR,Tvecs(:,rootIndx(j)))! r1
    r1norm=dot_product(r1,r1)
    r1norm=sqrt(r1norm)
    print*,"Norm of the residual r1 for vector:",rootIndx(j),r1norm
    if ((r1norm.gt.TOL).or.(itr.eq.1)) then ! aka vector not converged
      vecCount=vecCount+1
!      factor=1.0d0/(thetaOld(rootIndx(j))-HbarDiag(rootIndx(j)))
!      print*,'Factor:',thetaOld(rootIndx(j))-HbarDiag(rootIndx(j))

!      call dscal(Nsize,factor,r1,1) ! step 8
! ** copied logic from nextdav.f to normalize r1 by (Hbar-theta*I)^-1
      call vecdiv2(thetaOld(rootIndx(j)),r1,HbarDiag,r1,NSize) 
      deltanorm=dot_product(r1,r1)
      deltanorm=dsqrt(deltanorm)
      call dscal(Nsize,1.0d0/deltanorm,r1,1) !step9
#ifdef _DEBUG_LVLM
      call checksum("q1:",r1,NSize,s)
      print*,'Deposing on rec',totVecs+vecCount
#endif
      call putlst(r1,totVecs+vecCount,1,1,1,497)
      call blockQR(oldR,Iter_count,Nblocks,&
                    Nsize,IRREPX,totVecs+vecCount)!step 10

      call getlst(r1,totVecs+vecCount,1,1,1,497)
      deltapnorm=dot_product(r1,r1)
      print*,'norm of deltapnorm',deltapnorm
      ratio=deltapnorm/deltanorm
      if (ratio.gt.0.001) then
        call dscal(Nsize,1.0d0/sqrt(deltapnorm),r1,1)
        call putlst(r1,totVecs+vecCount,1,1,1,497)
      endif
    else
      if (itr.ne.1) then  
       print*,'** Dropping root:',rootIndx(j),'from the subspace expansion'
       unconvergedRoots(rootIndx(j))=.False.
!       do aaa=j+1,Nblocks
!         bkuprootIndx(aaa-1)=rootIndx(aaa)
!       enddo
!       print*,'** Remaining roots:', bkuprootIndx(1:tmproot-1)
       tmproot=tmproot-1
      endif
    endif


    print*
  enddo
  rootIndx=PACK([(zz, zz=1,Nblocks)], unconvergedRoots)!bkuprootIndx
  nroots=tmproot
  print*,'** number of unconverged roots for iteration:',itr,nroots

end subroutine
