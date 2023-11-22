subroutine partitionAlgo(CISmat,NSize,Nblock)
  integer,intent(in)::NSize,Nblock
  double precision,intent(in)::CISmat(NSize,NSize)

  double precision::Hbar(NSize,NSize),P(NSize,NBlock),Q(Nsize,NSize-NBlock)
  double precision::projP(NSize,NSize),projQ(NSize,NSize)
  logical::sVecAlgo,Converged,BWPT,RSPT
  integer::z,i,j,order,itr
 

  double precision::H0(NSize,NSize),V(NSize,NSize)
  double precision::R0(NSize,NSize)
  double precision:: tse(NBlock,NBlock),E0(NBlock,NBlock)
  double precision::eps0(NBlock,NBlock),tseOld(NBlock,NBlock)
  double precision::norm,pNew(NSize,2)
  double precision::fullSpace(NSize,100*Nblock)
  double precision::tmpProjP(NSize,NSize)!,evecs(NSize,NSize)
  double precision:: RHR(2,2),evecs(2,2),fatPhi(NSize,2)
  tseNew=0.0d0
  pNew=0.0d0
  tseOld=0.0d0
  Hbar=CISmat
  BWPT=.TRUE.
  RSPT=.FALSE.
  Converged=.False.
  sVecAlgo=.True.
  order=0 
  call definePQ(CISmat,NSize,1,P,Q,projP,projQ,sVecAlgo)
  itr=1
  do while (.not.Converged)
    print*
    print*,'************************'
    print*,'************************'
    print*,'Beginning itration:',itr
    if (itr.ne.1) then
!      projP=matmul(pNew,transpose(pNew(:,:Nblock)))
!      projP=matmul(fullSpace(:,:itr),transpose(fullSpace(:,:itr)))
!       projP=matmul(pNew,transpose(pNew(:,:1)))
!      tmpProjP=projP
!      call eig(tmpProjP,evecs,0,NSize,0)
!      print*,'eigvalues of projP:',tmpProjP
!      Hbar=matmul(matmul(projP,Hbar),projP)
!      call definePQ(Hbar,NSize,1,P,Q,projP,projQ,sVecAlgo)
    endif
    print*,'test on outerProj of H0 only'
    print*,matmul(matmul(projP,H0),projP)
    call defineH0(Hbar,NSize,Nblock,P,tse,H0)
    call defineV(Hbar,NSize,V)
    print*,'tse',tse
!    TSE=0.841916685500910
    E0=tse
    if (itr.eq.1) print*,'Initial eval guess:',tse
    print*,'Going into R0def'
    z=1
    call defineR0(R0,tse,NSize,Nblock,H0,q,RSPT,BWPT,sVecAlgo,z)
    call extPSpace(p,R0,V,NSize,NBlock,order,itr,fullSpace,sVecAlgo)
    call computeRHR(Hbar,tse,fullSpace,NSize,Nblock,order,itr,sVecAlgo,fatPhi)
    norm=dot_product(fatPhi(:,2),fatPhi(:,2))
    print*,'norm of fatPhi:',norm
!    norm=1.0d0/sqrt(norm)
!    fatPhi=norm*fatPhi
    print*,'Revised eigenvalue:',tse
    print*,'Correction:',tse-E0


    if (abs(abs(tseOld(1,1))-abs(tse(1,1))).gt.1.0D-9) then
      tseOld=tse
      print*,'fat phi:',fatPhi(:,1)
      print*,'fat phi2:',fatPhi(:,2)
      projP=matmul(fatPhi,transpose(fatPhi(:,:)))
      Hbar=matmul(matmul(projP,Hbar),projP)
      print*,'printing new Hbar:'
      call output(Hbar,1,NSize,1,NSize,NSize,NSize,1)
!      pNew(:,1)=fullSpace(:,1)+fullSpace(:,itr+1)
!      p(:,1)=fullSpace(:,itr+1)!p(:,1)+fullSpace(:,itr+1)
!      q=p
!      pNew(:,1)=fullSpace(:,2)!+fullSpace(:,2)
!      p=pNew
      exit
      print*,'test to first guess'
      print*,matmul(matmul(transpose(p(:,:)),Hbar),p)
    else
      print*,'************************'
      print*,'Converged on itration:',itr
      print*,'Converged eigenvalue:',tse
      Converged=.true.

    endif

    if (itr.gt.10) then
        print*,'over 10 itrations.... stopping'
        exit
    endif

    itr=itr+1
  enddo
  
end subroutine
