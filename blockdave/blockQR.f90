        subroutine blockQR(allRvecs,itr,Nblocks,Nsize,IRREPX,totVecs)
        integer,intent(in)::totVecs,Nblocks,Nsize,IRREPX,itr
        double precision,intent(in)::allRvecs(Nsize*totVecs)


        double precision:: tau(totVecs),work(totVecs)
        integer::lwork,info,j,k
        integer::tot
! ** VARIABLES: **
! newRvecs - the R vectors recently compiled after hbarxc
! allRvecs - scratch space meant to house the prior and current
! iterations R vectors
! Nblocks - Number of columns per block in our Davidson scheme
! Nsize - row dimension of R vectors
!
! ** PURPOSE **
! The routine loads/stores the old R vectors in 'allRvecs(:,:Nblocks)'
! and places the current iterations R vectors in
! 'allRvecs(:,Nblocks:2*Nblocks)'. Then, allRvecs is sent to BLAS QR
! routine. Upon exit, the new orthogonalized R vectors will be placed at
! allRvecs(:,Nblocks:2*Nblocks); copy this over to 'newRvecs' and exit. 

! Load in old R vectors first.
        print*,'**Block QR loading:',totVecs,'columns'
        print*
        k=1
        do j=1,totVecs
          call getlst(allRvecs(k),j,1,1,1,497)
          !call checksum("checkq1:",allRvecs(k),NSize,s)
          k=k+Nsize
        enddo

! Now copy over new R vectors -> allRvecs
!        call dcopy(Nsize*Nblocks,newRvecs,1,allRvecs(k),1)

! Now perform Lapack QR
!        lwork=Nblocks*itr
!        tot=Nblocks*itr
!        call dgeqrf(Nsize,tot,allRvecs,Nsize,tau,work,lwork,info)
!        if (info.eq.0) then
!          print*,'Optimal lwork QR:', work(1)
!        endif
!        if (info.lt.0) then
!          print*,'QR decomposition failed; check blockQR.f90'
!          stop
!        endif
!        call dorgqr(Nsize,tot,tot,allRvecs,Nsize,tau, &
!                        work,lwork,info)
!        if (info.lt.0) then
!          print*,'QR decomposition failed; check blockQR.f90'
!          stop
!        endif

        call GramSchmidt(allRvecs,Nsize,totVecs)
! Since the assumption is made that the previous Nblocks*(itr-1) blocks
! are already orthonormal, only write the final 
        k=1
        do j=1,totVecs
          call checksum("initVec:",allRvecs(k),NSize,s)
          print*,'QR overlap:',dot_product(allRvecs(k:k+Nsize-1),&
                                           allRvecs(k:k+Nsize-1))
          call putlst(allRvecs(k),j,1,1,1,497)
          !print*,'1stELE:', allRvecs(k)
          k=k+Nsize
        enddo

        end subroutine
