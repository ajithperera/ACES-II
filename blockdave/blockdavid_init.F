      Subroutine Block_david_driver(Irrepx,Work,Maxcor,Nsize,Iuhf,Iside,
     +                              Tol,Maxiter,Nblocks)

      Implicit Double Precision(A-H,O-Z)
      integer jj,nroots,totVecs,NUMSOL
      Dimension Work(Maxcor)
      Logical Converged,DONLS,CISflag
      integer filterIndx
      logical unconvergedRoots(Nblocks)
      double precision thetaOld(Nblocks),diff(Nblocks)
      double precision temp(NSize)
      integer rootIndx(Nblocks)
#include "machsp.com"
#include "syminf.com"
#include "sympop.com"
#include "flags.com"
#include "flags2.com"
      
      Common /Calcinfo/Nroot(8)

      Parameter(Maxblocks=20)

C This code assume that the TDA vectors are availble for the 
C requested number of roots via estate_sym (and less <= 10 for
C each symmtery block)
      
      print*,'Entered Block David for EOM'
      print*,Irrepx,Maxcor,Nsize,Iuhf,Iside,
     +             Tol,Maxiter,Nblocks
      print*
      unconvergedRoots=.True.
      Tol=0.0001d0!0.000000001d0 !0.0000001d0 
      nroots=Nblocks !initialize first iteration
      !Nblocks = Nroot(Irrepx)
      If (Nblocks .Gt. Maxblocks) Then
         Write(6,"(2a,i2,2a,i2)") " The requested number of roots",
     +                            " for irrep", Irrep, "greater",
     +                            " than the aximum allowed value", 
     +                               Maxblocks
         Call Errex
      Endif

C Pick the TDA eigenvectors for each root per symmetry block and
C carry out the multiplication and accmulate 
      !Tol=0.0001
      Do Ispin=3,3-2*Iuhf,-1
         Call Zerolist(Work,Maxcor,443+Ispin)
      Enddo 

C This is the start of the block Davidson loop

       Converged = .False.
       Iter_count = 0
       totVecs=0
       do i=1,NBlocks
        rootIndx(i)=i
       enddo
       CISfilterFlag=.True.
       Do while (.not.Converged) 
          Iter_count = Iter_count + 1
          Ioffc = 1
          print*,'********** Begin iteration:, ',Iter_count
          print*,'starting on vectors:',totVecs
          print*,'number of unconverged roots:', nroots
          print*,'list of unconverged roots:',rootIndx(1:nroots)
          print*,'*********************************'
!         ** Create list & put CIS vecs on oldRvec list 497
!         ** AFTER HBARXC::: go on 498
          if (Iter_count.eq.1) then
            Call Updmoi(Nblocks*Maxiter,nsize,Irrepx,497,0,0)
            Call Updmoi(Nblocks*Maxiter,nsize,Irrepx,498,0,0)
            Call Dzero(Work,Nsize)
            do k=1,Nblocks
                 Call Getlst(Work,k,1,1,Irrepx,94)
! ** filters out CIS vectors with single particle excitations outside
! QM1. Still adds this vector to the subspace so that higher energy
! vectors can make rough attempt to be orthogonal. The hope is to avoid
! subspace collapse
!                 if (CISfilterFlag) then
!                   filterIndx=0
!                 CALL filterRoots(thetaOld,Work(1),Nsize,1,NUMSOL,
!     +             CISflag,filterIndx)
!                   if (filterIndx.ne.0) then
!                     unconvergedRoots(filterIndx)=.False.
!                   endif
!                 endif
                 If (Iuhf .Eq. 0) Then
                    Fact = One/Dsqrt(Two)
                    Call Dscal(Irpdpd(Irrepx,9),Fact,Work,1)
                 Endif
           !      Work=0.0d0
           !      Work(k)=1.0d0
              call checksum("initVec:",Work,NSize,s)
              Call Putlst(Work,k,1,1,Irrepx,497)
            enddo
          endif

!! ** BEGIN PRIMARY INNER LOOP TO COMPUTE HBAR*C **
C It is assuemed that we need at least Nsize*Nblocks of memory to
C proceed. Note that this does not count any additional requirements 
C within blokdavison. Memory checks are essential since this is 
c very memory intensive process. 
#ifdef _DEBUG_LVLM
           print*,'Nsize is: ', Nsize
#endif
           Ioff_hc = 1

! *** CAN PRINT OUT FULL EOM HBAR MATRIX ***
           !call printHBAR(Work,Nsize,Maxcor,Iuhf,Irrepx)
           !stop

           Do Iblck = 1, nroots!Nblocks 
              I000 = Ione
              Iend = I000 + Nsize*nroots!Nblocks
              Memleft = Maxcor - Iend 
              If (Iend .Ge. Maxcor) Call Insmem("block_david_driver",
     +                                            Iend,Maxcor)

              Call Getlst(Work(Ioff_hc),Iblck+totVecs,1,0,Irrepx,497)

              Call LANCZOS_DUMP_VEC(Irrepx,Work(Ioff_hc),NSize,
     +                                 490,0,0,443,0,IUHF,.False.)

              Call Hbarxc(Work,Maxcor,Iuhf,1,Irrepx)

              Call Loadvec1(Irrepx,Work(Ioff_hc),Maxcor,Iuhf,490,2,460,
     +                      Nsize,.False.)

              Call Putlst(Work(Ioff_hc),Iblck+totVecs ,1,0,Irrepx,498)
#ifdef _DEBUG_LVLM
              call checksum("H*C:",Work(Ioff_hc),NSize,s)
#endif
              Ioff_hc = Ioff_hc + Nsize
           Enddo
          totVecs=totVecs+nroots 
#ifdef _DEBUG_LVLM
          print*,'total number of vectors in subspace:', totVecs
          print*,"outside primary loop!!"
#endif
         Ioff_hc = 1
         Ioff_c=Ioff_hc+NSize*totVecs !Nblocks*Iter_count
         Ioff_T=Ioff_c+NSize*totVecs !Nblocks*Iter_count
         Ioff_Tvec=Ioff_T+(totVecs)**2!(Nblocks*Iter_count)**2
         Ioff_Tval=Ioff_Tvec+(totVecs)**2!(Nblocks*Iter_count)**2
         Iend=Ioff_Tval+totVecs!Nblocks
              If (Iend .Ge. Maxcor) Call Insmem("block_david_driver",
     +                                            Iend,Maxcor)


! build the upper hessenberg matrix T= C^t Hbar C to diagonalize.
! Ritz eigenvalue approx to Hbar is eigenvalue of T, and corresponding vector is
! C*Tvec
        call constructT(Work(Ioff_hc),Work(Ioff_c),Nsize,Nblocks,
     +                          totVecs,
     +  Iter_count,Maxiter,Work(Ioff_T),Work(Ioff_Tvec),Work(Ioff_Tval))

#ifdef _DEBUG_LVLM
        print*,'evals (au):'
        print*,Work(Ioff_Tval:Ioff_Tval+Nblocks-1)
        print*,'evals (eV):'
        print*,Work(Ioff_Tval:Ioff_Tval+Nblocks-1)*27.2114
#endif
        thetaOld=Work(Ioff_Tval:Ioff_Tval+Nblocks-1)
! Use (Hbar*C - thetaOld*I*C)*Tvec == r1 to determine new
! guess vector for each unconverged root in the block
       call extrapBlock(Work(Ioff_hc),Work(Ioff_c),Work(Ioff_Tvec),
     +   thetaOld,TOL,Nblocks,NSize,Iter_count,totVecs,nroots,
     +   rootIndx,unconvergedRoots)
        print*,'evals AFTER extrapBlock:'
        print*,Work(Ioff_Tval:Ioff_Tval+Nblocks-1)
        print*,'evals (eV):'
        print*,Work(Ioff_Tval:Ioff_Tval+Nblocks-1)*27.2114

      print*,'*** new nroots after extrapBlock:',nroots
      diff=abs(thetaOld-Work(Ioff_Tval:Ioff_Tval+Nblocks-1))
      print*,'thetaOld',thetaOld
      print*,'newtheta',Work(Ioff_Tval:Ioff_Tval+Nblocks-1)
      print*,'diff:',diff
      thetaOld=Work(Ioff_Tval:Ioff_Tval+Nblocks-1)
      print*,'theta old',thetaOld
      print*,'TOL:',TOL
      print*,all(diff.lt.TOL)

       if (nroots.eq.0) then
        print*,"** Converged on iteration: ",Iter_count
        print*
        print*
        print*,'** List of converged roots: **'
        print*,thetaOld
        print*,'** in (eV) **'
        print*,thetaOld*27.2114
        Converged=.True.
!! Write final converged eVECS to 497, eVALS to 95 
        print*,'** Writing eVecs to disk **'
        print*,'number of blocks:',Nblocks
        call storeFinalVec(Work(Ioff_c),Work(Ioff_Tvec),NSize,
     +                  Nblocks,Iter_count,totVecs)
        do jj=1,Nblocks
          call putlst(thetaOld(jj),jj,1,1,1,95)
        enddo 
        exit
       endif


      if (Iter_count.gt.200) then
        stop
      endif
      Enddo 
!! ** Prioritizes the sequential storage of vectors having significant single particle
!excitations within QM1. The hope is to bias ACESII single vector
!Davidson toward the relevant roots. This plan still needs work.
!3/25/2021 ZWW
!      DONLS=.true.
!      if (DONLS) then
!        CISflag=.False.
!        CALL filterRoots(thetaOld,Work(Ioff_c),Nsize,Nblocks,NUMSOL,
!     +           CISflag,filterIndx)
!      endif
!! *** CALL THE NLS PROCEDURE TO GET ESTIMATE ROOTS
!! 
!! 
      print*,'**********************************************'
      print*,'**********************************************'
      print*,'**********************************************'
      print*,'**********************************************'
      print*,'** NOW PERFORMING NLS PORTION **'
      print*,'**********************************************'
      print*,'**********************************************'
!      Ioff_c=1
!      iroot=1
!      do i=1,Nblocks 
!        Work=0.0d0
!        call getlst(Work,i,1,0,Irrepx,497)
!        call NLSZEROPartTwo(Work,MAXCOR,irrepx,iuhf,iroot,nsize)
!        iroot=iroot+1
!        Ioff_c=Ioff_c+NSize
!      enddo
!      print*,'Starting Single Vector David*** '
!      call singledavid_init(Irrepx,Work,Maxcor,Nsize,Iuhf,Iside,
!     +                         Tol,Maxiter,Nblocks)
!      stop
      Return 
      End
   
