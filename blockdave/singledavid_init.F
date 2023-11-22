      subroutine singledavid_init(Irrepx,Work,Maxcor,Nsize,Iuhf,Iside,
     +                         Tol,Maxiter,Nblocks)
      Implicit Double Precision(A-H,O-Z)
      integer aa,jj
      Dimension Work(Maxcor)
      Logical Converged
      double precision thetaOld(1),diff(1)
      double precision temp(NSize)

      integer sdNblocks,ioff
      double precision eVecSCR(Nblocks*Nsize)



        Tol=1E-10
!! ** Read in guess vectors from block david scheme
          ioff=1
        do i=1,Nblocks
        call getlst(eVecSCR(ioff),i,1,1,1,497)
        call checksum("load initGuesses:",eVecSCR(ioff),NSize,s)
        ioff=ioff+NSize
        enddo


        sdNblocks=1
        sdioff=1


!!!!! ** IMPLEMENTED SINGLE VECTOR DAVIDSON TO BE USED ON BLOCK VECS
!!!!!
      do aa=1,Nblocks ! Iterate over all roots independently

        Converged = .False.
        Iter_count = 0
        print*
        print*
        print*,'Starting Single vector davidson for root', aa
        Do while (.not.Converged)
          Iter_count = Iter_count + 1
          Ioffc = 1
          print*,'********** Begin iteration:, ',Iter_count
          if (Iter_count.eq.1) then
            Call Dzero(Work,Nsize)
              Call putlst(eVecSCR(sdioff),1,1,1,Irrepx,497)
              call checksum("initVec:",eVecSCR(sdioff),NSize,s)
              sdioff=sdioff+NSize
          else ! dump prev iterations vectors to 498 AFTER QR DECOMP
            call blockQR(Work(Ioffc),Iter_count,sdNblocks,
     &                                          Nsize,IRREPX)
          endif

!! ** BEGIN PRIMARY INNER LOOP TO COMPUTE HBAR*C **
C It is assuemed that we need at least Nsize*sdNblocks of memory to
C proceed. Note that this does not count any additional requirements
C within blokdavison. Memory checks are essential since this is
c very memory intensive process.
           print*,'Nsize is: ', Nsize
           Ioff_hc = sdNblocks*NSize+1

           Do Iblck = 1,1! sdNblocks
              I000 = Ione
              Iend = I000 + Nsize*sdNblocks
              Memleft = Maxcor - Iend
              If (Iend .Ge. Maxcor) Call Insmem("block_david_driver",
     +                                            Iend,Maxcor)

              Call Getlst(Work(Ioff_hc),Iblck+sdNblocks*(Iter_count-1),
     +                                                  1,0,Irrepx,497)
              Call LANCZOS_DUMP_VEC(Irrepx,Work(Ioff_hc),NSize,
     +                                 490,0,0,443,0,IUHF,.False.)

              Call Hbarxc(Work,Maxcor,Iuhf,1,Irrepx)

              Call Loadvec1(Irrepx,Work(Ioff_hc),Maxcor,Iuhf,490,2,460,
     +                      Nsize,.False.)
              Call Putlst(Work(Ioff_hc),Iblck+sdNblocks*(Iter_count-1)
     +                                                  ,1,0,Irrepx,498)

              call checksum("H*C:",Work(Ioff_hc),NSize,s)
              Ioff_hc = Ioff_hc + Nsize

           Enddo
          print*,"outside primary loop!!"

         Ioff_hc = 1
         Ioff_c=Ioff_hc+NSize*sdNblocks*Iter_count
         Ioff_T=Ioff_c+NSize*sdNblocks*Iter_count
         Ioff_Tvec=Ioff_T+(sdNblocks*Iter_count)**2
         Ioff_Tval=Ioff_Tvec+(sdNblocks*Iter_count)**2
         Iend=Ioff_Tval+sdNblocks*Iter_count
              If (Iend .Ge. Maxcor) Call Insmem("block_david_driver",
     +                                            Iend,Maxcor)



        call constructT(Work(Ioff_hc),Work(Ioff_c),Nsize,
     + sdNblocks,Iter_count,Maxiter,Work(Ioff_T),Work(Ioff_Tvec),
     +                                  Work(Ioff_Tval))
        !stop
        print*,'evals:'
        print*,Work(Ioff_Tval:Iend-1)
        print*,'evals (eV):'
        print*,Work(Ioff_Tval:Iend-1)*27.2114
        print*,'before extrapBlock'
        call extrapBlockSD(Work(Ioff_hc),Work(Ioff_c),Work(Ioff_Tvec),
     +                Work(Ioff_Tval),sdNblocks,NSize,Iter_count,aa)

      if(Iter_count.eq.1) then
        thetaOld=1.0d0
      endif
      diff=abs(thetaOld-Work(Ioff_Tval:Ioff_Tval+sdNblocks-1))
      print*,'diff:',diff
      print*,'TOL:',TOL
      print*,all(diff.lt.TOL)

       if (all(diff.lt.TOL)) then
        print*,"** Converged on iteration: ",Iter_count
        print*
        print*
        print*,'** List of converged roots: **'
        print*,thetaOld
        print*,'** in (eV) **'
        print*,thetaOld*27.2114
        Converged=.True.
      endif
       thetaOld=Work(Ioff_Tval:Ioff_Tval+Iter_count*sdNblocks-1)
       print*,'theta old',thetaOld!*27.2114

      enddo
      enddo

      print*,'**** End of single vector Davidson'
      stop
      return
      end 
