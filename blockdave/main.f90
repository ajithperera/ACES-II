      subroutine blockdavid_tda(Hbar,Work,Nsize,Maxmem,Iuhf,cROOTStotal,&
                               MAXITER,NUMSOL,Tol,cVALS,cEVECS)
      implicit none
      integer:: i,j,Nsize,cROOTStotal,Ione,Maxmem,Iuhf
      integer i000,i010,MAXITER,Numsol
      double precision:: TOL
      double precision:: Hbar(nsize,nsize),cVALS(NUMSOL)
      double precision:: Work(Maxmem),cEVECS(Nsize,NUMSOL)

      Data Ione /1/

!! Initialize Guess C of size (Nsize x cROOTStotal), but entire C matrix
!! has dimension (Nsize x Nsize)
!! Can be later revised to size (Nsize x Nsize*MAX. DAVIDSON ITERATIONS)
!! to conserve memory
        print*,'Nsize',Nsize
        print*,'maxmem',Maxmem

!C Comment:
!C Change the meaning of the initvec and numsol to the following.
!C NUMSOL : the total number of roots requested
!C cROOTStotal : The total number of roots obtained (rename this variable 
!C to something else that is transparent (NUMSOL >= INitvec)

#ifdef _DEBUG_LVLM
        Write(6,"(a)") "The CIS matrix at entry"
!        call output(Hbar,1,Nsize,1,Nsize,Nsize,Nsize,1)
        print*,'Nsize:',Nsize
#endif
       I000 = Ione 
       I010 = I000 + Nsize
       call initGuess(Work(I000),Nsize,Nsize)

!! Perform block Davidson routine
        print*,'cROOTStotal:',cROOTStotal
       call davidsonDrive(Hbar,Work(I000),NSize,cROOTStotal,NUMSOL,&
                            MAXITER, cVALS,cEVECS,TOL)


!! ** IF USING NLMO FRAMEWORK, DISCARD THOSE EXCITATIONS HAVING
!!    ORIGIN OUTSIDE QM1
      ! call NLMOcorrect(cVals,cEVECS,Nsize,NUMSOL)

       return
       end
       
