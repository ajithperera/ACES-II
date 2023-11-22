c direct RHF/UHF/ROHF gradient with Gamess integrals

      implicit double precision (a-h,o-z)
#include "../gamess/gamess.par"
#include "mpif.h"
#include "parallel_aces.com"
#include "flags.com"
      common /GMSorder/ coeint(MXAO),iorder(MXAO),
     >               nAfrst(MXSH),nAstep(MXSH),sphrcl
      common // icore(1)
      common /imems/ i00,j00,i03,j03,i05,j05,i06
      common /nbas/ nao,nmo,nat
      COMMON /xGRAD  / DE(3,MXATM)
      COMMON /TIMEINFO/ TIMEIN, TIMENOW, TIMETOT, TIMENEW
      COMMON /MACHSP/ IINTLN, IFLTLN, IINTFP, IALONE, IBITWD
      logical sphrcl,rohf
      character*3 str
      dimension NOtab(2)

      call c_gtod(is0,ius0)
      call aces_init_rte
      call timer(0)

c ----------------------------------------------------------------------

c   o initialize the MPI environment
      call MPW_Init(iErr)
      if (iErr.ne.MPI_SUCCESS) then
         print *, '@SCFGRD: unable to initialize the MPI environment'
         call aces_exit(1)
      end if

c   o initialize the parallel_aces common block
      call aces_com_parallel_aces

#ifdef _HAVE_MPI
c   o change to the temp directory
      call cd_noderank(irank,iErr)
      if (iErr.ne.0) then
         print *, '@SCFGRD: unable to cd to the temp directory'
         print *, '         attempting to run in current dir'
      end if
#endif

c ----------------------------------------------------------------------

      call CRAPSI (icore(1),iuhf,-1)
      call nr2asc (irank,str)
      open (96,file='OUT.'//str)
      write (96,'(a,2i3)') 'nprocs,irank:',nprocs,irank
      
      rohf=.false.
c!      call getrec (1,'','UHFRHF  ',1,iuhf)
      iuhf=iflags(11)
      if (iuhf.eq.2) rohf=.true.
      call getrec (1,'','NOCCORB ',2,NOtab)
      na=NOtab(1)
      nb=na
      if (iuhf.ne.0) nb=NOtab(2)

      if (iuhf.eq.0) then
         write (96,'(a)') 'DIRECT RHF GRADIENT CALCULATION'
      else if (iuhf.eq.1) then
         write (96,'(a)') 'DIRECT UHF GRADIENT CALCULATION'
      else if (iuhf.eq.2) then
         write (96,'(a)') 'DIRECT ROHF GRADIENT CALCULATION'
      end if
      if (irank.eq.0) open (99,file='DIRGRD.LOG')
      open (98,file='MOL')
      call GETINF (98,.true.,99)
      call GRDINI (98,.true.,99, iuhf,na,nb, nat,MEMMAX,iintfp)
      close (98)

      call getrec (1,'','NAOBASFN',1,nao)
      write (96,*) 'no. of AOs:', nao
      call getrec (1,'','NBASTOT ',1,nmo)
      write (96,*) 'no. of MOs:', nmo
      write (96,'(a,2i4)') 'no. of occupied:', na,nb
      write (96,*)
      write (96,*) 'GRADIENT GAMESS initialization complete'

      mem1 = 0
      mem1 = mem1 +nmo         ! Eorba
      mem1 = mem1 +nmo         ! Eorbb
      mem1 = mem1 +nmo*nmo     ! coeSa
      mem1 = mem1 +nmo*nmo     ! coeSb
      mem1 = mem1 +NAO*nmo     ! coea 
      mem1 = mem1 +NAO*nmo     ! coeb
      mem1 = mem1 +NAO*nmo     ! coeGa 
      mem1 = mem1 +NAO*nmo     ! coeGb 
      mem1 = mem1 +NAO*nmo     ! Atab  
      mem1 = mem1 +NAO*NAO     ! densa
      mem1 = mem1 +NAO*NAO     ! densb
      mem1 = mem1 +NAO*NAO     ! Lagrange multipliers
      if (ROHF) then
      mem1 = mem1 +NAO*NAO     ! Focka 
      mem1 = mem1 +NAO*NAO     ! Fockb 
      mem1 = mem1 +NAO*NAO     ! FockGa
      mem1 = mem1 +NAO*NAO     ! FockGb
      end if

      write (96,'(a,i9,a)') 'main DIRGRD module fills',mem1,' words'
      if (mem1.gt.MEMMAX) then
         write (96,'(a)') 'error: limit exceeded'
         call MPIOFF
         stop
      end if
      memalo = mem1*IINTFP
      call ACES_MALLOC (memalo,icore,istart)
      if (icore(1).eq.0) then
         write (96,'(a)') 'unable to allocate...'
         call MPIOFF
         stop
      end if
      mem2lim=MEMMAX-mem1

      i00 = istart                    ! Eorba 
      j00 = i00 + nmo*IINTFP          ! Eorbb
      i01 = j00 + nmo*IINTFP          ! coeSa
      j01 = i01 + nmo*nmo*IINTFP      ! coeSb 
      i02 = j01 + nmo*nmo*IINTFP      ! coea 
      j02 = i02 + NAO*nmo*IINTFP      ! coeb
      i03 = j02 + NAO*nmo*IINTFP      ! coeGa
      j03 = i03 + NAO*nmo*IINTFP      ! coeGb
      i04 = j03 + NAO*nmo*IINTFP      ! Atab
      i05 = i04 + NAO*nmo*IINTFP      ! densa
      j05 = i05 + NAO*NAO*IINTFP      ! densb
      i06 = j05 + NAO*NAO*IINTFP      ! Lagrange multipliers
      i07 = i06 + NAO*NAO*IINTFP      ! end of core

      if (ROHF) then
      i07 = i06 + NAO*NAO*IINTFP      ! Focka
      j07 = i07 + NAO*NAO*IINTFP      ! Fockb
      i08 = j07 + NAO*NAO*IINTFP      ! FockGa 
      j08 = i08 + NAO*NAO*IINTFP      ! FockGb 
      i10 = j08 + NAO*NAO*IINTFP      ! end of core
      end if

      call getrec (1,'','SCFEVLA0',iintfp*nmo, icore(i00))
      call getrec (1,'','SCFEVCA0',iintfp*nmo*nmo, icore(i01))

      if (iuhf.ne.0) then
         call getrec (1,'','SCFEVLB0',iintfp*nmo, icore(j00))
         call getrec (1,'','SCFEVCB0',iintfp*nmo*nmo, icore(j01))
      end if

c transform MO's from SAO's to CAO's
      if (sphrcl) then
         call S2CAO (icore(i04),nao,nmo)
         call DGEMM ('N','N', nao,nmo,nmo, 1.d0,icore(i04),nao, 
     >          icore(i01),nmo, 0.d0,icore(i02),nao)
         if (iuhf.ne.0) 
     >   call DGEMM ('N','N', nao,nmo,nmo, 1.d0,icore(i04),nao, 
     >          icore(j01),nmo, 0.d0,icore(j02),nao)
      else
         call dcopy (nmo*nmo, icore(i01), 1, icore(i02), 1)
         if (iuhf.ne.0) 
     >   call dcopy (nmo*nmo, icore(j01), 1, icore(j02), 1)
      end if   

c transform MO's from ACES to GAMESS order:  coe --> coeG;
      call A2G (icore(i02),icore(i03),NAO,nmo,iorder)
      if (iuhf.ne.0)
     >call A2G (icore(j02),icore(j03),NAO,nmo,iorder)
 
c scale the eigenvectors to account for the differences between the 
c normalization of d,f,g functions in ACES and GAMESS;
c alternatively, both one- and two-electron gradient integrals could
c be scaled, but the call below is simpler
      call DFGscl (icore(i03),NAO,nmo)
      if (iuhf.ne.0) 
     >call DFGscl (icore(j03),NAO,nmo)

c calculating density matrix at i05 and energy-weighted density matrix
c (Lagrange multipliers) at i06
      call c2dens (iuhf,na,nb,icore(i00),icore(j00),
     > icore(i03),icore(j03), NAO,nmo,icore(i05),icore(j05),icore(i06))


c--------------------
c Lagrangian in the ROHF case requires the Fock matrices
      if (ROHF) then
      call getrec (1,'','FOCKA   ',iintfp*nmo*nmo, icore(i07))
      call getrec (1,'','FOCKB   ',iintfp*nmo*nmo, icore(j07))
      if (sphrcl) then
         write (96,'(a,a)') 'ROHF gradient is only supported in ',
     >                   'Cartesian basis set calculations'
         call MPIOFF
         stop
      end if
      call matA2G (icore(i07),icore(i08), NAO)
      call matA2G (icore(j07),icore(j08), NAO)
      call AFA (icore(i05),icore(i08),icore(i06),icore(i07),NAO,0.d0)
      call AFA (icore(j05),icore(j08),icore(i06),icore(j07),NAO,1.d0)
      call DSCAL (NAO*NAO,-1.d0,icore(i06),1)
      end if
c--------------------

      call full2t (icore(i05),NAO)
      if (iuhf.ne.0)
     >call full2t (icore(j05),NAO)
      call full2t (icore(i06),NAO)

c      memory=10000000
c      call ySETFM (memory)
      call ySETFM (mem2lim)

      call timer(1)
      write (96,'(a,f10.2)') 'init time:      ',timenew

      call timer(0)
      call xSTVDER
      call timer(1)
      write (96,'(a,f10.2)') '1-el time:      ',timenew

      call timer(0)
      call xJKDER
      call timer(1)
      write (96,'(a,f10.2)') '2-el time:      ',timenew
      write (96,*)
      if (irank.eq.0) then
      write (*,'(a)') 'SCF GRADIENT:'
      write (*,'(a,20x,a,17x,a,17x,a)') 'atom','x','y','z'
      do 7777 i=1,nat
 7777 write (*,'(i4,f29.9,2f18.9)') i,de(1,i),de(2,i),de(3,i)
      end if
      call putrec(1,'JOBARC','GRADIENT',iintfp*3*nat,de)
      call aces_ja_fin
      call MPIOFF
      call c_gtod(is1,ius1)
      print *, '@SCFGRD: Total wallclock time was ',
     &         (is1-is0)*1.d0+(ius1-ius0)*1.d-6,' seconds.'

      end
c--------------------------------------------------------------------------
      subroutine MPIOFF
#include "mpif.h"
      call MPW_Finalize(iErr)
      if (iErr.ne.MPI_SUCCESS) then
         print *, '@SCFGRD: Error in MPI_Finalize'
         call aces_exit(1)
      end if
      return
      end
c--------------------------------------------------------------------------
      subroutine S2CAO (A,n,nmo)
c compute transformation A from n Cartesian AO's to nmo spherical AO's;
c uses Gamess->Aces mappings established previously to avoid another analysis
c of the MOL file;
c the order of both CAO's and SAO's is the internal ACES order
      implicit double precision (a-h,o-z)
      dimension A(n,nmo)
#include "../gamess/gamess.par"

      COMMON /xNSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     >                CF(MXGTOT),CG(MXGTOT),
     >                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     >                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      common /GMSorder/ coeint(MXAO),iorder(MXAO),
     >               nAfrst(MXSH),nAstep(MXSH),sphrcl
      logical sphrcl
      dimension tabP(3,3),tabD(6,5),tabF(10,7),tabG(15,9)
      data tabP / 1.d0, 0.d0, 0.d0,
     >            0.d0, 1.d0, 0.d0,
     >            0.d0, 0.d0, 1.d0/
      data tabD /-1.d0, 0.d0, 0.d0,-1.d0, 0.d0, 2.d0, ! d0
     >            0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 0.d0, ! d-2           
     >            0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, ! d+1
     >            1.d0, 0.d0, 0.d0,-1.d0, 0.d0, 0.d0, ! d+2
     >            0.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0/ ! d-1
      data tabF /-1.d0, 0.d0, 0.d0,-1.d0, 0.d0, 4.d0, 
     >                        0.d0, 0.d0, 0.d0, 0.d0, ! f+1
     >            0.d0,-1.d0, 0.d0, 0.d0, 0.d0, 0.d0, 
     >                       -1.d0, 0.d0, 4.d0, 0.d0, ! f-1
     >            0.d0, 0.d0,-3.d0, 0.d0, 0.d0, 0.d0, 
     >                        0.d0,-3.d0, 0.d0, 2.d0, ! f0
     >            1.d0, 0.d0, 0.d0,-3.d0, 0.d0, 0.d0, 
     >                        0.d0, 0.d0, 0.d0, 0.d0, ! f+3
     >            0.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 
     >                        0.d0, 0.d0, 0.d0, 0.d0, ! f-2
     >            0.d0, 3.d0, 0.d0, 0.d0, 0.d0, 0.d0,  
     >                       -1.d0, 0.d0, 0.d0, 0.d0, ! f-3
     >            0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 
     >                        0.d0,-1.d0, 0.d0, 0.d0/ ! f+2
      data tabG / 3.d0, 0.d0, 0.d0, 6.d0, 0.d0,-24.d0,0.d0, 0.d0,
     >                  0.d0, 0.d0, 3.d0, 0.d0,-24.d0,0.d0, 8.d0, ! g0    
     >            0.d0,-1.d0, 0.d0, 0.d0, 0.d0, 0.d0,-1.d0, 0.d0,
     >                  6.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, ! g-2
     >            0.d0, 0.d0,-3.d0, 0.d0, 0.d0, 0.d0, 0.d0,-3.d0,
     >                  0.d0, 4.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, ! g+1
     >            1.d0, 0.d0, 0.d0,-6.d0, 0.d0, 0.d0, 0.d0, 0.d0,
     >                  0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 0.d0, ! g+4
     >            0.d0, 0.d0, 0.d0, 0.d0, 3.d0, 0.d0, 0.d0, 0.d0,
     >                  0.d0, 0.d0, 0.d0,-1.d0, 0.d0, 0.d0, 0.d0, ! g-3
     >           -1.d0, 0.d0, 0.d0, 0.d0, 0.d0, 6.d0, 0.d0, 0.d0,
     >                  0.d0, 0.d0, 1.d0, 0.d0,-6.d0, 0.d0, 0.d0, ! g+2
     >            0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 0.d0,-1.d0, 0.d0,
     >                  0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, ! g-4
     >            0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 0.d0,-3.d0,
     >                  0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, ! g+3
     >            0.d0, 0.d0, 0.d0, 0.d0,-3.d0, 0.d0, 0.d0, 0.d0,
     >                  0.d0, 0.d0, 0.d0,-3.d0, 0.d0, 4.d0, 0.d0/ ! g-1

      call zero (A,n*nmo)

      do 10 ish=1,NSHELL
         L = ktype(ish)-1
         nSinL = 2*L+1
         nCinL = ((L+1)*(L+2))/2
         nCAO1 = iorder( KLOC(ish) )
         do 20 i=1,nSinL
            nSAO = nAfrst(ish) +(i-1)*nAstep(ish)
            do 30 j=1,nCinL
               nCAO = nCAO1 + (j-1)*nAstep(ish)
               if (L.eq.0) then
                  A(nCAO,nSAO) = 1.d0
               else if (L.eq.1) then   
                  A(nCAO,nSAO) = tabP(j,i)
               else if (L.eq.2) then   
                  A(nCAO,nSAO) = tabD(j,i)
               else if (L.eq.3) then   
                  A(nCAO,nSAO) = tabF(j,i)
               else if (L.eq.4) then
                  A(nCAO,nSAO) = tabG(j,i)
               else
                  write (96,'(a,i6)') 'S2CAO: unsupported function:',L
                  call MPIOFF
                  stop
               end if
 30         continue
 20      continue
 10   continue
      end
c--------------------------------------------------------------------------
      subroutine A2G (coe,coeG,N,nmo,iorder)
      implicit double precision (a-h,o-z)
      dimension coe(*),coeG(*),iorder(*)

      do 5 iG=1,N
      iA=iorder(iG)
      ind=iG
      do 5 j=1,nmo
      coeG(ind)=coe(iA)
      iA=iA+N
      ind=ind+N
 5    continue
      end
c--------------------------------------------------------------------------
      subroutine c2dens (iuhf,na,nb,Eorba,Eorbb,ca,cb,n,nmo,
     >                   densa,densb,ewdm)
      implicit double precision (a-h,o-z)
      dimension Eorba(*),Eorbb(*),ca(n,nmo),cb(n,nmo),
     >          densa(n,n),densb(n,n),ewdm(n,n)
      data zero,two /0.d0,2.d0/

      do 10 j=1,na
      do 10 i=1,n
 10   densa(i,j)=-ca(i,j)*Eorba(j)
      call DGEMM ('N','T', n,n,na, 1.d0,ca,n, densa,n, 0.d0, ewdm,n)
      call DGEMM ('N','T', n,n,na, 1.d0,ca,n,    ca,n, 0.d0,densa,n)

      if (iuhf.eq.0) then
         l3=n*n
         call DSCAL (l3,two,densa,1)
         return
      end if

      do 20 j=1,nb
      do 20 i=1,n
 20   densb(i,j)=-cb(i,j)*Eorbb(j)
      call DGEMM ('N','T', n,n,nb, 1.d0,cb,n, densb,n, 1.d0, ewdm,n)
      call DGEMM ('N','T', n,n,nb, 1.d0,cb,n,    cb,n, 0.d0,densb,n)


c      do 10 i=1,n
c      do 10 j=1,i
c      sum=zero
c      sumE=zero
c      do 20 k=1,na
c         sum=sum +ca(i,k)*ca(j,k)
c         sumE=sumE -Eorba(k)*ca(i,k)*ca(j,k)
c 20   continue
c      densa(i,j)=sum
c      densa(j,i)=sum
c      ewdm(i,j)=sumE
c      ewdm(j,i)=sumE
c 10   continue
      
c      do 11 i=1,n
c      do 11 j=1,i
c      sum=zero
c      sumE=zero
c      do 21 k=1,nb
c         sum=sum +cb(i,k)*cb(j,k)
c         sumE=sumE -Eorbb(k)*cb(i,k)*cb(j,k)
c 21   continue
c      densb(i,j)=sum
c      densb(j,i)=sum
c      ewdm(i,j)=ewdm(i,j)+sumE
c      ewdm(j,i)=ewdm(i,j)
c 11   continue

      end
c--------------------------------------------------------------------------
      subroutine DFGscl (c,n,nmo)
      implicit double precision (a-h,o-z)
      dimension c(n,nmo)
#include "../gamess/gamess.par"
      common /GMSorder/ coeint(MXAO),iorder(MXAO),
     >               nAfrst(MXSH),nAstep(MXSH),sphrcl
      logical sphrcl

      do 10 j=1,nmo
      do 10 i=1,n
 10   c(i,j)=c(i,j)*coeint(i)

      end
c--------------------------------------------------------------------------
      subroutine matA2G (A,G, n)
c calculates the N*N matrix G in the Gamess Cartesian AO basis
c given the N*N matrix A in the ACES Cartesian AO basis
      implicit double precision (a-h,o-z)
      dimension A(n,n),G(n,n),W(n,n)
#include "../gamess/gamess.par"
      common /GMSorder/ coeint(MXAO),iorder(MXAO),
     >               nAfrst(MXSH),nAstep(MXSH),sphrcl
      logical sphrcl

      do 10 j=1,n
         jx=iorder(j)
         cjx=coeint(j)
      do 10 i=1,n
         ix=iorder(i)
         cix=coeint(i)
         G(i,j)=A(ix,jx) / (cix*cjx)
 10   continue

      end
c--------------------------------------------------------------------------
      subroutine AFA (A,F,G,W, n,x)
c calculates G(n,n)=A(n,n)*F(n,n)*A(n,n) + x*G(n,n) using the work array W(n,n)
      implicit double precision (a-h,o-z)
      dimension A(n,n),F(n,n),G(n,n),W(n,n)

      call DGEMM ('N','N',n,n,n, 1.d0,A,n, F,n, 0.d0,W,n)
      call DGEMM ('N','N',n,n,n, 1.d0,W,n, A,n, x,G,n)

      end
c--------------------------------------------------------------------------
      subroutine full2t (A,n)
c packs the full triangle of a symmetric matrix A(n,n) into the first
c (n*(n+1))/2 elements of A  
      implicit double precision (a-h,o-z)
      dimension A(*)

      ind=0
      do 10 j=1,n
      do 10 i=1,j
      ind=ind+1
      a(ind)=a((j-1)*n+i)
 10   continue
      end
