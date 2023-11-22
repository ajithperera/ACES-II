      subroutine heffdiag(iuhf,icycle)
      implicit integer (i-n)
      implicit real*8  (a-h,o-z)
#include "bwcc.com"
cjp
cjp this routine diagonalizes the effective hamiltonian from common bwcc
cjp selects the root wanted and stores the energy in common bwcc
cjp
c
      logical o
      dimension heff2(maxref,maxref)
      dimension work(8*maxref)
      dimension ip(maxref)
      dimension wr(maxref),wi(maxref),vl(maxref,maxref),
     &          vr(maxref,maxref)
cjp tolerance how large the imaginary part is neglected 
      parameter(tolimag=1e-5)

cjp parallelization
cjp broadcast columns of Heff between processors
cjp in a parallel run
       o=.true.
c
cjp continue diagonalization on all processors
        if(o)write(6,*)
        if(o)write(6,*)
cjp test if heff is a normal matrix
      tmp=0d0
      do i=1,nref
      do j=1,nref
      heff2(i,j)=0d0
      do k=1,nref
        heff2(i,j)=heff2(i,j)+heff(i,k)*heff(j,k)-heff(k,i)*heff(k,j)
      enddo
      if(dabs(heff2(i,j)).gt.tmp) tmp=dabs(heff2(i,j))
      enddo
      enddo
        if(o)write(6,*) 'departure of Heff from normal matrix',tmp
        if(o)write(6,*)
      do i=1,nref
       if(o)write(6,123)  (heff2(i,j),j=1,nref)
      enddo
        if(o)write(6,*)
        if(o)write(6,*)


cjp dgeev is destructive on input matrix
      if(nref.eq.1) then
cjp trivial
      heffevalr(1)=heff(1,1)
      heffevali(1)=0d0
      heffevecl(1,1)=1d0
      heffevecr(1,1)=1d0
      else
      do i=1,nref
      do j=1,nref
       heff2(j,i)=heff(j,i)
      enddo
      enddo
cjp
cjp NOTICE!!! dgeev normalizes eigenvectors in a strange way
cjp not to have L**H . R = 1, but euclidean norm of individual eigenvectors to be one
cjp this is to be considered if they should be used in any later computation!!!
cjp
      call xgeev('V','V', nref, heff2, maxref, wr, wi, vl, maxref, vr,
     &            maxref,work, 8*maxref, info )
      if (info.ne.0) stop 'ERROR in Heff diagonalization'
cjp permute the eigenstuff in increasing order of real part of eigenvalue
      do i=1,nref
      ip(i)=i
      heffevalr(i)=wr(i)
      enddo
      call sort2(nref,heffevalr,ip)
      

cjp store permuted results
      do i=1,nref
      ii=ip(i)
      heffevalr(i)=wr(ii)
      heffevali(i)=wi(ii)
      do j=1,nref
      heffevecl(j,i)=vl(j,ii)
      heffevecr(j,i)=vr(j,ii)
      enddo
      enddo

cjp
cjp renormalize eigenvectors to become bi-unitary
cjp this calculation is valid only if they are real-valued!
cjp do only if it turns out they will be used later
cjp

cjp end of the nontrivial case
      endif

cjp print results:
      if(o)write(6,*) 
      if(o)write(6,*) 'MR-BW-CC EFFECTIVE HAMILTONIAN'
        if(o)write(6,*)

      do i=1,nref
c      if(o)write(6,*) (heff(i,j),j=1,nref)
       if(o)write(6,123)  (heff(i,j),j=1,nref)
      enddo
        if(o)write(6,*)
      if(o)write(6,*) ' results of diagonalization (for complex 
     & eivectors for storage cf. man dgeev)'
123   format(10f16.10)
       if(o)write(6,*)
      if(o)write(6,*) 'eigenvalues (real and imaginary)'
        if(o)write(6,*)
        if(o)write(6,123) ( heffevalr(i),i=1,nref)
        if(o)write(6,123) ( heffevali(i),i=1,nref)
       if(o)write(6,*)
      if(o)write(6,*) 'right eigenvectors'
        if(o)write(6,*)
        do j=1,nref
        if(o)write(6,123) ( heffevecr(j,i),i=1,nref)
        enddo
       if(o)write(6,*)
      if(o)write(6,*) 'left eigenvectors'
        if(o)write(6,*)
        do j=1,nref
        if(o)write(6,123) ( heffevecl(j,i),i=1,nref)
        enddo
        if(o)write(6,*)
      if(iuhf.gt.0) then
        if(o)write(6,*)
      if(o)write(6,*) 'NOTE for spin coupling of open shell references'
      if(o)write(6,*) 'references like ...ab, ...ba actually corr. to'
      if(o)write(6,*) '                ...hl, ...hl'
      if(o)write(6,*)
     +    'determinants |...,h-alpha,l-beta| and |...,l-alpha,h-beta|'
      if(o)write(6,*) 'not |...,h-beta,l-alpha|.Therefore triplet state'
      if(o)write(6,*) 'corr. to coupling of ...ab, ...ba with - sign'
        if(o)write(6,*)
      endif
        if(o)write(6,*)

cjp compute some simple criteria for tworeference problems
      if(nref.eq.2) then
      if(o)write(6,*)'|H01/H00| = ',abs(heff(1,2)/heff(1,1))
      if(o)write(6,*)'2*SQRT(H01*H10)/|H00-H11| = ',
     +   2d0*sqrt(heff(1,2)*heff(2,1))/abs(heff(1,1)-heff(2,2))
      if(o)write(6,*)'E0-H00 = ',heffevalr(1)-heff(1,1)
      endif

cjp test for swapped eigenvectors (only for those we are interested in)
      if(icycle.gt.1) then
      do i=1,nbwstates
        k=ibwstate(i)
        s=0.
        do j=1,nref
        s=s+heffevecrold(j,k)*heffevecl(j,k)
        enddo
        if(abs(s).lt.0.5) then
         if(o)write(6,*)'@HEFFDIAG-W: Eivec ',k,' swapped!!! icycle= ',
     &             icycle
         bwwarning(3)=.true.
         bwwarntext(3)='@HEFFDIAG-W: Eigenvector swap occured'
        endif
      enddo
      endif
cjp store eigenvectors for this test in next iteration
      do i=1,nref
      do j=1,nref
      heffevecrold(j,i)=heffevecr(j,i)
      enddo
      enddo


cjp test for complex roots
      iflag=0
      do i=1,nref
        if(abs(heffevali(i)).ge.tolimag) iflag=iflag+1
      enddo
      if(iflag.ne.0) then
          if(o)write(6,*)'@HEFFDIAG-W: ',iflag,
     &       ' complex Heff eigenvalues detected'
          bwwarning(4)=.true.
          bwwarntext(4)='@HEFFDIAG-W: complex Heff eigenvalue detected'
      endif

cjp   write the energy for later use
      epsilon0=0.0
      do i=1,nbwstates
cjp problem if the requested eigenvalue is complex
          if(abs(heffevali(ibwstate(i))).ge.tolimag) then
           if(o)write(6,*) '@HEFFDIAG-E: MR-BW-CC energy is complex'
           stop '@HEFFDIAG-E: MR-BW-CC energy is complex'
          endif
          epsilon0=epsilon0+heffevalr(ibwstate(i))*cbwstate(i)
      enddo

      if(bwgossip.and.o)write(6,*)'CORRELATION ENERGY FROM MR-BW-CC = ',
     &     epsilon0
       if(o)write(6,*)
      return
      end
