      subroutine exchange(buf,ibuf,dens,dens1,
     &    cdens,imap,ehfkin,ehfnatr,ehfx,
     &    ehfcoul,ehar,escf,iflg)

      implicit none
#include <aces.h>
#include <iuhf.com>
#include <vmol.com>
c TODO: these need proper gfname lookup calls
      integer iiiiio, ijijio, iijjio, ijklio, vpoutio
      character*(*) iiiifil, ijijfil, iijjfil, ijklfil, vpoutfil
      parameter (iiiiio=10)
      parameter (ijijio=21)
      parameter (iijjio=22)
      parameter (ijklio=23)
      parameter (vpoutio=30)
      parameter (iiiifil='IIII')
      parameter (ijijfil='IJIJ')
      parameter (iijjfil='IIJJ')
      parameter (ijklfil='IJKL')
      parameter (vpoutfil='VPOUT')
#include <intgrtflags.com>
#include <machsp.com>
#include <flags.com>
#include <mol.com>
#include <symm2_ks.com>
#include <dftfunc.com>

      integer ibuf(ilnbuf),imap(nbastot),iflg
      M_REAL
     &    buf(ilnbuf),dens(nbastot,nbastot,iuhf+1),
     &    cdens(nbastot,nbastot,iuhf+1),
     &    dens1(nbastot*(nbastot+1)/2,iuhf+1),
     &    ehfkin,ehfnatr,ehfx,ehfcoul,ehar,escf,ierr

      integer
     &    iupki,iupkj,iupkk,iupkl,i,icnt,roff,isym,j,int,
     &    nut,k,l,i23,i24,i25,ind,index,ispin,alpha,beta,unit(3),
     &    iunit,nunit
      M_REAL
     &    ecoulab,ecoulaa,ecoulbb,eexchaa,eexchbb,x,tmp1a,tmp2a,tmp2b,
     &    tmp3a,tmp3b,fnuc,field,dipnuc,ehf1e,factor

      integer indxt,irrep,ibot,iposp,length,nbasir,indx1,indx2,iposf

      character*8 label
      character*8 dum1,dum2,dum3,title
      character*32 junk
      character*20 files(3)
      logical exist

      indxt(i,j)=i+(j*(j-1))/2

c Statement functions for packing and unpacking indices.
      iupki(int)=iand(int,ialone)
      iupkj(int)=iand(ishft(int,-ibitwd),ialone)
      iupkk(int)=iand(ishft(int,-2*ibitwd),ialone)
      iupkl(int)=iand(ishft(int,-3*ibitwd),ialone)

      call callstack_push('EXCHANGE')
      factor=1.d0
      alpha=1
      beta=2
      if (iuhf.eq.0) then
        factor=2.d0
        beta=1
      endif

c Make sure that all integral files are present
      nunit=3
      if (c1symmet.eq.1) nunit=1
      unit(1)=iiiiio
      unit(2)=ijijio
      unit(3)=iijjio
      files(1)=iiiifil
      files(2)=ijijfil
      files(3)=iijjfil
      do 10 iunit=1,nunit
        inquire(file=files(iunit),exist=exist)
        if (.not.exist) then
          write (*,9200) files(iunit)
          call errex
        endif
   10 continue
 9200 format('@EXCHANGE-F, file ',a12,' not found')

c fill alpha density array for contraction with one electron integrals
      icnt=0
      roff=0
      do 710 isym=1,nirrep
        do 712 i=1,numbasir(isym)
          do 715 j=1,i
            icnt=icnt+1
            do 711 ispin=1,iuhf+1
              dens1(icnt,ispin)=dens(roff+j,roff+i,ispin)
              if (i.ne.j) dens1(icnt,ispin)=2.d0*dens1(icnt,ispin)
  711       continue
  715     continue
  712   continue
        roff=roff+numbasir(isym)
  710 continue

c multiply off diagonal elements of the alpha density matrix by two
c for coulomb contraction with two electron integrals

      do 55 ispin=1,iuhf+1
        do 50 i=1,nbastot
          do 40 j=1,nbastot
            cdens(j,i,ispin)=2.d0*dens(j,i,ispin)
   40     continue
          cdens(i,i,ispin)=dens(i,i,ispin)
   50   continue
   55 continue

c read integrals and contract them with the appropriate density
c matrix for one-electron, coulomb and exchange energies

      open(iiiiio,file=iiiifil,status='old',form='unformatted',
     &    access='sequential')

c one electron energy

      ehf1e=0.d0
      call locate(iiiiio,'ONEHAMIL')
   60 read(iiiiio)buf,ibuf,nut
      do 80 ispin=1,iuhf+1
        do 70 int=1,nut
          ehf1e=ehf1e+dens1(ibuf(int),ispin)*buf(int)
   70   continue
   80 continue
      if (nut.ne.-1) goto 60
      ehf1e=ehf1e*factor

c kinetic energy

      ehfkin=0.d0
      call locate(iiiiio,'KINETINT')
   90 read(iiiiio)buf,ibuf,nut
      do 110 ispin=1,iuhf+1
        do 100 int=1,nut
          ehfkin=ehfkin+dens1(ibuf(int),ispin)*buf(int)
  100   continue
  110 continue
      if (nut.ne.-1) goto 90
      ehfkin=ehfkin*factor

      close(iiiiio,status='keep')

c two electron energy

      ecoulab=0.d0
      ecoulaa=0.d0
      ecoulbb=0.d0
      eexchaa=0.d0
      eexchbb=0.d0

      do 400 iunit=1,nunit

        open(unit(iunit),file=files(iunit),status='old',
     &      form='unformatted',access='sequential')
        call locate(unit(iunit),'TWOELSUP')
    1   read(unit(iunit))buf,ibuf,nut
        do 7 int=1,nut
          i=iupki(ibuf(int))
          j=iupkj(ibuf(int))
          k=iupkk(ibuf(int))
          l=iupkl(ibuf(int))
          x=buf(int)

          tmp1a=0.d0
          tmp2a=0.d0
          tmp2b=0.d0
          tmp3a=0.d0
          tmp3b=0.d0

c coulomb energy
          if (iunit.eq.1 .or. iunit.eq.3) then
c           alpha beta
            tmp1a=x*cdens(j,i,alpha)*cdens(l,k,beta)
            if (i.ne.k .or. j.ne.l)
     &          tmp1a=tmp1a+x*cdens(l,k,alpha)*cdens(j,i,beta)
c           alpha alpha
            tmp2a=x*cdens(j,i,alpha)*cdens(l,k,alpha)
            tmp2b=x*cdens(j,i,beta)*cdens(l,k,beta)

c           off diagonal elements
            if (i.ne.k .or. j.ne.l) then
              tmp2a=2.d0*tmp2a
              tmp2b=2.d0*tmp2b
            endif
          endif

c exchange energy
          if (iunit.eq.1 .or. iunit.eq.2) then
c           alpha alpha
            tmp3a=x*dens(l,i,alpha)*dens(k,j,alpha)
            tmp3b=x*dens(l,i,beta)*dens(k,j,beta)

            if (i.ne.j .or. k.ne.l) then
              tmp3a=tmp3a+x*dens(k,i,alpha)*dens(l,j,alpha)
              tmp3b=tmp3b+x*dens(k,i,beta)*dens(l,j,beta)
            endif

            if (i.ne.j .and. k.ne.l) then
              tmp3a=2.d0*tmp3a
              tmp3b=2.d0*tmp3b
            endif

c           off diagonal elements
            if (i.ne.k .or. j.ne.l) then
              tmp3a=2.d0*tmp3a
              tmp3b=2.d0*tmp3b
            endif
          endif

          ecoulab=ecoulab+tmp1a
          ecoulaa=ecoulaa+tmp2a
          ecoulbb=ecoulbb+tmp2b
          eexchaa=eexchaa+tmp3a
          eexchbb=eexchbb+tmp3b

    7   continue
        if (nut.ne.-1) goto 1
        close(unit(iunit),status='keep')
  400 continue

      ehfcoul=ecoulab+0.5d0*(ecoulaa+ecoulbb)
      ehfx=-0.5d0*(eexchaa+eexchbb)

c if finite field then nuclear dipole correction
      fnuc=0.d0
      if(iflags(23).ne.0.or.iflags(24).ne.0.or.iflags(25).ne.0) then
        open(vpoutio,file=vpoutfil,form='unformatted',status='old')
        rewind(vpoutio)
c ----------------------------------------------------------------------
c makmap()
        length=(nbastot*(nbastot+1))/2
        call izero(imap,length)
        ibot=1
        iposp=0
        do irrep=1,nirrep
           nbasir=numbasir(irrep)
           do indx1=ibot,ibot+nbasir-1
              do indx2=ibot,indx1
                 iposp=iposp+1
                 iposf=indxt(indx2,indx1)
                 imap(iposf)=iposp
              end do
           end do
           ibot=ibot+nbasir
        end do
c ----------------------------------------------------------------------

  312   i23=iflags(23)
        i24=iflags(24)
        i25=iflags(25)
        if(iflags(23).ne.0)then
          label='     X  '
          ind=iflags(23)
          iflags(23)=0
        elseif(iflags(24).ne.0)then
          label='     Y  '
          ind=iflags(24)
          iflags(24)=0
        elseif(iflags(25).ne.0)then
          label='     Z  '
          ind=iflags(25)
          iflags(25)=0
        else
          goto 23
        endif
        field=dble(ind)*1.d-6
c ----------------------------------------------------------------------
c seeklb()
        rewind(vpoutio)
        title=' '
        do while (title.ne.label)
           read(vpoutio,end=120,err=121)dum1,dum2,dum3,title
           goto 122
  120      write(*,*) '@EXCHANGE: end-of-file on unit ',vpoutio
           call errex
  121      write(*,*) '@EXCHANGE: error reading from unit ',vpoutio
           call errex
  122      continue
        end do
c ----------------------------------------------------------------------
        backspace(vpoutio)
        read(vpoutio)junk,dipnuc
        fnuc=fnuc+field*dipnuc

  123   read(vpoutio)buf,ibuf,nut
        if (nut.eq.-1) goto 312
        do 313 index=1,nut
          roff=imap(ibuf(index))
          if (roff.ne.0)
     &        ehf1e=ehf1e+(dens1(roff,alpha)+
     &        dens1(roff,beta))*buf(index)*field
  313   continue
        goto 123

   23   close(vpoutio,status='keep')
        iflags(23)=i23
        iflags(24)=i24
        iflags(25)=i25
      endif

      ehfnatr=ehf1e-ehfkin
      escf=ehf1e+ehfcoul+ehfx+nucrep+fnuc
      if (int_ks) then

c     escf is the expectation value of H w.r.t. tke KS
c          single determinant.

           call putrec(1,'JOBARC','KSSCFENG',iintfp,escf)
      end if

      ehar=escf-ehfx

      if (int_printscf.gt.int_printlev) then
        if(iflg.eq.0)then
          write(*,9000) ehfnatr
          write(*,9010) ehfkin
          write(*,9030) ehfcoul
          write(*,9060) ehfx
          write(*,9020) ehf1e
          write(*,9080) escf
        else
          write(*,9100) ehfnatr
          write(*,9110) ehfkin
          write(*,9120) ehf1e
          write(*,9130) ehfcoul
          write(*,9160) ehfx
          write(*,9180) escf
        endif
      endif

 9000 format(/'The SCF nuclear-electron attraction energy =',f20.12)
 9010 format('The SCF kinetic energy                     =',f20.12)
 9020 format('The SCF one electron energy                =',f20.12)
 9030 format('The SCF coulomb energy                     =',f20.12)
 9060 format('The SCF exchange energy                    =',f20.12)
 9080 format('The SCF total energy                       =',f20.12)

 9100 format(/'The correlated nuclear-electron attraction energy =',
     &    f20.12)
 9110 format('The correlated kinetic energy                     =',
     &    f20.12)
 9120 format('The correlated one electron energy                =',
     &    f20.12)
 9130 format('The correlated coulomb energy                     =',
     &    f20.12)
 9160 format('The correlated exchange energy                    =',
     &    f20.12)
 9180 format('The correlated density total SCF energy           =',
     &    f20.12)

      call callstack_pop
      return
      end
