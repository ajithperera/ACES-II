      subroutine A3_Ecp_basis_grinit(Natoms, Nshells, Namat, Iqmstr,
     &                       Jcostr,
     &                       Nucstr, Nrcstr, Jstrt, Jrs, Expa, Nalpha,
     &                       Cont, Npcoef, Cord, Charge, Nharm)
c  
c-----------------------------------------------------------------------      
c     subroutine doit fills the common-blocks for subroutine paloop
c     input:
c     namat : names of the atoms 
c     iqmstr: max. l-quantum-number+1 of a shell
c     jcostr: number of shells
c     nucstr: number of the uncontracted functions of the shells > NUCO
c             (exponents)
c     nrcstr: number of blocks(general contracted) > NRCO
c             (contractions)
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
C-----------------------------------------------------------------------
C#include "ecp.par"
C#include "ecp.com"
#include "ecplims.par"
#include "cvwecp.h"
C-----------------------------------------------------------------------
      PARAMETER (LUCMD = 5, LUPRI = 6)
      PARAMETER (MXCENT = 200, MXCOOR = 3*MXCENT)
      PARAMETER (MXQN=5, MXAQN=MXQN*(MXQN+1)/2, MXAQNS=MXAQN*MXAQN)
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, D2 = 2.0D0, D4 = 4.0D0,
     *           DP25 = 0.25D0, DP5 = 0.5D0, DP75 = 0.75D0,
     *           PI = 3.14159 26535 89793 D0,
     *           BOHRAD = 0.52917 706 D0,
     *           THRMIN = 1.D-15)
C-----------------------------------------------------------------------
      character*4 namat(Max_centers)
      Logical Nharm
C

      Dimension iqmstr(Max_centers),jcostr(Max_centers,Maxjco),
     &          nucstr(Max_centers,Maxang,Maxjco),
     &          nrcstr(Max_centers,Maxang,Maxjco),
     &          Jstrt(Max_shells),Jrs(Max_shells), Expa(Nalpha),
     &          Cont(Npcoef),Cord(3,Max_shells),
     &          Charge(Max_centers)
C
      Dimension
     &          ipeoff(Max_centers,Maxang,Maxjco),
     &          iccoff(Max_centers,Maxang,Maxjco),
     &          newoff(Max_centers,Maxang,Maxjco),
     &          numpri(Max_centers,Maxang,Maxjco),
     &          imin(Maxang),imax(Maxang),nfcts(Maxang),
     &          ifeld(Max_centers)
C
c-----------------------------------------------------------------------
c     namsml              : sml = small
c     gout                : intbatch
c     ipeoff              : offset for primitive exponents
c     iccoff              : offset for contraction coefficients
c     newoff              : offsets for primitives
c     numpri              : number of primitive functions of a shell
c     imin,imax           : lower(upper) offset of a cartesian component
c     nfcts               : number of types of a cartesian component 
c     ifeld               : checkarray 
c-----------------------------------------------------------------------
      common /nshel / expnt(ndi13),contr(ndi13,ndico),numcon(ndi13),
     &                katom(ndi9),ktype(ndi9),
     1                kprim(ndi9),kbfn(ndi9),kmini(ndi9),kmaxi(ndi9),
     2                nprims(ndi9),ndegen(ndi9),nshell,nbf
      common /modez / zetm(ndi9)

#ifdef _DEBUG_LVL0
      Write(6,*) 
      Write(6,*) (iqmstr(i), i=1, nucind)
      Write(6,*) ((jcostr(i,j),j=1,iqmstr(i)),i=1,nucind)
      Write(6,*) (jstrt(i), i=1, 10)
      Write(6,*) (jrs(i), i=1, 10)
#endif
      icount=0
      call izero(ifeld,ndi10)
      do 1200 iat=1,NUCIND
        do 1250 ilq=1,iqmstr(iat)
          do 1300 ijco=1,jcostr(iat,ilq)
            icount=icount+1
            if (ifeld(iat).ne.1) then
              ipeoff(iat,ilq,ijco)=jstrt(icount)
              iccoff(iat,ilq,ijco)=jrs(icount)

#ifdef _DEBUG_LVL0
      write(*,1) iat,ilq,ijco,ipeoff(iat,ilq,ijco)
   1  format('iat: ',i3,'ilq: ',i3,' ijco: ',i3,' ipeoff: ',i3)
      write(*,2) iat,ilq,ijco,iccoff(iat,ilq,ijco)
   2  format('iat: ',i3,'ilq: ',i3,' ijco: ',i3,' iccoff: ',i3)
#endif 
            endif
 1300     continue
 1250   continue
        ifeld(iat)=1 
 1200 continue 

      call zero(expnt,ndi13)
      call zero(contr,ndi13)
      call izero(numpri,ndi10*ndilmx*maxjco)
      icount=0
      kcount=0
      ncount=0
      do 1500 iat=1,NUCIND
        do 1550 ilq=1,iqmstr(iat)
          do 1575 ijco=1,jcostr(iat,ilq)
            ncount=ncount+1
            do 1600 inuco=1,nucstr(iat,ilq,ijco)
              icount=icount+1
              do 1650 inrco=1,nrcstr(iat,ilq,ijco)

CSSS                contr(icount,inrco)=PRICCF(iccoff(iat,ilq,ijco)+inuco+   
CSSS     &          nucstr(iat,ilq,ijco)*(inrco-1),inrco)

                contr(icount,inrco)=PRICCF(ipeoff(iat,ilq,ijco)+inuco,
     &                                      inrco)
#ifdef _DEBUG_LVL0
      write(*,3) icount,inrco,contr(icount,inrco)
   3  format('icount: ',i3,' inrco: ',i3,' contr: ',f12.8)
#endif     
 1650         continue

              expnt(icount)=PRIEXP(ipeoff(iat,ilq,ijco)+inuco)
              kcount=kcount+1
#ifdef _DEBUG_LVL0
      write(*,4) icount,expnt(icount)
   4  format('icount: ',i3,' expnt: ',f12.8)
#endif  

 1600       continue

            numpri(iat,ilq,ijco)=kcount
            kcount=0
            numcon(ncount)=nrcstr(iat,ilq,ijco)

#ifdef _DEBUG_LVL0
      write(*,6) ncount,numcon(ncount)
   6  format('ncount: ',i3,' numcon: ',i3)
#endif


 1575     continue
 1550   continue
 1500 continue

c     create new offsets of primitive exponents

      call izero(newoff,ndi10*ndilmx*maxjco)
      khelp=1
      do 3000 iat=1,NUCIND
        do 3100 ilq=1,iqmstr(iat)
          do 3200 ijco=1,jcostr(iat,ilq)
            newoff(iat,ilq,ijco)=khelp
            khelp=khelp+numpri(iat,ilq,ijco)
 3200     continue
 3100   continue
 3000 continue

c-----------------------------------------------------------------------
c     fill in common infoa
c-----------------------------------------------------------------------
#ifdef _NOSKIP
      call zero(xyz,3*ndi10)
c --- store only the nonsymmetryredundant coordinates to xyz
      natom=NUCIND
      do 7000 i=1,NUCIND
        do 7100 j=1,3
          xyz(j,i)=CORD(i,j)
 7100   continue
 7000 continue
#endif
c-----------------------------------------------------------------------
c     fill in common /typinf/
c-----------------------------------------------------------------------

      ltmax=1
      do 10 iat=1,NUCIND
        if (ltmax.le.iqmstr(iat)) then 
          ltmax=iqmstr(iat)
        endif  
   10 continue      
cch      nftdrv=ltmax*(ltmax+1)*(ltmax+2)/6
      ltdrv=ltmax+1

c-----------------------------------------------------------------------
c     fill in common /nshell/
c-----------------------------------------------------------------------

c     doing the k's from common nshel
      
      if (NHARM.eq.1) then
c     spherical harmonics are used!
        do 5500 l=1,ndilmx
          imin(l)=(l-1)*l*(l+1)/6+1
          imax(l)=l*(l+1)*(l+2)/6 
          nfcts(l)=2*l-1
 5500   continue
      else 
        do 5000 l=1,ndilmx
          imin(l)=(l-1)*l*(l+1)/6+1
          imax(l)=l*(l+1)*(l+2)/6 
          nfcts(l)=l*(l+1)/2
 5000   continue
      endif

      nbf=0
      nshll=0
      indsh=0
c     loop over atoms
      do 600 iat=1,NUCIND
c       loop over l-quantum-numbers
        do 800 ilq=1,iqmstr(iat)
c         loop over contracted basisfunctions
          do 900 ijco=1,jcostr(iat,ilq)
            indsh=indsh+1
CSSS            mulilq=ISTBAO(indsh)
              nshll=nshll+1
              katom(nshll)=iat
              ktype(nshll)=ilq
              kmini(nshll)=imin(ilq)
              kmaxi(nshll)=imax(ilq)
              ndegen(nshll)=nfcts(ilq)             
              kbfn(nshll)=nbf+1
CSSS              mulsh(nshll)=mulilq
              nprims(nshll)=numpri(iat,ilq,ijco)
              kprim(nshll)=newoff(iat,ilq,ijco)
              nbf=nbf+ndegen(nshll)*fmult(mulilq)
cch1000       continue
  900     continue
  800   continue
  600 continue

      nshell=nshll

#ifdef _DEBUG_LVL0
      write(*,*)
      do  i=1,nshell
         write(*,*) 'nshell=',i
         write(*,3913) i,katom(i),i,ktype(i)
3913     format('katom<',i2,'>  ',i2,' ktype<',i2,'>  ',i2)
         write(*,3914) i,kprim(i),i,kbfn(i)
3914     format('kprim<',i2,'>  ',i2,' kbfn<',i2,'>   ',i2)
         write(*,3915) i,nprims(i),i,ndegen(i)
3915     format('nprims<',i2,'> ',i2,' ndegen<',i2,'> ',i2)
         write(*,3916) i,kmini(i),i,kmaxi(i)
3916     format('kmini<',i2,'>  ',i2,' kmaxi<',i2,'>  ',i2)
         write(*,3917) i, mulsh(i)
3917     format('mulsh<',i2,'>  ',i2)
         write(*,*)
      enddo
      write(*,*) 'nbf=',nbf
      write(*,*)
#endif

c-----------------------------------------------------------------------
c     fill in common modez
c-----------------------------------------------------------------------
      do 4000 ii=1,nshell
        zetm(ii)=expnt(kprim(ii))
        do 4100 ipr=kprim(ii)+1,kprim(ii)+nprims(ii)-1
          zetm(ii)=min(expnt(ipr),zetm(ii))
 4100   continue
 4000 continue

c-----------------------------------------------------------------------
c     fill in common pairij
c-----------------------------------------------------------------------

      call setipq(ipq,ndi4)

c-----------------------------------------------------------------------
c     fill in common pairij
c-----------------------------------------------------------------------

      call gustav(ltdrv,jx,jy,jz,iadr,imasc,nftmax,ndilmx,imin,imax,
     &            nftdrv)

c-----------------------------------------------------------------------
      nshdim=nshell
      TLA=1.0d-12
      toll=-log10(TLA)+8.0
       toll=8.0    
c-----------------------------------------------------------------------
c     Ende 
c-----------------------------------------------------------------------
      return
      end 
