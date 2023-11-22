c build the Fock matrix (-ces) by calling GAMINT for each possible shell
c quadruplet and combining the integrals with the density matrix (-ces)

      subroutine SCFINT 
     >      (rhf,focka,fockb,da,db,tab,imap,nbast,nsize,xints,dsh)
      implicit double precision (a-h,o-z)
#include "gamess.par"
      PARAMETER (MXDD=16*MXGSH*MXGSH)  

      common /GMSfast/ ghondo,i,j,k,l,
     >    Lish,Ljsh,Lksh,Llsh, iandj,kandl,same, ijind,klind
      
      dimension focka(*),da(*),imap(nbast,*)
      dimension fockb(*),db(*)
      dimension xints(*),dsh(*)

#undef _PARALLEL_ACES_COM_
#include "parallel_aces.com"
      COMMON /xNSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     >                CF(MXGTOT),CG(MXGTOT),
     >                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     >                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      COMMON /xSHLT  / TOL,CUTOFF,OUT,ICOUNT
      common /GMSorder/ coeint(MXAO),iorder(MXAO),
     >               nAfrst(MXSH),nAstep(MXSH),sphrcl
      common /GMSdirect/ dirscf,fdiff,schwrz,pople2,i011,i012,i013
      dimension ghondo(50625),ddij(MXDD)
      dimension ijind(225),klind(225) 
      dimension nbai(15),nbaj(15),nbak(15),nbal(15)
      dimension nbaiA(15),nbajA(15),nbakA(15),nbalA(15)
      logical OUT
      logical rhf,schskp,dirscf,pople2,schwrz,fdiff
      logical skipa,skipb,skipc,done,pople,iandj,kandl,same
      logical first,sphrcl
      data fourth,half /0.25d0,0.5d0/
      data first /.true./

      if (irank.eq.0) write (99,*) 'SCFINT begin'
      CUTSV  = CUTOFF
      CUTOFF = MIN(CUTOFF,1.0D-10)
      IF(FDIFF) CUTOFF=CUTOFF*half

      cutint=CUTOFF
      nquad=0
      nschwz=0
      nnotdn=0

      nloop=-1
      do 10 ii=1,nshell
      do 20 jj=1,ii
      do 30 kk=1,jj
      do 40 ll=1,kk

         nloop=nloop+1
         nloop=mod(nloop,nprocs)
         if (nloop.ne.irank) goto 40

         skipa=jj.eq.kk
         skipb=(ii.eq.kk).or.(jj.eq.ll)
         skipc=(ii.eq.jj).or.(kk.eq.ll)

         iexch=1
         ish=ii
         jsh=jj
         ksh=kk
         lsh=ll
         goto 780
 740     if (skipa) goto 760
         iexch=2
         ish=ii
         jsh=kk
         ksh=jj
         lsh=ll
         goto 780
 760     if (skipb.or.skipc) goto 840
         iexch=3
         ish=ii
         jsh=ll
         ksh=jj
         lsh=kk
 780     continue

         qq4=1.d0
         pople=pople2
         nquad=nquad+1

         schskp=.false.
      IF(SCHWRZ) THEN
         IJIJ = (ISH*ISH-ISH)/2 + JSH
         KLKL = (KSH*KSH-KSH)/2 + LSH
         TEST = QQ4*XINTS(IJIJ)*XINTS(KLKL)
         IF(DIRSCF) THEN
            DENMAX = xSCHWDN(DSH,ISH,JSH,KSH,LSH)
            TEST = TEST*DENMAX
         END IF
         SCHSKP = TEST.LT.CUTOFF
c         write (99,'(a,4i3)') 'ijkl:',ish,jsh,ksh,lsh
c         write (99,'(2g20.12,l2)') denmax,test,schskp
         IF(SCHSKP) NSCHWZ = NSCHWZ + 1
      END IF
      if (schskp) goto 999

c         write (99,'(a,4i3,i6)') 'calling GAMINT:',ish,jsh,ksh,lsh,iexch
         call GAMINT (ish,jsh,ksh,lsh,iexch,qq4, done,
     >                ghondo,ddij,pople, .false.,.false.,
     >                iandj,kandl,same,nofi,nofj,nofk,nofl,
     >                ijind,klind,nbai,nbaj,nbak,nbal)
         if (.not.done) then
            nnotdn=nnotdn+1
            goto 999
         end if
c         write (99,'(a,4i3)')    'after GAMINT:  ',ish,jsh,ksh,lsh

         if (SPHRCL) then
            call GAMdfg (ghondo,nofi,nofj,nofk,nofl,
     >           nbai,nbaj,nbak,nbal, iandj,kandl,same, ijind,klind)
            call shftbs (ish,jsh,ksh,lsh, nofi,nofj,nofk,nofl,
     >           nbaiA,nbajA,nbakA,nbalA)
            Lish=(nofi-1)/2
            Ljsh=(nofj-1)/2
            Lksh=(nofk-1)/2
            Llsh=(nofl-1)/2
         end if

c------------ loops over basis functions in four shells ---------------
         ijn=0
         jmax=nofj
         do 360 i=1,nofi
           if (iandj) jmax=i
           do 340 j=1,jmax
             ijn=ijn+1
             if (pople) then
               print *, '@SCFINT: POPLE type not allowed'
               call aces_exit(1)
             else
               n1=ijind(ijn)
             end if

             kln=0
             lmax=nofl
             do 320 k=1,nofk
                if (kandl) lmax=k
                do 300 l=1,lmax
                   kln=kln+1
                   if (same.and.kln.gt.ijn) goto 340
c       if (ish.eq.4.and.jsh.eq.4.and.ksh.eq.2.and.lsh.eq.4) then
c          write (99,'(a,4i3,i6,g25.15)') 'orig.indices,ind,val',
c     >    i,j,k,l,n1+klind(kln),ghondo(n1+klind(kln))
c       end if   

                   if (SPHRCL) then

c                      call GAMsph (ghondo, i,j,k,l, val,
c     >             Lish,Ljsh,Lksh,Llsh, iandj,kandl,same, ijind,klind)
                      call GAMsph (val)

                      ibas=nbaiA(i)
                      jbas=nbajA(j)
                      kbas=nbakA(k)
                      lbas=nbalA(l)
c	write (99,'(4i3,g25.15)') ibas,jbas,kbas,lbas,val
                   else
                      if (pople) then
                         print *, '@SCFINT: POPLE type not allowed'
                         call aces_exit(1)
                      else
                         nn=n1+klind(kln)
                         val=ghondo(nn)
                      end if
c basis function numbers (according to GAMESS)
                      ibasG=nbai(i)
                      jbasG=nbaj(j)
                      kbasG=nbak(k)
                      lbasG=nbal(l)
c some components of d,f,and g must be multiplied by factors <> 1
c to account for differences between their ACES and GAMESS normalizations
        val=val*coeint(ibasG)*coeint(jbasG)*coeint(kbasG)*coeint(lbasG)

c basis function numbers (according to ACES)
                      ibas=iorder(ibasG)
                      jbas=iorder(jbasG)
                      kbas=iorder(kbasG)
                      lbas=iorder(lbasG)
                   end if

                   if (abs(val).lt.cutint) goto 300

                  if(ibas.eq.jbas) val = val*half
                  if(kbas.eq.lbas) val = val*half
                  if(ibas.eq.kbas .and. jbas.eq.lbas) val = val*half

                  nij = imap(ibas,jbas)
                  nkl = imap(kbas,lbas)
                  nik = imap(ibas,kbas)
                  nil = imap(ibas,lbas)
                  njk = imap(jbas,kbas)
                  njl = imap(jbas,lbas)

                  if (rhf) then
                     focka(nij) = focka(nij)+ val*da(nkl)
                     focka(nkl) = focka(nkl)+ val*da(nij)
                     focka(nik) = focka(nik)- val*da(njl)*fourth
                     focka(nil) = focka(nil)- val*da(njk)*fourth
                     focka(njk) = focka(njk)- val*da(nil)*fourth
                     focka(njl) = focka(njl)- val*da(nik)*fourth
                  else
                     dum=val*(da(nkl)+db(nkl))
                     focka(nij) = focka(nij)+ dum
                     fockb(nij) = fockb(nij)+ dum
                     dum=val*(da(nij)+db(nij))
                     focka(nkl) = focka(nkl)+ dum
                     fockb(nkl) = fockb(nkl)+ dum
                     focka(nik) = focka(nik)- val*da(njl)*half
                     fockb(nik) = fockb(nik)- val*db(njl)*half
                     focka(nil) = focka(nil)- val*da(njk)*half
                     fockb(nil) = fockb(nil)- val*db(njk)*half
                     focka(njk) = focka(njk)- val*da(nil)*half
                     fockb(njk) = fockb(njk)- val*db(nil)*half
                     focka(njl) = focka(njl)- val*da(nik)*half
                     fockb(njl) = fockb(njl)- val*db(nik)*half
                  end if
 300  continue       
 320  continue       
 340  continue       
 360  continue       
c------------ loops over basis functions in four shells (END)


 999     continue
         goto (740,760,840), iexch
 840     continue
 40   continue
 30   continue
 20   continue
 10   continue

      call MPIAR (focka,tab,nsize)
      if (.not.rhf) call MPIAR (fockb,tab,nsize)
      call MPIAR1i (nschwz,nsctot)
      call MPIAR1i (nnotdn,nnotot)
      if (first) then
         first=.false.
         call MPIAR1i (nquad,nqutot)
         if (irank.eq.0) 
     >     write (99,'(a,i10)') 'Total shell quadruplets:',nqutot
      end if   
      if (irank.eq.0) then
         write (99,'(a,i10)') 'Skipped by Schwarz:     ',nsctot
         write (99,'(a,i10)') 'Skipped in GAMINT:      ',nnotot
         write (99,*) 'SCFINT end'
      end if   
      CUTOFF = CUTSV
      end
c--------------------------------------------------------------------------
c info transfer ACES -> GAMES;
c ----------------------------
c NATS, ATPOS(3,nats): no. and coordinates of atoms
c NPRIMS,EXPS(nprims),COEx(nprims): no.of primitives, their exponents 
c   and contraction coeffs in x-subshell;
c   all components resulting from rotations (e.g. px,py,pz) are counted
c   as ONE primitive;
c NSHLS: no. of shells;
c KFRST(i), KLAST(i): numbers of first and last primitive in the i-th shell;
c KATOMS(i): atom where the i-th shell is placed;
c KANGU(i): 0=s, 1=p, 2=d, 3=f, 4=g, -1=sp, -2=spd;
c KORDER(i)=j means that GAMESS's i-th basis function is j-th for ACES
c CUTOFF: integrals with abs(value) < cutoff will be skipped
c IOUT: I/O no. for log messages

      subroutine GAMINI (sphrc,nats,atpos,  nprims,exps,
     >           coes,coep,coed,coef,coeg,
     >           nshls,kfrst,klast,katoms,kangu,
     >           korder,doeint,nfrst,nstep,cutoff,iout )
      implicit double precision (a-h,o-z)
      dimension atpos(3,*)
      dimension exps(*),coes(*),coep(*),coed(*),coef(*),coeg(*)
      dimension kfrst(*),klast(*),katoms(*),kangu(*),korder(*)
      dimension doeint(*),nfrst(*),nstep(*)

#undef _PARALLEL_ACES_COM_
#include "parallel_aces.com"
c the commons below contains variables needed by the GAMESS integrals
#include "gamess.par"
      COMMON /xINFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     >                ZAN(MXATM),C(3,MXATM)
      COMMON /xIOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /xNSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     >                CF(MXGTOT),CG(MXGTOT),
     >                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     >                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      COMMON /xOUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK
      COMMON /xPAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /xPKFIL / PK,PANDK,BLOCK
      COMMON /xRESTAR/ TIMLIM,IREST,NREC,INTLOC,IST,JST,KST,LST
      COMMON /xSYMTRY/ MAPSHL(MXSH,48),MAPCTR(MXATM,48),
     >                T(432),INVT(48),NT
      common /GMSorder/ coeint(MXAO),iorder(MXAO),
     >               nAfrst(MXSH),nAstep(MXSH),sphrcl
      common /zmins/ zminx(mxsh)
      logical sphrc,sphrcl
      LOGICAL GOPARR,DSKWRK,MASWRK
      LOGICAL PK,PANDK,BLOCK

      dimension dum(1),idum(1)
c different quantities used in integrals routines, depending on the
c shells' type (angular momentum)
      dimension ktypes(-2:4),kmins(-2:4),kmaxs(-2:4),nbfs(-2:4)
      data ktypes /3,2, 1,2,3,4,5/
      data kmins /1,1, 1,2,5,11,21/
      data kmaxs /10,4, 1,4,10,20,35/
      data nbfs /10,4, 1,3,6,10,15/

      sphrcl=sphrc
      NAT=nats
      do 10 i=1,nats
      do 10 k=1,3
      C(k,i)=atpos(k,i)
 10   continue

      if (irank.eq.0) then
         write (99,*) 'primitives in GAMINI'
         write (99,'(1x,a,1x,a,1x,5(3x,a,3x))') 
     >     'No.','exponent','s-coef','p-coef','d-coef','f-coef','g-coef'
      end if
      do 20 i=1,nprims
      if (irank.eq.0) write (99,'(i4,6g12.4)') 
     >      i,exps(i),coes(i),coep(i),coed(i),coef(i),coeg(i)
      EX(i)=exps(i)
      CS(i)=coes(i)
      CP(i)=coep(i)
      CD(i)=coed(i)
      CF(i)=coef(i)
      CG(i)=coeg(i)
 20   continue
      if (irank.eq.0) write (99,*) 

      nbasis=0
      if (irank.eq.0) write (99,*) 'shells in GAMINI'
      do 30 i=1,nshls
      zmini=1.d30
      do 31 ig=kfrst(i),klast(i)
 31   zmini=min(zmini,ex(ig))
      zminx(i)=zmini
      KSTART(i)=kfrst(i)
      KNG(i)=klast(i)-kfrst(i)+1
      KATOM(i)=katoms(i)
      k=kangu(i)
      KMIN(i)=kmins(k)
      KMAX(i)=kmaxs(k)
      KTYPE(i)=ktypes(k)
      KLOC(i)=nbasis+1
      MAPSHL(i,1)=i
      nbasis=nbasis+nbfs(k)
      if (irank.eq.0) write (99,'(i3,a,i2,a,i1,a,i4,a,i4,a,i5,a,i5)') 
     > i,' katom=',katom(i),'  ktype=',ktype(i),'  orbs:',
     > nbasis-nbfs(k)+1,' -',nbasis, '  prims:',kfrst(i),' -',klast(i)
 30   continue
      if (irank.eq.0) write (99,*)
      
      if (irank.eq.0) write (99,'(a)') ' GAMESS->ACES orbital mapping'
      do 40 i=1,nbasis
      if (irank.eq.0) write (99,'(a,i4,a,i4,a,g23.16)') 
     >  'GAMESS:',i,'  ACES:',korder(i), '  coef:',doeint(i)
      coeint(i)=doeint(i)
 40   iorder(i)=korder(i)

      if (irank.eq.0) then
         write (99,*)
         if (sphrcl) then
            write (99,*) ' Spherical harmonics will be used'
         else
            write (99,*) ' Cartesian Gaussians will be used'
         end if
      end if

      do 50 i=1,nshls
         nAfrst(i)=nfrst(i)
         nAstep(i)=nstep(i)
 50   continue

      NSHELL=nshls
      NUM=nbasis
      IW=iout
      TIMLIM=1.d99 ! no time limit
      IREST=0 ! normal run, no restart
      NT=1  ! so far, only C1 symmetry
      ME=0
      NPROC=1
      IBTYP=0  ! nonzero means dynamic load balancing
      MASWRK=.true. ! I am the master node
      GOPARR=.false. ! no parallelization
      PK=.false.
      NPRINT=7
      ITOL=20  ! primitive products with exp less than 10**(-ITOL) skipped
      ICUT=-log10(cutoff)
      NORMF=0 ! normalize the basis functions
      NORMP=0 ! contraction coeffs are wrt normalized primitives
      INTG76=1 ! always call xgamgen from xdebut (not always needed
c                but doesn't require an additional parameter to gamini)
      call xDEBUT(INTG76,.true.,dum,dum,idum,idum,idum,.false.)
      end
c--------------------------------------------------------------------------
c establish mapping between the standard order of orbitals in GAMESS and ACES
c and correct the difference between norm.constants of d,f,g orbitals
c in ACES and GAMESS
c (for a fragment of basis set: single atom, one value of L)
c
c NBAS1: No. of first orbital in this fragment
c LCURR: current value of the angular momentum L (0..4)
c NSH: No. of shells in this fragment
c IORDER(i)=j means that GAMESS's i-th basis function is j-th for ACES
c (i,j number the whole basis set, not only current fragment)
c COEINT(i)=x: GAMESS's i-th basis function must be multiplied by x
c so that two-electron integrals are consistent with ACES's one-electron
c
c Spherical harmonic basis set uses a different orbital mapping:
c Orbital No. i of the shell ISH (see routine GAMsph to see standard
c ACES orbital order within a P,D,F, and G shell) has the global number
c NAFRST(ISH) + (i-1)*NASTEP(ISH) 

      subroutine SETORD (nbas1,Lcurr,nsh,iorder,coeint,nAfrst,nAstep)
      implicit double precision (a-h,o-z)
      save nshls,nskip
      parameter (zero=0.d0, one=1.d0)
      parameter (s3=1.732050807568877d0, s9=3.d0)
      parameter (s15=3.872983346207417d0, s105=10.24695076595960d0)
      dimension iorder(*),coeint(*),nAfrst(*),nAstep(*)

c Lorder(i,L)=j means that in shells with ang.mom.=L the i-th orbital
c for GAMESS is the j-th orbital for ACES
      dimension Lorder(15,0:4)
      data Lorder /1,14*0,  1,2,3,12*0,  1,4,6,2,3,5,9*0,
     > 1,7,10,2,3,4,8,6,9,5,5*0,  1,11,15,2,3,7,12,10,14,4,6,13,5,8,9/
      data nshls,nskip / 0,0 /

      double precision Lcoef(15,0:4)
      data Lcoef /one,14*zero, 3*one,12*zero, 3*s3,3*one,9*zero,
     > 3*s15,6*s3,one,5*zero, 3*s105,6*s15,3*s9,3*s3/

      NinL=((Lcurr+1)*(Lcurr+2))/2
      NinLs=2*Lcurr+1
      iG=nbas1
      do 10 ish=1,nsh
         ishtot=nshls+ish
         do 20 N=1,NinL
            iA=(nbas1-1) +(Lorder(N,Lcurr)-1)*nsh +ish
            iorder(iG)=iA
            coeint(iG)=Lcoef(N,Lcurr)
            iG=iG+1
            if (N.eq.1) then
               nAfrst(ishtot)=iA-nskip
               nAstep(ishtot)=nsh
            end if
 20      continue
 10   continue

      nshls=nshls +nsh
      nskip=nskip +nsh*(NinL-NinLs)
c nskip keeps track of the number by which the basis set was reduced
c because of the Cartesian--->spherical remapping
      end
c--------------------------------------------------------------------------
c compress density to shell pairs, for spherical harmonic shells;
c this version works for S,P,D,F,G shells (ktype=1..5), not for L or M
c (which are obsolete anyway);
c called from GAMESS's routine SHLDEN if SPHRCL=.true.

      subroutine SHDsph (alpha,beta,DA,DB,DSH,imap,nbas)
      implicit double precision (a-h,o-z)
#include "gamess.par"
      COMMON /xNSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     *                CF(MXGTOT),CG(MXGTOT),
     *                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     *                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      logical alpha,beta
      dimension da(*),db(*),dsh(*),imap(nbas,nbas)
      common /GMSorder/ coeint(MXAO),iorder(MXAO),
     >               nAfrst(MXSH),nAstep(MXSH),sphrcl
      logical sphrcl
      data zero / 0.d0 /

      ijsh=0
      do 240 ish=1,nshell
         i1=nAfrst(ish)
         istep=nAstep(ish)
         Li=ktype(ish)-1
         in=i1 +2*Li*istep
         do 230 jsh=1,ish
            j1=nAfrst(jsh)
            jstep=nAstep(jsh)
            Lj=ktype(jsh)-1
            jn=j1 +2*Lj*jstep
            ijsh=ijsh+1
            dmax=zero
            do 130 i=i1,in,istep
               if (ish.eq.jsh) jn=i
               do 120 j=j1,jn,jstep
                  ij=imap(i,j)
                  if(alpha) then
                     if(abs(da(ij)).gt.dmax) dmax = abs(da(ij))
                  end if
                  if(beta) then
                     if(abs(db(ij)).gt.dmax) dmax = abs(db(ij))
                  end if
 120           continue
 130        continue
            dsh(ijsh) = dmax
c            write (99,'(2i3,g20.12)') ish,jsh,dmax
 230     continue
 240  continue
      end
c-----------------------------------------------------------------
c calculate exchange integrals (ij|ij) for later use in Schwarz
c estimations; this is the routine EXCHNG from GAMESS, modified to use
c the integral front-end GAMINT and generalized for sph.harm. basis sets

      subroutine GAMexc (xints,intg76)
      implicit double precision(a-h,o-z)
#include "gamess.par"
      PARAMETER (MXDD=16*MXGSH*MXGSH)
      dimension xints(*),ghondo(50625),ddij(MXDD)

      common /GMSfast/ ghondo,i,j,k,l,
     >    Lish,Ljsh,Lksh,Llsh, iandj,kandl,same, ijind,klind

      common /GMSorder/ coeint(MXAO),iorder(MXAO),
     >               nAfrst(MXSH),nAstep(MXSH),sphrcl
      logical sphrcl,done,iandj,kandl,same,pople,OUT
      dimension ijind(225),klind(225)
      dimension nbai(15),nbaj(15),nbak(15),nbal(15)
      dimension nbaiA(15),nbajA(15),nbakA(15),nbalA(15)

      COMMON /xNSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     >                CF(MXGTOT),CG(MXGTOT),
     >                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     >                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      COMMON /xSHLT  / TOL,CUTOFF,OUT,ICOUNT

      data zero,one / 0.d0,1.d0 /

      nsh2=(nshell*(nshell+1))/2
      tolsv = tol
      tol = 75.0d+00
      iexch = 1
      qq4   = one
      nint  = 0
      ijij = 0
      do 600 ishx = 1,nshell
         do 500 jshx = 1,ishx
c	write (99,'(a,2i6)') 'ish,jsh:',ishx,jshx
            ish=ishx
            ksh=ishx
            jsh=jshx
            lsh=jshx
            ijij = ijij+1
c            POPLE=.TRUE.
c            IF(INTG76.EQ.0) POPLE=.FALSE.
c            IF(KTYPE(ISH).GT.2) POPLE=.FALSE.
c            IF(KTYPE(JSH).GT.2) POPLE=.FALSE.
            call GAMINT (ish,jsh,ksh,lsh,iexch,qq4, done,
     >           ghondo,ddij,pople, .false.,.false.,
     >           iandj,kandl,same,nofi,nofj,nofk,nofl,
     >           ijind,klind,nbai,nbaj,nbak,nbal)
         if (SPHRCL) then
            call GAMdfg (ghondo,nofi,nofj,nofk,nofl,
     >           nbai,nbaj,nbak,nbal, iandj,kandl,same, ijind,klind)
            call shftbs (ish,jsh,ksh,lsh, nofi,nofj,nofk,nofl,
     >           nbaiA,nbajA,nbakA,nbalA)
            Lish=(nofi-1)/2
            Ljsh=(nofj-1)/2
         end if

            vmax = zero

            ijn = 0
            jmax=nofj
            do 300 i=1,nofi
               if(iandj) jmax = i
               do 200 j=1,jmax
                     ijn = ijn+1
                   if (SPHRCL) then

                      k=i
                      l=j
                      Lksh=Lish
                      Llsh=Ljsh
c                      call GAMsph (ghondo, i,j,i,j, val,
c     >             Lish,Ljsh,Lish,Ljsh, iandj,kandl,same, ijind,klind)
                      call GAMsph (val)

                   else
                      nn = ijind(ijn) + klind(ijn)
                      val = ghondo(nn)
c	write (99,'(a,2i6,g25.15)') 'i,j,val:',i,j,val
                      ibasG=nbai(i)
                      jbasG=nbaj(j)
c     account for differences in norm.constants (see subr.GAMINT)
      val=val*coeint(ibasG)*coeint(ibasG)*coeint(jbasG)*coeint(jbasG)
c      write (99,'(a,2i6,g25.15)') 'ibas,jbas,newval',ibasG,jbasG,val
                   end if

                  IF(VAL.GT.ZERO) NINT=NINT+1
                  IF(VAL.GT.VMAX) VMAX=VAL
 200           CONTINUE
 300        CONTINUE
            XINTS(IJIJ)=VMAX
 500     CONTINUE
 600  CONTINUE
      DO 700 I=1,NSH2
         XINTS(I) = SQRT(XINTS(I))
 700  CONTINUE
      TOL = TOLSV

      RETURN
      END
c--------------------------------------------------------------------
      subroutine mkimap (nirrep,nbast,nbas,imap)
      dimension nbas(8),imap(nbast,nbast)
      INDX(I,J)=J+(I*(I-1))/2

         CALL IZERO(IMAP,NBAST*NBAST)
         IOFF=0
         IOFFT=0
         DO IRREP=1,NIRREP
            DO I=1,NBAS(IRREP)
               DO J=1,NBAS(IRREP)
                  INDS=INDX(MAX(I,J),MIN(I,J))+IOFFT
                  IMAP(I+IOFF,J+IOFF)=INDS
               END DO
            END DO
            IOFF=IOFF+NBAS(IRREP)
            IOFFT=IOFFT+NBAS(IRREP)*(NBAS(IRREP)+1)/2
         END DO
         end
c----------------------------------------------------------------------
      subroutine matadd(from,cf,to,ct,n)
      implicit double precision (a-h,o-z)
      dimension from(*),to(*)
      do 10 i=1,n
 10   to(i)=ct*to(i) +cf*from(i)
      end
c---------------------------------------------------------------------
      subroutine shlask (nsh)
      implicit double precision (a-h,o-z)
#include "gamess.par"
      COMMON /xNSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     >                CF(MXGTOT),CG(MXGTOT),
     >                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     >                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      nsh=NSHELL
      end
c---------------------------------------------------------------------
      subroutine dump (a)
      real*8 a(*)
#undef _PARALLEL_ACES_COM_
#include "parallel_aces.com"
      if (irank.eq.0) write (*,'(3g20.10)') (a(i),i=1,3)
      end
c---------------------------------------------------------------------------
      subroutine getnp (np)
      integer np
#undef _PARALLEL_ACES_COM_
#include "parallel_aces.com"
      np=nprocs
      end
c---------------------------------------------------------------------------
      subroutine nr2asc (nr,str)
      character*3 str
      character*1 digits(0:9)

      data digits /'0','1','2','3','4','5','6','7','8','9'/

      if (nr.gt.999) then
         str='NaN'
         return
      end if
      i=nr/100
      j=(nr-100*i)/10
      k=nr-100*i-10*j

      str(1:1)=digits(i)
      str(2:2)=digits(j)
      str(3:3)=digits(k)

c      str=digits(k)
c      if (i.ne.0.or.j.ne.0) str=digits(j)//str
c      if (i.ne.0) str=digits(i)//str

      end
