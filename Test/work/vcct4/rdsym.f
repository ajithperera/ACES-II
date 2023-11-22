      subroutine grdsym(n1,n2,n3,n4,nsq1,nsq2,t2,t2s,s1,s2,lb1,lb2,nlst)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      character*8 lb1,lb2
      integer s1,s2
      dimension t2(1),t2s(1),s1(nsq1),s2(nsq2)
      call getrec(20,'JOBARC',lb1,nsq1,s1)
      call getrec(20,'JOBARC',lb2,nsq2,s2)
      call rdsym(n1,n2,n3,n4,nsq1,nsq2,t2,t2s,s1,s2,nlst)
      return
      end
      subroutinegrdsym3(n1,n2,n3,n4,nsq1,nsq2,t2,t2s,s1,s2,lb1,lb2,nlst)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      character*8 lb1,lb2
      integer s1,s2
      dimension t2(1),t2s(1),s1(nsq1),s2(nsq2)
      call getrec(20,'JOBARC',lb1,nsq1,s1)
      call getrec(20,'JOBARC',lb2,nsq2,s2)
      call rdsym3(n1,n2,n3,n4,nsq1,nsq2,t2,t2s,s1,s2,nlst)
      return
      end
      SUBROUTINE rdsym(n1,n2,n3,n4,nsq1,nsq2,t2,t2s,
     *SYVEC1,SYVEC2,nlist)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION SYVEC1(1),SYVEC2(1),t2(n1,n2,n3,n4),t2s(1)
      integer syvec1,syvec2,dirprd
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),N(18)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),DIRPRD(8,8)
      no=nocco(1)
      nu=nvrto(1)
      kkk=0
      kij=0
      DO 10 IRREP=1,NIRREP
         NSYDIS=IRPDPD(IRREP,ISYTYP(2,nLIST))
         NSYDSZ=IRPDPD(IRREP,ISYTYP(1,nLIST))
         NSIZ=NSIZ+NSYDIS*NSYDSZ
         do 9 inum=1,nsydis
         kij=kij+1
         call getlst(t2s,inum,1,1,irrep,nlist)
         do 8 is=1,nsydsz
            kkk=kkk+1
            call pckindsk(kkk,3,syvec1,syvec2,nsq1,nsq2,i,j,ia,ib,nlist)
            t2(ia,ib,i,j)=t2s(is)
 8       continue
 9    CONTINUE
 10   CONTINUE
      RETURN
      END
      SUBROUTINE rdsym3(n1,n2,n3,n4,nsq1,nsq2,t2,t2s,
     *SYVEC1,SYVEC2,nlist)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION SYVEC1(1),SYVEC2(1),t2(n1,n2,n3,n4),t2s(1)
      integer syvec1,syvec2,dirprd
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),N(18)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),DIRPRD(8,8)
      no=nocco(1)
      nu=nvrto(1)
      kkk=0
      kij=0
      DO 10 IRREP=1,NIRREP
         NSYDIS=IRPDPD(IRREP,ISYTYP(2,nLIST))
         NSYDSZ=IRPDPD(IRREP,ISYTYP(1,nLIST))
         NSIZ=NSIZ+NSYDIS*NSYDSZ
         do 9 inum=1,nsydis
         kij=kij+1
         call getlst(t2s,inum,1,1,irrep,nlist)
         do 8 is=1,nsydsz
            kkk=kkk+1
         call pckindsk(kkk,3,syvec1,syvec2,nsq1,nsq2,i,j,ia,ib,nlist)
            t2(ia,ib,j,i)=t2s(is)
 8       continue
 9    CONTINUE
 10   CONTINUE
      RETURN
      END
      SUBROUTINE rdsymv(nu,t2,t2s,SYVEC1,syvec2,ifr,inr,nlist)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION SYVEC1(1),syvec2(1),t2(nu,nu),t2s(1),ifr(1),inr(1)
      integer syvec1,dirprd,a,b,c,d,pop,vrt,spop,svrt,syvec2
      logical smlv
      character*8 lb
      common/countv/irec
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),N(18)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),DIRPRD(8,8)
      common/sym/pop(8,2),vrt(8,2),ntaa(6)
      common/sumpop/spop(8),svrt(8)
      common/flags/iflags(100)
      equivalence(icllvl,iflags(2))
      smlv=icllvl.eq.11.or.icllvl.eq.13.or.icllvl.eq.14.or.
c     *icllvl.eq.15.or.icllvl.eq.16.or.icllvl.eq.17.or.icllvl.eq.19.or.
     *icllvl.eq.15.or.icllvl.eq.19.or.
     *icllvl.eq.21.or.icllvl.eq.22
      lb='SVAVB2  '
      no=nocco(1)
      nu=nvrto(1)
      nu2=nu*nu
      call getrec(20,'JOBARC',lb,nu2,syvec1)
      nfl=35
      if(smlv)then
         call mksyv(nu,syvec2)
         nfl=36
      else
         call veccop(nu2,syvec2,syvec1)
      endif
      kkk=0
      nrec=nu2/1024+1
      irec=2*nrec
      DO 10 IRREP=1,NIRREP
         NSYDIS=IRPDPD(IRREP,ISYTYP(2,nLIST))
         NSYDSZ=IRPDPD(IRREP,ISYTYP(1,nLIST))
         NSIZ=NSIZ+NSYDIS*NSYDSZ
         do 9 inum=1,nsydis
            call zeroma(t2,1,nu2)
            kij=kij+1
            call getlst(t2s,inum,1,1,irrep,nlist)
            do 8 is=1,nsydsz
               kkk=kkk+1
               call pckindsk(kkk,3,syvec2,syvec1,nu2,nu2,a,b,c,d,nlist)
               t2(c,d)=t2s(is)
 8          continue
            call vpakv(nfl,nu,a,b,ifr,inr,t2)
 9       CONTINUE
 10   CONTINUE
      RETURN
      END
      SUBROUTINE PCKINDsk(IVALUE,ISPIN,SYVEC1,SYVEC2,NSMSZ1,NSMSZ2,
     &                  I,J,A,B,list)
C
C THIS ROUTINE RETURNS I,J,A AND B FOR IVALUE, WHICH IS THE ABSOLUTE
C  OFFSET IN THE SYMMETRY PACKED LIST.  T
C  IF ISPIN IS 1 OR 2, THE LIST STRUCTURE IS ASSUMED TO BE A<B;I<J AND
C  a<b;i,j RESPECTIVELY.  IF ISPIN IS 3, IT IS ASSUMED TO BE A,b;I,j.
C  NDISSZ IS THE SIZE OF THE DISTRIBUTIONS.
C
CEND
C
      IMPLICIT INTEGER (A-Z)
      DIMENSION SYVEC1(NSMSZ1),SYVEC2(NSMSZ2)
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),N(18)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
      ILRG(IX)=INT(0.5D0*(1.D0+DSQRT(8.d0*IX-7)))+1
      NNM1O2(IX)=(IX*(IX-1))/2
      IGETI(LSTELM,NUMA)=1+(LSTELM-1)/NUMA
      IGETA(LSTELM,NUMA)=LSTELM-(IGETI(LSTELM,NUMA)-1)*NUMA
      IF(IVALUE.Eq.0)THEN
       I=0
       J=0
       A=0
       B=0
       RETURN
      ENDIF
c      IF(ISPIN.EQ.3)LIST=46
c      IF(ISPIN.EQ.2)LIST=45
c      IF(ISPIN.EQ.1)LIST=44
C
C DETERMINE THE DPD IRREP FOR IVALUE.
C
      NSIZ=0
      IBOTL=0
      IBOTR=0
      DO 10 IRREP=1,NIRREP
       NSYDIS=IRPDPD(IRREP,ISYTYP(2,LIST))
       NSYDSZ=IRPDPD(IRREP,ISYTYP(1,LIST))
       NSIZ=NSIZ+NSYDIS*NSYDSZ
       IF(IVALUE.LE.NSIZ)GOTO 20
       IBOTR=IBOTR+NSYDIS
       IBOTL=IBOTL+NSYDSZ
c       write(6,*)'ibotr,ibotl:',ibotr,ibotl
10    CONTINUE
20    NOFF=NSIZ-NSYDIS*NSYDSZ
C
C THE DPD IRREP IS KNOWN, SO DETERMINE THE OFFSET WITHIN THIS IRREP.
C
      IOFF=IVALUE-NOFF
      NSYDIS=IRPDPD(IRREP,ISYTYP(2,LIST))
      NSYDSZ=IRPDPD(IRREP,ISYTYP(1,LIST))
c      write(6,*)'irrep,ioff,nsydis,nsydsz:',irrep,ioff,nsydis,nsydsz
C
C DETERMINE OFFSETS IN BOTH ELEMENTS OF THE DPD LIST FOR THIS IRREP.
C
      if(ioff.eq.0)then
        idisp=0
        ielemp=0
       else
      IDISP=IGETI(IOFF,NSYDSZ)
      IELEMP=IGETA(IOFF,NSYDSZ)
      endif
c      write(6,*)'idisp,ielemp:',idisp,ielemp
C
C DETERMINE OFFSETS IN UNPACKED LIST.
C
c      write(6,*)'ibotl+ielemp,ibotr+idisp',(ibotl+ielemp),(ibotr+idisp)
      IELEMU=SYVEC1(IBOTL+IELEMP)
      IDISU=SYVEC2(IBOTR+IDISP)
      if(list.eq.13)then
      J=1+(IDISU-1)/NOCCO(1)    !ijkL
      I=IDISU-(J-1)*NOCCO(1)    !ijKl
      B=1+(IELEMU-1)/NOCCO(1)   !iJkl
      A=IELEMU-(B-1)*NOCCO(1)   !Ijkl
      return
      endif
      if(list.eq.10)then
      J=1+(IDISU-1)/NOCCO(1)    !ijkA
      I=IDISU-(J-1)*NOCCO(1)    !ijKa
      B=1+(IELEMU-1)/NOCCO(1)   !iJka
      A=IELEMU-(B-1)*NOCCO(1)   !Ijka
      return
      endif
      if(list.eq.16)then
      J=1+(IDISU-1)/NOCCO(1)    !abiJ
      I=IDISU-(J-1)*NOCCO(1)    !abIj
      B=1+(IELEMU-1)/NVRTO(1)   !aBij
      A=IELEMU-(B-1)*NVRTO(1)   !Abij
      return
      endif
      if(list.eq.25)then
      J=1+(IDISU-1)/NVRTO(1)    !biaJ
      I=IDISU-(J-1)*NVRTO(1)    !biAj
      B=1+(IELEMU-1)/NVRTO(1)   !bIaj
      A=IELEMU-(B-1)*NVRTO(1)   !Biaj
      return
      endif
      if(list.eq.30)then
      J=1+(IDISU-1)/NVRTO(1)    !abcI
      I=IDISU-(J-1)*NVRTO(1)    !abCi
      B=1+(IELEMU-1)/NVRTO(1)   !aBci
      A=IELEMU-(B-1)*NVRTO(1)   !Abci
      endif
      if(list.eq.46)then        
      J=1+(IDISU-1)/NOCCO(1)    !abiJ
      I=IDISU-(J-1)*NOCCO(1)    !abIj
      B=1+(IELEMU-1)/NVRTO(1)   !aBij
      A=IELEMU-(B-1)*NVRTO(1)   !Abij
      return
      endif
      if(list.eq.233)then
      J=1+(IDISU-1)/NVRTO(1)    !abcD
      I=IDISU-(J-1)*NVRTO(1)    !abCd
      B=1+(IELEMU-1)/NVRTO(1)   !aBcd
      A=IELEMU-(B-1)*NVRTO(1)   !Abcd
      endif
      RETURN
      END
      SUBROUTINE PCKINDsk22(IVALUE,ISPIN,SYVEC1,SYVEC2,NSMSZ1,NSMSZ2,
     &                  I,J,A,B,list)
C
C THIS ROUTINE RETURNS I,J,A AND B FOR IVALUE, WHICH IS THE ABSOLUTE
C  OFFSET IN THE SYMMETRY PACKED LIST.  T
C  IF ISPIN IS 1 OR 2, THE LIST STRUCTURE IS ASSUMED TO BE A<B;I<J AND
C  a<b;i,j RESPECTIVELY.  IF ISPIN IS 3, IT IS ASSUMED TO BE A,b;I,j.
C  NDISSZ IS THE SIZE OF THE DISTRIBUTIONS.
C
CEND
C
      IMPLICIT INTEGER (A-Z)
      logical smlv
      DIMENSION SYVEC1(NSMSZ1),SYVEC2(NSMSZ2)
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),N(18)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
      common/flags/iflags(100)
      equivalence(icllvl,iflags(2))
      ILRG(IX)=INT(0.5D0*(1.D0+DSQRT(8.d0*IX-7)))+1
      NNM1O2(IX)=(IX*(IX-1))/2
      IGETI(LSTELM,NUMA)=1+(LSTELM-1)/NUMA
      IGETA(LSTELM,NUMA)=LSTELM-(IGETI(LSTELM,NUMA)-1)*NUMA
      smlv=icllvl.eq.11.or.icllvl.eq.13.or.icllvl.eq.14.or.
c     *icllvl.eq.15.or.icllvl.eq.16.or.icllvl.eq.17.or.icllvl.eq.19.or.
     *icllvl.eq.15.or.icllvl.eq.19.or.
     *icllvl.eq.21.or.icllvl.eq.22
      IF(IVALUE.Eq.0)THEN
       I=0
       J=0
       A=0
       B=0
       RETURN
      ENDIF
c      IF(ISPIN.EQ.3)LIST=46
c      IF(ISPIN.EQ.2)LIST=45
c      IF(ISPIN.EQ.1)LIST=44
C
C DETERMINE THE DPD IRREP FOR IVALUE.
C
      NSIZ=0
      IBOTL=0
      IBOTR=0
      DO 10 IRREP=1,NIRREP
       NSYDIS=IRPDPD(IRREP,ISYTYP(2,LIST))
       NSYDSZ=IRPDPD(IRREP,ISYTYP(1,LIST))
       if(list.eq.233.and.(smlv))nsydsz=nsydis
       NSIZ=NSIZ+NSYDIS*NSYDSZ
       IF(IVALUE.LE.NSIZ)GOTO 20
       IBOTR=IBOTR+NSYDIS
       IBOTL=IBOTL+NSYDSZ
c       write(6,*)'ibotr,ibotl:',ibotr,ibotl
10    CONTINUE
20    NOFF=NSIZ-NSYDIS*NSYDSZ
C
C THE DPD IRREP IS KNOWN, SO DETERMINE THE OFFSET WITHIN THIS IRREP.
C
      IOFF=IVALUE-NOFF
      NSYDIS=IRPDPD(IRREP,ISYTYP(2,LIST))
      NSYDSZ=IRPDPD(IRREP,ISYTYP(1,LIST))
       if(list.eq.233.and.(smlv))nsydsz=nsydis
c      write(6,*)'irrep,ioff,nsydis,nsydsz:',irrep,ioff,nsydis,nsydsz
C
C DETERMINE OFFSETS IN BOTH ELEMENTS OF THE DPD LIST FOR THIS IRREP.
C
      if(ioff.eq.0)then
        idisp=0
        ielemp=0
       else
      IDISP=IGETI(IOFF,NSYDSZ)
      IELEMP=IGETA(IOFF,NSYDSZ)
      endif
c      write(6,*)'idisp,ielemp:',idisp,ielemp
C
C DETERMINE OFFSETS IN UNPACKED LIST.
C
c      write(6,*)'ibotl+ielemp,ibotr+idisp',(ibotl+ielemp),(ibotr+idisp)
      IELEMU=SYVEC1(IBOTL+IELEMP)
      IDISU=SYVEC2(IBOTR+IDISP)
      if(list.eq.13)then
      J=1+(IDISU-1)/NOCCO(1)    !ijkL
      I=IDISU-(J-1)*NOCCO(1)    !ijKl
      B=1+(IELEMU-1)/NOCCO(1)   !iJkl
      A=IELEMU-(B-1)*NOCCO(1)   !Ijkl
      return
      endif
      if(list.eq.10)then
      J=1+(IDISU-1)/NOCCO(1)    !ijkA
      I=IDISU-(J-1)*NOCCO(1)    !ijKa
      B=1+(IELEMU-1)/NOCCO(1)   !iJka
      A=IELEMU-(B-1)*NOCCO(1)   !Ijka
      return
      endif
      if(list.eq.16)then
      J=1+(IDISU-1)/NOCCO(1)    !abiJ
      I=IDISU-(J-1)*NOCCO(1)    !abIj
      B=1+(IELEMU-1)/NVRTO(1)   !aBij
      A=IELEMU-(B-1)*NVRTO(1)   !Abij
      return
      endif
      if(list.eq.25)then
      J=1+(IDISU-1)/NVRTO(1)    !biaJ
      I=IDISU-(J-1)*NVRTO(1)    !biAj
      B=1+(IELEMU-1)/NVRTO(1)   !bIaj
      A=IELEMU-(B-1)*NVRTO(1)   !Biaj
      return
      endif
      if(list.eq.30)then
      J=1+(IDISU-1)/NVRTO(1)    !abcI
      I=IDISU-(J-1)*NVRTO(1)    !abCi
      B=1+(IELEMU-1)/NVRTO(1)   !aBci
      A=IELEMU-(B-1)*NVRTO(1)   !Abci
      endif
      if(list.eq.46)then        
      J=1+(IDISU-1)/NOCCO(1)    !abiJ
      I=IDISU-(J-1)*NOCCO(1)    !abIj
      B=1+(IELEMU-1)/NVRTO(1)   !aBij
      A=IELEMU-(B-1)*NVRTO(1)   !Abij
      return
      endif
      if(list.eq.233)then
      J=1+(IDISU-1)/NVRTO(1)    !abcD
      I=IDISU-(J-1)*NVRTO(1)    !abCd
      B=1+(IELEMU-1)/NVRTO(1)   !aBcd
      A=IELEMU-(B-1)*NVRTO(1)   !Abcd
      endif
      RETURN
      END
