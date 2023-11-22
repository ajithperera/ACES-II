      subroutine zer(scr,max)
      implicit double precision(a-h,o-z)
      dimension scr(max/2)
      integer dirprd
      common/syminf/nstart,nirrep,irrepa(255,2),dirprd(8,8)
      common/sympop/irpdpd(8,22),ISYTYP(2,500),id(18)
      common/dsym/irrepx
C
      DO 1000 irrepr=1,nirrep
       irrepl=dirprd(irrepx,irrepr)
       numdis=irpdpd(irrepr,isytyp(2,10))
       numdsz=irpdpd(irrepl,isytyp(1,10))
       call getlst(scr,1,numdis,1,irrepr,310)
       ind=0
       do 1001 I=1,NUMDIS
       call checksum(' 310   ',scr(ind+1),numdsz)
        do 1001 j=1,NUMDsz
        ind=ind+1
        write(*,*)'310', i,j,scr(ind)
1001   continue
       go to 1000
       numdis=irpdpd(irrepr,isytyp(2,10))
       numdsz=irpdpd(irrepl,isytyp(1,10))
       call getlst(scr,1,numdis,1,irrepr,310)
       ind=0
       do 1002 I=1,NUMDIS
        do 1002 j=1,NUMDsz
        ind=ind+1
        write(*,*) '310',i,j,scr(ind)
1002   continue
       numdis=irpdpd(irrepr,isytyp(2,16))
       numdsz=irpdpd(irrepl,isytyp(1,16))
       call getlst(scr,1,numdis,1,irrepr,316)
       ind=0
       do 1003 I=1,NUMDIS
        do 1003 j=1,NUMDsz
        ind=ind+1
        write(*,*) '316',i,j,scr(ind)
1003   continue
       numdis=irpdpd(irrepr,isytyp(2,30))
       numdsz=irpdpd(irrepl,isytyp(1,30))
       call getlst(scr,1,numdis,1,irrepr,330)
       ind=0
       do 1004 I=1,NUMDIS
        do 1004 j=1,NUMDsz
        ind=ind+1
        write(*,*) '330',i,j,scr(ind)
1004   continue
       numdis=irpdpd(irrepr,isytyp(2,21))
       numdsz=irpdpd(irrepl,isytyp(1,21))
       call getlst(scr,1,numdis,1,irrepr,321)
       ind=0
       do 1005 I=1,NUMDIS
        do 1005 j=1,NUMDsz
        ind=ind+1
        write(*,*) '321',i,j,scr(ind)
1005   continue
       numdis=irpdpd(irrepr,isytyp(2,25))
       numdsz=irpdpd(irrepl,isytyp(1,25))
       call getlst(scr,1,numdis,1,irrepr,325)
       ind=0
       do 1006 I=1,NUMDIS
        do 1006 j=1,NUMDsz
        ind=ind+1
        write(*,*) '325',i,j,scr(ind)
1006   continue
1000  continue
      return 
      end
