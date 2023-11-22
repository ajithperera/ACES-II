      subroutine check(zabab,zaaaa,zbbbb,scr,irrepx)
      implicit double precision (a-h,o-z)
      integer dirprd,pop,vrt,disabab,disaaaa
      dimension zabab(*),zaaaa(*),zbbbb(*),scr(*)
      common/machsp/iintln,ifltln,iintfp,ialone,ibitwd
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON/SYMINF/NSTART,NIRREP,IRREPY(255,2),DIRPRD(8,8)
c
c antisymmetrize abab list explicitly
c
      lenabab=idsymsz(irrepx,isytyp(1,46),isytyp(2,46))
      ioffabab=1
      ioffaaaa=1
      do 10 irrepr=1,nirrep
       irrepl=dirprd(irrepr,irrepx)
       numabab=irpdpd(irrepr,isytyp(2,46))
       disabab=irpdpd(irrepl,isytyp(1,46))
       numaaaa=irpdpd(irrepr,isytyp(2,44))
       disaaaa=irpdpd(irrepl,isytyp(1,44))
       call scopy(disabab*numabab,zabab(ioffabab),
     &    1,scr,1)
       call assym2(irrepr,pop(1,1),disabab,scr)
       call sqsym (irrepl,vrt(1,1),disaaaa,disabab,numaaaa,
     &             scr(numabab*disabab+1),scr)
       write(6,*)' antisymmetrized abab list '
       call prvecr(scr(numabab*disabab+1),numaaaa*disaaaa)
       write(6,*)' aaaa list passed in '
       call prvecr(zaaaa(ioffaaaa),numaaaa*disaaaa) 
       write(6,*)' bbbb list passed in '
       call prvecr(zbbbb(ioffaaaa),numaaaa*disaaaa) 
       ioffabab=ioffabab+numabab*disabab
       ioffaaaa=ioffaaaa+numaaaa*disaaaa
10    continue
c
      return
      end
