

      subroutine newamp(icore,maxcor,iuhf,r0)
c
c overwrites zeta lists with zeta + r0*L
c
      implicit integer (a-z)
      double precision r0
      dimension icore(maxcor)
      common/machsp/iintln,ifltln,iintfp,ialone,ibitwd
      common/sympop/irpdpd(8,22),isytyp(2,500),id(18)
      common/sym/pop(8,2),vrt(8,2),nt(2),nfmi(2),nfea(2)
      common/syminf/nstart,nirrep,irrepy(255,2),dirprd(8,8)
c
      do 10 ispin=3,3-2*iuhf,-1
       listl=443+ispin
       listz=143+ispin
       do 20 irrep=1,nirrep
        numdis=irpdpd(irrep,isytyp(2,listl))
        dissiz=irpdpd(irrep,isytyp(1,listl))
        i000=1
        i010=i000+iintfp*numdis*dissiz
        call getlst(icore(i000),1,numdis,1,irrep,listz)
        call getlst(icore(i010),1,numdis,1,irrep,listl)
        call saxpy (numdis*dissiz,r0,icore(i010),1,icore(i000),1)
        call putlst(icore(i000),1,numdis,1,irrep,listz)
20     continue
10    continue
c
      do 30 ispin=1,1+iuhf
       length=nt(ispin)
       i000=1
       i010=i000+iintfp*length
       call getlst(icore(i000),1,1,1,ispin,190)
       call getlst(icore(i010),1,1,1,ispin,490)
       call saxpy (length,r0,icore(i010),1,icore(i000),1)
       call putlst(icore(i000),1,1,1,ispin,190)
30    continue
      return
      end
