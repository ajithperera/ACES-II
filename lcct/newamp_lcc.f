

      subroutine newamp_lcc(icore,maxcor,iuhf)
c
c overwrites lambda lists with lambda+zeta (pi)
c
      implicit integer (a-z)
      dimension icore(maxcor)
      common/machsp/iintln,ifltln,iintfp,ialone,ibitwd
      common/sympop/irpdpd(8,22),isytyp(2,500),id(18)
      common/sym/pop(8,2),vrt(8,2),nt(2),nfmi(2),nfea(2)
      common/syminf/nstart,nirrep,irrepy(255,2),dirprd(8,8)
c
      do 10 ispin=3,3-2*iuhf,-1
       listl=143+ispin
       listz=343+ispin
       do 20 irrep=1,nirrep
        numdis=irpdpd(irrep,isytyp(2,listl))
        dissiz=irpdpd(irrep,isytyp(1,listl))
        i000=1
        i010=i000+iintfp*numdis*dissiz
        call getlst(icore(i000),1,numdis,1,irrep,listl)
        call getlst(icore(i010),1,numdis,1,irrep,listz)
        call saxpy (numdis*dissiz,1.0d0,icore(i010),1,icore(i000),1)
        call putlst(icore(i000),1,numdis,1,irrep,listl)
20     continue
10    continue
c
      do 30 ispin=1,1+iuhf
       length=nt(ispin)
       i000=1
       i010=i000+iintfp*length
       call getlst(icore(i000),1,1,1,ispin,190)
       call getlst(icore(i010),1,1,1,ispin,390)
       call saxpy (length,1.0d0,icore(i010),1,icore(i000),1)
       call putlst(icore(i000),1,1,1,ispin,190)
30    continue
      return
      end
