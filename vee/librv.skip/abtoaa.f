
      subroutine abtoaa(icore,maxcor,iuhf)
c
c generate aaaa spin case list for rhf
c
      implicit integer (a-z)
      dimension icore(maxcor)
      common/machsp/iintln,ifltln,iintfp,ialone,ibitwd
      common/sympop/irpdpd(8,22),isytyp(2,500),id(18)
      common/sym/pop(8,2),vrt(8,2),nt(2),nfmi(2),nfea(2)
      common/syminf/nstart,nirrep,irrepy(255,2),dirprd(8,8)
      if(iuhf.eq.0)then
       DO 200 IRREP=1,NIRREP
        listab=146
        listaa=144
        NUMAB=IRPDPD(IRREP,ISYTYP(2,listab))
        DSZAB=IRPDPD(IRREP,ISYTYP(1,listab))
        NUMAA=IRPDPD(IRREP,ISYTYP(2,listaa))
        DSZAA=IRPDPD(IRREP,ISYTYP(1,listaa))
        NSIZ1=DSZAB*NUMAA
        NSIZ2=DSZAB*NUMAB
        I000=1
        I010=I000+IINTFP*NSIZ1
        I020=I010+IINTFP*NSIZ2
        CALL GETLST(ICORE(I010),1,NUMAB,1,IRREP,LISTAB)
        CALL ASSYM(IRREP,POP,DSZAB,DSZAB,ICORE(I000),ICORE(I010))
        CALL SQSYM(IRREP,VRT,DSZAA,DSZAB,NUMAA,ICORE(I010),ICORE(I000))
        CALL PUTLST(ICORE(I010),1,NUMAA,1,IRREP,LISTAA)
200    CONTINUE 
      endif
      return
      end
