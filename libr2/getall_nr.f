        subroutine getall_nr(target,scr,ScrIDim,list,irrep)
c
c wrapper around getall to get the redundant W/T/L(ij,ab)
c PR
        implicit integer (a-z)
        integer list, reflist,maxcor,iscrsizi,irrepsiz
        integer irrep,tarsiz,tardis,refsiz,refdis
        integer ibgn,i001,i002,i003
        logical bSST002, bSST003, bSSTRNG
        double precision target(*), scr(*), fact
        character *4 szCase, szS003
        COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
        COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
c        COMMON/SYMPOP2/IRPDPD(8,22)
c        COMMON/SYMPOP/IRP_DM(8,22),ISYTYP(2,500),NTOT(18)
        COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),NJUNK(18)
        COMMON /SYM2/ POP(8,2),VRT(8,2),NT(2),NFEA(2),NFMI(2)
        COMMON /INFO/ NOCCO(2),NVRTO(2)

c   o the units of scratch are INTEGERS, but the array is cast in DOUBLES
      if (ScrIDim.lt.1) then
         print *, '@GETALL_NR: Assertion failed.'
         print *, '            scratch dimension is ',ScrIDim
         call aces_exit(1)
      end if
      maxcor = ScrIDim/iintfp

        bSST002 = .false.
        bSST003 = .false.
        bSSTRNG = .false.
        iofflist = 0
        fact = 1.0d0
        szS003 = "AIBJ"
        mylist = list

c        print *,"GETALL_NR: list ",list

c for L(ijab), avoid repeated conditional and cast 'list' into T space
        if(mylist.gt.133 .and. mylist.lt.140) then
          iofflist = 100
          mylist=mylist-100
          szS003 = "AJBI"
        endif
c integrals
        if(mylist.eq.17) then
           reflist = 16
           szCase="BBAA"
           bSST002=.true.
        else if(mylist.eq.18) then
           reflist = 16
           szCase="AABB"
           bSST002=.true.
        else if(mylist.eq.19) then
           reflist = 14
           szCase="AAAA"
           bSST003=.true.
        else if(mylist.eq.20) then
           reflist = 15
           szCase="BBBB"
           bSST003=.true.
        else if(mylist.eq.21) then
           reflist = 16
           szCase="BBAA"
           bSST002=.true.
           bSSTRNG=.true.
        else if(mylist.eq.22) then
           reflist = 16
           szCase="AABB"
           bSST002=.true.
           bSSTRNG=.true.
c amplitudes
        else if(mylist.eq.34) then
           reflist = 44+iofflist
           szCase="AAAA"
           bSST003=.true.
           szS003 = "AJBI"
        else if(mylist.eq.35) then
           reflist = 45+iofflist
           szCase="BBBB"
           bSST003=.true.
           szS003 = "AJBI"
        else if(mylist.eq.36) then
           reflist = 46+iofflist
           szCase="BBAA"
           bSST002=.true.
        else if(mylist.eq.37) then
           reflist = 46+iofflist
           szCase="AABB"
           bSST002=.true.
        else if(mylist.eq.38) then
           reflist = 46+iofflist
           szCase="BBAA"
           bSST002=.true.
           bSSTRNG=.true.
        else if(mylist.eq.39) then
           reflist = 46+iofflist
           szCase="AABB"
           bSST002=.true.
           bSSTRNG=.true.
        else
           print *,"@GETALL_NR: WRONG LIST ARGUMENT"
c           print *,"IRREP, LIST ",irrep, list
           call aces_exit(-1)
        endif

        tarsiz = isymsz(isytyp(1,mylist+iofflist),
     &                 isytyp(2,mylist+iofflist))
        refsiz = isymsz(isytyp(1,reflist),isytyp(2,reflist))

c        print *,"tarsiz, refsiz: ", tarsiz, refsiz
c        print *,"list, reflist: ",list, reflist

        nao = max(nocco(1),nocco(2)) + max(nvrto(1),nvrto(2))
        iscrsiz = nao*nao
C
        IBGN = 1
        I001 = IBGN + MAX(tarsiz,refsiz)
        I002 = I001 + MAX(tarsiz,refsiz)
        I003 = I002 + iscrsiz
        IF (I003.GT.MAXCOR) THEN
           call insmem('GETALL_NR',(i003*iintfp)-1,ScrIDim)
        END IF

        call getall(scr(I001), MAX(tarsiz,refsiz), 1, reflist)
c
c choose the correct transformation
c
        if(bSST002) then
           call sst002(scr(I001), scr(IBGN), refsiz, tarsiz,
     &                 scr(I002), szCase)
        endif
        if(bSST003) then
           call sst003(scr(I001), scr(IBGN), refsiz, tarsiz,
     &                 scr(I002), szCase, szS003)
        endif
        if(bSSTRNG) then
           call sstrng(scr(IBGN), scr(I001), tarsiz, tarsiz,
     &                 scr(I002), szCase)
           call dcopy(tarsiz, scr(I001), 1, scr(IBGN),1)
        endif

        call dcopy(tarsiz,scr(IBGN),1,target,1)

        return
        end

