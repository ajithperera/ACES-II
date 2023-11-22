
      subroutine formpot(icore,maxcor,uhf,listv,singles,include_fai,
     &                   cycle,maxiter)
      implicit none
C Common blocks
      integer pop(8,2),vrt(8,2),nt(2),nd1(2),nd2(2)
      common/sym/pop,vrt,nt,nd1,nd2
      integer iintln,ifltln,iintfp
      common/machsp/iintln,ifltln,iintfp
C Input variables
      integer maxcor,uhf,listv,cycle,maxiter
      logical singles,include_fai
      double precision iterations(maxiter+1),errmax
C Pre-allocated Local variables
      integer icore(maxcor)
C Local variables
      integer spin,ivai,ifai,idoo,idvv,iend,itmp

csss      Write(6,*) "varibales singles and fai", singles, include_fai
      call getrec(-1,'JOBARC','OOMBITER',(maxiter+1)*iintfp,iterations)
      do spin = 1,uhf+1
        itmp = 0
        call aces_list_touch(1,nt(spin),spin,listv,itmp)
        call aces_list_touch(1,nt(spin),2+spin,93,itmp)
        if (.not. include_fai) then
          call zero(icore,nt(spin))
        else
          call getlst(icore,1,1,1,2+spin,93)
        endif
        call putlst(icore,1,1,1,spin,listv)
C Vai  = SUM(bj) Tbj <ab||ij>
        if (singles) call makvme(icore,maxcor,uhf,spin,listv)
#ifdef _DEBUG_LVLM1
        Write(6,*) "Tbj<ab||ij>"
        call print_vai(icore, nt(1), 1, listv)
#endif
C Vai += Fai(Dii-Daa)
        ivai = 1
        ifai = ivai+nt(spin)*iintfp
        idoo = ifai+nt(spin)*iintfp
        idvv = idoo+nd1(spin)*iintfp
        iend = idvv+nd2(spin)*iintfp
        if (iend .gt. maxcor) call insmem('FORMPOT',iend,maxcor)
        
        errmax = iterations(cycle)
        Write(6,"(a,F10.5)") "@-Formpot errmax:", errmax
       
        if (errmax .lt. 1.0D4) then
         call addfai(icore(ivai),icore(ifai),icore(idoo),icore(idvv),
     &               nt(spin),nd1(spin),nd2(spin),spin,listv)
        endif 
#ifdef _DEBUG_LVLM1
         Write(6,*) "After addfai"
         call print_vai(icore, nt(1), 1, listv)
#endif
C Vai += SUM(bj) Tbj <aj||ib>
        if (singles) call t1inv1(icore,maxcor,uhf,spin,listv)
#ifdef _DEBUG_LVL0
         Write(6,*) "After Tbj<aj|ib>"
         call print_vai(icore, nt(1), 1, listv)
#endif

      end do

C Vai += SUM(j) Djj <aj||ij> + SUM(b) Dbb <ab||ib>
      if (errmax .lt. 1.0D4) call dwinv1(icore,maxcor,uhf,listv)
#ifdef _DEBUG_LVLM1
      Write(6,*) " Djj + Dii"
      call print_vai(icore, nt(1), 1, listv)
#endif
C Vai += 1/2 SUM(jbc) Tbicj <aj||bc> - 1/2 SUM(jkb) Tajbk <jk||ib>
      call t2inv1(icore,maxcor,uhf,listv)
#ifdef _DEBUG_LVLM1
      Write(6,*) "Tbicj<aj||bc>"
      call print_vai(icore, nt(1), 1, listv)
#endif
      return
      end

      subroutine print_vai(vai,len,spin,listv)
      implicit none
      integer pop(8,2),vrt(8,2),nt(2),nd1(2),nd2(2)
      common/sym/pop,vrt,nt,nd1,nd2
      integer nstart,nirrep
      common/syminf/nstart,nirrep
      integer len,spin,listv
      double precision vai(len)
      integer irrep,occ,virt,ivai

      call getlst(vai,1,1,1,spin,listv)
      ivai = 1
      Write(6,*) "Vai in formpot"
      do irrep = 1,nirrep
        write(6,*) 'irrep=',irrep
        occ = pop(irrep,spin)
        virt = vrt(irrep,spin)
        call output(vai(ivai), 1, occ, 1, virt, occ, virt, 1)
        ivai = ivai + occ*virt
      end do

      return
      end

