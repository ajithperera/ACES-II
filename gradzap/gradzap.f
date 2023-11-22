      program gradzap
c---------------------------------------------------------------------------
c   Used to restore GRADIENT entry in JOBARC at the beginning of each
c   iteration of the gradient calculations.
c--------------------------------------------------------------------------

      implicit none
      double precision dGradient(3*2048)
      integer nInts
      do nInts = 1, 3*2048
         dGradient(nInts) = 0.d0
      end do
      call aces_init_rte
      call aces_ja_init
      call getrec(0,'JOBARC','GRADIENT',nInts,dGradient)
      if (nints .gt. 0) 
     *   call putrec(1,'JOBARC','GRADIENT',nInts,dGradient)
      call aces_ja_fin
      end
