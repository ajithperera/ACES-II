
      PROGRAM OOMBPT2
      implicit none
C     Common blocks
      integer icore
      common//icore(1)
C Local variables
      integer i0,maxcor,uhf,listv,cycle,maxiter
      logical conv, singles, include_fai

C Initialize everything
      call aces_init(icore,i0,maxcor,uhf,.true.)
C Setup calculation
      call setmet(singles,include_fai,listv,cycle)
C Step 3: Form density matrices
      call onerdm(icore(i0),maxcor,uhf,singles)
C Generalization: Need the density matrices from an arbitrary method
C in particular, want the response density matrices. Need to check
C Step 4: Form potential
      call getrec(20,'JOBARC','OOMBPTCY',1,maxiter)
      call formpot(icore(i0),maxcor,uhf,listv,singles,include_fai,
     &             cycle,maxiter)
C Step 5: Check for convergence
	  call conv_check_wrap(icore(i0),maxcor,uhf,conv,listv,
     &                         cycle)
      if (conv) then
        call conv_set(icore(i0),maxcor,uhf,cycle)
      else
C Step 5: Form new fock operator and diagonalize for new orbs
C F' = F + <i|v_s|a> + <a|v_s|i>
c        call formfock_wrap(icore(i0),maxcor,uhf,listv)
        call oombpt_step_wrap(icore(i0),maxcor,uhf,listv,cycle)
      endif
      call aces_fin
      stop
      end

