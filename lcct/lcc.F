      Subroutine lcc(work, maxcor, uhf )
      implicit none
      integer maxcor,work(maxcor),uhf,listl1,listl1off,listl2,listl2rs,
     &        listr2,listr2rs,listgrl,listgtl,listgrlof,listgtlof,
     &        listtmp,listtmpof,listt1,listt1off,listr1,listr1off,
     &        listz1,listz1off,listz2
      logical nonhf,anderiv
      double precision edummy,one
      data one /1.0d0/
C The LCC is developed to compute Lambda based CCSD(T) energies
C and gradients. The L-CCSD(T) has a simple structure. Instead
C of the convered CCSD T2 amplitudes, the converged CCSD L2
C amplitudes are used to compute the two triple correction diagrams.
C
#include "flags.h"
#include "flags2.com"
      anderiv = iflags2(h_IFLAGS2_grad_calc) .eq. 1
      Call LTrps(work,maxcor,uhf,edummy,"CCSD(T)")

      if (anderiv) then
C     The sigma equations are of the form H{bar}Sigma = <S,D|WT3|0>=<S,D|d3WT2|0>
C     This is formed and stored within Ltrps
C     Setup some common block variables
        call setderiv(listl1,listl1off,listl2,listl2rs,listr1,listr1off,
     &                listr2,listr2rs,listgrl,listgtl,listgrlof,
     &                listgtlof,listtmp,listtmpof,listt1,listt1off,
     &                listz1,listz1off,listz2,uhf,nonhf)
C Form Hbar
        call w5t1abcd(work,maxcor,uhf)
        call formw1(work,maxcor,uhf,.true.)
        call hbarabcd(work,maxcor,uhf,.true.)
        call hbrijka0(work,maxcor,uhf)
        call fixfbar(work,maxcor,uhf)
        call modaibc(work,maxcor,uhf,one)

        call modf(work,maxcor,uhf,1)
        if (uhf .eq. 0) call makess(work,maxcor,uhf)

        write(6,*) 'Forming Sigma Matrix'
C Call the solver for Hbar from vcceh to form sigma
        call formsig(work,maxcor,uhf,listr1,listr1off,listr2,listr2rs)

C     At this point the sigma term for L-CCSD(T) has been constructed
C     Construct RHS and solve for the second intermediate vector.
        write(6,*) 'Forming Pi Matrix'
        call formpi(work,maxcor,uhf,listl1,listl1off,listl2,listl2rs,
     &              listr1,listr1off,listr2,listr2rs,listz1,listz1off,
     &              listz2)

        if (uhf .eq. 0) call resetss(work,maxcor,uhf)
        call modf(work,maxcor,uhf,-1)
C Unform Hbar
        call modaibc(work,maxcor,uhf,-one)
        call modiajk(work,maxcor,uhf,-one)
        call hbarabcd(work,maxcor,uhf,.false.)
        call formw1(work,maxcor,uhf,.false.)
        write(6,*) 'Forming Lambda-Sigma Density Elements'
        call form_dens_lamsig(work,maxcor,uhf,nonhf,listl1,
     &                        listl1off,listl2,listl2rs,listr1,
     &                        listr1off,listr2,listr2rs,listgrl,listgtl,
     &                        listgrlof,listgtlof,listtmp,listtmpof,
     &                        listt1,listt1off)

        call newamp_lcc(work,maxcor,uhf)
        if (uhf .eq. 0) call abtoaa(work,maxcor,uhf,listl2)
        call rnabij(work,maxcor,uhf,'L')
        call gformg(1,1,144,44,100,work,maxcor,0,one,uhf)
        if (uhf .eq. 1) then
          call newtyp2(1,56,9,10,.TRUE.)
          call newtyp2(1,57,10,9,.TRUE.)
        endif
      endif

      return
      End


