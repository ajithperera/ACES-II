      subroutine hfperbond(core,fock,nbas,nocc,energy,
     & energyn,enuc)
      implicit none

      integer iii
      
      integer nbas, nocc

      double precision core(nbas,nbas), fock(nbas,nbas),
     & energy(nocc), energyn(nocc), enuc, total, totaln

      total = 0.0D0
      do iii = 1, nocc
         energy(iii) = core(iii,iii) + fock(iii,iii)
         write(*,10) iii, energy(iii)
         total = total + energy(iii)
      end do
      write(*,20) 
      write(*,30) total
      write(*,*)
      
      totaln = 0.0D0
      do iii = 1, nocc
         energyn(iii) = (enuc*energy(iii))/total
         write(*,40) iii, energyn(iii)
         totaln = totaln + energyn(iii)
      end do
      write(*,20) 
      write(*,50) totaln
      write(*,*)
     
      total = 0.0D0 
      do iii = 1, nocc
         energy(iii) = energy(iii) + energyn(iii)
         write(*,60) iii, energy(iii)
         total = total + energy(iii)
      end do
      write(*,20) 
      write(*,70) total
      write(*,*)

 10   format('Orbital ',I6,' has electronic energy ',F20.15)
 40   format('Orbital ',I6,' has nuclear energy ',F20.15)
 60   format('Orbital ',I6,' has total energy ',F20.15)
 20   format('----------------------------------------')
 30   format('Electronic energy is ',F25.15)
 50   format('Nuclear energy is ',F25.15)
 70   format('Total energy is ',F25.15)

      return
      end

