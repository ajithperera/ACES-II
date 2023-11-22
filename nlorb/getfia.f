      subroutine getfia(nbastot,nocc,nvir,fock,
     & sumoff)
      implicit none

      integer iii, aaa

      integer nocc, nvir, nbastot

      double precision sumoff, cutoff

      double precision fock(nbastot,nbastot)

      parameter (
     & cutoff = 1.0D-6
     & )

      sumoff = 0.0D0
      do iii = 1, nocc
         do aaa = 1, nvir
            sumoff = sumoff + abs(fock(iii,aaa+nocc))
         end do
      end do

      sumoff = sumoff/(dble(nocc*nvir))
      write(*,*)
      write(*,10) sumoff
      write(*,*)

      if (sumoff.ge.cutoff) then
         write(*,11)
C         stop
      end if 

 10   format('absolute magnitude of occ-vir block of fock',F20.15) 
 11   format('final set of orbitals are far from HF orbitals',/,
     & 'meaning that the NLMOs are far from NBOs - restart',/,
     & 'with larger bondsize - program will exit')

      return
      end 

