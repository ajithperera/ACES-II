      subroutine doitran(row,trans,scr,centbf,mombf,centbfo,mombfo,
     & nbas,nbastot,angmax)
      implicit none

CCC hughes

C      integer iii

CCC hughes
 
      integer nbas, nbastot, angmax, row(nbas-nbastot),
     & trans(nbastot,nbas), scr(nbas,nbas), centbf(nbas),
     & mombf(nbas), centbfo(nbastot), mombfo(nbastot)

CCC hughes

C      do iii = 1, nbas
C         write(*,*) iii, centbf(iii)
C      end do

CCC hughes

      call mkrow(row,nbas-nbastot,mombf,angmax,nbas,nbastot)
      call mkitrans(trans,nbas,nbastot,scr,row,nbas-nbastot)
      call imult(trans,nbastot,nbas,centbf,nbas,1,centbfo,nbastot,1)
      call imult(trans,nbastot,nbas,mombf,nbas,1,mombfo,nbastot,1)

      call izerosec(row,nbas-nbastot)
      call izerosec(trans,nbastot*nbas)
      call izerosec(scr,nbas*nbas)

CCC hughes

C      do iii = 1, nbastot
C         write(*,*) iii, centbfo(iii)
C      end do

CCC hughes

      return
      end
