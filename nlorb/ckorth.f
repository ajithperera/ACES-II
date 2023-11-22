      subroutine ckorth(ovrlpao,ovrlpmo,coef,nbas)
      implicit none

      integer iii, jjj

      integer nbas

      double precision ovrlpao(nbas,nbas), ovrlpmo(nbas,nbas),
     & coef(nbas,nbas)

      call orthock(ovrlpmo,ovrlpao,coef,nbas)

      return
      end
