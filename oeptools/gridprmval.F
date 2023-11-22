      subroutine gridprmval (
     &     natm, nprm, maxshl, maxang, 
     &     nshlatm, nprmshl, angfct, alpha, 
     &     dist, prmgrd)
c
c The subroutine evaluates the values of primitives at
c a given grid point.
c
c Igor Schweigert, Oct 2002
c $Id: gridprmval.FPP,v 1.1.1.1 2008/05/29 19:35:40 taube Exp $
c 
      implicit none
c
c     .. input parameters ..
c
      integer
     &     natm, nprm, maxshl, maxang,
     &     nshlatm (natm), nprmshl (maxshl, natm), angfct (maxang,3)
c
      double precision
     &     dist (natm, 3), alpha (nprm), prmgrd (nprm)
c
c     .. local variables ..
c
      integer
     &     iprm, iatm, ishl, i
c
      double precision
     &     dst, pre
c
      iprm = 1
      do iatm = 1, natm
c
c        .. distance from the primitive center to the grid point ..
c
         dst = 
     &        dist (iatm,1)**2 + dist (iatm,2)**2 + dist (iatm,3)**2
c
         do ishl = 1, nshlatm (iatm)
c
c        .. the angular moment dependent prefactor ..
c
            pre =
     &           dist (iatm, 1)**angfct (ishl, 1)*
     &           dist (iatm, 2)**angfct (ishl, 2)*
     &           dist (iatm, 3)**angfct (ishl, 3)
c
            do i = 1, nprmshl (ishl, iatm)
c
c              .. the value of the primitive at the grid  point ..
c
               prmgrd (iprm) = pre * dexp (-1.d0 * alpha (iprm) * dst)
               iprm = iprm + 1
c
            enddo
         enddo
      enddo
c     
      return
      end
