      subroutine gridprmhes (
     &     natm, nprm, maxshl, 
     &     nshlatm, nprmshl, angfct,
     &     atmcrd, alpha, pntcrd, prmhes)
c
c The subroutine evaluates the values of primitive hessians at a given
c grid point. PRMHES is sorted as XX, XY, YY, XZ, YZ, ZZ
c
c Igor Schweigert, Mar 2005
c $Id: gridprmhes.FPP,v 1.1.1.1 2008/05/29 19:35:40 taube Exp $
c 
      implicit none
c
c     Arguments
c
      integer
     &     natm, nprm, maxshl, 
     &     nshlatm (natm), nprmshl (maxshl, natm), angfct (3, maxshl)
c
      double precision
     &     atmcrd (3, natm), alpha (nprm), pntcrd (3), prmhes (nprm,6)
c
c     Local variables
c
      integer
     &     iprm, iatm, ishl, ixyz, i, j, fct,
     &     pow0 (3,6), pow1 (3,6), pow2 (3,6)
c
      double precision
     &     dist (3), dst, e, pre0 (6), pre1 (6), pre2 (6)
c
      data 
     &     pow0 /
     &     -2,  0,  0,
     &     -1, -1,  0,
     &      0, -2,  0,
     &     -1,  0, -1,
     &      0, -1, -1,
     &      0,  0, -2/,
     &     pow1 /
     &      0,  0,  0,
     &      1, -1,  0,
     &      0,  0,  0,
     &      1,  0, -1,
     &      0,  1, -1,
     &      0,  0,  0/,
     &     pow2 /
     &      2,  0,  0,
     &      1,  1,  0,
     &      0,  2,  0,
     &      1,  0,  1,
     &      0,  1,  1,
     &      0,  0,  2/
c
c     Calculate the hessians of the primitives: <> Loop over atoms. <>
c     Calculate the distance from the current atom to the grid point. <>
c     Loop over primitive shells on the atom. <> Calculate the
c     exponent-independent prefactors. <> Loop over primitives in the
c     current shell and calculate the hessians.
c
      iprm = 1
      do iatm = 1, natm
c
         dist (1) = pntcrd (1) - atmcrd (1, iatm)
         dist (2) = pntcrd (2) - atmcrd (2, iatm)
         dist (3) = pntcrd (3) - atmcrd (3, iatm)
         dst = dist (1)**2 + dist (2)**2 + dist (3)**2
c
         do ishl = 1, nshlatm (iatm)
c
            do j = 1, 3
               do i = 1, j
                  ixyz = j*(j-1)/2 + i
c
                  pre0 (ixyz) = 0d0
                  fct = angfct (i,ishl) * angfct (j,ishl)
                  if (i.eq.j) fct = fct - angfct (i,ishl)
                  if (fct.ne.0)
     &                 pre0 (ixyz) = pre0 (ixyz) + dble (fct) * 
     &                 dist (1) ** (angfct(1,ishl) + pow0(1,ixyz))*
     &                 dist (2) ** (angfct(2,ishl) + pow0(2,ixyz))*
     &                 dist (3) ** (angfct(3,ishl) + pow0(3,ixyz))
c     
                  pre1 (ixyz) = 0d0
                  fct = -2 * angfct(j,ishl)
                  if (i.eq.j) fct = fct - 1
                  if (fct.ne.0) 
     &                 pre1 (ixyz) = pre1 (ixyz) + dble (fct) * 
     &                 dist (1) ** (angfct(1,ishl) + pow1(1,ixyz))*
     &                 dist (2) ** (angfct(2,ishl) + pow1(2,ixyz))*
     &                 dist (3) ** (angfct(3,ishl) + pow1(3,ixyz))
                  fct = -2 * angfct (i,ishl)
                  if (i.eq.j) fct = fct - 1
                  if (fct.ne.0) 
     &                 pre1 (ixyz) = pre1 (ixyz) + dble (fct) * 
     &                 dist (1) ** (angfct(1,ishl) - pow1(1,ixyz))*
     &                 dist (2) ** (angfct(2,ishl) - pow1(2,ixyz))*
     &                 dist (3) ** (angfct(3,ishl) - pow1(3,ixyz))
c     
                  pre2 (ixyz) = 4d0 * 
     &                 dist (1) ** (angfct(1,ishl) + pow2(1,ixyz))*
     &                 dist (2) ** (angfct(2,ishl) + pow2(2,ixyz))*
     &                 dist (3) ** (angfct(3,ishl) + pow2(3,ixyz))
               enddo
            enddo
c     
            do i = 1, nprmshl (ishl, iatm)
               e = dexp ( - alpha (iprm) * dst)
               do ixyz = 1, 6
                  prmhes (iprm, ixyz) = e * (
     &                 pre0 (ixyz) +
     &                 pre1 (ixyz) * alpha (iprm) +
     &                 pre2 (ixyz) * alpha (iprm) **2 )
               enddo
               iprm = iprm + 1
            enddo
c     
         enddo
      enddo
c     
      return
      end
