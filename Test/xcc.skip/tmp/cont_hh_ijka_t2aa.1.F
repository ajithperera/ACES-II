
c This routine contracts <MC||jk> with T2(Mi,AB) for all A, B, C, and M
c under a given i, j, and k.

c NOTE:
c    Orbitals i, j, and k do not necessarily follow the convention of
c the calling routine, which means the orbitals are arbitrarily ordered.
c This matters when the T2 amplitudes are masked with a particular sign
c to ensure proper anti-symmetry.
c    The notation of comments is not entirely consistent. Capital letters
c denote all orbitals of a particular irrep while lowercase letters
c represent a single orbital in that irrep.

c NOTE:
c    If you are unlucky enough to have to debug this, then there are a
c few things which you need to know...
c 1) Pay CLOSE attention to the scope of each variable. itmp, iadd,
c    irp_x, irp_y, irp_z, etc. are used all over the place as temporary
c    variables. If this were C (instead of Fortran), then each segment
c    needing a temp integer would be surrounded by braces.
c 2) [more to come as the code gets more complex]

c INPUT

c ----------------------------------------------------------------------

#define _PRINT_TERMS
#ifdef _UNICOS
#  define B_GEMM sgemm
#else
#  define B_GEMM dgemm
#endif /* _UNICOS */

      subroutine cont_hh_ijka_t2aa(
     &              dcore,    dscr,           dim_dscr,
     &              ispin,    dist_altb,      dist_iltj,      dist_ai,
     &              i, irp_i, max_i, irp_ij,  irp_i_eq_irp_j,
     &              j, irp_j, max_j, irp_jk,  irp_j_eq_irp_k,
     &              k, irp_k, max_k, irp_ijk, flip_sign,
     &              add_wjk,  start_t3,       ndx_t3vrt,      off_t2occ
     &           )
      implicit none

c ARGUMENT LIST
      integer dim_dscr
      double precision dcore(*), dscr(dim_dscr)
      integer ispin, dist_altb, dist_iltj, dist_ai
      integer i, irp_i, max_i, irp_ij
      logical                           irp_i_eq_irp_j
      integer j, irp_j, max_j, irp_jk
      logical                           irp_j_eq_irp_k
      integer k, irp_k, max_k, irp_ijk
      logical                           flip_sign
      integer add_wjk, start_t3, ndx_t3vrt(36), off_t2occ(36)

c INTERNAL VARIABLES
      integer pq, irp_pq, max_pq
      integer     irp_a,  max_a
      integer     irp_b,  max_b
      integer     irp_c,  max_c
      integer m,  irp_m,  max_m
      integer p,  irp_p,  max_p, from_p, to_p
      integer q,  irp_q,  max_q, from_q, to_q
      integer r,  irp_r,  max_r, from_r, to_r
      integer     irp_x,  max_x
      integer     irp_y,  max_y,               irp_yz
      integer     irp_z,  max_z
      logical do_irp_r,       do_irp_q,       do_t3
      logical irp_m_eq_irp_i, irp_p_eq_irp_q
      integer irrep, scrsize, ioff, itmp
      integer start_w, start_t2, iadd, ndx

c COMMON BLOCKS
#include "xcc.com"
#include "sym.com"

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      write(*,*) '@XCC: Entered file ',__FILE__
#endif

#ifdef _PRINT_TERMS
      WRITE(*,*) ispin,': T2(M',irp_i,'[',i+1,'],PQ) * ',
     &          '<MR||',irp_j,'[',j+1,']',irp_k,'[',k+1,']>'
#endif

c ----------------------------------------------------------------------

c The plan is:
c   o loop over irrep R of <MR||jk>
c      - create the rectangular matrix T2(Mi,PQ) for all M and P<Q
c      o loop over r orbitals
c         - gemm T2(Mi,PQ) with <Mr||jk>
c         o loop over irrep Q
c            - add in T3(ijk,PQr)

      start_w = add_wjk

c   o loop over R irreps of <MR||jk>
      do irp_r = 1, nirrep
c      o check for integrals
         irp_m = dirprd(irp_r,irp_jk)
         max_m = pop(irp_m,ispin)
         max_r = vrt(irp_r,ispin)
         do_irp_r = ((max_r.ne.0).and.(max_m.ne.0))
c      o check for T2s
         if (do_irp_r) then
            irp_pq = dirprd(irp_m,irp_i)
            max_pq = irpdpd(irp_pq,dist_altb)
            if (irp_pq.eq.1) then
               itmp = max_m
            else
               itmp = max_m+1
            end if
            do_irp_r = ((max_pq.ne.0).and.(itmp.ne.1))
            if (dirprd(irp_r,irp_pq).ne.irp_ijk) then
c            o Do we need to check this? I have never seen it stop.
               write(*,*) 'BROKEN SYMMETRY: RxPQ != IJK'
               call crapso
               stop 1
            end if
         end if
      if (do_irp_r) then

c      o save space for the product vector T2(Mi,PQ)*<Mr||jk> at dscr(1)
         start_t2 = 1 + max_pq

c      o create the block T2(Mi,PQ) for all M and P<Q
         ndx = start_t2
         if (irp_m.eq.irp_i) then
            irp_m_eq_irp_i = .true.
            itmp = rshift(irp_i*(irp_i-1),1) + irp_i
            ioff = off_t2occ(itmp) + rshift(i*(i-1),1)
            iadd = ndx_t2(irp_pq,ispin) + ( ioff * max_pq )
            if (i.ne.0) then
c            o load all m<i
               itmp = i * max_pq
               do pq = 0, itmp-1
                  dscr(ndx+pq) = dcore(iadd+pq)
               end do
               iadd = iadd + lshift(itmp,1)
               ndx  = ndx  + itmp
            end if
c         o skip T2(ii,pq)
            do pq = 0, max_pq-1
               dscr(ndx+pq) = 0.0d0
            end do
            if (i.ne.(max_i-1)) then
c            o load all m>i
               do m = 1, max_i-(i+1)
                  ndx = ndx + max_pq
                  do pq = 0, max_pq-1
                     dscr(ndx+pq) = -dcore(iadd+pq)
                  end do
                  iadd = iadd + ( (m+i) * max_pq )
               end do
            end if
c        else (irp_m.ne.irp_i)
         else
            irp_m_eq_irp_i = .false.
            if (irp_m.lt.irp_i) then
               itmp = rshift(irp_i*(irp_i-1),1) + irp_m
               ioff = off_t2occ(itmp) + ( i * max_m )
               iadd = ndx_t2(irp_pq,ispin) + ( ioff * max_pq )
               itmp = max_m * max_pq
               do pq = 0, itmp-1
                  dscr(ndx+pq) = dcore(iadd+pq)
               end do
            else
               itmp = rshift(irp_m*(irp_i-m),1) + irp_i
               ioff = off_t2occ(itmp) + i
               iadd = ndx_t2(irp_pq,ispin) + ( ioff * max_pq )
               itmp = max_i * max_pq
               do m = 1, max_m
                  do pq = 0, max_pq-1
                     dscr(ndx+pq) = -dcore(iadd+pq)
                  end do
                  ndx  = ndx  + max_pq
                  iadd = iadd + itmp
               end do
            end if
         end if
c      o flip the sign if j<i<k
         if (flip_sign) then
            itmp = max_pq * max_m
            do pq = 0, itmp-1
               dscr(start_t2+pq) = -dscr(start_t2+pq)
            end do
         end if
c      WRITE(*,*) 'I: ',irp_i,'[',i+1,'], M:',irp_m,'[1-',max_m,']'
c      do pq = 0, max_pq-1
c         WRITE(*,*) pq, (dscr(start_t2+pq+m*max_pq),m=0,max_m-1)
c      end do

c      o loop over virtual orbitals of <MR||jk>
         do r = 0, max_r-1

c         o GEMM T2(Mi,PQ) with <Mr||jk>
            if (max_m.eq.1) then
               do pq = 0, max_pq-1
                  dscr(1+pq) = dscr(start_t2+pq) * dcore(start_w)
               end do
            else
               call B_GEMM('N','N',max_pq,1,max_m,
     &                     1.0d0,dscr(start_t2),max_pq,
     &                           dcore(start_w),max_m,
     &                     0.0d0,dscr(1),       max_pq)
            end if
            start_w = start_w + max_m

c ----------------------------------------------------------------------

c   o At this point, we don't know anything about the relationship between
c     irp_r and irp_pq. This means we need to be very careful about analyzing
c     cases where irp_r=irp_p and irp_q=irp_r, since we may not have enough
c     virtual orbitals in irp_r to get a non-zero T3. The address variable
c     'ndx' will step sequentially through the list of T3 pieces.
c     BEWARE: This variable has a scope much larger than any other temporary
c     variable.

      ndx=1

c   o loop over virtual irreps of T3(ijk,PQr)
      do irp_q = 1, nirrep
         irp_p = dirprd(irp_q,irp_pq)
         if (irp_p.eq.irp_q) then
            max_p = vrt(irp_p,ispin)
            max_q = max_p
            do_irp_q = (max_q.gt.1)
         else
            if (irp_p.lt.irp_q) then
               max_p = vrt(irp_p,ispin)
               max_q = vrt(irp_q,ispin)
               do_irp_q = ((max_p.ne.0).and.(max_q.ne.0))
            else
               do_irp_q = .false.
            end if
         end if
      if (do_irp_q) then

c      o skip irp_q if there aren't enough virtual orbitals for T3
c        (itmp is the number of p<q pairs to skip if we fail do_t3)
         if (irp_p.eq.irp_q) then
            irp_p_eq_irp_q = .true.
            if (irp_r.eq.irp_p) then
               itmp  = 1
               do_t3 = (max_q.ne.2)
            else
               do_t3 = .true.
            end if
         else
            irp_p_eq_irp_q = .false.
            if (irp_r.eq.irp_p) then
                  itmp  = max_q
                  do_t3 = (max_p.ne.1)
            else
               if (irp_r.eq.irp_q) then
                  itmp  = max_p
                  do_t3 = (max_q.ne.1)
               else
                  do_t3 = .true.
               end if
            end if
         end if
         if (do_t3) then

c         o point to the irreps whose T3s are in the buffer
            irp_c = max(irp_p,irp_q,irp_r)
            irp_a = min(irp_p,irp_q,irp_r)
            irp_b = irp_p + irp_q + irp_r - irp_a - irp_c
            itmp  = rshift(irp_c*(irp_c-1),1) + irp_b
            iadd  = ndx_t3vrt(itmp)

c         o Now, we need to consider how p and q run across r.

            if (irp_r.lt.irp_p) then
               if (irp_p_eq_irp_q) then
                  itmp = rshift(max_q*(max_q-1),1)
               else
                  itmp = ( max_p * max_q )
               end if
               iadd = iadd + r
               do pq = 0, itmp-1
                  dcore(iadd) = dcore(iadd) + dscr(ndx+pq)
                  iadd = iadd + max_r
               end do
               ndx = ndx + itmp
            else

               if (irp_q.lt.irp_r) then
                  if (irp_p_eq_irp_q) then
                     itmp = rshift(max_q*(max_q-1),1)
                  else
                     itmp = ( max_q * max_p )
                  end if
                  iadd = iadd + ( r * itmp )
                  do pq = 0, itmp-1
                     dcore(iadd+pq) = dcore(iadd+pq) + dscr(ndx+pq)
                  end do
                  ndx = ndx + itmp
               else

                  if (irp_p_eq_irp_q) then
                     to_q = max_q-1
                     if (r.eq.0) then
                        ndx = ndx + 1
                        do q = 2, to_q
                           ndx = ndx + 1
                           do p = 1, q-1
                              dcore(iadd) = dcore(iadd) + dscr(ndx)
                              iadd = iadd + p
                              ndx  = ndx  + 1
                           end do
                        end do
                     else
                        ioff = rshift(r*(r-1),1)
                        if (r.gt.1) then
c                           iadd = iadd + ( r * (r-1) * (r-2) ) / 6
                           iadd = iadd + ( (r*r) + (-3*r) + 2 ) * r / 6
                           do pq = 0, ioff-1
                              dcore(iadd+pq) =   dcore(iadd+pq)
     &                                         + dscr(ndx+pq)
                           end do
                           iadd = iadd + ioff
                           ndx  = ndx  + ioff + r
                        else
                           ndx  = ndx  + 1
                        end if
                        if (r.ne.to_q) then
                           iadd = iadd + ioff
                           do p = 0, r-1
                              dcore(iadd+p) =   dcore(iadd+p)
     &                                        - dscr(ndx+p)
                           end do
                           iadd = iadd + r
                           ndx  = ndx  + r + 1
                           if (r.ne.(to_q-1)) then
                              do q = r+2, to_q
                                 iadd = iadd + ioff
                                 do p = 0, r-1
                                    dcore(iadd+p) =   dcore(iadd+p)
     &                                              - dscr(ndx+p)
                                 end do
                                 ndx  = ndx  + r + 1
                                 iadd = iadd + r
                                 do p = 1, q-(r+1)
                                    iadd = iadd + r
                                    dcore(iadd) =   dcore(iadd)
     &                                            + dscr(ndx)
                                    ndx  = ndx  + 1
                                    iadd = iadd + p
                                 end do
                              end do
c                          end if (r.ne.(to_q-1))
                           end if
c                       end if (r.ne.to_q)
                        end if
c                    end if (r.eq.0)
                     end if
                  else

c                 else (irp_p.ne.irp_q)
                     if (irp_p.eq.irp_r) then
                        to_p = max_p-1
                        ioff = rshift(r*(r-1),1)
                        do q = 0, max_q-1
                           if (r.ne.0) then
                              iadd = iadd + ioff
                              do p = 0, r-1
                                 dcore(iadd+r) =   dcore(iadd+r)
     &                                           - dscr(ndx+r)
                              end do
                              ndx  = ndx  + r + 1
                              iadd = iadd + r
                           end if
                           if (r.ne.to_p) then
                              do p = 1, to_p-r
                                 iadd = iadd + r
                                 dcore(iadd) = dcore(iadd) + dscr(ndx)
                                 ndx  = ndx  + 1
                                 iadd = iadd + p
                              end do
                           end if
                        end do
                     else

                        if (irp_r.eq.irp_q) then
                           if (r.ne.0) then
                              iadd =   iadd
     &                               + ( max_p * rshift(r*(r-1),1) )
                              itmp = max_p * r
                              do pq = 0, itmp-1
                                 dcore(iadd+pq) =   dcore(iadd+pq)
     &                                            + dscr(ndx+pq)
                              end do
                              ndx  = ndx  + itmp + max_p
                              iadd = iadd + lshift(itmp,1)
                           end if
                           to_q = max_q-1
                           if (r.ne.to_q) then
                              to_p = max_p-1
                              do q = 1, to_q-r
                                 do p = 0, to_p
                                    dcore(iadd+p) =   dcore(iadd+p)
     &                                              - dscr(ndx+p)
                                 end do
                                 iadd = iadd + ( (q+r) * max_p )
                                 ndx  = ndx  + max_p
                              end do
                           end if
                        else

c                       else ((irp_p.lt.irp_r).and.(irp_r.lt.irp_q)) FINALLY
                           ioff = max_p * max_r
                           iadd = iadd + ( r * max_p )
                           to_p = max_p-1
                           do q = 1, max_q
                              do p = 0, to_p
                                 dcore(iadd+p) =   dcore(iadd+p)
     &                                           - dscr(ndx+p)
                              end do
                              ndx  = ndx  + max_p
                              iadd = iadd + ioff
                           end do
c                       end if (irp_r.eq.irp_q)
                        end if

c                    end if (irp_p.eq.irp_r)
                     end if
c                 end if (irp_p_eq_irp_q)
                  end if
c              end if (irp_q.lt.irp_r)
               end if
c           end if (irp_r.lt.irp_p)
            end if

c      o close the T3 conditional
         else
c         o itmp is the number of p<q pairs to skip
            ndx = ndx + itmp
c        end if (do_t3)
         end if

c   o close the loop over virtual irreps
c     end if (do_irp_q)
      end if
c     end do irp_q = 1, nirrep
      end do

c ----------------------------------------------------------------------

c      o get the next M block and close the loop over virtual orbitals
c        end do r = 0, max_r-1
         end do

c   o close the loop over R irreps of <MR||jk>
c     else (.not.do_irp_r)
c      o increment start_w even if there are no T2s
         start_w = start_w + ( max_m * max_r )
c     end if (do_irp_r)
      end if
c     end do irp_r = 1, nirrep
      end do

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      write(*,*) '@XCC: Leaving file ',__FILE__
#endif

      return
      end

