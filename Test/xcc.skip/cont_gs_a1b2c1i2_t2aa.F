
c This routine contracts <Bc|Ek> with T2(IJ,EA) for all A, B, and c
c under a given I, J, and k.

c NOTE:
c    Orbital I must be less than orbital J for the antisymmetry to
c wash out; therefore, unless you know what you are doing, t3_factor
c should be positive.

c INPUT
c    [yeah, right]

c ----------------------------------------------------------------------

c#define _DEBUG_CONT_GS_A1B2C1I2_T2AA

#ifdef _UNICOS
#  define B_GEMM sgemm
#else
#  define B_GEMM dgemm
#endif /* _UNICOS */

      subroutine cont_gs_a1b2c1i2_t2aa(
     &              dcore,  dscr,   dim_dscr,
     &              ispin1, ispin2,
     &              list_a1b2c1i2, dist_a1b2,
     &              add_t2ij, off_a1ltb1,
     &              i, irp_i, max_i, irp_ij,
     &              j, irp_j, max_j, irp_jk,
     &              k, irp_k, max_k, irp_ijk,
     &              t3_factor, start_t3, ndx_t3vrt
     &           )
      implicit none

c ARGUMENT LIST
      integer dim_dscr
      double precision dcore(*), dscr(dim_dscr)
      integer ispin1, ispin2
      integer list_a1b2c1i2, dist_a1b2
      integer add_t2ij, off_a1ltb1
      integer i, irp_i, max_i, irp_ij
      integer j, irp_j, max_j, irp_jk
      integer k, irp_k, max_k, irp_ijk
      double precision t3_factor
      integer start_t3, ndx_t3vrt(*)

c INTERNAL VARIABLES
      integer pq, irp_pq, max_pq
      integer e,  irp_e,  max_e, from_e,       last_e, num_e
      integer p,  irp_p,  max_p,         to_p
      integer q,  irp_q,  max_q
      integer r,  irp_r,  max_r, from_r, to_r, last_r, num_r
      logical do_irp_q, do_irp_r, do_t3
      integer size_scr, ioff, itmp
      integer ndx_ek, start_t2ijer, ndx, iadd, max_dists
      double precision dtmp_t3_factor

c COMMON BLOCKS
#include "xcc.com"
#include "sym.com"
#include "machsp.com"

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      print *, '@XCC: Entered file ',__FILE__
#endif

c ----------------------------------------------------------------------

#ifdef _ASSERT
      itmp = 0
c   o assert dim_dscr > 0
      if (dim_dscr.lt.1) then
         print *, '@CONT_GS_A1B2C1I2_T2AA: Assertion failed.'
         print *, '   dim_dscr = ',dim_dscr
         itmp = -1
      end if
c   o assert ispin1 and ispin2 are 1 or 2 and not equal
      if (.not.(((ispin1.eq.1).and.(ispin2.eq.2)).or.
     &          ((ispin1.eq.2).and.(ispin2.eq.1))    )) then
         print *, '@CONT_GS_A1B2C1I2_T2AA: Assertion failed.'
         print *, '   ispin1 = ',ispin1
         print *, '   ispin2 = ',ispin2
         itmp = -1
      end if
c   o assert i, j, and k are well-behaved
      if ((max_i.le.i).or.(i.lt.0).or.
     &    (max_j.le.j).or.(j.lt.0).or.
     &    (max_k.le.k).or.(k.lt.0)    ) then
         print *, '@CONT_GS_A1B2C1I2_T2AA: Assertion failed.'
         print *, '   i = ',i,', max_i = ',max_i
         print *, '   j = ',j,', max_j = ',max_j
         print *, '   k = ',k,', max_k = ',max_k
         itmp = -1
      end if
c   o assert i<j
      if ((irp_i.gt.irp_j).or.((irp_i.eq.irp_j).and.(i.gt.j))) then
         print *, '@CONT_GS_A1B2C1I2_T2AA: Assertion failed.'
         print *, '   i = ',1+i,'[',irp_i,']'
         print *, '   j = ',1+j,'[',irp_j,']'
         itmp = -1
      end if
c   o assert irp_i, irp_j, and irp_k are in [1,nirrep]
      if ((irp_i  .lt.1).or.(nirrep.lt.irp_i  ).or.
     &    (irp_j  .lt.1).or.(nirrep.lt.irp_j  ).or.
     &    (irp_k  .lt.1).or.(nirrep.lt.irp_k  ).or.
     &    (irp_ij .lt.1).or.(nirrep.lt.irp_ij ).or.
     &    (irp_jk .lt.1).or.(nirrep.lt.irp_jk ).or.
     &    (irp_ijk.lt.1).or.(nirrep.lt.irp_ijk)    ) then
         print *, '@CONT_GS_A1B2C1I2_T2AA: Assertion failed.'
         print *, '   irp_i   = ',irp_i
         print *, '   irp_j   = ',irp_j
         print *, '   irp_k   = ',irp_k
         print *, '   irp_ij  = ',irp_ij
         print *, '   irp_jk  = ',irp_jk
         print *, '   irp_ijk = ',irp_ijk
         itmp = -1
      end if
c   o assert dcore pointers are in [1,inf)
      if ((add_t2ij.lt.1).or.(start_t3.lt.1)) then
         print *, '@CONT_GS_A1B2C1I2_T2AA: Assertion failed.'
         print *, '   add_t2ij = ',add_t2ij
         print *, '   start_t3 = ',start_t3
         itmp = -1
      end if
      if (itmp.ne.0) then
         call xcc_errex
      end if
c   o warn if t3_factor = 0
      if (t3_factor.eq.0.0d0) then
         print *, '@CONT_GS_A1B2C1I2_T2AA: Warning.'
         print *, '   (t3_factor == 0.0d0) -> do nothing'
      end if
#endif /* _ASSERT */

      if (t3_factor.ne.0.0d0) then
         dtmp_t3_factor = -t3_factor
      else
#ifdef _TRACE_XCC
         print *, '@XCC: Leaving file ',__FILE__
#endif
         return
      end if

c The plan is:
c   o loop over irreps of R in T2(IJ,ER)
c      o loop over orbital E
c         - load as many <Pq|Ek> lists as possible (all Pq for each E)
c         - if (ispin1.eq.2) then <Qp|Ke> -> <pQ|eK>
c         o loop over orbital R in T2(IJ,ER)
c            - create T2(IJ,ER) (all E for each R)
c            - gemm <Pq|Ek> with T2(IJ,ER)
c            o loop over irreps of q in T3(IJk,PqR)
c               - add in the T3(IJk,PqR) contributions

c   o loop over irreps of R in T2(IJ,ER)
      do irp_r = 1, nirrep
c      o check for T2 amplitudes
         irp_e = dirprd(irp_r,irp_ij)
         max_r = vrt(irp_r,ispin1)
         if (irp_e.ne.irp_r) then
            max_e = vrt(irp_e,ispin1)
            do_irp_r = ((max_e.ne.0).and.(max_r.ne.0))
         else
            max_e = max_r
            do_irp_r = (max_e.gt.1)
         end if
         if (do_irp_r) then
c         o check for integrals
            irp_pq = dirprd(irp_e,irp_k)
            max_pq = irpdpd(irp_pq,dist_a1b2)
            do_irp_r = (max_pq.ne.0)
         end if
      if (do_irp_r) then

c      o initialize values for these irreps
         last_e = max_e-1
         last_r = max_r-1
c        size_scr : scratch space minus at least one T3(IJk,PqR) product
c                   matrix and T2(IJ,ER) factor vector
         size_scr = dim_dscr - ( max_pq + max_e )
         itmp = mod(size_scr,max_pq)
         max_dists = ( size_scr - itmp ) / max_pq
         if (max_dists.eq.0) then
            print *, '@CONT_GS_A1B2C1I2_T2AA: Not enough memory.'
            print *, '                        At least ',
     &               lshift((max_pq-size_scr),3),
     &               ' more bytes are needed.'
            call xcc_errex
         end if

c      o point to the first |Ek> (or |Ke>) list
         if (ispin1.eq.1) then
            itmp = nirrep*(irp_k-1) + irp_e
            ndx_ek = 1 + off_pq(itmp,h_Offsets_Ai)+ ( max_e * k )
         else
            itmp = nirrep*(irp_e-1) + irp_k
            ndx_ek = 1 + off_pq(itmp,h_Offsets_Ia) + k
         end if

c      o loop over orbital E in <Pq|Ek>
         from_e = 0
         do while (from_e.lt.max_e)

c         o get as many <Pq|Ek> lists as possible
            num_e = max_e-from_e
            num_e = min(num_e,max_dists)
            if (ispin1.eq.1) then
               call getlst(dscr(1),ndx_ek,num_e,2,irp_pq,list_a1b2c1i2)
               ndx_ek = ndx_ek + num_e
            else
               do e = 0, num_e-1
                  call getlst(dscr(1+(e*max_pq)),ndx_ek,1,2,irp_pq,
     &                        list_a1b2c1i2)
                  ndx_ek = ndx_ek + max_k
               end do
               call transp_dist(dscr(1),vrt(1,1),vrt(1,2),
     &                          irp_pq,num_e,
     &                          .false.)
            end if

c         o reset the T2(IJ,ER) pointer
c           This is the scratch address of the temporary factor.
            start_t2ijer = 1 + ( num_e * max_pq )

c         o loop over orbital R in T2(IJ,ER)
            from_r = 0
            do while (from_r.lt.max_r)

c            o We know there is enough memory to process at least one R orbital
c              (from size_scr), so we do not need to check for memory.
c              ndx : (temp) a T2, T3 column pair (factor, product)
               ndx   = max_e + max_pq
               itmp  = 1 + dim_dscr - start_t2ijer
               ioff  = mod(itmp,ndx)
c              to_r : (temp) the number of column pairs that can fit in dscr
               to_r  = ( itmp - ioff ) / ndx
               num_r = min(to_r,(max_r-from_r))
               to_r  = from_r + num_r - 1

c            o create T2(IJ,ER) (all E for each R)
               call load_t2_virts(
     &                 dcore(add_t2ij), dscr(start_t2ijer),
     &                 irp_e, max_e,               ispin1,
     &                 irp_r, max_r, from_r, to_r, ispin1, off_a1ltb1
     &              )

c            o gemm <Pq|Ek> with T2(IJ,ER)
               iadd = start_t2ijer + from_e
               ndx  = start_t2ijer + ( max_e * num_r )
               call B_GEMM('N','N',max_pq,num_r,num_e,
     &                     dtmp_t3_factor,dscr(1),   max_pq,
     &                                    dscr(iadd),max_e,
     &                     0.0d0,         dscr(ndx), max_pq)

c            o loop over T3(IJk,PqR) contributions
               do r = from_r, to_r

c ----------------------------------------------------------------------

c   o add in the T3 contributions
#include "put_t3aba.F.src"

c ----------------------------------------------------------------------

c            o close the loop over T3(IJk,PqR) contributions
c              end do r = from_r, to_r
               end do

c         o increment and close the loop over orbital R in T2(IJ,ER)
               from_r = from_r + num_r
c           end do while (from_r.lt.max_r)
            end do

c      o increment and close the loop over orbital E in <Pq|Ek>
            from_e = from_e + num_e
c        end do while (from_e.lt.max_e)
         end do

c   o close the loop over irreps of R in T2(IJ,ER)
c     end if (do_irp_r)
      end if
c     end do irp_r = 1, nirrep
      end do

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      print *, '@XCC: Leaving file ',__FILE__
#endif

      return
      end

