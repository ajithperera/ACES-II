
c This routine contracts <AB||Ei> with T2(jk,CE) for all A, B, C, and E
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
c    double precision dcore(*) : core space for doubles
c    double precision dscr(*)  : first free scratch address
c    integer          dim_dscr : the dimension of the scratch array

c    integer ispin     : the spin case to contract
c    integer list_abci : the handle for <AB||EI>
c    integer dist_altb : the handle for irpdpd type a<b

c    [orbital info]
c    logical flip_sign : TRUE if j<i<k

c    integer add_t2jk      : the starting address of T2(jk,*)
c    integer start_t3      : the starting address of T3(ijk,*)
c    integer ndx_t3vrt(36) : the starting addresses of T3(ijk,ABC)
c                            in upper-triangular packed form

c ----------------------------------------------------------------------

#define _PRINT_TERMS

      subroutine cont_hh_abci_t2aa(
     &              dcore,    dscr,           dim_dscr,
     &              ispin,    list_abci,      dist_altb,
     &              i, irp_i, max_i, irp_ij,  irp_i_eq_irp_j,
     &              j, irp_j, max_j, irp_jk,  irp_j_eq_irp_k,
     &              k, irp_k, max_k, irp_ijk, flip_sign,
     &              add_t2jk, start_t3,       ndx_t3vrt
     &           )
      implicit none

c ARGUMENT LIST
      integer dim_dscr
      double precision dcore(*), dscr(dim_dscr)
      integer ispin, list_abci, dist_altb
      integer i, irp_i, max_i, irp_ij
      logical                           irp_i_eq_irp_j
      integer j, irp_j, max_j, irp_jk
      logical                           irp_j_eq_irp_k
      integer k, irp_k, max_k, irp_ijk
      logical                           flip_sign
      integer add_t2jk, start_t3, ndx_t3vrt(36)

c INTERNAL VARIABLES
      integer    irp_pq, max_pq
      integer    irp_a,  max_a
      integer    irp_b,  max_b
      integer    irp_c,  max_c
      integer e, irp_e,  max_e,         to_e, num_ei
      integer p, irp_p,  max_p,         to_p
      integer q, irp_q,  max_q, from_q, to_q
      integer r, irp_r,  max_r, from_r, to_r
      integer    irp_x,  max_x
      integer    irp_y,  max_y,               irp_yz
      integer    irp_z,  max_z
      logical do_irp_q,       do_irp_e,       do_t3
      logical irp_p_eq_irp_q, irp_r_eq_irp_e
      integer irrep, size_scr, ioff, itmp
      integer add_w, ndx_ei, start_t2jkre, ndx, iadd

c PARAMETERS
      double precision thresh_w
      parameter (thresh_w = 1.0d-15)

c COMMON BLOCKS
#include "xcc.com"
#include "sym.com"

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      write(*,*) '@XCC: Entered file ',__FILE__
#endif

#ifdef _PRINT_TERMS
      WRITE(*,*) ispin,': <PQ||E',irp_i,'[',i+1,']> * T2(',
     &           irp_j,'[',j+1,'] ',irp_k,'[',k+1,'],*)'
#endif

c ----------------------------------------------------------------------

c The plan is:
c   o loop over the common irrep E
c      - locate <PQ||Ei> and T2(jk,RE)
c      o loop over e orbitals
c         - load as many <PQ||ei> lists as possible
c         o loop over p<q orbitals
c            - contracting each integral with T2(jk,Re)
c            - add in T3(ijk,pqR)

c   o loop over E irreps
      do irp_e = 1, nirrep
c      o check for T2 amplitudes
         irp_r = dirprd(irp_e,irp_jk)
         max_r = vrt(irp_r,ispin)
         if (irp_r.eq.irp_e) then
            max_e = max_r
            do_irp_e = (max_e.gt.1)
         else
            max_e = vrt(irp_e,ispin)
            do_irp_e = ((max_r.ne.0).and.(max_e.ne.0))
         end if
c      o check for integrals
         if (do_irp_e) then
            irp_pq = dirprd(irp_e,irp_i)
            max_pq = irpdpd(irp_pq,dist_altb)
            do_irp_e = (max_pq.ne.0)
            if (dirprd(irp_r,irp_pq).ne.irp_ijk) then
c            o Do we need to check this? I have never seen it stop.
               write(*,*) 'BROKEN SYMMETRY: RxPQ != IJK'
               call crapso
               stop 1
            end if
         end if
      if (do_irp_e) then

c      o point to the first T2(jk,RE) amplitude
         start_t2jkre = add_t2jk
         if (irp_r.eq.irp_e) then
            irp_r_eq_irp_e = .true.
            if (irp_e.gt.1) then
               do irrep = 1, irp_e-1
                  itmp = vrt(irrep,ispin)
                  start_t2jkre = start_t2jkre + rshift(itmp*(itmp-1),1)
               end do
            end if
         else
            irp_r_eq_irp_e = .false.
            itmp = max(irp_r,irp_e)
            do irp_b = 1, itmp-1
               irp_a = dirprd(irp_b,irp_jk)
               if (irp_a.lt.irp_b) then
                  max_a = vrt(irp_a,ispin)
                  max_b = vrt(irp_b,ispin)
                  start_t2jkre = start_t2jkre + ( max_a * max_b )
               end if
            end do
         end if

c      o point to the first |Ei> list
         ndx_ei = 1 + ( i * max_e )
         if (irp_i.gt.1) then
            do irp_y = 1, irp_i-1
               irp_x = dirprd(irp_y,irp_pq)
               max_y = pop(irp_y,ispin)
               max_x = vrt(irp_x,ispin)
               ndx_ei = ndx_ei + ( max_x * max_y )
            end do
         end if

c      o initialize values for this irrep
c        num_ei : the number of |Ei> lists in memory
         num_ei = 0
c        size_scr : scratch space minus 2 T2(jk,Re) vectors
         size_scr = dim_dscr - lshift(max_r,1)

c      o loop over virtual orbitals
         to_e = max_e-1
         do e = 0, to_e

c         o get the next |ei> distribution (could be another thread)
            if (num_ei.eq.0) then
c            o get as many <PQ||Ei> vectors as possible
               itmp = mod(size_scr,max_pq)
               itmp = ( size_scr - itmp ) / max_pq
               num_ei = min((max_e-e),itmp)
               if (num_ei.eq.0) then
                  write(*,*) '@MAKE_T3AAA_UHF: There is not enough ',
     &                       'memory for another <PQ||ei> distribution.'
                  write(*,*) '                 At least ',
     &                       (max_pq-size_scr),' floats are needed.'
                  call crapso
                  stop 1
               end if
c            o reset the current <pq||ei> address
               add_w = lshift(max_r,1) + 1
c            o load the integrals
               call getlst(dscr(add_w),ndx_ei,num_ei,2,irp_pq,list_abci)
c            o point to the next |ei> list
               ndx_ei = ndx_ei + num_ei
            end if
            num_ei = num_ei - 1

c         o put T2(jk,Re) into the second T2 vector
            iadd         = start_t2jkre
            start_t2jkre = start_t2jkre + max_r
            if (irp_r_eq_irp_e) then
               ndx = max_r + 1
               if (e.ne.0) then
c               o load all r<e
                  iadd = iadd + rshift(e*(e-1),1)
                  do r = 0, e-1
                     dscr(ndx+r) = dcore(iadd+r)
                  end do
                  iadd = iadd + lshift(e,1)
                  ndx  = ndx  + e
               end if
               dscr(ndx) = 0.0d0
               if (e.ne.to_e) then
c               o load all r>e
                  do r = 1, to_e-e
                     dscr(ndx+r) = -dcore(iadd)
                     iadd = iadd + r + e
                  end do
               end if
c           else (irp_r.ne.irp_e)
            else
               if (irp_r.lt.irp_e) then
c               o jump to column e+1 and load all r<e
                  iadd = iadd + ( e * max_r ) - 1
                  do r = 1, max_r
                     dscr(max_r+r) = dcore(iadd+r)
                  end do
               else
c               o drop to row e+1 and load all r>e
                  iadd = iadd + e
                  do r = 1, max_r
                     dscr(max_r+r) = -dcore(iadd)
                     iadd = iadd + max_e
                  end do
               end if
c           end if (irp_r_eq_irp_e)
            end if
c         o flip the sign if j<i<k
            if (flip_sign) then
               do r = 1, max_r
                  dscr(max_r+r) = -dscr(max_r+r)
               end do
            end if
c      do r = 1, max_r
c         WRITE(*,*) 'T2(',irp_j,'[',j+1,'] ',irp_k,'[',k+1,'] , ',
c     &                    irp_r,'[',r  ,'] ',irp_e,'[',e+1,']) :',
c     &              dscr(max_r+r)
c      end do
c      WRITE(*,*)
c         o done retrieving T2(jk,Re)

c         o loop over virtual irreps of <PQ|
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

c            o skip irp_q if there aren't enough virtual orbitals for T3
c              (itmp is the number of p<q pairs to skip if we fail do_t3)
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

c ----------------------------------------------------------------------

c   o loop over p<q virtual orbitals but skip small integrals
      to_p = max_p-1
      to_q = max_q-1
      if (irp_p_eq_irp_q) then
         from_q = 1
      else
         from_q = 0
      end if
      do q = from_q, to_q
         if (irp_p_eq_irp_q) to_p = (q-1)
      do p = 0, to_p
      if (thresh_w.lt.abs(dscr(add_w))) then

c      o scale the T2(jk,Re) amplitudes by <pq||ei>
         do r = 1, max_r
            dscr(r) = dscr(max_r+r) * dscr(add_w)
         end do

c      o point to the irreps whose T3s are in the buffer
         irp_c = max(irp_p,irp_q,irp_r)
         irp_a = min(irp_p,irp_q,irp_r)
         irp_b = irp_p + irp_q + irp_r - irp_a - irp_c
         itmp  = rshift(irp_c*(irp_c-1),1) + irp_b
         iadd  = ndx_t3vrt(itmp)

c      o Now, we need to consider how r runs across p and q.

         if (irp_r.lt.irp_p) then
            if (irp_p_eq_irp_q) then
               itmp = rshift(q*(q-1),1) + p
            else
               itmp = ( q * max_p ) + p
            end if
            iadd = iadd + ( itmp * max_r )
            do r = 0, max_r-1
c      WRITE(*,*) 'updating T3(ijk, ',irp_r,'[',r+1,'] ',
c     &                               irp_p,'[',p+1,'] ',
c     &                               irp_q,'[',q+1,']) at ',iadd+r,
c     &           ' with ',1+r
               dcore(iadd+r) = dcore(iadd+r) + dscr(1+r)
            end do
         else

            if (irp_q.lt.irp_r) then
               if (irp_p_eq_irp_q) then
                  ioff = rshift(max_q*(max_q-1),1)
                  iadd = iadd + rshift(q*(q-1),1) + p
               else
                  ioff = max_q * max_p
                  iadd = iadd + ( q * max_p ) + p
               end if
               do r = 1, max_r
c      WRITE(*,*) 'updating T3(ijk, ',irp_p,'[',p+1,'] ',
c     &                               irp_q,'[',q+1,'] ',
c     &                               irp_r,'[',r  ,']) at ',iadd,
c     &           ' with ',r
                  dcore(iadd) = dcore(iadd) + dscr(r)
                  iadd = iadd + ioff
               end do
            else

               if (irp_p_eq_irp_q) then
c                  iadd = iadd + ( q * (q-1) * (q-2) ) / 6
                  iadd = iadd + ( (q*q) + (-3*q) + 2 ) * q / 6
                  if (p.ne.0) then
                     iadd = iadd + rshift(p*(p-1),1)
                     do r = 0, p-1
c      WRITE(*,*) 'updating T3(ijk, ',irp_q,'[',r+1,'] ',
c     &                               irp_q,'[',p+1,'] ',
c     &                               irp_q,'[',q+1,']) at ',iadd+r,
c     &           ' with ',r+1
                        dcore(iadd+r) = dcore(iadd+r) + dscr(1+r)
                     end do
                     iadd = iadd + p
                  end if
                  if (p.ne.to_p) then
                     ndx = p+1
                     do r = 1, to_p-p
                        iadd        = iadd        + p
c      WRITE(*,*) 'updating T3(ijk, ',irp_q,'[',p+1,'] ',
c     &                               irp_q,'[',ndx+r,'] ',
c     &                               irp_q,'[',q+1,']) at ',iadd,
c     &           ' with ',ndx+r
                        dcore(iadd) = dcore(iadd) - dscr(ndx+r)
                        iadd        = iadd        + r
                     end do
                  end if
                  if (q.ne.to_q) then
                     ndx  = q+1
                     iadd = iadd + rshift(q*(q-1),1) + p
                     do r = 1, to_q-q
c      WRITE(*,*) 'updating T3(ijk, ',irp_q,'[',p+1,'] ',
c     &                               irp_q,'[',q+1,'] ',
c     &                               irp_q,'[',ndx+r,']) at ',iadd,
c     &           ' with ',ndx+r
                        dcore(iadd) = dcore(iadd) + dscr(ndx+r)
                        itmp        = (q+r)*(q+r-1)
                        iadd        = iadd + rshift(itmp,1)
                     end do
                  end if
               else

c              else (irp_p.ne.irp_q)
                  if (irp_p.eq.irp_r) then
                     iadd = iadd + ( q * rshift(max_p*(max_p-1),1) )
                     if (p.ne.0) then
                        iadd = iadd + rshift(p*(p-1),1)
                        do r = 0, p-1
c      WRITE(*,*) 'updating T3(ijk, ',irp_p,'[',r+1,'] ',
c     &                               irp_p,'[',p+1,'] ',
c     &                               irp_q,'[',q+1,']) at ',iadd+r,
c     &           ' with ',1+r
                           dcore(iadd+r) = dcore(iadd+r) + dscr(1+r)
                        end do
                        iadd = iadd + lshift(p,1)
                     end if
                     if (p.ne.to_p) then
                        ndx  = p+1
                        do r = 1, max_p-p
c      WRITE(*,*) 'updating T3(ijk, ',irp_p,'[',p+1,'] ',
c     &                               irp_p,'[',ndx+r,'] ',
c     &                               irp_q,'[',q+1,']) at ',iadd,
c     &           ' with ',ndx+r
                           dcore(iadd) = dcore(iadd) - dscr(ndx+r)
                           iadd = iadd + r + p
                        end do
                     end if
                  else

                     if (irp_r.eq.irp_q) then
                        iadd = iadd + p
                        if (q.ne.0) then
                           iadd = iadd + ( max_p * rshift(q*(q-1),1) )
                           do r = 1, q
c      WRITE(*,*) 'updating T3(ijk, ',irp_p,'[',p+1,'] ',
c     &                               irp_q,'[',r  ,'] ',
c     &                               irp_q,'[',q+1,']) at ',iadd,
c     &           ' with ',r
                              dcore(iadd) = dcore(iadd) - dscr(r)
                              iadd        = iadd + max_p
                           end do
                           iadd = iadd + ( q * max_p )
                        end if
                        if (q.ne.to_q) then
                           ndx  = q+1
                           do r = 1, max_q-q
c      WRITE(*,*) 'updating T3(ijk, ',irp_p,'[',p+1,'] ',
c     &                               irp_q,'[',q+1,'] ',
c     &                               irp_q,'[',ndx+r,']) at ',iadd,
c     &           ' with ',ndx+r
                              dcore(iadd) = dcore(iadd) + dscr(ndx+r)
                              iadd = iadd + ( (r+q) * max_p )
                           end do
                        end if
                     else

c                    else ((irp_p.lt.irp_r).and.(irp_r.lt.irp_q)) FINALLY
                        iadd = iadd + ( q * max_p * max_r ) + p
                        do r = 1, max_r
c      WRITE(*,*) 'updating T3(ijk, ',irp_p,'[',p+1,'] ',
c     &                               irp_r,'[',r  ,'] ',
c     &                               irp_q,'[',q+1,']) at ',iadd,
c     &           ' with ',r
                           dcore(iadd) = dcore(iadd) - dscr(r)
                           iadd        = iadd + max_p
                        end do
c                    end if (irp_r.eq.irp_q)
                     end if

c                 end if (irp_p.eq.irp_r)
                  end if
c              end if (irp_p_eq_irp_q)
               end if
c           end if (irp_q.lt.irp_r)
            end if
c        end if (irp_r.lt.irp_p)
         end if

c   o get the next integral and close the p<q flow-control loops
c     end if (thresh_w.lt.abs(dscr(add_w)))
      end if
         add_w = add_w + 1
c     end do p = from_p, to_p
      end do
c     end do q = from_q, to_q
      end do

c ----------------------------------------------------------------------

c            o close the T3 conditional
               else
c               o itmp is the number of p<q pairs to skip
                  add_w = add_w + itmp
c              end if (do_t3)
               end if

c         o close the <PQ| irrep flow-control loops
c           end if (do_irp_q)
            end if
c           end do irp_q = 1, nirrep
            end do
 
c      o close the |Ei> flow-control loop
c        end do e = 0, to_e
         end do

c   o close the E irrep flow-control loops
c     end if (do_irp_e)
      end if
c     end do irp_e = 1, nirrep
      end do

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      write(*,*) '@XCC: Leaving file ',__FILE__
#endif

      return
      end

