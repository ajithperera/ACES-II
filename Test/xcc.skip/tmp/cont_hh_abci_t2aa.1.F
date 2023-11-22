
c This routine contracts <AB||EI> with T2(JK,CE) for all A, B, C, and E
c under a given I, J, and K.

c NOTE:
c    Orbitals i, j, and k do not necessarily follow the convention of
c the calling routine, which means the orbitals are arbitrarily ordered.
c The only time this matters is when the T2 amplitudes are masked with
c a particular sign to ensure proper anti-symmetry.

c INPUT
c    double precision dcore(*) : core space for doubles
c    double precision dscr(*)  : first free scratch address
c    integer          dscrdim  : the dimension of the scratch array

c    integer ispin           : the spin case to contract
c    integer list_abci       : the handle for <AB||EI>
c    integer dist_altb       : the handle for irpdpd type a<b

c    integer irpdpd_t3vrt(8) : the number of a<b<c amplitudes per irrep

c ----------------------------------------------------------------------

c#define _PRINT_TERMS

      subroutine cont_hh_abci_t2aa(dcore,dscr,dscrdim,
     &                             ispin,list_abci,dist_altb,
     &                             irpdpd_t3vrt)
      implicit none

c ARGUMENT LIST
      integer dscrdim
      double precision dcore(*), dscr(dscrdim)
      integer ispin, list_abci, dist_altb, irpdpd_t3vrt(8)

c INTERNAL VARIABLES
c      integer irp_ij, irp_ik, irp_jk, irp_ijk
c      integer pq, irp_pq, max_pq
c      integer i,  irp_i,  max_i, from_i, to_i
c      integer j,  irp_j,  max_j, from_j, to_j
c      integer k,  irp_k,  max_k, from_k, to_k
c      integer m,  irp_m,  max_m,               num_m
c      integer a,  irp_a,  max_a, from_a, to_a
c      integer b,  irp_b,  max_b, from_b, to_b
c      integer c,  irp_c,  max_c, from_c, to_c
c      integer e,  irp_e,  max_e,               num_e
c      integer p,  irp_p,  max_p, from_p, to_p
c      integer q,  irp_q,  max_q, from_q, to_q
c      integer r,  irp_r,  max_r, from_r, to_r
c      integer x,  irp_x,  max_x, from_x, to_x
c      integer y,  irp_y,  max_y, from_y, to_y
c      logical do_irp_j, do_irp_q, do_irp_e, do_irp_m
c      logical irp_i_eq_irp_j, irp_j_eq_irp_k
c      logical irp_r_eq_irp_e, irp_p_eq_irp_q
c      integer irrep, scrsize, ioff, itmp, count, irpdpd_t3occ(8)
c      integer w_add, t2_start, t2jk_start, t3_start, ndx, iadd
c      double precision dtmp

c PARAMETERS
c      double precision w_thresh
c      parameter (w_thresh = 1.0d-15)

c EXTERNAL FUNCTIONS
c      double precision ddot

c COMMON BLOCKS
#include "xcc.com"
#include "method.com"
#include "sym.com"
#include "info.com"
c#include "machsp.com"

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      write(*,*) '@XCC: Entered file ',__FILE__
#endif

c ----------------------------------------------------------------------

c The plan is:
c   o load as many <PQ||EI> integrals as possible
c   o walk through the integrals contracting each one with all T2(JK,RE)
c   o immediately add the T3 contribution back into T3(IJK,*)

c   o find the T2(JK,*) amplitudes
      t2jk_start = ndx_t2(irp_jk,ispin)
      ioff = 0
      if (irp_j_eq_irp_k) then
         if (irp_k.gt.1) then
c         o jump to irp_k
            do irrep = 1, irp_k-1
               itmp = pop(irrep,ispin)
               ioff = ioff + rshift(itmp*(itmp-1),1)
            end do
         end if
         ioff = ioff + rshift(k*(k-1),1) + j
      else
         do irp_y = 1, irp_k-1
            irp_x = dirprd(irp_y,irp_jk)
            if (irp_x.lt.irp_y) then
               max_y = pop(irp_y,ispin)
               max_x = pop(irp_x,ispin)
               ioff  = ioff + ( max_x * max_y )
            end if
         end do
         ioff = ioff + ( k * max_j ) + j
      end if
      itmp       = irpdpd(irp_jk,dist_altb)
      t2jk_start = t2jk_start + ( itmp * ioff )
      WRITE(*,*) 'SKIPPING ',ioff,' A<B T2 DISTRIBUTIONS'

c   o loop over E irreps of <PQ||EI>
      do irp_e = 1, nirrep
         max_e = vrt(irp_e,ispin)
         irp_r = dirprd(irp_e,irp_jk)
         max_r = vrt(irp_r,ispin)
         do_irp_e = .true.
c         o check for T2 amplitudes
         if (irp_r.eq.irp_e) then
            if (max_r.lt.2) do_irp_e = .false.
            irp_r_eq_irp_e = .true.
         else
            if ((max_r.eq.0).or.(max_e.eq.0)) do_irp_e = .false.
            irp_r_eq_irp_e = .false.
         end if
         if (do_irp_e) then
c         o check for integrals
            irp_pq = dirprd(irp_e,irp_i)
            max_pq = irpdpd(irp_pq,dist_altb)
            if (max_pq.eq.0) do_irp_e = .false.
c            o Do we need to check this? I have never seen it stop.
            if (dirprd(irp_r,irp_pq).ne.irp_ijk) then
c               do_irp_e = .false.
               WRITE(*,*) 'BROKEN SYMMETRY: RxPQ != IJK'
               call crapso
               stop 1
            end if
         end if
      if (do_irp_e) then

c      o initialize values for this irrep
c         o the number of |EI> lists in memory
         num_e = 0
c         o scratch space minus 2 T2(JK,RE) vectors
         scrsize = dscrdim - lshift(max_r,1)
c         o where the masked T2 amplitudes start
         t2_start = max_r + 1
c         o loop over |EI> distributions
         do e = 0, max_e-1

c         o get T2(JK,RE) for all r
            iadd = t2jk_start
            if (irp_r_eq_irp_e) then
               if (irp_e.gt.1) then
c               o jump to irp_e
                  do irrep = 1, irp_e-1
                     itmp = vrt(irrep,ispin)
                     iadd = iadd + rshift(itmp*(itmp-1),1)
                  end do
               end if
               ndx = 1
               if (e.ne.0) then
c               o load all r<e
                  iadd = iadd + rshift(e*(e-1),1)
                  do r = 0, e-1
                     dscr(ndx+r) = dcore(iadd+r)
                  end do
                  ndx  = ndx  + e
                  iadd = iadd + lshift(e,1)
               end if
c               o skip T2(JK,EE)
               dscr(ndx) = 0.0d0
               ndx = ndx + 1
               if (e.ne.(max_e-1)) then
c               o load all r>e
                  do r = 0, max_e-(e+1)
                     dscr(ndx+r) = -dcore(iadd)
                     iadd = iadd + (e+1) + r
                  end do
               end if
c           else (irp_r.ne.irp_e)
            else
c            o jump to max(irp_r,irp_e)
               do irp_b = 1, (max(irp_r,irp_e)-1)
                  irp_a = dirprd(irp_b,irp_jk)
                  if (irp_a.lt.irp_b) then
c                  o who cares about 0 virtuals for this op
                     max_b = vrt(irp_b,ispin)
                     max_a = vrt(irp_a,ispin)
                     iadd  = iadd + ( max_a * max_b )
                  end if
               end do
               if (irp_r.lt.irp_e) then
c               o jump to column e+1 and load all r<e
                  iadd = iadd + ( e * max_r )
                  do r = 0, max_r-1
                     dscr(1+r) = dcore(iadd+r)
                  end do
               else
c               o drop to row e+1 and load all r>e
                  iadd = iadd + e
                  do r = 1, max_r
                     dscr(r) = dcore(iadd)
                     iadd = iadd + max_e
                  end do
                  do r = 1, max_r
                     dscr(r) = -dscr(r)
                  end do
               end if
c           end if (irp_r_eq_irp_e)
            end if
c         o done retrieving T2(JK,RE) for all r

c         o get the next |EI> distribution (could be another thread)
            if (num_e.eq.0) then
c            o get as many <PQ||EI> vectors as possible
               num_e = min((max_e-e),mod(scrsize,max_pq))
               if (num_e.eq.0) then
                  write(*,*) '@MAKE_T3AAA_UHF: There is not enough ',
     &                       'memory for another <PQ||EI> distribution.'
                  call crapso
                  stop 1
               end if
c               o reset the current W address
               w_add = lshift(max_r,1) + 1
               WRITE(*,*) 'Getting ',num_e,' <ab||ci> lists for ',
     &                    irp_e,'(',e+1,')',' under ',irp_i,'(',i+1,')',
     &                    ' <> GETLST :',e+1,num_e,2,irp_pq,list_abci
               call getlst(dscr(w_add),(e+1),num_e,2,irp_pq,list_abci)
            end if
            num_e = num_e - 1

c ----------------------------------------------------------------------

c   o loop over virtual irreps
      do irp_q = 1, nirrep
         max_q = vrt(irp_q,ispin)
         irp_p = dirprd(irp_q,irp_pq)
         max_p = vrt(irp_p,ispin)
         do_irp_q = .false.
         if (irp_p.eq.irp_q) then
            if (max_p.gt.1) do_irp_q = .true.
            irp_p_eq_irp_q = .true.
         else
            if ((irp_p.lt.irp_q).and.
     &          (max_p.ne.0).and.(max_q.ne.0) ) do_irp_q = .true.
            irp_p_eq_irp_q = .false.
         end if
      if (do_irp_q) then

c      o loop over p<q virtual orbitals
         to_p   = max_p-1
         from_p = 0
         to_q   = max_q-1
         from_q = 0
         if (irp_p_eq_irp_q) from_q = 1
         do q = from_q, to_q
            if (irp_p_eq_irp_q) to_p = (q-1)
         do p = from_p, to_p
c         o do work IFF the integral is large enough
         if (w_thresh.lt.abs(dscr(w_add))) then

c         o copy the T2(JK,RE) amplitudes
            do r = 0, max_r-1
               dscr(t2_start+r) = dscr(1+r)
            end do

c         o mask the sign and eliminate R=P and R=Q
            if ((irp_p.le.irp_r).and.(irp_r.le.irp_q)) then
               to_r   = max_r-1
               from_r = 0
               if (irp_r.eq.irp_p) then
                  dscr(t2_start+p) = 0.0d0
                  from_r = p
               end if
               if (irp_r.eq.irp_q) then
                  dscr(t2_start+q) = 0.0d0
                  to_r = q
               end if
               do r = from_r, to_r
                  dscr(t2_start+r) = -dscr(t2_start+r)
               end do
            end if
            if ((irp_j.le.irp_i).and.(irp_i.le.irp_k)) then
               if ((j.lt.i).and.(i.lt.k)) then
c               o flip the sign if j<i<k
                  do r = 0, max_r-1
                     dscr(t2_start+r) = -dscr(t2_start+r)
                  end do
               end if
            end if

c         o scale by <PQ||EI>
            do r = 0, max_r-1
               dscr(t2_start+r) =   dscr(t2_start+r)
     &                            * dscr(w_add)
            end do

c         o add these terms to the proper T3(IJK,ABC) amplitude
            do r = 0, max_r-1
               dtmp = dtmp + dscr(t2_start+r)
            end do
c      do r = 0, max_r-1
c      WRITE(*,'(A,6(I1,A,I2,A),F16.12)') 
c     &           'T3( ',irp_i,'[',i+1,'] ',
c     &                  irp_j,'[',j+1,'] ',
c     &                  irp_k,'[',k+1,'] , ',
c     &                  irp_p,'[',p+1,'] ',
c     &                  irp_r,'[',r+1,'] ',
c     &                  irp_q,'[',q+1,'] ) += ',
c     &           dscr(t2_start+r)
c      end do

c      o get the next integral and close the p<q flow-control loops
         else
c            WRITE(*,*) 'SKIPPING W(',mod(w_add-(max_r*2),max_pq),
c     &      ' of ',max_pq,') = ',abs(dscr(w_add)),' < ',w_thresh
c        end if (w_thresh.lt.abs(dscr(w_add)))
         end if
         w_add = w_add + 1
c        end do p = from_p, to_p
         end do
c        end do q = from_q, to_q
         end do

c   o close the <PQ| irrep flow-control loops
c     end if (do_irp_q)
      end if
c     end do irp_q = 1, nirrep
      end do
 
c ----------------------------------------------------------------------

c      o close the |EI> flow-control loop
c        end do e = 0, max_e-1
         end do

c   o close the E irrep flow-control loops for <PQ||EI>
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

