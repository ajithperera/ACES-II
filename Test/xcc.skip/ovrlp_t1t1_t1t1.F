
c This routine adds the quantity <0|adj(T1)^2*T1^2|0>/4 to the overlap.

c INPUT
c    double precision dcore(*) : address space for xcc.com
c    double precision dscr(1)  : first free scratch address

c#define _DEBUG_OVRLP_T1T1_T1T1

      subroutine ovrlp_t1t1_t1t1(dcore,dscr)
      implicit none

c ARGUMENT LIST
      double precision dcore(*), dscr(*)

c INTERNAL VARIABLES
      integer nocc, nvrt, to_occ, to_vrt, num_ab
      integer ispin, irrep, j, i, a, b, ndx, t1_start
      double precision ovrlp_tmp, dtmp

c COMMON BLOCKS
#include "xcc.com"
#include "method.com"
#include "sym.com"

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      print *, '@XCC: Entered file ',__FILE__
#endif

c ----------------------------------------------------------------------

c The plan is:
c   o create the quantity
c            (  A B A B )   1     (  A B A B )
c        SUM ( t t t t  ) + - SUM ( t t t t  )
c        I<J (  I I J J )   2  I  (  I I I I )
c        A,B                  A,B
c   o subtract it from 0.5*<0|adj(T1)*T1|0>^2

      ovrlp_tmp = 0.0d0

c   o loop over spin
      do ispin = 1, 1+iuhf

c      o loop over the non-empty T1 irreps
         do irrep = 1, nirrep
            nocc  = pop(irrep,ispin)
            nvrt  = vrt(irrep,ispin)
         if ((nocc*nvrt).ne.0) then

c         o maximum values
            to_occ = nocc-1
            to_vrt = nvrt-1
            num_ab = nvrt*nvrt

c ----------------------------------------------------------------------

c [ T(I,A)*T(I,B)*T(I,B)*T(I,A) ] / 2

      do i = 0, to_occ

c      o load T1(I,A)*T1(I,B) into the scratch vector
         t1_start = ndx_t1(irrep,ispin) + nvrt*i
         ndx = 1
         do b = 0, to_vrt
            dtmp = dcore(t1_start+b)
            do a = 0, to_vrt
               dscr(ndx) = dcore(t1_start+a) * dtmp
               ndx = ndx + 1
            end do
         end do

c      o add [ T1(I,A)*T1(I,B)/4 ] ^ 2 to the temp overlap
         dtmp = 0.0d0
         do ndx = 1, num_ab
            dtmp = dtmp + dscr(ndx)*dscr(ndx)
         end do
         ovrlp_tmp = ovrlp_tmp + ( dtmp * 0.50d0 )

c     end do i = 0, to_occ
      end do

c ----------------------------------------------------------------------

c [ T(I,A)*T(I,B)*T(J,A)*T(J,B) ] ; I<J

c   o do work only if the appropriate T1 amplitudes exist
      if (nocc.gt.1) then

         do j = 1, to_occ

c         o load T1(J,A)*T1(J,B) into the scratch vector
            t1_start = ndx_t1(irrep,ispin) + nvrt*j
            ndx = 1
            do b = 0, to_vrt
               dtmp = dcore(t1_start+b)
               do a = 0, to_vrt
                  dscr(ndx) = dcore(t1_start+a) * dtmp
                  ndx = ndx + 1
               end do
            end do

            do i = 0, j-1

c            o load T1(I,A)*T1(I,B) into the scratch vector
               t1_start = ndx_t1(irrep,ispin) + nvrt*i
               ndx = num_ab + 1
               do b = 0, to_vrt
                  dtmp = dcore(t1_start+b)
                  do a = 0, to_vrt
                     dscr(ndx) = dcore(t1_start+a) * dtmp
                     ndx = ndx + 1
                  end do
               end do

c            o add the dot product to the temp overlap
               do ndx = 1, num_ab
                  ovrlp_tmp =   ovrlp_tmp
     &                        + ( dscr(ndx) * dscr(num_ab+ndx) )
               end do

c           end do i = 0, j-1
            end do

c        end do j = 1, to_occ
         end do

c     end if (nocc.gt.1)
      end if

c ----------------------------------------------------------------------

c      o close the loop over T1 irreps
c        end if ((nocc*nvrt).ne.0)
         end if
c        end do irrep = 1, nirrep
         end do

         if (iuhf.eq.0) ovrlp_tmp = ovrlp_tmp * 2.0d0

c   o close the loop over spin
c     end do ispin = 1, (1+iuhf)
      end do

c   o The first part of this term was already calculated (sorta). Simply
c     take half of <0|adj(T1)*T1|0>^2.
c        NOTE: Since this includes the contribution from (i=j,a=b),
c     the (i=j,a=b) cross terms had to be included in the steps above
c     to correct the overlap.
      dtmp = 0.50d0 * t_sqr(1) * t_sqr(1)
      ovrlp_tmp = -ovrlp_tmp

#ifdef _DEBUG_OVRLP_T1T1_T1T1
      print *, 'The T1^2*T1^2 cross term is ~',ovrlp_tmp,' a.u.'
      print *, 'The T1^2*T1^2 disconnected term is ~',dtmp,' a.u.'
#endif

      ovrlp_tmp = ovrlp_tmp + dtmp

      overlap = overlap + ovrlp_tmp
      print *, '<0|adj(T1)^2*T1^2|0>/4 = ',ovrlp_tmp

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      print *, '@XCC: Leaving file ',__FILE__
#endif

      return
      end

