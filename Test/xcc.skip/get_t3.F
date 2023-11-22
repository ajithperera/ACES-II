
c This routine loads the T3 amplitudes and calculates <0|adj(T3)*T3|0>.

c INPUT
c    double precision dcore(*) : address space for xcc.com

#include "xcc.h"

#define _DO_T3AAB
#define _DO_T3BBA /* requires _DO_T3AAB */
#define _DO_T3AAA
#define _DO_T3BBB /* requires _DO_T3AAA */
c#define _PRINT_T3AAA_AMPLITUDES
c#define _PRINT_T3AAB_AMPLITUDES

      subroutine get_t3(dcore)
      implicit none

#include "mpif.h"

c ARGUMENT LIST
      double precision dcore(*)

c INTERNAL VARIABLES
      integer ispin, irrep, isize
      integer t3occ(8), t3vrt(8)
      double precision e4taaa,  e4tbbb,  e4taab,  e4tbba
      double precision e5staaa, e5stbbb, e5staab, e5stbba
      double precision dtmp, dtmp2, tmp_ovrlp
      integer i
c#ifdef PRINT_TIMES
      double precision cpu_i, cpu_f, cs_d
      integer year_i, mon_i, mday_i, hour_i, min_i, sec_i
      integer year_f, mon_f, mday_f, hour_f, min_f, sec_f
      integer ch_d, cm_d, cpu, wh_d, wm_d, ws_d, wall
      integer rtime_sec, rtime_usec
c#endif /* PRINT_TIMES */
      double precision dSendBuf(8), dRecvBuf(8)

c PARAMETERS
      double precision zero
      parameter (zero=0.0d0)

c COMMON BLOCKS
#include "xcc.com"
#include "method.com"
#include "sym.com"
#include "machsp.com"
#include "flags.com"
#include "timeinfo.com" /* for timenow */

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      print *, '@XCC: Entered file ',__FILE__
#endif

c ----------------------------------------------------------------------

      t_sqr(3) = zero

c   o At this point, we don't do anything with the correlation energy.
c     Just calculate it and print it.
      e4taaa   = zero
      e5staaa  = zero
      e4tbbb   = zero
      e5stbbb  = zero
      e4taab   = zero
      e5staab  = zero
      e4tbba   = zero
      e5stbba  = zero

c NOTE: ispin for these mixed spin cases refers to the FIRST spin.
c       E.G., T3AAB would have ispin=1 and start at dcore(ndx_t3(1,3)).

#ifdef _DO_T3AAB
c   o T3(IJk,ABc); A<B,c, I<J,k
#ifdef ONLY_ONE_T3_DIST
      call irpdpd_pltpq(vrt(1,1),vrt(1,2),t3vrt)
      isize = 0
      do irrep = 1, nirrep
         ndx_t3(irrep,3) = freecore
         isize = max(isize,t3vrt(irrep))
      end do
c     freecore will recover this space right after make_t3aab...
      freecore = freecore + isize
#else
      call irpdpd_pltpq(pop(1,1),pop(1,2),t3occ)
      call irpdpd_pltpq(vrt(1,1),vrt(1,2),t3vrt)
      do irrep = 1, nirrep
         ndx_t3(irrep,3) = freecore
         freecore = freecore + ( t3occ(irrep) * t3vrt(irrep) )
      end do
      isize = freecore - ndx_t3(1,3)
#endif /* ONLY_ONE_T3_DIST */
      if (freecore.gt.1+dcoresize) then
         print *, '@GET_T3: Not enough memory.'
         call xcc_errex
      end if
      tmp_ovrlp = zero
c     dtmp is carried along far below this line if the reference is RHF
      dtmp = zero
           call wallclock(year_i,mon_i,mday_i,hour_i,min_i,sec_i)
           print *, 'make_t3aab in  : ',
     &                    year_i,mon_i,mday_i,hour_i,min_i,sec_i
           call c_gtod(rtime_sec,rtime_usec)
           cpu_i = (rtime_usec*1.d-6) + rtime_sec
      call make_t3aab(dcore,1,e4taab,e5staab,tmp_ovrlp,
     &                        e4taaa,e5staaa,dtmp)
           call c_gtod(rtime_sec,rtime_usec)
           cpu_f = (rtime_usec*1.d-6) + rtime_sec
           call wallclock(year_f,mon_f,mday_f,hour_f,min_f,sec_f)
           wall = cpu_f-cpu_i
           print *, 'make_t3aab out : ',
     &                    year_f,mon_f,mday_f,hour_f,min_f,sec_f,
     &              ' (',wall,')'
#ifdef PRINT_TIMES
c      print '(a,a,i4,a,i2,a,f5.2,a,i8,a,a,i4,a,i2,a,i2,a,i8,a)',
c     &      'T3AAB TIMING:  ',
c     &      'cpu ', ch_d,':',cm_d,':',cs_d,' (', cpu,') ;  ',
c     &      'real ',wh_d,':',wm_d,':',ws_d,' (',wall,')'
#endif /* PRINT_TIMES */
#ifdef ONLY_ONE_T3_DIST
      freecore = ndx_t3(1,3)
#endif /* ONLY_ONE_T3_DIST */
      t_sqr(3) = t_sqr(3) + tmp_ovrlp
#ifdef _DEBUG_XCC
      print *, '<0|adj(T3AAB)*T3AAB|0> = ', tmp_ovrlp
#endif /* _DEBUG_XCC */

#ifdef _DO_T3BBA
c   o T3(ijK,abC); a<b,C, i<j,K
      if (iuhf.eq.0) then
         do irrep = 1, nirrep
            ndx_t3(irrep,4) = ndx_t3(irrep,3)
         end do
         t_sqr(3) = t_sqr(3) + tmp_ovrlp
         e4tbba  = e4taab
         e5stbba = e5staab
      else
#ifdef ONLY_ONE_T3_DIST
         call irpdpd_pltpq(vrt(1,2),vrt(1,1),t3vrt)
         isize = 0
         do irrep = 1, nirrep
            ndx_t3(irrep,4) = freecore
            isize = max(isize,t3vrt(irrep))
         end do
c        freecore will recover this space right after make_t3aab...
         freecore = freecore + isize
#else
         call irpdpd_pltpq(pop(1,2),pop(1,1),t3occ)
         call irpdpd_pltpq(vrt(1,2),vrt(1,1),t3vrt)
         do irrep = 1, nirrep
            ndx_t3(irrep,4) = freecore
            freecore = freecore + ( t3occ(irrep) * t3vrt(irrep) )
         end do
         isize = freecore - ndx_t3(1,4)
#endif /* ONLY_ONE_T3_DIST */
         if (freecore.gt.1+dcoresize) then
            print *, '@GET_T3: Not enough memory.'
            call xcc_errex
         end if
         tmp_ovrlp = zero
              call wallclock(year_i,mon_i,mday_i,hour_i,min_i,sec_i)
              print *, 'make_t3bba in  : ',
     &                       year_i,mon_i,mday_i,hour_i,min_i,sec_i
              call c_gtod(rtime_sec,rtime_usec)
              cpu_i = (rtime_usec*1.d-6) + rtime_sec
         call make_t3aab(dcore,2,e4tbba,e5stbba,tmp_ovrlp,
     &                           e4tbbb,e5stbbb,dtmp)
              call c_gtod(rtime_sec,rtime_usec)
              cpu_f = (rtime_usec*1.d-6) + rtime_sec
              call wallclock(year_f,mon_f,mday_f,hour_f,min_f,sec_f)
              wall = cpu_f-cpu_i
              print *, 'make_t3bba out : ',
     &                       year_f,mon_f,mday_f,hour_f,min_f,sec_f,
     &                 ' (',wall,')'
#ifdef PRINT_TIMES
c      print '(a,a,i4,a,i2,a,f5.2,a,i8,a,a,i4,a,i2,a,i2,a,i8,a)',
c     &      'T3BBA TIMING:  ',
c     &      'cpu ', ch_d,':',cm_d,':',cs_d,' (', cpu,') ;  ',
c     &      'real ',wh_d,':',wm_d,':',ws_d,' (',wall,')'
#endif /* PRINT_TIMES */
#ifdef ONLY_ONE_T3_DIST
         freecore = ndx_t3(1,4)
#endif /* ONLY_ONE_T3_DIST */
         t_sqr(3) = t_sqr(3) + tmp_ovrlp
      end if
#ifdef _DEBUG_XCC
      print *, '<0|adj(T3BBA)*T3BBA|0> = ', tmp_ovrlp
#endif /* _DEBUG_XCC */
#endif /* _DO_T3BBA */
#endif /* _DO_T3AAB */

#ifdef _DO_T3AAA
c   o T3(IJK,ABC); A<B<C, I<J<K
      if (iuhf.eq.0) then
c         call irpdpd_pltpltp(pop(1,1),t3occ)
c         call irpdpd_pltpltp(vrt(1,1),t3vrt)
c         do irrep = 1, nirrep
c            ndx_t3(irrep,1) = freecore
c            freecore = freecore + ( t3occ(irrep) * t3vrt(irrep) )
c         end do
c         do irrep = 1, nirrep
c            ndx_t3(irrep,2) = ndx_t3(irrep,1)
c         end do
c         isize = freecore - ndx_t3(1,1)
c         call make_t3aaa_rhf(dcore,e4taaa)
         e4tbbb  = e4taaa
         e5stbbb = e5staaa
         t_sqr(3) = t_sqr(3) + ( dtmp * 2.0d0 )
#ifdef _DEBUG_XCC
         print *, '<0|adj(T3AAA)*T3AAA|0> = ', dtmp
         print *, '<0|adj(T3BBB)*T3BBB|0> = ', dtmp
#endif /* _DEBUG_XCC */
      else
#ifdef _DO_T3BBB
         do ispin = 1, 2
#else
         do ispin = 1, 1
#endif /* _DO_T3BBB */
#ifdef ONLY_ONE_T3_DIST
            call irpdpd_pltpltp(vrt(1,ispin),t3vrt)
            isize = 0
            do irrep = 1, nirrep
               ndx_t3(irrep,ispin) = freecore
               isize = max(isize,t3vrt(irrep))
            end do
c           freecore will recover this space right after make_t3aaa_uhf...
            freecore = freecore + isize
#else
            call irpdpd_pltpltp(pop(1,ispin),t3occ)
            call irpdpd_pltpltp(vrt(1,ispin),t3vrt)
            do irrep = 1, nirrep
               ndx_t3(irrep,ispin) = freecore
               freecore = freecore + ( t3occ(irrep) * t3vrt(irrep) )
            end do
            isize = freecore - ndx_t3(1,ispin)
#endif /* ONLY_ONE_T3_DIST */
            if (freecore.gt.1+dcoresize) then
               print *, '@GET_T3: Not enough memory.'
               call xcc_errex
            end if
            tmp_ovrlp = zero
            dtmp      = zero
            dtmp2     = zero
                 call wallclock(year_i,mon_i,mday_i,hour_i,min_i,sec_i)
                 print *, 'make_t3aaa in  : ',
     &                          year_i,mon_i,mday_i,hour_i,min_i,sec_i
                 call c_gtod(rtime_sec,rtime_usec)
                 cpu_i = (rtime_usec*1.d-6) + rtime_sec
            call make_t3aaa_uhf(dcore,ispin,dtmp,dtmp2,tmp_ovrlp)
                 call c_gtod(rtime_sec,rtime_usec)
                 cpu_f = (rtime_usec*1.d-6) + rtime_sec
                 call wallclock(year_f,mon_f,mday_f,hour_f,min_f,sec_f)
                 wall = cpu_f-cpu_i
                 print *, 'make_t3aaa out : ',
     &                          year_f,mon_f,mday_f,hour_f,min_f,sec_f,
     &                    ' (',wall,')'
#ifdef PRINT_TIMES
c      if (ispin.eq.1) then
c      print '(a,a,i4,a,i2,a,f5.2,a,i8,a,a,i4,a,i2,a,i2,a,i8,a)',
c     &      'T3AAA TIMING:  ',
c     &      'cpu ', ch_d,':',cm_d,':',cs_d,' (', cpu,') ;  ',
c     &      'real ',wh_d,':',wm_d,':',ws_d,' (',wall,')'
c      else
c      print '(a,a,i4,a,i2,a,f5.2,a,i8,a,a,i4,a,i2,a,i2,a,i8,a)',
c     &      'T3BBB TIMING:  ',
c     &      'cpu ', ch_d,':',cm_d,':',cs_d,' (', cpu,') ;  ',
c     &      'real ',wh_d,':',wm_d,':',ws_d,' (',wall,')'
c      end if
#endif /* PRINT_TIMES */
#ifdef ONLY_ONE_T3_DIST
            freecore = ndx_t3(1,ispin)
#endif /* ONLY_ONE_T3_DIST */
            if (ispin.eq.1) then
               e4taaa  = dtmp
               e5staaa = dtmp2
            else
               e4tbbb  = dtmp
               e5stbbb = dtmp2
            end if
            t_sqr(3) = t_sqr(3) + tmp_ovrlp
#ifdef _DEBUG_XCC
            if (ispin.eq.1) then
               print *, '<0|adj(T3AAA)*T3AAA|0> = ', tmp_ovrlp
            else
               print *, '<0|adj(T3BBB)*T3BBB|0> = ', tmp_ovrlp
            end if
#endif /* _DEBUG_XCC */
c        end do ispin = 1, ?
         end do
      end if
#endif /* _DO_T3AAA */

      overlap = overlap + t_sqr(3)
      print *, '<0|adj(T3)*T3|0> = ', t_sqr(3)

#ifndef ONLY_ONE_T3_DIST
#ifdef _PRINT_T3AAA_AMPLITUDES
      call irpdpd_pltpltp(pop(1,1),t3occ)
      call irpdpd_pltpltp(vrt(1,1),t3vrt)
      do irrep = 1, nirrep
         isize = t3occ(irrep) * t3vrt(irrep)
         print *, 'There are ',isize,' T3AAA amplitudes in irrep ',
     &            irrep,'.'
         do i = 0, isize-1
            print '(i5,a,i1,a,f16.10)',
     &         1+i,'[',irrep,'] : ',dcore(ndx_t3(irrep,1)+i)
         end do
      end do
      if (iuhf.eq.1) then
         call irpdpd_pltpltp(pop(1,2),t3occ)
         call irpdpd_pltpltp(vrt(1,2),t3vrt)
         do irrep = 1, nirrep
            isize = t3occ(irrep) * t3vrt(irrep)
            print *, 'There are ',isize,' T3BBB amplitudes in irrep ',
     &               irrep,'.'
            do i = 0, isize-1
               print '(i5,a,i1,a,f16.10)',
     &            1+i,'[',irrep,'] : ',dcore(ndx_t3(irrep,2)+i)
            end do
         end do
      end if
#endif /* _PRINT_T3AAA_AMPLITUDES */
#ifdef _PRINT_T3AAB_AMPLITUDES
      call irpdpd_pltpq(pop(1,1),pop(1,2),t3occ)
      call irpdpd_pltpq(vrt(1,1),vrt(1,2),t3vrt)
      do irrep = 1, nirrep
         isize = t3occ(irrep) * t3vrt(irrep)
         print *, 'There are ',isize,' T3AAB amplitudes in irrep ',
     &            irrep,'.'
         do i = 0, isize-1
            print '(i5,a,i1,a,f16.10)',
     &         1+i,'[',irrep,'] : ',dcore(ndx_t3(irrep,3)+i)
         end do
      end do
      if (iuhf.eq.1) then
         call irpdpd_pltpq(pop(1,2),pop(1,1),t3occ)
         call irpdpd_pltpq(vrt(1,2),vrt(1,1),t3vrt)
         do irrep = 1, nirrep
            isize = t3occ(irrep) * t3vrt(irrep)
            print *, 'There are ',isize,' T3BBA amplitudes in irrep ',
     &               irrep,'.'
            do i = 0, isize-1
               print '(i5,a,i1,a,f16.10)',
     &            1+i,'[',irrep,'] : ',dcore(ndx_t3(irrep,4)+i)
            end do
         end do
      end if
#endif /* _PRINT_T3AAB_AMPLITUDES */
#endif /* ONLY_ONE_T3_DIST */

c ----------------------------------------------------------------------

c   o prep the receive buffer (might be unnecessary)
      do i = 1, 8
         dRecvBuf(i) = 0.d0
      end do

c   o prep the send buffer
      dSendBuf(1) = e4taaa
      dSendBuf(2) = e4tbbb
      dSendBuf(3) = e4taab
      dSendBuf(4) = e4tbba
      dSendBuf(5) = e5staaa
      dSendBuf(6) = e5stbbb
      dSendBuf(7) = e5staab
      dSendBuf(8) = e5stbba

c   o reduce the data
      call MPI_REDUCE(dSendBuf,dRecvBuf,8,
     &                MPI_DOUBLE_PRECISION,MPI_SUM,
     &                0,MPI_COMM_WORLD,i)
      if (i.ne.MPI_SUCCESS) call mpi_errex('MPI_REDUCE',i)

c   o get out of here early
      if (CNID.ne.0) return

c   o extract the data
      e4taaa  = dRecvBuf(1)
      e4tbbb  = dRecvBuf(2)
      e4taab  = dRecvBuf(3)
      e4tbba  = dRecvBuf(4)
      e5staaa = dRecvBuf(5)
      e5stbbb = dRecvBuf(6)
      e5staab = dRecvBuf(7)
      e5stbba = dRecvBuf(8)

      print *, 'E[4]TAAA  = ', e4taaa
      print *, 'E[4]TBBB  = ', e4tbbb
      print *, 'E[4]TAAB  = ', e4taab
      print *, 'E[4]TBBA  = ', e4tbba
      print *, 'E[5]STAAA = ', e5staaa
      print *, 'E[5]STBBB = ', e5stbbb
      print *, 'E[5]STAAB = ', e5staab
      print *, 'E[5]STBBA = ', e5stbba
c      if (iuhf.eq.0) e5qftaaa = e5qftaaa * 2.0d0
c      print *, 'E[5]QfTAAA = ', e5qftaaa
c      if (iuhf.eq.0) e5qftaab = e5qftaab * 2.0d0
c      print *, 'E[5]QfTAAB = ', e5qftaab

c   o This is only legit with CCSD or CCSD(T) calculations.
      call getrec(1,'JOBARC','TOTENERG',iintfp,dtmp)

      if (iflags(2).eq.10) then
      print '(a,f20.12)', 'E CCSD         = ', dtmp
      dtmp = dtmp + (e4taaa+e4tbbb+e4taab+e4tbba)
      print '(a,f20.12)', 'E CCSD+T(CCSD) = ', dtmp
      dtmp = dtmp + (e5staaa+e5stbbb+e5staab+e5stbba)
      print '(a,f20.12)', 'E CCSD(T)      = ', dtmp
      dtmp = dtmp + 0.5d0*(e5qftaaa+e5qftaab)
c      print '(a,f20.12)', 'E CCSD(TQf)    = ', dtmp

      else
      if (iflags(2).eq.22) then
      print '(a,f20.12)', 'E CCSD(T)      = ', dtmp
c      print '(a,f20.12)', 'E CCSD(TQf)    = ',
c     &       dtmp + 0.5d0*(e5qftaaa+e5qftaab)
      dtmp = dtmp - (e5staaa+e5stbbb+e5staab+e5stbba)
      print '(a,f20.12)', 'E CCSD+T(CCSD) = ', dtmp
      dtmp = dtmp - (e4taaa+e4tbbb+e4taab+e4tbba)
      print '(a,f20.12)', 'E CCSD         = ', dtmp

      else
      dtmp =        (e4taaa+e4tbbb+e4taab+e4tbba)
      print '(a,f20.12)', 'deltaE T(CCSD) = ', dtmp
      dtmp = dtmp + (e5staaa+e5stbbb+e5staab+e5stbba)
      print '(a,f20.12)', 'deltaE (T)     = ', dtmp
      dtmp = dtmp + 0.5d0*(e5qftaaa+e5qftaab)
      print '(a,f20.12)', 'deltaE (TQf)   = ', dtmp

c     end if (iflags(2).eq.22)
      end if
c     end if (iflags(2).eq.10)
      end if

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      print *, '@XCC: Leaving file ',__FILE__
#endif

      return
      end

