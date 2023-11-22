#ifndef _XCC_COM_
#define _XCC_COM_
c __FILE__ : begin

      double precision overlap, t_sqr(3), e5qftaaa, e5qftaab
      integer dcoresize, freecore
      integer ndx_eval(8,2,2), ndx_t1(8,2), ndx_t2(8,3), ndx_t2p1(8,3),
     &        ndx_t3(8,4), ndx_s1vo(8,2)

#define h_Offsets_IltJ 1
#define h_Offsets_AltB 2
#define h_Offsets_iltj 3
#define h_Offsets_altb 4
      integer off_pltp(36,4)

#define h_Offsets_IJ 1
#define h_Offsets_Ij 2
#define h_Offsets_IA 3
#define h_Offsets_Ia 4
#define h_Offsets_iJ 5
#define h_Offsets_ij 6
#define h_Offsets_iA 7
#define h_Offsets_ia 8
#define h_Offsets_AI 9
#define h_Offsets_Ai 10
#define h_Offsets_AB 11
#define h_Offsets_Ab 12
#define h_Offsets_aI 13
#define h_Offsets_ai 14
#define h_Offsets_aB 15
#define h_Offsets_ab 16
      integer off_pq(64,16)

      common /xcc_com/
     &                 overlap, t_sqr, e5qftaaa, e5qftaab,
     &                 dcoresize, freecore,
     &                 ndx_eval, ndx_t1, ndx_t2, ndx_t2p1, ndx_t3,
     &                 ndx_s1vo,
     &                 off_pltp, off_pq
      save   /xcc_com/

      integer            CNTOT, CNID
      common /mpi_stats/ CNTOT, CNID
      save   /mpi_stats/

#include "syminf.com"
#include "sympop.com"
c __FILE__ : end
#endif /* _XCC_COM_ */
