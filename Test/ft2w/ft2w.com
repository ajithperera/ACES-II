#ifndef _FT2W_COM_
#define _FT2W_COM_
c __FILE__ : begin

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

      integer dcoresize, freecore
      integer ndx_t1(8,2), ndx_t2(8,3), ndx_tau2(8,3)

      common /ft2w_com/
     &                 dcoresize, freecore,
     &                 ndx_t1, ndx_t2, ndx_tau2,
     &                 off_pltp, off_pq
      save   /ft2w_com/

c __FILE__ : end
#endif /* _FT2W_COM_ */
