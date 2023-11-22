      subroutine mn_align_geom(natoms, atmass, q, refq, scr, mxcor)
c
      implicit none     
c
      integer mxcor, natoms
      double precision scr(mxcor)
      integer i, j, idummy, ixyz_tran, iscr
      double precision q(3*natoms), 
     $     refq(3*natoms), atmass(natoms)
      logical print
      print = .true.
C
      call mn_align_COM(natoms, Q, atmass)
c
c align rotation as best as possible with reference geometry
c
         idummy = 1
         ixyz_tran = idummy + 3*natoms
         iscr = ixyz_tran + 3*3
c
         call mn_xyztran(natoms, Q, refQ, scr(idummy),
     $        atmass, scr(ixyz_tran), scr(iscr), mxcor-iscr+1)
c
c the new coordinates are in idummy. Copy to Q
c
         call DCOPY(3*natoms, scr(idummy), 1, Q, 1)
c
         return
         end
