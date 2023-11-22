
#ifndef _INTGRTFLAGS_COM_
#define _INTGRTFLAGS_COM_

c This contains flags that are set in the INTPACK namelist.  See the
c file initpack.F for a description of each of them.
c
      integer int_numradpts,int_radtyp,int_partpoly,int_radscal,
     &    int_parttyp,int_fuzzyiter,int_defanglgrid,int_defgridtype,
     &    int_overlp

      M_REAL
     &    int_radlimit

      common /intgrtflags/  int_numradpts,int_radtyp,int_partpoly,
     &    int_radscal,int_parttyp,int_fuzzyiter,int_defanglgrid,
     &    int_defgridtype,int_overlp

      common /intgrtflagsd/ int_radlimit

      save /intgrtflags/
c      save /intgrtflagsl/
      save /intgrtflagsd/

c The following are parameters used in the namelist

      integer int_radtyp_handy,int_radtyp_gl
      parameter (int_radtyp_handy=1)
      parameter (int_radtyp_gl   =2)

      integer int_partpoly_equal,int_partpoly_bsrad,
     &    int_partpoly_dynamic
      parameter (int_partpoly_equal  =1)
      parameter (int_partpoly_bsrad  =2)
      parameter (int_partpoly_dynamic=3)

      integer int_radscal_none,int_radscal_slater
      parameter (int_radscal_none  =1)
      parameter (int_radscal_slater=2)

      integer int_parttyp_rigid,int_parttyp_fuzzy
      parameter (int_parttyp_rigid=1)
      parameter (int_parttyp_fuzzy=2)

      integer int_gridtype_leb
      parameter (int_gridtype_leb=1)

#endif /* _INTGRTFLAGS_COM_ */

