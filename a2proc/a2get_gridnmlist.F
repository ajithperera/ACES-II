      subroutine A2get_gridnmlist(Stride,Derivatives,Grid_file,
     &                            Npoints)
C
      implicit double precision (a-h,o-z)
      logical printdef,Grid_file
      character*80 wrk
      character*12 Derivatives 
   

C Call to nl_init find the *PES_SCAN namelist if it is in the ZMAT

      printdef = .True.
      call nl_init('PES_GRID', ierr, printdef)

C Step size control the size of the individual steps.

      call nl_real('STEP_SIZE', 0.3D0, Stride)

C Numerical differentiation type (single or double) sided.

      call nl_str('DERIVATIVE', "DOUBLE", Wrk)

      call nl_log('GRID_FILE', .FALSE., Grid_file)

      call nl_int('HOW_MANY', 10, Npoints)

      call nl_term

      If (Wrk(1:6) .EQ. "SINGLE") Derivatives="Single_sided"
      If (Wrk(1:6) .EQ. "DOUBLE") Derivatives="Double_sided"

      RETURN
      END

    


      
