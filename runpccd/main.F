#include "flags.h"
      Program Runpccd 

      integer t0, t1
      character*80 Fname 


c   o gettimeofday
      call c_gtod(t0,itmp)
      dTimeStart = t0 + 1.d-6*itmp

c   o print the ACES banner

      call title
      call flushout
      call runit('rm -f FILES')
      call runit("xjoda")

c   o initialize the ACES environment and turn off the job archive subsystem
      call crapsi(ijunk,iuhf,-1)
      call aces_ja_fin

      call Job_control(Iuhf)

      call c_gtod(t1,itmp)
      print '(/)'
      print *, '@ACES2: ',
     &         'The ACES2 program has completed successfully in ',
     &         t1+1-t0,' seconds.'
      print '(/)'

      call c_exit(0)
      end

