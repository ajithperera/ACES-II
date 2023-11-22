#include "flags.h"
      Program Fo_aces2 

      integer t0, t1
      logical Occnums_file
      character*80 Fname 


c   o gettimeofday
      call c_gtod(t0,itmp)
      dTimeStart = t0 + 1.d-6*itmp

c   o print the ACES banner

      Call aces_init_rte
      Call aces_ja_init 
      Call Getrec(20,"JOBARC",'NDROPGEO',1,Idrop)
      If (Idrop .Gt. 0) Then
         Write(6,"(a)") " Frozen core adjustments are being made"
         Call Pccd_runchg()
         Call aces_ja_fin
         Write(6,*)
      Else 
         call title
         call flushout
         call runit('rm -f FILES')
         call runit("xjoda")

c   o initialize the ACES environment and turn off the job archive subsystem
         call crapsi(ijunk,iuhf,-1)
         call aces_ja_fin
   
         Call Gfname("OCCNUMS",Fname,Ilength)
         Inquire(File=Fname(1:7),Exist=Occnums_file)

         If (Occnums_file) Then
             Write(6,"(2a)") " The occupation number file exists.",
     +                       " creating a duplicate copy OCCNUMS_DUP"   
             Write(6,"(a)")  " and deleting the original."
             call runit("mv OCCNUMS OCCNUMS_DUP")
         Endif 

         call Job_control(Occnums_file)

         call c_gtod(t1,itmp)
         print '(/)'
         print *, '@ACES2: ',
     &         'The ACES2 program has completed successfully in ',
     &         t1+1-t0,' seconds.'
         print '(/)'

      Endif 

      call c_exit(0)
      end

