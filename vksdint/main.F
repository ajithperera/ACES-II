        Program main

      implicit none
      integer iuhf,ihfdftgrad
     
#include"icore.com"
#include"istart.com"
         call aces_init(icore,i0,icrsiz,iuhf,.true.)
         
         call getrec(20,'JOBARC','HFDFTGRA',1,ihfdftgrad)
           if (ihfdftgrad.eq.1) then
             call vhfksdint(icore(i0),icrsiz,iuhf)
             call vksdint
           else
             call vksdint
           end if
        call aces_fin
        end             
