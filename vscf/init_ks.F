      Subroutine init_ks
C
C Initilize the /SYMM_KS/ common block and allocate memory 
C to be used exclusively for DFT related work. The amonut
C that is being allocated is read from the user input (in
C ks namelist) and is independent from the MEM_SIZE specified
C in ACESII namelist.
C
      
      Implicit none
C
#include <aces.h>
#include <flags.com>
#include <flags2.com>
#include <iuhf.com>
#include <mol.com>
#include <kscore.com>
#include <symm2_ks.com>
#include <timeinfo.com>
c#include <../intgrt/include/int.com>

#include<intgrtflags.com>

      Integer ks_maxmem, kstmp 
C
      call callstack_push('init_ks')
C
      call timer(0)
C
C No longer needed; in fact can be problematic, 11/07, Ajith Perera
C
CSSS      call getrec(1,'JOBARC','IFLAGS  ',100,iflags)
CSSS      call getrec(1,'JOBARC','IFLAGS2 ',500,iflags2)
C
      iuhf=0
      if (iflags(11).ge.1) iuhf=1
      call getrec(-1,'JOBARC','COMPNIRR',1,nirrep)
      call sb_com_symm2_ks
C
C Prkash: set up KS-namelist input for KS_maxmem

C---------------------------------------------
      ks_maxmem = int_ksmem
      maxicore  =500000

      do memknown=-1,-1
#ifdef _DEBUG_INTGRT
         write(*,*) '@INTGRT: Starting main loop.'
         write(*,*) '         memknown is ', memknown
#endif /* _DEBUG_INTGRT */
         kstmp = 1
         if (memknown.ne.0) then
            kscore(1) =  0
            do while ((kscore(1).eq.0).and.(KS_maxmem.gt.0))
               call aces_malloc(ks_maxmem,kscore,kstmp)
               if (kscore(1).eq.0) ks_maxmem = ks_maxmem - 1024*1024
            end do
         end if 
         call ks_mem_init(kscore(kstmp),ks_maxmem)
#ifdef _DEBUG
         Print*, "Entering sb_com_mol"
#endif
         call sb_com_mol
C
      end do
C   
      
      call callstack_pop
C
      return
      end
