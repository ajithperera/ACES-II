      program intpack

c Numerical integration used for calculating overlap among
c three or four basis functions. Modifications by S. Ivanov based 
c on intgrt. It is to be used to calculate elements needed for TDKS. 

c Written by Nevin Oliphant as hfdft, Quantum Theory Project,
c University of Florida, 1993.

c ReWritten by Nevin Oliphant and Sullivan Beck as potdens, 
c Quantum Theory Project, University of Florida, 1994.

c ReWritten by Nevin Oliphant and Sullivan Beck as intgrt, 
c Quantum Theory Project, University of Florida, 1995.
#include <aces.h>
#include <flags.com>
#include <flags2.com>
#include <iuhf.com>
#include <mol.com>
#include <kscore.com>
#include <symm2_ks.com>
#include <timeinfo.com>
#include <int.com>
#include <intgrtflags.com>
          integer islaterON
c INITIALIZATION
c Initializes constants used throughout the program.

c      call initpack
      int_ks = .false.
      call initintgrt(kshf,intnumradpts,islaterON)

c IVS Jan 2006: The "memknown=0,1" loop no longer works. Neither can we
c use the MEMORY ASV, because to get its value we need to go through 
c libra3/crapsi, which need (in case of memknown=-1) memory already 
c initialized. Let me set maxicore to ACES's default value
c of 15MW, and if it's not enough, well, someone is SOL.

c      maxicore = 15000000
c      memknown = -1
c      needbas  = 1
c      call crapsi(1)
c
      call timer(0) 
      call aces_init_rte
      call aces_ja_init
      call getrec(1,'JOBARC','IFLAGS  ',100,iflags)
      call getrec(1,'JOBARC','IFLAGS2 ',500,iflags2)
      iuhf=0
      if (iflags(11).ge.1) iuhf=1
      call getrec(-1,'JOBARC','COMPNIRR',1,nirrep)
      call sb_com_symm2_ks
      maxmem=iflags(36)
      maxicore=500000
      do memknown=-1,-1
#ifdef _DEBUG_INTGRT
         write(6,*) '@INTPCK: Starting main loop.'
         write(6,*) '         memknown is ', memknown
#endif /* _DEBUG_INTPCK */
         i0tmp = 1
         if (memknown.ne.0) then
            kscore(1) = 0
            do while ((kscore(1).eq.0).and.(maxmem.gt.0))
               call aces_malloc(maxmem,kscore,i0tmp)
               if (kscore(1).eq.0) maxmem = maxmem - 1024*1024
            end do
         end if
         call ks_mem_init(kscore(i0tmp),maxmem)
         call sb_com_mol

c Set up grid for numerical integration

      call setupgrd(max_angpts)
c Start here for either SCF or Relaxed correlated density integrations
      idns=0

c Set up the actual numerical grid for each symmetry unique atom and
c integrate over the symmetry unique octants of each symmetry unique
c atom

      call numint

      call sb_mem_fin

#ifdef _DEBUG_INTGRT
         call final_check
#endif /* _DEBUG_INTGRT */

c     end do memknown=0,1
      end do

      call aces_ja_fin
      call timer(1)
      write(6,9020) timetot 
 9020 format(t3,'@CHECKOUT-I, Total execution time : ',f12.4,
     &       ' seconds.')

      call callstack_term

      stop
      end
