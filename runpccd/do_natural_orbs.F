
      subroutine do_natural_orbs
      implicit integer (a-z)
      integer ndropa,icalc,iabcd,igamma_abcd
      logical NOS_Exist
     
#include "flags.h"
#include "flags.com"
#include "flags2.com"
      
      ndrgeo1=0
      call aces_ja_init
      call getrec(20,'JOBARC','IFLAGS  ',100,iflags)
      icalc=iflags(h_IFLAGS_calc)
      iflags(h_IFLAGS_calc)=1
      iflags(h_IFLAGS_density)=1
      iflags(h_IFLAGS_scf_maxcyc)=0
      iflags(h_IFLAGS_guess)=4
      iflags2(h_IFLAGS2_grad_calc)=1
      iabcd=iflags(h_IFLAGS_abcdtype)
      iflags(h_IFLAGS_abcdtype)=0
      ivtran=iflags(h_IFLAGS_vtran)
      iflags(h_IFLAGS_vtran)=0
      igamma_abcd=iflags(h_IFLAGS_gamma_abcd)
      iflags(h_IFLAGS_gamma_abcd)=0
      call putrec(20,'JOBARC','IFLAGS  ',100,iflags)
      call putrec(20,'JOBARC','IFLAGS2 ',500,iflags2)
      call getrec(20,'JOBARC','NUMDROPA',1,ndropa)
      call putrec(20,'JOBARC','FNOFREEZ',1,ndropa)
      call putrec(20,'JOBARC','NUMDROPA',1,0)
      if (iflags(11).ne.0) then
         call putrec(20,'JOBARC','NUMDROPB',1,0)
      end if
      if ((geom_opt.or.vib_specs).and.analytical_gradient) then
         print *,'@RUNFNO FNO gradients not yet supported.'
         call aces_exit(1)
         call putrec(20,'JOBARC','NDROPGEO',1,0)
      endif
      call aces_ja_fin
      call runit('xvtran')
      call runit('xintprc')
      call runit('xvcc')
      call runit('xdens')

      call aces_ja_init
      call getrec(20,'JOBARC','IFLAGS  ',100,iflags)
      call getrec(20,'JOBARC','IFLAGS2 ',500,iflags2)
      iflags(h_IFLAGS_calc)=icalc
      iflags(h_IFLAGS_abcdtype)=iabcd
      iflags(h_IFLAGS_vtran)=ivtran
      iflags(h_IFLAGS_gamma_abcd)=igamma_abcd
      call putrec(20,'JOBARC','IFLAGS  ',100,iflags)
      call putrec(20,'JOBARC','IFLAGS2 ',500,iflags2)
      call aces_ja_fin
      inquire(file='NATORBS',exist=NOS_Exist)
      if (NOS_Exist) Then
         call runit('xvmol')
         call runit('xvmol2ja')
         call runit('mv NATORBS OLDMOS')
         call runit('xvscf')
	 call runit('xvtran')
         call runit('xintprc')
         call runit('xvcc')
      else
         write(6,"(2a)") " The pCCD natural orbital optimization", 
     &                   "requires NATORBS file and it does not exist"
         call errex
      endif 
      return
      end

