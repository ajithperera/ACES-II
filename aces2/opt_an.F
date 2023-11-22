
      subroutine opt_an

#include "aces2.h"

      integer icore, is, ius
      logical do_derint

      ndrgeo=0
      call a2getrec(20,'JOBARC','NDROPGEO',1,NDRGEO)

c   o generate an initial approximate Hessian. 12/04, Ajith Perera
      if (mopac_guess) call runit('xmopac')

      icycle=0
1000  continue
      icycle=icycle+1
#ifdef _DEBUG
      write(*,*) '@ACES2: Beginning main opt_an loop. Iteration ',
     &           icycle
#endif

      call scf

      ndrgeo=0
      call a2getrec(20,'JOBARC','NDROPGEO',1,NDRGEO)

      if (plain_scf) then
         do_derint=.TRUE.
c        Does hessian have to be recalculated?
         if (iflags(55).ne.0) then
            ijunk=mod(icycle-1,iflags(55))
            if (ijunk.eq.0) then
               do_derint=.FALSE.
               itmp=iflags(54)
               itmp2=iflags(3)
               iflags(3)=2
               iflags(54)=1
               call a2putrec(20,'JOBARC','IFLAGS  ',100,iflags)
               call runit('xvtran')
               call runit('xintprc')
c AP - McKinley goes here
               call runit('xvdint')
               call runit('xcphf')
               iflags(3)=itmp2
               iflags(54)=itmp
               call a2putrec(20,'JOBARC','IFLAGS  ',100,iflags)
            end if
         end if
         if (do_derint) call runit(der_integral_package)
         call c_gtod(is,ius)
         print '(a,f10.1,a)', 'ACES2: Total elapsed time is ',
     &                        is+1.d-6*ius-dTimeStart,' seconds'
         call a2putrec(1,'JOBARC','DIRTYFLG',1,0)
         call runit('xjoda')
         call a2getrec(1,'JOBARC','JODADONE',1,istat)
         if (istat.ne.0) return
         call rmfiles
         goto 1000
      else
         if (ndrgeo.ne.0)  then
            call runit('cp JOBARC JOBARC_AM')
            call runit('cp JAINDX JAINDX_AM')
         end if
         call runit('xvtran')
         if (iflags2(152).eq.1) call runit('cp HF2 HF2N')
         call runit('xintprc')
         if (iflags2(152).eq.1) call runit('mv HF2N HF2')
         if (iflags(55).ne.0) then
            ijunk=mod(icycle-1,iflags(55))
            if (ijunk.eq.0) then
c NO - lines commented out are to make this work
c               istat=ishell('cp JOBARC jobarc')
c               istat=ishell('cp JAINDX jaindx')
               itmp0=iflags(2)
               itmp=iflags(54)
               itmp2=iflags(3)
               iflags(2)=0
               iflags(3)=2
               iflags(54)=1
               call a2putrec(20,'JOBARC','IFLAGS  ',100,iflags)
c AP - McKinley goes here
               call runit('xvdint')
               call runit('xcphf')
               iflags(2)=itmp0
               iflags(3)=itmp2
               iflags(54)=itmp
               call a2putrec(20,'JOBARC','IFLAGS  ',100,iflags)
c               call aces_ja_init
c               call aces_io_init(icore,1,0,.false.)
c               call aces_io_remove(51,'GAMLAM')
c               call aces_io_fin
c               call aces_ja_fin
c               istat=ishell('mv jobarc JOBARC')
c               istat=ishell('mv jaindx JAINDX')
            end if
         end if
         if (iflags(87).eq.0) then
            call runvcc
            call runit('xlambda')
         else if (iflags(87).eq.1) then
            call runit('xvee')
         else
            call runvcc
            call runit('xlambda')
            call runit('xvee')
         end if
         if (ndrgeo.ne.0) call runchg
         call runit('xdens')
         call runit('xanti')
         call runit('xbcktrn')
         call runit(der_integral_package)

c Modify (extrapolate) the gradients externally.
c 01/2006. Ajith Perera.
         if (bGExtrap.or.bCExtrap) then
            call runit('xa2proc extrap gradient')
         end if

         call c_gtod(is,ius)
         print '(a,f10.1,a)', 'ACES2: Total elapsed time is ',
     &                        is+1.d-6*ius-dTimeStart,' seconds'
         call a2putrec(1,'JOBARC','DIRTYFLG',1,0)
         call runit('xjoda')
         call a2getrec(1,'JOBARC','JODADONE',1,istat)
         if (istat.ne.0) return
         call rmfiles
         goto 1000
      end if

      return
      end

