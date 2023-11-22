      Subroutine Scf_eneg

      integer ia100(100), iguess

#include "aces2.h"
#include "flags.h"

#ifdef _TRACE_ACES2
      write(*,*) '@ACES2: Entered ', __FILE__
#endif

c This calls to  integrals do not fit in here, perhaps 
c we should create a small function call integral and 
c move all the calls to integral packages outside. 

      if (.not.direct) call runit(integral_package)
      call runit('xvmol2ja')
      If (iflags(h_IFLAGS_xfield).NE.0 .OR.
     &    iflags(h_IFLAGS_yfield).NE.0 .OR.
     &    iflags(h_IFLAGS_zfield).NE.0 ) call runit('xvprops')
C
C Extension of DKH optimizations added on 05/2012. Ajith Perera.
C
         If (Iflags2(h_IFLAGS2_dkh_order) .Gt. 0 .and. 
     &      .not. contract) Call runit('xvprops')

      if (hf_scf) then
         if (nddo_guess) call runit('xnddo')
         If (dkh .and. contract) Call Runit("xdkh")
         call runit('xvscf')

         if (iflags(74).ne.0) then

C The following can be accomplished with stability analysis
C 03/2003, Ajith Perera.
C
C 1. HF stability is ON, it will perform stability analysis
C    irrespective of the calculation level.
C
C 2. HF stablity is FOLLOW, proceed to generate appropriate
C    HF orbitls and do the user specified calculation.

c         o cache any DROPMO info
            call a2getrec(20,'JOBARC','NUMDROPA',1,ndropa)
            call a2putrec(20,'JOBARC','NUMDROPA',1,0)
            if (iflags(11).ne.0) then
            call a2getrec(20,'JOBARC','NUMDROPB',1,ndropb)
            call a2putrec(20,'JOBARC','NUMDROPB',1,0)
            end if

c         o do the instability analysis
            call runit('xvtran')
            call runit('xintprc')
            call a2getrec(-1,'JOBARC','SCFKICK',1,itest)

c         o warn the user 
            if (itest.eq.1) then
               if (iflags(74).ne.2) then
c               o there are valid instabilities, but we cannot follow
                  if (.not.plain_scf.or.geom_opt.or.
     &               iflags(18).ne.0.or.iflags(54).ne.0) then
                     print '()'
                     print *, '@ACES2: WARNING - There are valid ',
     &                        'instabilities in irrep 1, but'
                     print *, '        HFSTABILITY != FOLLOW. Results ',
     &                        'should be interpreted with caution!'
                     print '()'
                  end if
               end if
            else
C
C This block is to handle the cases where there are no instablities
C but follow is requested. If we do not rerun scf (converge in onje
C iteraltions) the data file are incorrect for subsequent calculations
C (especially gradients). 09/04, Ajith Perera.
C
               if (iflags(74).eq.2) then
                  If (dkh .and. contract) Call Runit("xdkh")
                  If (.not.plain_scf.or.geom_opt.or. iflags(18).ne.0
     &                .or.iflags(54).ne.0) call runit('xvscf')
               endif
c           end if (itest.eq.1)
            endif 

c         o restore any DROPMO info
            if (ndropa.ne.0) then
               call a2putrec(20,'JOBARC','NUMDROPA',1,ndropa)
            end if
            if (ndropb.ne.0.and.iflags(11).ne.0) then
               call a2putrec(20,'JOBARC','NUMDROPB',1,ndropb)
            end if
c
c         o re-run vscf to follow an instability OR to setup dropped MOs.
c           Also, Set GUESS=MOREAD so that scf will converge quickly.
c
            If (ndropa.ne.0.or.(itest.eq.1.and.iflags(74).eq.2)) then
            call a2getrec(20,'JOBARC','IFLAGS',100,ia100)
              iguess=ia100(45)
              ia100(45)=0
              call a2putrec(20,'JOBARC','IFLAGS',100,ia100)
              If (dkh .and. contract) Call Runit("xdkh")
              call runit('xvscf')
C         o Reset guess parameter
              call a2getrec(20,'JOBARC','IFLAGS',100,ia100)
              ia100(45)=iguess
              call a2putrec(20,'JOBARC','IFLAGS',100,ia100)
            Endif
c
c         o switch off HFSTAB for post-scf calculations
            call a2getrec(20,'JOBARC','IFLAGS',100,ia100)
c         o save the original value of hfstability for posterity. 
            call a2putrec(20,'JOABRC',"STAB_TAB",1,iflags(74))
            ia100(74)=0
            call a2putrec(20,'JOBARC','IFLAGS',100,ia100)


c        end if (iflags(74).ne.0)
         end if

      else if(ks_scf) then
         If (dkh .and. contract) Call Runit("xdkh")
         call runit('xvscf')
         call runit('xintgrt')
      else if(hf_dft) then
         If (dkh .and. contract) Call Runit("xdkh")
         call runit('xvscf')
CSSS         If (.Not. (Iflags(h_IFLAGS_calc)  .gt. 0 .and. 
CSSS     &              Iflags(h_IFLAGS_props) .eq. 1))
             call runit('xintgrt')
      else if (sl_oep) then
         call runit('xintpack')
         call runit('xintpack2vmol')
         If (dkh .and. contract) Call Runit("xdkh")
         call runit('xvscf')
      end if
c
c Also we can do TDA single point excitation energy calculations. This can be
c done via EOM-CC code or Hideo Sekino's TDHF code (there are some restrictions in
c Hideo's code so let's use it only of property calculations.
c
      If (Iflags(h_IFLAGS_excite) .EQ. 1) Then
         If (.not.(Iflags(h_IFLAGS_calc) .eq. 48 .or.
     +             Iflags(h_IFLAGS_calc) .eq. 49)) Then
            Call Change_orbital_space
            Call Runit('xvtran')
            Call Runit('xintprc')
            Call Runit('xvee')
         Else

C If rCCD or drCCD set along with cis; we could simply run the ringcc
C and ringcl, but this route is ill advised. 09/2018, Ajith Perera

            Call Change_orbital_space
            Call Runit('xvtran')
            Call Runit('xintprc')
            Call Runit('xringcc')
            Call Runit('xringcl')
         Endif 
      Endif
C
C This is to run the rpa code developed recently. Ajith Perera, 08/2016.
C
      If (Iflags(h_IFLAGS_excite) .EQ. 2) Then
         Call Change_orbital_space
         Call Runit('xvtran')
         Call Runit('xintprc')
         Call Runit('xrpa')
      Endif
C
      Return
      End

