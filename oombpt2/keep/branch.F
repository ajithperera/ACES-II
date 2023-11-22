
c     Modifications for "Release 2" :
c
c     NO  3/29/94  HF-DFT single points
c     JDW 3/31/94  5th order energy corrections, tdhf, EOM-EA single points
c     JDW 9/25/95  KKB's' dropped core derivatives
c     SG  3/4/96   P-EOM-MBPT(2) gradients
c     AP  12/07/99 KS-DFT
c     AP  07/2000  seward, alaska, and mckinley

#include "flags.h"

      subroutine branch

#include "aces2.h"
#include "parallel_aces.com"

      integer icore, is, ius
      character*79 szHyper, szMolden
      logical bHyper, bMolden, blCCSDT, bCCSD2T
#if !defined(NO_EXTERNAL) && defined(_GFORTRAN)
      call dummyproc(aces_bd_parallel_aces)
#endif
      szHyper = 'xa2proc hyperchem'
      szMolden = 'xa2proc molden'
      bHyper  = (iflags2(h_IFLAGS2_external).eq.1)
      bMolden = (iflags2(h_IFLAGS2_external).eq.2)
      blCCSDT = (iflags(h_IFLAGS_calc).eq.42)
      bCCSD2T = (iflags(h_IFLAGS_calc).eq.51)
      bGExtrap = (iflags2(h_IFLAGS2_extrap).eq.1)
      bEExtrap = (iflags2(h_IFLAGS2_extrap).eq.2)
      bCExtrap = (iflags2(h_IFLAGS2_extrap).eq.3)

      if (geomopt) then
         if (analytical_gradient) then
            call opt_an
         else
            call opt_num
         end if
c      o generate HyperChem input file after optimizing the geometry
         if (bHyper) call runit(szHyper)
         return
      end if

      call a2getrec(20,'JOBARC','NDROPGEO',1,NDRGEO)

      if (plain_scf) then

         if (iflags(h_IFLAGS_vib).eq.0) then

            call scf

            if (iflags(h_IFLAGS_props).eq.1) then
               call runit('xprops')
            else if (iflags(h_IFLAGS_props).eq.2) then
               call runit('xvtran')
               call runit('xintprc')
               call runit('xvdint')
               call runit('xcphf')
            else if (iflags(h_IFLAGS_props).eq.3) then
               call runit('xvtran')
               call runit('xintprc')
               call runit('xvdint')
               call runit('xcphf')
               call runit('xjoda')
            end if

c        else if (VIB != 0)
         else

            if (iflags(h_IFLAGS_vib).eq.1) then

               call scf
               call runit('xvtran')
               call runit('xintprc')
               call runit('xvdint')
               call runit('xcphf')
               call runit('xjoda')

            else if (iflags(h_IFLAGS_vib).eq.3) then

c This shell script (with minor mods) can perform a restartable SCF vib freq
c rmjunk; restore; xjoda || exit 1
c analytical_gradient=$(some_test)
c while true
c do xvmol; xvmol2ja; xvscf
c    test analytical_gradient && xvdint
c    xa2proc clrdirty
c    backup; restore; xjoda || exit 1
c    test $(xa2proc jareq FNDFDONE i 1 | tail -1) -ne 0 && break
c done
c clean; rmbackup
c exit

               istat=0
               do while (istat.eq.0)
                  call scf
                  if (raman) then
c                  o JDW & AP 09/98 SCF Raman intensites
                     call runit('xvtran')
                     call runit('xintprc')
                     call runit(der_integral_package)
                     call runit('xcphf')
                  else
                     if (analytical_gradient) then
                        call runit(der_integral_package)
                     end if
                  end if
                  call a2getrec(1,'JOBARC','LASTGEOM',1,istat)
                  if (istat.eq.0) then
                     call c_gtod(is,ius)
                     print '(a,f10.1,a)',
     &                     'ACES2: Total elapsed time is ',
     &                     is+1.d-6*ius-dTimeStart,' seconds'
                     call a2putrec(1,'JOBARC','DIRTYFLG',1,0)
                     call runit(xjoda)
                  else
#ifdef _HAVE_MPI
                     if (nprocs.ne.1) call sync_fd
#endif
                     if (irank.eq.0) then
                        call a2putrec(1,'JOBARC','DIRTYFLG',1,0)
                        call runit(xjoda)
                     end if
                  end if
                  call rmfiles
               end do

            end if

c        end if (VIB == 0)
         end if

c ----------------------------------------------------------------------

c     else if (.not.plain_scf) then
      else

         if (iflags(h_IFLAGS_vib).eq.0) then

            call scf
            if (dirmp2) then
               call runit('xdirmp2')
            else
               call precc(ndrgeo)
            end if

            if (iflags(h_IFLAGS_props).eq.0) then

c            o single point energies
               if (iflags(h_IFLAGS_excite).eq.1) then
                  call runit('xvee')
               else
                  if (blCCSDT) then
                     iflags(h_IFLAGS_calc)=10
                     iflags(h_IFLAGS_deriv_lev)=1
                     call a2putrec(20,'JOBARC','IFLAGS  ',100,iflags)
                  end if
                  if (bCCSD2T) then
                     iflags(h_IFLAGS_calc)=10
                     iflags(h_IFLAGS_hbar)=1
                     iflags(h_IFLAGS_estate_prop)=1
                     call a2putrec(20,'JOBARC','IFLAGS  ',100,iflags)
                  endif
                  if ((iflags(h_IFLAGS_calc) .eq. 48) .or.
     &                (iflags(h_IFLAGS_calc) .eq. 49)) then
                    call runoombpt
                  else if (.not.dirmp2) then
                    call runvcc
                  endif
               end if
               if (iflags(h_IFLAGS_hbar).eq.1)
     &            call runit('xlambda')
               if (iflags(h_IFLAGS_excite).gt.1) call runit('xvee')
               if (iflags2(h_IFLAGS2_ea_calc).eq.5) call runit('xvea')
               if (blCCSDT) then
                  call runit('xlambda')
                  iflags(h_IFLAGS_calc)=42
                  iflags(h_IFLAGS_deriv_lev)=0
                  call a2putrec(20,'JOBARC','IFLAGS  ',100,iflags)
                  call runit('xlcct')
               end if
               if (bCCSD2T) then
                  iflags(h_IFLAGS_calc)=51
                  iflags(h_IFLAGS_hbar)=0
                  iflags(h_IFLAGS_estate_prop)=0
                  call a2putrec(20,'JOBARC','IFLAGS  ',100,iflags)
                  call runit('xlcct')
               endif

c           else if (VIB == 0 && PROPS != 0) then
            else

c            o single point properties
               if (iflags(h_IFLAGS_excite).eq.1) then
                  call runit('xvee')
               else
                  if (blCCSDT) then
                     iflags(h_IFLAGS_calc)=10
                     call a2putrec(20,'JOBARC','IFLAGS  ',100,iflags)
                  endif
                  call runvcc
                  call runit('xlambda')
                  if (blCCSDT) then
                     iflags(h_IFLAGS_calc)=42
                     call a2putrec(20,'JOBARC','IFLAGS  ',100,iflags)
                     call runit('xlcct')
                  endif
                  if (iflags(h_IFLAGS_excite).gt.1) call runit('xvee')
               end if

               if ((iflags(h_IFLAGS_props).eq.8 ).or.
     &             (iflags(h_IFLAGS_props).eq.9 ).or.
     &             (iflags(h_IFLAGS_props).eq.10).or.
     &             (iflags(h_IFLAGS_props).eq.11).or.
     &             (iflags(h_IFLAGS_props).eq.13).and.(.not.raman)) then
c               o AP 10/97 EOM-based second order properties
c                          polarizabilities, hyperpolarizabilities,
c                          NMR spin-spin coupling constants
                  call runit('xvcceh')
                  return
               else if (iflags(h_IFLAGS_props).gt.100) then
c               o AP 10/97 At the moment, I have no idea what this is for.
c                 YAU 5/04 This looks like an OLD hack that shifted PROPS by 100
                  call runit('xvdint')
                  call runit('xvcceh')
                  call runit('xjoda')
                  return
               end if

               if (ndrgeo.ne.0) call runchg
               call runit('xdens')
               if (iflags(h_IFLAGS_props).eq.1) then
                  call runit('xprops')
               else
                  call runit('xvdint')
                  call runit('xcphf')
                  if (iflags(h_IFLAGS_props).eq.2) then
                     print *, '@BRANCH: ',
     &                        'exact 2nd order CC properties are broken'
                     return
                     call runit('xsquint')
                     call runit('xdtran')
                  else if (iflags(h_IFLAGS_props).eq.3) then
                     call runit('xnmr')
                     if (iflags(h_IFLAGS_treat_pert).ne.0) then
c                     o sequential treatment of B-field components requires
c                       reevaluation of GIAO integrals
                        call runit('rm -f I*X')
                        call runit('xvdint')
                        call runit('xnmr')
                        call runit('rm -f I*Y')
                        call runit('xvdint')
                        call runit('xnmr')
                        call runit('rm -f I*Z')
                     end if
                     call runit('xjoda')
                  end if
               end if

c           end if (VIB == 0 && PROPS==0)
            end if

c        else if (VIB != 0) then
         else

            if (iflags(h_IFLAGS_vib).eq.1) then

               print *, '@BRANCH: exact post-SCF frequencies are broken'
               return
               call scf
               call precc(ndrgeo)
               call runvcc
               if (ndrgeo.ne.0) call runchg
               call runit('xdens')
c              AP - McKinley goes here
               call runit('xvdint')
               call runit('xcphf')
               call runit('xsquint')
               call runit('xdtran')
               call runit('xbcktrn')
               call runit('xvdint')
               call runit('xjoda')

            else if (iflags(h_IFLAGS_vib).eq.3) then

               istat=0
               do while (istat.eq.0)
                  call scf
                  if (dirmp2) then
                     call runit('xdirmp2')
                  else
                     call precc(ndrgeo)
                     if (blCCSDT) then
   	                   iflags(2)=10
                       call a2putrec(20,'JOBARC','IFLAGS  ',100,iflags)
                     end if
                     call runvcc
                  end if
                  if (analytical_gradient) then
                     if (raman) then
c                     o JDW & AP 8/98 EOM-CCSD/MBPT Raman Intensities
                        if (iflags(h_IFLAGS_props).eq.11)  then
                           call runit('xlambda')
                           call runit('xvcceh')
                           call runit('xdens')
                           call runit('xanti')
                           call runit('xbcktrn')
                           call runit(der_integral_package)
                           call aces_ja_init
                           call aces_io_init(icore,1,0,.false.)
                           call aces_io_remove(54,'DERGAM')
                           call aces_io_fin
                           call aces_ja_fin
                        else if (iflags(h_IFLAGS_props).eq.2) then
                           print *, '@BRANCH: exact 2nd order ',
     &                              'properties are broken'
                           return
c YAU - xlambda will be in or out depending on who codes MBPT(2) an. Hessians
c                           call runit('xlambda')
                           call runit('xdens')
c                          AP - McKinley goes here
                           call runit('xvdint')
                           call runit('xcphf')
                           call runit('xsquint')
                           call runit('xanti')
                           call runit('xbcktrn')
                           call runit('xvdint')
                        end if
c                    else if (.not.raman) then
                     else
                        call runit('xlambda')
                        if (blCCSDT) then
                           iflags(h_IFLAGS_calc)=42
                           call a2putrec(20,'JOBARC','IFLAGS  ',100,
     &                                   iflags)
                           call runit('xlcct')
                        end if
                        if (iflags(h_IFLAGS_excite).gt.1) then
                           call runit('xvee')
                        end if
                        if (ndrgeo.ne.0) call runchg
                        call runit('xdens')
                        call runit('xanti')
                        call runit('xbcktrn')
                        call runit(der_integral_package)

c Modify (extrapolate) the gradients (not for SCF calcs)
c 01/2006, Ajith Perera. When POST-HF exact frequencies are availble
c the extrpolation has to be done differently.
                        if (bGExtrap.or.bCExtrap) then
                           call runit('xa2proc extrap gradient')
                        end if

c                    end if (raman)
                     end if
c                 else if (.not.analytical_gradient) then
                  else
                     if (iflags(h_IFLAGS_excite).gt.1) then
                        call runit('xlambda')
                        call runit('xvee')
                     end if
c                 end if (analytical_gradient)
                  end if
                  call a2getrec(1,'JOBARC','LASTGEOM',1,istat)
                  if (istat.eq.0) then
                     call c_gtod(is,ius)
                     print '(a,f10.1,a)',
     &                     'ACES2: Total elapsed time is ',
     &                     is+1.d-6*ius-dTimeStart,' seconds'
                     call a2putrec(1,'JOBARC','DIRTYFLG',1,0)
                     call runit(xjoda)
                  else
#ifdef _HAVE_MPI
                     if (nprocs.ne.1) call sync_fd
#endif
                     if (irank.eq.0) then
                        call a2putrec(1,'JOBARC','DIRTYFLG',1,0)
                        call runit(xjoda)
                     end if
                  end if
                  call rmfiles
               end do

c           end if (VIB == 1)
            end if

c        end if (VIB == 0)
         end if

c     end if (plain_scf)
      end if

c   o generate HyperChem input file after calculating vibrational frequencies
      if (bHyper) call runit(szHyper)
      if (bMolden) call runit(szMolden)

      return
c     end subroutine branch
      end

