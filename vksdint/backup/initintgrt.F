      subroutine initintgrt

c This program initializes parameters from the INTGRT namelist.

#include <aces.h>
#include <io.par>
#include <const.par>
#include <intgrtflags.com>
#include <dftfunc.com>

      integer err,ptr,num,i,cutoff

      character*20 list(100),ele
      logical printdef, streq
      M_REAL
     &    coef,thresh_power

      printdef=.true.
      if (int_ks .and. .not.int_ks_finaliter) printdef=.false.

      call callstack_push('INITINTGRT')
      call nl_init('INTGRT',err,printdef)

c Namelist entries:
c --------------------------------------------------------------------------
c PRINT_NL   : whether to print the values read in in the namelist
c              Default: true
c --------------------------------------------------------------------------
c Print flags
c
c              All print flags have one of the following values:
c                never   - never print out the property
c                dft     - only print out the property during a DFT calculation
c                ks      - print it out durint the final iteration of KS
c                          or during a DFT calculation
c                always  - print it out during all KS iterations as well
c                          as during a DFT calculation
c
c PRINT_SCF  : whether or not to print a summary of SCF (or correlated in
c              the case of a relaxed density) results
c              Default: dft
c PRINT_INT  : whether to print out a summary of integration parameters
c              Default: dft
c PRINT_SIZE : whether to print relative size information
c              Default: never
c PRINT_ATOM : whether to print out integration information for each atom
c              Default: never
c
c PRINT_MOS  : whether to print the MOs from the SCF or correlated calculation
c              Default: never
c PRINT_OCC  : whether to print the occupancy of each MO
c              Default: never
c --------------------------------------------------------------------------
c Integration flags
c
c POTRADPTS  : the number of radial points to use in the integration around
c              each atomic center in KS calc (default: 50)
c ENERADPTS  : the number of radial points to use in the integration around
c              each atomic center in post KS calc (default: 50)
c RADTYP     : the type of radial integration available are:
c                Handy           : (default)
c                Gauss-Legendre  :
c PARTPOLY   : partition polyhedra
c                equal     : equally regardless of atomic center
c                bsrad     : by Bragg-Slater radii (default)
c                dynamic   : dynamic (minimum in density) sizing
c RADSCAL    : how to do radial scaling
c                none      : no radial scaling
c                slater    : use Slater's rules for radial scaling (default)
c PARTTYP    : type of partitioning
c                rigid     : rigid partitioning
c                fuzzy     : fuzzy partitioning (default)
c FUZZYITER  : the number of iterations to use in fuzzy partitioning
c              (default: 4)
c RADLIMIT   : the length of the radial integration (default 3.00)
c
c ENEGRID    : the angular grid to use by default at each radial point
c              in an energy calculation (default: 4)
c POTGRID    : the angular grid to use by default at each radial point
c              in a Kohn-Sham calculation of the potential (default: 4)
c ENETYPE    : the angular grid type to use by default at each radial point
c              in an energy calculation
c                lebedev    : default
c POTTYPE    : the angular grid type to use by default at each radial point
c              in a Kohn-Sham calculation of the potential
c                lebedev    : default
c --------------------------------------------------------------------------
c Functionals
c
c FUNC       : The type of functional (and potential) to use in any
c              calculation other than simple HFDFT.  In the case of HFDFT,
c              all energy functionals are calculated.  In the case of
c              Kohn-Sham, only one potential will be calculated.  In the
c              case of an energy search, only one energy is stored in the
c              JOBARC file.
c
c              This keyword is of one of the forms:
c                FUNC=func
c                FUNC=exch,corr
c                FUNC=func,coef,func,coef [,func,coef]*
c              where func is the name of one of the functionals below and
c              coef is the coefficient to use.  In the final form, at least
c              two functionals must be included in the list.
c
c              Recognized correlation functionals are:
c                vwn     : Vosko-Wilk-Nusair, used as VWN5 in Gaussian94 
c                lyp     : Lee-Yang-Parr
c              pbe_cor   : Perdew-Burke-Ernzerhof
c              pw91_cor  : Perdew-Wang 91
c
c                wl      : Wilson-Levy
c                wi      : Wilson-Ivanov
c                wl and wi do not have their corresponding potentials         
c                wi has form for open-shell species only
c
c              Exchange functionals are:
c                lda     : LDA (Slater) exchange
c                becke   : Becke exchange
c               pbe_ex   : Perdew-Burke-Ernzerhof
c               pw91_ex  : Perdew-Wang 91
c                hf      : Hartree-Fock-like exchange
c
c              Hybrid functionals are:
c                b3lyp   : Becke III LYP
c
c              A hybrid functional may only be entered as FUNC=hybrid.
c
c              The default functional is none.
c========================================================================
c                         More Flags                                    =
c========================================================================
c
c     numacc         Numerical accelerator. Default: on
c
c     cutoff         Cutoff value for the numerical integration.
c                    Default: 12 Densities < 10**(-12) are assumed 0.
c     exact_ex       Exact LOCAL EXCHANGE. Default: false
c     tdks           Time-dependnet KS is requested. Default: false 

      list(1)='never'
      list(2)='dft'
      list(3)='ks'
      list(4)='always'
      call nl_ele('print_scf' ,list,4,int_prt_dft  ,int_printscf)
      call nl_ele('print_int' ,list,4,int_prt_dft  ,int_printint)
      call nl_ele('print_atom',list,4,int_prt_never,int_printatom)
      call nl_ele('print_size',list,4,int_prt_never,int_printsize)
      call nl_ele('print_mos' ,list,4,int_prt_never,int_printmos)
      call nl_ele('print_occ' ,list,4,int_prt_never,int_printocc)

      call nl_int('eneradpts',50,int_numradpts)
      call nl_int('potradpts',50,potradpts)

      list(1)='Handy'
      list(2)='Gauss-Legendre'
      call nl_ele('radtyp',list,2,int_radtyp_handy,int_radtyp)

      list(1)='equal'
      list(2)='bsrad'
      list(3)='dynamic'
      call nl_ele('partpoly',list,3,int_partpoly_bsrad,int_partpoly)

      list(1)='none'
      list(2)='Slater'
      call nl_ele('radscal',list,2,int_radscal_slater,int_radscal)

      list(1)='rigid'
      list(2)='fuzzy'
      call nl_ele('parttyp',list,2,int_parttyp_fuzzy,int_parttyp)

      call nl_int('fuzzyiter',4,int_fuzzyiter)

      call nl_real('radlimit',three,int_radlimit)
      call nl_log('numacc',.true.,num_acc_ks)
      call nl_log('exact_ex',.false.,ks_exact_ex)
      call nl_log('tdks',.false.,int_tdks)

      call nl_int('enegrid',4,int_defenegrid)
      call nl_int('potgrid',1,int_defpotgrid)
      list(1)='lebedev'
      call nl_ele('enetype',list,1,int_gridtype_leb,int_defenetype)
      call nl_ele('pottype',list,1,int_gridtype_leb,int_defpottype)

c Parse the functional to use:

      do i=1,fun_num_exch
        coef_exch(i)=zero
      end do
      do i=1,fun_num_corr
        coef_corr(i)=zero
      end do

      nam_exch(fun_exch_lda)  = 'LDA (Slater, Xalpha)'
      nam_exch(fun_exch_becke)= 'Becke'
      nam_exch(fun_exch_pbe)= 'Perdew-Burke-Ernzerhof'
      nam_exch(fun_exch_pw91)= 'Perdew-Wang 91'
      nam_exch(fun_exch_hf)   = 'Exact Nonlocal Exchange'

      nam_corr(fun_corr_vwn)  = 'Vosko-Wilk-Nusair'
      nam_corr(fun_corr_lyp)  = 'Lee-Yang-Parr'
      nam_corr(fun_corr_pbe)  = 'Perdew-Burke-Ernzerhof'
      nam_corr(fun_corr_pw91) = 'Perdew-Wang 91'
      nam_corr(fun_corr_wl)   = 'Wilson-Levy'
      nam_corr(fun_corr_wi)   = 'Wilson-Ivanov'
 
      nam_hyb(fun_hyb_b3lyp)  = 'Becke III LYP'

      abb_exch(fun_exch_lda)  = 'lda'
      abb_exch(fun_exch_becke)= 'becke'
      abb_exch(fun_exch_pbe)  = 'pbe_ex'
      abb_exch(fun_exch_pw91) = 'pw91_ex'
      abb_exch(fun_exch_hf)   = 'hf'

      abb_corr(fun_corr_vwn)  = 'vwn'
      abb_corr(fun_corr_lyp)  = 'lyp'
      abb_corr(fun_corr_pbe)  = 'pbe_cor'
      abb_corr(fun_corr_pw91) = 'pw91_cor'
      abb_corr(fun_corr_wl)   = 'wl'
      abb_corr(fun_corr_wi)   = 'wi'

      abb_hyb(fun_hyb_b3lyp) = 'b3lyp'

      call nl_str('func','none',nam_func)

      ptr=0
      call list_next_ele(nam_func,ptr,ele,err)

      if (err.ne.0) then
        write (stdout,910)
        call errex
        stop
      endif

      if (streq(ele,'none',.true.)) then
         int_dft_fun=fun_dft_none
         go to 100
      else
         int_dft_fun=fun_special
      end if  

      call list_next_num(nam_func,ptr,coef,err)

  910 format ('ERROR: invalid functional.')
  920 format ('ERROR: functional ',A10,'name used twice.')
  921 format ('ERROR: functional name unrecognized.')
  922 format ('ERROR: invalid functional coefficient.')
  923 format ('ERROR: 2 or more functionals required.')

c If we're working with the form: FUNC=func,coef,func,coef [,func,coef]*

      if (err.eq.0) then
        num=0
   10   continue

c       At this point, ele contains the name of the functional and coef
c       contains it's coefficient.

c       Check exchange functionals...

        do i=1,fun_num_exch
          if (streq(ele,abb_exch(i),.true.)) then
            if (coef_exch(i).ne.zero) then
              write (stdout,920) abb_exch(i)
              call errex
              stop
            endif
            coef_exch(i)=coef
            num=num+1
            goto 20
          endif
        end do

c       Check correlation functional

        do i=1,fun_num_corr
          if (streq(ele,abb_corr(i),.true.)) then
            if (coef_corr(i).ne.zero) then
              write (stdout,920) abb_corr(i)
              call errex
              stop
            endif
            coef_corr(i)=coef
            num=num+1
            goto 20
          endif
        end do

c       If neither, die.

        write (stdout,921)
        call errex
        stop

   20   continue

c       Get the next functional and coefficient.  If they don't exist,
c       make sure there have been at least 2.

        call list_next_ele(nam_func,ptr,ele,err)

        if (err.eq.0) then
          call list_next_num(nam_func,ptr,coef,err)
          if (err.eq.0) goto 10
          write (stdout,922)
          call errex
          stop
        else if (num.lt.2) then
          write (stdout,923)
          call errex
          stop
        endif

c We're working with the form FUNC=func or FUNC=corr,exch or FUNC=hybrid

      else
        num=0
c       num=0  : nothing entered
c           1  : exchange entered
c           2  : correlation entered

        do i=1,fun_num_exch
          if (streq(ele,abb_exch(i),.true.)) then
            coef_exch(i)=one
            num=num+1
            goto 60
          endif
        end do

        do i=1,fun_num_corr
          if (streq(ele,abb_corr(i),.true.)) then
            coef_corr(i)=one
            num=num+2
            goto 60
          endif
        end do

   60   continue
        if (num.eq.0) go to 70
         call list_next_ele(nam_func,ptr,ele,err)
          if (err.ne.zero) go to 100

        do i=1,fun_num_exch
          if (streq(ele,abb_exch(i),.true.)) then
             if (num.eq.2) then
               coef_exch(i)=one
               goto 100
             else
               write (stdout,910)
               call errex
               stop
             end if 
          end if
        end do


        do i=1,fun_num_corr
          if (streq(ele,abb_corr(i),.true.)) then
             if (num.eq.1) then
               coef_corr(i)=one
               goto 100
             else
               write (stdout,910)
               call errex
               stop
             end if
          end if
        end do

           write (stdout,910)
           call errex
           stop

   70   continue

c       Hybrid functionals...

        if (streq(ele,abb_hyb(fun_hyb_b3lyp),.true.)) then
           int_dft_fun=fun_hyb_b3lyp
           coef_exch(fun_exch_hf)=b3lypa
           coef_exch(fun_exch_lda)=(one-b3lypa-b3lypb)
           coef_exch(fun_exch_becke)=b3lypb
           coef_corr(fun_corr_vwn)=one-b3lypc
           coef_corr(fun_corr_lyp)=b3lypc
        else
          write (stdout,910)
          call errex
          stop
        endif
       endif
  100 continue



c KSPOT      : Availiable potentials : 
c                vwn     : Vosko-Wilk-Nusair 
c                lyp     : Lee-Yang-Parr
c
c                lda     : LDA (Slater) exchange 
c                becke   : Becke exchange
c                hf      : Nonlocal HF-like (default)
c
c          The format and the sructure are the same as those 
c          for the key word "func".

c     Initialize the coefficients needed in hybrid schemes.
c     Note that hybrid scheme is a more general definition and the KS
c     method with multiplicative potentials is one of the possible choices.

      do i=1,fun_num_exch
        coef_pot_exch(i)=zero
      end do
      do i=1,fun_num_corr
        coef_pot_corr(i)=zero
      end do
 
      call nl_str('kspot','lda',nam_kspot)
      
      ptr=0
      call list_next_ele(nam_kspot,ptr,ele,err)
      if (err.ne.0) then
        write (stdout,1010)
        call errex
        stop
      endif
      int_kspot=fun_special
      call list_next_num(nam_kspot,ptr,coef,err)

 1010 format ('ERROR: invalid potential.')
 1020 format ('ERROR: potential ', A7, ' name used twice.')
 1021 format ('ERROR: potential name unrecognized.')
 1022 format ('ERROR: invalid potential coefficient.')
 1023 format ('ERROR: 2 or more potentials required.')
 1024 format ('ERROR: potential not implemented.')
c If we're working with the form: kspot=pot,coef,pot,coef [,pot,coef]*
c We could have any combination of exch and corr. 
c Each pot appears only once.

      if (err.eq.0) then
        num=0
  110   continue

c       At this point, ele contains the name of the potential and coef
c       contains it's coefficient.

c       Check exchange potentials...

        do i=1,fun_num_exch
          if (streq(ele,abb_exch(i),.true.)) then
            if (coef_pot_exch(i).ne.zero) then
              write (stdout,1020) abb_exch(i)
              call errex
              stop
            endif
            coef_pot_exch(i)=coef
            num=num+1
            goto 120
          endif
        end do

c       Check correlation potentials...

       do i=1,fun_num_corr
          if (streq(ele,abb_corr(i),.true.)) then
            if (coef_pot_corr(i).ne.zero) then
              write (stdout,1020) abb_corr(i)
              call errex
              stop
            endif
            coef_pot_corr(i)=coef
            num=num+1
            goto 120
          endif
        end do

c       If neither, die.

        write (stdout,1021)
        call errex
        stop

  120   continue

c       Get the next potential and coefficient.  If they don't exist,
c       make sure there have been at least 2.

        call list_next_ele(nam_kspot,ptr,ele,err)

        if (err.eq.0) then
          call list_next_num(nam_kspot,ptr,coef,err)
          if (err.eq.0) goto 110
          write (stdout,1022)
          call errex
          stop
        else if (num.lt.2) then
          write (stdout,1023)
          call errex
          stop
        endif

c We're working with the form kspot=pot, (exch or corr),
c or kspot=exch,corr (or corr,exch), or kspot=standard hybrid

      else
        num=0

c       num=0  : nothing entered
c           1  : exchange entered
c           2  : correlation entered

        do i=1,fun_num_exch
          if (streq(ele,abb_exch(i),.true.)) then
            coef_pot_exch(i)=one
            num=num+1
            goto 160
          endif
        end do

        do i=1,fun_num_corr
          if (streq(ele,abb_corr(i),.true.)) then
            coef_pot_corr(i)=one
            num=num+2
            goto 160
          endif
        end do

  160   continue
        if (num.eq.0) go to 170
         call list_next_ele(nam_kspot,ptr,ele,err)
         if (err.ne.0) go to 200

        do i=1,fun_num_exch
          if (streq(ele,abb_exch(i),.true.)) then
             if (num.eq.2) then
               coef_pot_exch(i)=one
               goto 200
             else
               write (stdout,1010)
               call errex
               stop
             end if
          end if
        end do


        do i=1,fun_num_corr
          if (streq(ele,abb_corr(i),.true.)) then
             if (num.eq.1) then
               coef_pot_corr(i)=one
               goto 200
             else
               write (stdout,1010)
               call errex
               stop
             end if
          end if
        end do

        write (stdout,1010)
        call errex
        stop 

  170   continue

c       Hybrid functionals...

        if (streq(ele,abb_hyb(fun_hyb_b3lyp),.true.)) then
           int_kspot=fun_hyb_b3lyp
           coef_pot_exch(fun_exch_hf)=b3lypa
           coef_pot_exch(fun_exch_lda)=(one-b3lypa-b3lypb)
           coef_pot_exch(fun_exch_becke)=b3lypb
           coef_pot_corr(fun_corr_vwn)=one-b3lypc
           coef_pot_corr(fun_corr_lyp)=b3lypc
        else
          write (stdout,1010)
          call errex
          stop
        endif
       endif

  200 continue

c Save the value of coef_pot_exch(fun_exch_hf) also as
c coef_pot_nonlocal. Needed for Vscf.

        if (ks_exact_ex) then 
           coef_pot_nonlocal=one
        else
           coef_pot_nonlocal=coef_pot_exch(fun_exch_hf)
        end if

c If kspot=wl or wi 
        if ((coef_pot_corr(fun_corr_wl).ne.zero).or.
     &    (coef_pot_corr(fun_corr_wi).ne.zero)) then
          write (stdout,1024)
          call errex
          stop
        endif

      call nl_int('cutoff',12,cutoff)
      thresh_power=dble(cutoff)
      thresh_power=dabs(thresh_power)

       if (thresh_power.lt.8) then
         write (stdout,5010)
      end if

 5010 format(' !!! WARNING !!! LOW NUMERICAL ACCURACY! ')

      thresh=1.0d0/(10.0d0**thresh_power) 
      call nl_term
      call callstack_pop
      return
      end

c Local Variables: c
c mode: fortran c
c fortran-do-indent: 2 c
c fortran-if-indent: 2 c
c fortran-continuation-indent: 4 c
c fortran-comment-indent-style: nil c
c End: c
