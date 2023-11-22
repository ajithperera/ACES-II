      subroutine prep_mrcc(nbas, nbasx, eigvec, eigval, scr, mxcor,
     $    nirrep, nocc, nbfirr, iintfp)
c     
c     in this routine the various types of active orbitals are
c     determined that are needed in a MRCC calculation.
c     In the first turn they will be determined by MRCC, but
c     from then on they are determined by overlap criteria. 
c     This will allow the symmetry of a molecule to be lowered in
c     vibrational frequency calculations.
c     
c     The SCF orbitals are permuted such that the ordering
c     is appropriate (i.e. sequential) for use in MRCC.
c     Also the records CCSYM_V, CCSYM_O, IPSYM_A and EASYM_A
c     are written to JOBARC
c     
c     the information to determine active orbitals is obtained from
c     various projectors that are written by previous runs of MRCC.
c
c     mxcor = 3*nbasx*nbasx + 2*nbas    (or   nbas*nbas + 3*nbas)
c     
      implicit none
      integer maxirrep
      parameter (maxirrep=8)
c
      integer nbas, nbasx, mxcor, nirrep, nocc(maxirrep),
     $    nbfirr(maxirrep), iintfp, iunit, nrow, ncol, idone,
     $    perline
      double precision eigvec(nbas, nbas), eigval(nbas), scr(mxcor),
     $    xnorm, one, zilch, sdot
      logical print, yesno, skip
c
      integer nscr, im, ij, iscr, iact1(maxirrep), iact2(maxirrep),
     $    ioff, n, ip1, ip2, irp, iscr2, ioff0, ie, ib, ierror,
     $    nvrt(maxirrep), i1, i2, ione, i, iscr1, iscr3, nbas2
c
c First read information from file ACT_MRCC made in xmrcc
c The procedure is a little stupid for historical reasons
c It was not realized that xsymcor destroys the JOBARC file.
c
      one = 1.0d0
      ione = 1
      zilch = 0.0d0
      skip = .true.
      print = .false.
c
C     Open file ACT_MRCC for writing MOs in AO basis.
C
      inquire(file='ACT_MRCC',exist=yesno)
      iunit = 69
      if(yesno)then
        open(iunit,file='ACT_MRCC',status='old',access='sequential',
     &      form='formatted')
c
        perline = 4
c
c read active orbital energies from file ACT_MRCC and put on JOBARC
c        
        read(iunit, '(I12)') idone
        call putrec(-1, 'JOBARC',   'NUMACT_M', ione, idone)
        if (idone .ne. 0) then
          call get_matrix(scr, 1, idone, perline, iunit)
          call putrec(-1, 'JOBARC', 'EVAL_M  ',idone*iintfp, scr)
          call get_matrix(scr, 1, idone, perline, iunit)
          call putrec(-1, 'JOBARC', 'PROP_M  ',idone*iintfp, scr)
        endif
c
        read(iunit, '(I12)') idone
        call putrec(-1, 'JOBARC', 'NUMACT_J', ione, idone)
        if (idone .ne. 0) then
          call get_matrix(scr, 1, idone, perline, iunit)
          call putrec(-1, 'JOBARC', 'EVAL_J  ',idone*iintfp, scr)
          call get_matrix(scr, 1, idone, perline, iunit)
          call putrec(-1, 'JOBARC', 'PROP_J  ',idone*iintfp, scr)
        endif
c
        read(iunit, '(I12)') idone
        call putrec(-1, 'JOBARC', 'NUMACT_E', ione, idone)
        if (idone .ne. 0) then
          call get_matrix(scr, 1, idone, perline, iunit)
          call putrec(-1, 'JOBARC', 'EVAL_E  ',idone*iintfp, scr)
          call get_matrix(scr, 1, idone, perline, iunit)
          call putrec(-1, 'JOBARC', 'PROP_E  ',idone*iintfp, scr)
        endif
c
        read(iunit, '(I12)') idone
        call putrec(-1, 'JOBARC', 'NUMACT_B', ione, idone)
        if (idone .ne. 0) then
          call get_matrix(scr, 1, idone, perline, iunit)
          call putrec(-1, 'JOBARC', 'EVAL_B  ',idone*iintfp, scr)
          call get_matrix(scr, 1, idone, perline, iunit)
          call putrec(-1, 'JOBARC', 'PROP_B  ',idone*iintfp, scr)
        endif
c 
        if (.not. skip) then
c
c proces the four types of projector
c
          read(iunit, '(3I12)') idone, nrow, ncol
          call putrec(-1, 'JOBARC', 'OCCACT_M', ione, idone)
          if (idone .eq. 1) then
            call get_matrix(scr, nbasx, nbasx, perline, iunit)
            call putrec(-1, 'JOBARC', 'ZAO_MM  ',nbasx*nbasx*iintfp,
     $          scr)
          endif
c
          read(iunit, '(3I12)') idone, nrow, ncol
          call putrec(-1, 'JOBARC', 'OCCACT_J', ione, idone)
          if (idone .eq. 1) then
            call get_matrix(scr, nbasx, nbasx, perline, iunit)
            call putrec(-1, 'JOBARC', 'ZAO_JJ  ',nbasx*nbasx*iintfp,
     $          scr)
          endif
c
          read(iunit, '(3I12)') idone, nrow, ncol
          call putrec(-1, 'JOBARC', 'VRTACT_E', ione, idone)
          if (idone .eq. 1) then
            call get_matrix(scr, nbasx, nbasx, perline, iunit)
            call putrec(-1, 'JOBARC', 'ZAO_EE  ',nbasx*nbasx*iintfp,
     $          scr)
          endif
c
          read(iunit, '(3I12)') idone, nrow, ncol
          call putrec(-1, 'JOBARC', 'VRTACT_B', ione, idone)
          if (idone .eq. 1) then
            call get_matrix(scr, nbasx, nbasx, perline, iunit)
            call putrec(-1,'JOBARC', 'ZAO_BB  ',nbasx*nbasx*iintfp,
     $          scr)
          endif
c
        endif
c
        close(iunit, status='keep')
c
c also calculate expectation value of x2 + y2 + z2 for SCF MO eigenvectors
c and put on JOBARC
c
        iscr1 = 1
        iscr2 = iscr1 + nbas*nbas
        iscr3 = iscr2 + nbas*nbas
        nbas2 = (nbas * (nbas+1)) /2
c
        call getrec(20, 'JOBARC', '2NDMO_XX', nbas2*iintfp,
     &   scr(iscr1))
        call getrec(20, 'JOBARC', '2NDMO_YY', nbas2*iintfp,
     &   scr(iscr2))
        call SAXPY(nbas2, one, scr(iscr2), 1, scr(iscr1),1)
        call getrec(20, 'JOBARC', '2NDMO_ZZ', nbas2*iintfp,
     &   scr(iscr2))
        call SAXPY(nbas2, one, scr(iscr2), 1, scr(iscr1),1)
        call expnd2(scr(iscr1), scr(iscr2), nbas)
c
c transform to MO basis
c
        call getrec(20, 'JOBARC', 'SCFEVCA0', nbas*nbas*iintfp,
     $      scr(iscr1))
        call xgemm('N', 'N', nbas, nbas, nbas, one, scr(iscr2), nbas,
     $      scr(iscr1), nbas, zilch, scr(iscr3), nbas)
        call xgemm('T', 'N', nbas, nbas, nbas, one, scr(iscr1), nbas,
     $      scr(iscr3), nbas, zilch, scr(iscr2), nbas)
c
        ioff = iscr2
        do i = 1, nbas
          scr(i) = scr(ioff)
          ioff = ioff + nbas + 1
        enddo
c
        call putrec(20, 'JOBARC', 'SCFPROP0', nbas*iintfp, scr)
c
      else
        write(6,*) ' File MRCC_ACT does not exist'
        idone = 0
        call putrec(-1, 'JOBARC', 'NUMACT_M', ione, idone)
        call putrec(-1, 'JOBARC', 'NUMACT_J', ione, idone)
        call putrec(-1, 'JOBARC', 'NUMACT_E', ione, idone)
        call putrec(-1, 'JOBARC', 'NUMACT_B', ione, idone)
c
        if (.not. skip) then
          call putrec(-1, 'JOBARC', 'OCCACT_M', ione, idone)
          call putrec(-1, 'JOBARC', 'OCCACT_J', ione, idone)
          call putrec(-1, 'JOBARC', 'VRTACT_E', ione, idone)
          call putrec(-1, 'JOBARC', 'VRTACT_B', ione, idone)
        endif
      endif
c
      call getrec(20, 'JOBARC', 'SCFEVCA0', nbas*nbas*iintfp, eigvec)
      call getrec(20, 'JOBARC', 'SCFEVLA0', nbas*iintfp, eigval)
c
      if (print) then
        write(6,*) ' @prep_mrcc : original eigenvalues'
        call output(eigval, 1, 1, 1, nbas, 1, nbas, 1)
c
c check orthogonality of eigenvectors
c
        iscr1 = 1
        iscr2 = iscr1 + nbas*nbas
        call getrec(20, 'JOBARC', 'AOOVRLAP', nbas*nbas*iintfp,
     $      scr(iscr1))
        call xgemm('N', 'N', nbas, nbas, nbas, one, scr(iscr1), nbas,
     $      eigvec, nbas, zilch, scr(iscr2), nbas)
        call xgemm('T', 'N', nbas, nbas, nbas, one, eigvec, nbas,
     $      scr(iscr2), nbas, zilch, scr(iscr1), nbas)
        call zero(scr(iscr2), nbas*nbas)
        ioff = iscr2
        do i = 1, nbas
          scr(ioff) = one
          ioff = ioff + nbas + 1
        enddo
        call saxpy(nbas*nbas, -one, scr(iscr1), 1, scr(iscr2), 1)
        xnorm = sdot(nbas*nbas, scr(iscr2), 1, scr(iscr2), 1)
        write(6,*) ' Difference from orthonormality ', xnorm
        if (xnorm .gt. 1.0d-5) then
          call output(scr(iscr1), 1, nbas, 1, nbas, nbas, nbas, 1)
        endif
      endif
c
      if (yesno) then
        do irp = 1, nirrep
          nvrt(irp) = nbfirr(irp) - nocc(irp)
        enddo
c
c     first treat the occupied orbitals.
c     
        ione = 1
cmn      call getrec(-1, 'JOBARC','OCCACT_J', ione, i1)
cmn      call getrec(-1, 'JOBARC','OCCACT_M', ione, i2)
        call getrec(-1, 'JOBARC','NUMACT_M', ione, i1)
        call getrec(-1, 'JOBARC','NUMACT_J', ione, i2)
        if (i1 .ne. 0 .or. i2 .ne. 0) then
c
          nscr = 3*nbasx*nbasx
cmn        call mrccact(scr, nscr, 'mm', nbas, nbasx)
          call mrcc_eval(scr, nscr, 'mm', nbas, iintfp)
          im = 1
          ij = im + nbas
          iscr = ij + nbas
          call mrcc_eval(scr(ij), nscr, 'jj', nbas, iintfp)
cmn        call mrccact(scr(ij), nscr, 'jj', nbas, nbasx)
c
          if (print) then
            write(6,*) ' overall occupied projection'
            write(6,*) ' nbfirr'
            write(6,990) (nbfirr(i), i=1,nirrep)
 990        format(8i4)
            write(6,*) ' active m'
            call output(scr(im), 1, 1, 1, 1, nbas, 1, 1)
            write(6,*) ' active j'
            call output(scr(ij), 1, 1, 1, nbas, 1, nbas, 1)
          endif
c     
c     the projected norms are contained in scr(im) and scr(ij)
c     
          call izero(maxirrep, iact1)
          call izero(maxirrep, iact2)
          ioff = 0
          do irp = 1, nirrep
            n = nocc(irp)
            ip1 = iscr + nbas*n
            ip2 = ip1 + n
            iscr2 = ip2 + n
            call ordermrcc(scr(im+ioff), scr(ij+ioff), n, nbas, 
     $          eigvec(1,1+ioff), scr(iscr), scr(ip1),
     $          scr(ip2), .true., ierror, iact1(irp), iact2(irp),
     $          eigval(1+ioff), scr(iscr2))
            if (ierror .eq. 1) then
              write(6,*) '@prep_mrcc: something wrong occupieds'
              call output(scr(ip1), 1, n, 1, n, n, n, 1)
              call output(scr(ip2), 1, n, 1, n, n, n, 1)
              call errex
            endif
            ioff = ioff + nocc(irp)+nvrt(irp)
          enddo
c
c occupied orbitals are done. Write iact1 and iact2 to JOBARC
c 
          write(6,*) '  @prep_mrcc: occupied orbs'
          write(6, 999) (iact1(i), i=1,nirrep)
          write(6,998) (iact2(i), i=1,nirrep)
 999      format('  IPSYM_A: ', 8i4)
 998      format('  CCSYM_O: ', 8i4)
          call putrec(20, 'JOBARC', 'IPSYM_A ', maxirrep, iact1)
          call putrec(20, 'JOBARC', 'CCSYM_O ', maxirrep, iact2)
c
        else
          write(6,*) ' @prep_mrcc: occupied orbitals not treated'
        endif
c
c Now treat virtuals in similar fashion.
c
cmn      call getrec(-1, 'JOBARC','VRTACT_E', ione, i1)
cmn      call getrec(-1, 'JOBARC','VRTACT_B', ione, i2)
c
        call getrec(-1, 'JOBARC','NUMACT_E', ione, i1)
        call getrec(-1, 'JOBARC','NUMACT_B', ione, i2)
        if (i1 .ne. 0 .or. i2 .ne. 0) then
c
          nscr = 3*nbasx*nbasx
cmn        call mrccact(scr, nscr, 'ee', nbas, nbasx)
          call mrcc_eval(scr, nscr, 'ee', nbas, iintfp)
          ie = 1
          ib = ie + nbas
          iscr = ib + nbas
          call mrcc_eval(scr(ib), nscr, 'bb', nbas, iintfp)
cmn        call mrccact(scr(ib), nscr, 'bb', nbas, nbasx)
c
          if (print) then
            write(6,*) ' overall virtual projection'
            write(6,990) (nbfirr(i), i=1,nirrep)
            write(6,*) ' active e'
            call output(scr(ie), 1, 1, 1, nbas, 1, nbas, 1)
            write(6,*) ' active b'
            call output(scr(ib), 1, 1, 1, nbas, 1, nbas, 1)
          endif
c     
c     the projected norms are contained in scr(ie) and scr(ib)
c     
          call izero(maxirrep, iact1)
          call izero(maxirrep, iact2)
          ioff = 0
          do irp = 1, nirrep
            n = nvrt(irp)
            ip1 = iscr + nbas*n
            ip2 = ip1 + n
            iscr2 = ip2+n
            ioff = ioff + nocc(irp)
            call ordermrcc(scr(ie+ioff), scr(ib+ioff), n, nbas, 
     $          eigvec(1,1+ioff), scr(iscr), scr(ip1),
     $          scr(ip2), .false., ierror, iact1(irp), iact2(irp),
     $          eigval(1+ioff), scr(iscr2))
            if (ierror .eq. 1) then
              write(6,*) '@perp_mrcc: something wrong occupieds'
              call output(scr(ip1), 1, n, 1, n, n, n, 1)
              call output(scr(ip2), 1, n, 1, n, n, n, 1)
              call errex
            endif
            ioff = ioff + n
          enddo
c
c virtual orbitals are done. Write iact1 and iact2 to JOBARC
c 
          write(6,*) '  @prep_mrcc: virtual orbs'
          write(6, 997) (iact1(i), i=1,nirrep)
          write(6,996) (iact2(i), i=1,nirrep)
 997      format('  EASYM_A: ', 8i4)
 996      format('  CCSYM_V: ', 8i4)
          call putrec(20, 'JOBARC', 'EASYM_A ', maxirrep, iact1)
          call putrec(20, 'JOBARC', 'CCSYM_V ', maxirrep, iact2)
c
        else
c
          write(6,*) ' @prep_mrcc: virtual orbitals not treated'
c
        endif
c
        if (print) then
          write(6,*) ' @prep_mrcc : reordered eigenvalues'
          call output(eigval, 1, 1, 1, nbas, 1, nbas, 1)
c
c check orthogonality of eigenvectors
c
          iscr1 = 1
          iscr2 = iscr1 + nbas*nbas
          call getrec(20, 'JOBARC', 'AOOVRLAP', nbas*nbas*iintfp,
     $        scr(iscr1))
          call xgemm('N', 'N', nbas, nbas, nbas, one, scr(iscr1), nbas,
     $        eigvec, nbas, zilch, scr(iscr2), nbas)
          call xgemm('T', 'N', nbas, nbas, nbas, one, eigvec, nbas,
     $        scr(iscr2), nbas, zilch, scr(iscr1), nbas)
          call zero(scr(iscr2), nbas*nbas)
          ioff = iscr2
          do i = 1, nbas
            scr(ioff) = one
            ioff = ioff + nbas + 1
          enddo
          call saxpy(nbas*nbas, -one, scr(iscr1), 1, scr(iscr2), 1)
          xnorm = sdot(nbas*nbas, scr(iscr2), 1, scr(iscr2), 1)
          write(6,*) ' Difference from orthonormality ', xnorm
          if (xnorm .gt. 1.0d-5) then
            call output(scr(iscr1), 1, nbas, 1, nbas, nbas, nbas, 1)
          endif
        endif
c
c also write orbitals to JOBARC
c
        call putrec(20, 'JOBARC', 'SCFEVCA0', nbas*nbas*iintfp,
     $      eigvec)
        call putrec(20, 'JOBARC', 'SCFEVCB0', nbas*nbas*iintfp,
     $      eigvec)
c
c and the eigenvalues...
c
        call putrec(20, 'JOBARC', 'SCFEVLA0', nbas*iintfp,
     $      eigval)
        call putrec(20, 'JOBARC', 'SCFEVLB0', nbas*iintfp,
     $      eigval)
c
      endif
      return
      end

c     
c Local Variables: c
c mode: fortran c
c fortran-do-indent: 2 c
c fortran-if-indent: 2 c
c fortran-continuation-indent: 4 c
c fortran-comment-indent-style: nil c
c End: c
